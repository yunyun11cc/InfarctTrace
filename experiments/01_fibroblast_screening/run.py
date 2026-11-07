"""Exploratory D7-vs-uninjured fibroblast screening; see the frozen protocol."""
import csv
import gzip
import io
import json
import re
import tarfile
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, ttest_ind
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/01_fibroblast_screening'
OUT.mkdir(parents=True, exist_ok=True)


def save_table(frame, name):
    frame.to_csv(OUT / name, sep='\t', index=True, index_label=frame.index.name or 'gene_id')


def annotation():
    genes, transcripts = {}, {}
    source = ROOT / 'datasets/annotation/Mus_musculus.GRCm38.87.gtf.gz'
    with gzip.open(source, 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            if fields[2] not in ('gene', 'transcript'):
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', fields[8]))
            gid = attrs['gene_id'].split('.')[0]
            if fields[2] == 'gene':
                genes[gid] = dict(symbol=attrs.get('gene_name', gid),
                                  biotype=attrs.get('gene_biotype', 'unknown'))
            else:
                tid = attrs['transcript_id'].split('.')[0]
                if tid in transcripts and transcripts[tid] != gid:
                    raise ValueError(f'Ambiguous annotation transcript: {tid}')
                transcripts[tid] = gid
    gene_info = pd.DataFrame.from_dict(genes, orient='index')
    save_table(gene_info, 'gene_annotation.tsv.gz')
    return gene_info, transcripts


def load_fb():
    meta = pd.read_csv(ROOT / 'datasets/GSE252027/metadata/samples_FB.csv').set_index('gsm')
    matrices = {name: {} for name in ('TPM', 'NumReads', 'EffectiveLength')}
    qc, excluded = [], None
    with tarfile.open(ROOT / 'datasets/GSE252027/original/GSE252027_RAW.tar') as archive:
        entries = {m.name.split('_')[0]: m for m in archive.getmembers() if m.isfile()}
        for gsm, row in meta.iterrows():
            with gzip.open(archive.extractfile(entries[gsm]), 'rt') as stream:
                frame = pd.read_csv(stream, sep='\t').set_index('Name')
            if not np.isfinite(frame.to_numpy()).all() or (frame.to_numpy() < 0).any():
                raise ValueError(f'Invalid quantification: {gsm}')
            gene_rows = frame.index.str.match(r'^ENSMUSG\d+(\.\d+)?$')
            if excluded is None:
                excluded = frame.loc[~gene_rows].copy()
                excluded.to_csv(OUT / 'excluded_non_gene_rows_example.tsv', sep='\t')
                frame.head(12).to_csv(OUT / 'raw_quant_first12rows.tsv', sep='\t')
            selected = frame.loc[gene_rows].copy()
            selected.index = selected.index.str.replace(r'\.\d+$', '', regex=True)
            if not selected.index.is_unique:
                raise ValueError('Version stripping produced duplicate gene IDs')
            if (selected.EffectiveLength <= 0).any():
                raise ValueError('Nonpositive effective length')
            for name in matrices:
                matrices[name][gsm] = selected[name]
            qc.append(dict(gsm=gsm, day=int(row.day), replicate=int(row.replicate),
                           feature_rows=len(frame), gene_rows=int(gene_rows.sum()),
                           excluded_non_gene_rows=int((~gene_rows).sum()),
                           estimated_reads_all=float(frame.NumReads.sum()),
                           estimated_reads_retained=float(selected.NumReads.sum()),
                           tpm_sum_all=float(frame.TPM.sum()),
                           genes_tpm_ge1=int((selected.TPM >= 1).sum())))
    matrices = {key: pd.DataFrame(value) for key, value in matrices.items()}
    if any(frame.isna().any().any() for frame in matrices.values()):
        raise ValueError('Gene sets differ between FB libraries')
    meta.to_csv(OUT / 'discovery_samples.csv')
    pd.DataFrame(qc).to_csv(OUT / 'discovery_qc.csv', index=False)
    for name, frame in matrices.items():
        save_table(frame, f'FB_{name}.tsv.gz')
    return meta, matrices


def length_scaled_counts(tpm, effective_length, estimated_reads):
    """The formula in tximport::makeCountsFromAbundance, lengthScaledTPM."""
    unscaled = tpm.mul(effective_length.mean(axis=1), axis=0)
    counts = unscaled.mul(estimated_reads.sum(axis=0) / unscaled.sum(axis=0), axis=1)
    if not np.allclose(counts.sum(), estimated_reads.sum(), rtol=1e-10):
        raise ValueError('Length scaling changed library totals')
    return counts.round().astype('int64')


def discovery(meta, matrices, gene_info):
    selected = meta.index[meta.day.isin([0, 7])]
    counts = length_scaled_counts(*(matrices[k][selected] for k in
                                   ('TPM', 'EffectiveLength', 'NumReads')))
    save_table(counts, 'D0_D7_lengthScaledTPM_counts_all_genes.tsv.gz')
    keep = (counts >= 10).sum(axis=1) >= 4
    design = pd.DataFrame({'condition': np.where(meta.loc[selected, 'day'] == 0, 'D0', 'D7')}, index=selected)
    inference = DefaultInference(n_cpus=2)
    dds = DeseqDataSet(counts=counts.loc[keep].T, metadata=design, design='~condition',
                       refit_cooks=True, inference=inference, quiet=False)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=['condition', 'D7', 'D0'], alpha=0.05,
                      inference=inference, quiet=True)
    stats.summary()
    result = stats.results_df.join(gene_info)
    result['discovery_pass'] = (result.padj < 0.05) & (result.log2FoldChange.abs() >= 1)
    save_table(result, 'discovery_D7_vs_D0.tsv')
    normalized = pd.DataFrame(dds.layers['normed_counts'].T, index=counts.index[keep], columns=selected)
    save_table(normalized, 'D0_D7_normalized_counts.tsv.gz')
    return result, int(keep.sum())


def validation(tx2gene):
    annotations, audit = {}, dict(platform_rows=0, excluded_category=0,
                                  unmapped=0, ambiguous_multigene=0, mapped_probe_clusters=0)
    with gzip.open(ROOT / 'datasets/GSE111059/metadata/GSE111059_family.soft.gz', 'rt') as stream:
        for line in stream:
            if line.startswith('!platform_table_begin'):
                header = next(stream).rstrip('\n').split('\t')
                break
        for line in stream:
            if line.startswith('!platform_table_end'):
                break
            values = line.rstrip('\n').split('\t')
            record = dict(zip(header, values))
            audit['platform_rows'] += 1
            if record['category'] != 'main':
                audit['excluded_category'] += 1
                continue
            tids = set(re.findall(r'ENSMUST\d+', record['mrna_assignment']))
            gids = {tx2gene[t] for t in tids if t in tx2gene}
            if len(gids) == 1:
                annotations[record['ID']] = next(iter(gids))
                audit['mapped_probe_clusters'] += 1
            elif len(gids) > 1:
                audit['ambiguous_multigene'] += 1
            else:
                audit['unmapped'] += 1
    with gzip.open(ROOT / 'datasets/GSE111059/original/GSE111059_series_matrix.txt.gz', 'rt') as stream:
        for line in stream:
            if line.startswith('!series_matrix_table_begin'):
                break
        lines = []
        for line in stream:
            if line.startswith('!series_matrix_table_end'):
                break
            lines.append(line)
    probe_matrix = pd.read_csv(io.StringIO(''.join(lines)), sep='\t', index_col=0)
    probe_matrix.index = probe_matrix.index.astype(str)
    pd.Series(annotations, name='gene_id').to_csv(OUT / 'validation_probe_mapping.tsv', sep='\t', index_label='probe_id')
    mapped = probe_matrix.loc[probe_matrix.index.isin(annotations)].copy()
    mapped['gene_id'] = mapped.index.map(annotations)
    matrix = mapped.groupby('gene_id', sort=True).median()
    meta = pd.read_csv(ROOT / 'datasets/GSE111059/metadata/samples.csv').set_index('geo_accession')
    if set(matrix.columns) != set(meta.index):
        raise ValueError('Validation sample mismatch')
    matrix = matrix[meta.index]
    if not np.isfinite(matrix.to_numpy()).all():
        raise ValueError('Invalid validation expression')
    save_table(matrix, 'validation_gene_log2_expression.tsv.gz')
    case = matrix.loc[:, meta.days_after_MI == 7].to_numpy()
    control = matrix.loc[:, meta.days_after_MI == 0].to_numpy()
    test = ttest_ind(case, control, axis=1, equal_var=False)
    pvals = test.pvalue
    finite = np.isfinite(pvals)
    padj = np.full(len(matrix), np.nan)
    padj[finite] = false_discovery_control(pvals[finite], method='bh')
    result = pd.DataFrame({'validation_log2FC': case.mean(axis=1)-control.mean(axis=1),
                            'validation_pvalue': pvals, 'validation_padj': padj}, index=matrix.index)
    save_table(result, 'validation_D7_vs_D0_Welch.tsv')
    audit.update(matrix_probe_rows=len(probe_matrix), mapped_gene_rows=len(matrix),
                 finite_gene_tests=int(finite.sum()), samples=len(meta),
                 expression_min=float(matrix.min().min()), expression_max=float(matrix.max().max()))
    return result, matrix, meta, audit


def plots(meta, tpm, merged, validation_matrix, validation_meta):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    logtpm = np.log2(tpm + 1)
    variable = logtpm.var(axis=1).nlargest(1000).index
    pca = PCA(n_components=2)
    coordinates = pca.fit_transform(logtpm.loc[variable].T)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout='constrained')
    for day in sorted(meta.day.unique()):
        mask = meta.day.to_numpy() == day
        axes[0].scatter(coordinates[mask, 0], coordinates[mask, 1], label=f'D{day}', s=45)
    for i, (_, row) in enumerate(meta.iterrows()):
        axes[0].annotate(f'D{row.day}r{row.replicate}', coordinates[i], fontsize=7, xytext=(3, 3), textcoords='offset points')
    axes[0].set(xlabel=f'PC1 ({pca.explained_variance_ratio_[0]:.1%})',
                ylabel=f'PC2 ({pca.explained_variance_ratio_[1]:.1%})',
                title='All 24 submitted FB libraries: descriptive PCA')
    axes[0].legend(ncol=2, fontsize=8)
    data = merged.dropna(subset=['pvalue'])
    y = -np.log10(np.maximum(data.pvalue, 1e-300))
    axes[1].scatter(data.log2FoldChange, y, s=3, color='#b6bec6', alpha=0.4)
    selected = data.cross_cohort_supported
    axes[1].scatter(data.loc[selected, 'log2FoldChange'], y[selected], s=8, color='#d05535', label='Cross-cohort supported')
    axes[1].axvline(-1, color='gray', linestyle='--', linewidth=.7)
    axes[1].axvline(1, color='gray', linestyle='--', linewidth=.7)
    axes[1].set(xlabel='Discovery log2 fold change: D7 / D0', ylabel='-log10(raw p-value)',
                title='Discovery effect sizes (selection uses BH-adjusted p)')
    axes[1].legend(fontsize=8)
    fig.savefig(OUT / '01_qc_and_screen.png', dpi=170)
    plt.close(fig)
    chosen = merged.loc[merged.cross_cohort_supported].head(6)
    if chosen.empty:
        return
    fig, axes = plt.subplots(len(chosen), 2, figsize=(12, 2.5 * len(chosen)), layout='constrained', squeeze=False)
    for row, (gid, record) in enumerate(chosen.iterrows()):
        for col, (matrix, samples, days, ylabel) in enumerate([
            (logtpm, meta.index, meta.day, 'log2(TPM + 1)'),
            (validation_matrix, validation_meta.index, validation_meta.days_after_MI, 'Processed log2 expression'),
        ]):
            ax = axes[row, col]
            xs = sorted(days.unique())
            means = []
            for k, day in enumerate(xs):
                vals = matrix.loc[gid, samples[days == day]].to_numpy()
                ax.scatter(k + np.linspace(-.09, .09, len(vals)), vals, s=18, alpha=.8)
                means.append(vals.mean())
            ax.plot(range(len(xs)), means, color='black', linewidth=1)
            ax.set(xticks=range(len(xs)), xticklabels=[f'D{d}' for d in xs], ylabel=ylabel,
                   title=f"{record.symbol} | {'Discovery RNA-seq' if col == 0 else 'Independent microarray'}")
    fig.suptitle('Display examples ranked by max adjusted p; lines connect group means, not animals', fontsize=11)
    fig.savefig(OUT / '02_candidate_time_profiles.png', dpi=160)
    plt.close(fig)


def main():
    print('Reading annotation and FB quantifications...', flush=True)
    genes, tx2gene = annotation()
    meta, matrices = load_fb()
    print('Fitting discovery D7-vs-D0 model...', flush=True)
    result, tested = discovery(meta, matrices, genes)
    print('Mapping and testing independent microarray...', flush=True)
    external, matrix, external_meta, mapping = validation(tx2gene)
    merged = result.join(external)
    merged['same_direction'] = merged.log2FoldChange * merged.validation_log2FC > 0
    merged['cross_cohort_supported'] = (merged.discovery_pass & merged.same_direction &
                                        (merged.validation_padj < .05) & (merged.validation_log2FC.abs() >= .5))
    merged['display_rank_max_padj'] = merged[['padj', 'validation_padj']].max(axis=1, skipna=False)
    merged = merged.sort_values(['cross_cohort_supported', 'display_rank_max_padj'], ascending=[False, True])
    save_table(merged, 'all_evidence.tsv')
    save_table(merged.loc[merged.cross_cohort_supported], 'cross_cohort_supported_candidates.tsv')
    summary = dict(discovery_submitted_FB_libraries=len(meta), gene_rows=len(matrices['TPM']),
                   discarded_non_gene_rows=369, discovery_D0_n=4, discovery_D7_n=4,
                   discovery_genes_after_count_filter=tested,
                   discovery_finite_padj=int(result.padj.notna().sum()),
                   discovery_pass=int(result.discovery_pass.sum()),
                   discovery_candidates_mappable_to_validation=int((merged.discovery_pass & merged.validation_padj.notna()).sum()),
                   candidates_same_direction=int((merged.discovery_pass & merged.same_direction).sum()),
                   cross_cohort_supported=int(merged.cross_cohort_supported.sum()),
                   supported_up=int((merged.cross_cohort_supported & (merged.log2FoldChange > 0)).sum()),
                   supported_down=int((merged.cross_cohort_supported & (merged.log2FoldChange < 0)).sum()),
                   annotation='Ensembl GRCm38 release 87', validation_mapping=mapping,
                   interpretation='Exploratory expression association; not causal, drug-target, or clinical validation')
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    plots(meta, matrices['TPM'], merged, matrix, external_meta)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
