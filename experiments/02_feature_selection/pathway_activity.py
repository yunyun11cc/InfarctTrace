"""Score the 22 libraries with ssGSEA and correlate the scores with the top eight genes."""

import argparse
import json
from pathlib import Path
from textwrap import wrap

import gseapy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, rankdata

from select_features import MAIN_LIBRARY_COUNTS, SEED, STAGE_ONE, UNION_OUTPUT


HERE = Path(__file__).resolve().parent
N_PERMUTATIONS = 19_999
MIN_SIZE = 5  # The mouse Reactome NLRP3 set contains only ten source genes.
MAX_SIZE = 500


def score_pathways(vst, gene_sets):
    result = gseapy.ssgsea(data=vst, gene_sets=gene_sets, organism='mouse', outdir=None,
                          sample_norm_method='rank', correl_norm_type='rank',
                          min_size=MIN_SIZE, max_size=MAX_SIZE, permutation_num=0,
                          weight=0.25, threads=1, seed=SEED, no_plot=True)
    scores = result.res2d.pivot(index='Name', columns='Term', values='NES')
    return scores.loc[vst.columns, list(gene_sets)].astype(float)


def standardized_ranks(values):
    ranks = rankdata(values, axis=0)
    ranks = ranks - ranks.mean(axis=0)
    norms = np.linalg.norm(ranks, axis=0)
    if (norms == 0).any():
        raise ValueError('核心基因表达或通路得分恒定，无法计算 Spearman 相关')
    return ranks / norms


def spearman_permutation(expression, scores, permutations=N_PERMUTATIONS):
    """Two-sided marginal tests with shared random permutations of the library pairing."""
    if not expression.index.equals(scores.index):
        raise ValueError('核心基因表达与 ssGSEA 得分的文库顺序不一致')
    gene = standardized_ranks(expression.to_numpy())
    pathway = standardized_ranks(scores.to_numpy())
    observed = gene.T @ pathway
    rng = np.random.default_rng(SEED)
    order = rng.permuted(np.tile(np.arange(len(expression)), (permutations, 1)), axis=1)
    null = np.einsum('pna,nb->pab', gene[order], pathway)
    p = (1 + (np.abs(null) >= np.abs(observed) - 1e-12).sum(axis=0)) / (permutations + 1)
    return observed, p


def plot_correlations(genes, pathways, rho, p, figure_path, library_count):
    fig, ax = plt.subplots(figsize=(10.6, 6.0))
    fig.subplots_adjust(left=0.135, right=0.90, bottom=0.13, top=0.69)
    fig.text(0.135, 0.95, 'Core gene expression and pathway enrichment', fontsize=17, fontweight='bold')
    fig.text(0.135, 0.905, f'{library_count} libraries  |  Spearman rank correlation with ssGSEA scores',
             fontsize=10, color='#4A565F')
    image = ax.imshow(rho, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    labels = ['\n'.join(wrap(name, width=15, break_long_words=False)) for name in pathways['Pathway']]
    labels = [label.replace('Cytokine–receptor\ninteraction', 'Cytokine–\nreceptor interaction') for label in labels]
    ax.set_xticks(range(len(pathways)), labels)
    ax.set_yticks(range(len(genes)), genes)
    ax.xaxis.tick_top()
    ax.tick_params(axis='both', length=0, labelsize=9, pad=8)
    for tick in ax.get_yticklabels():
        tick.set_fontstyle('italic')
    for group, color in [('Inflammation', '#35566B'), ('Fibrosis', '#855743')]:
        columns = np.flatnonzero(pathways['Group'].eq(group))
        left, right = columns[0] - 0.44, columns[-1] + 0.44
        ax.plot([left, right], [1.18, 1.18], transform=ax.get_xaxis_transform(), color=color,
                linewidth=2, clip_on=False)
        ax.text((left + right) / 2, 1.205, group, transform=ax.get_xaxis_transform(),
                ha='center', va='bottom', color=color, fontsize=10, fontweight='bold')
    ax.axvline(pathways['Group'].eq('Inflammation').sum() - 0.5, color='white', linewidth=2)
    for row, column in np.ndindex(rho.shape):
        stars = '**' if p[row, column] < 0.01 else '*' if p[row, column] < 0.05 else ''
        value = rho[row, column]
        ax.text(column, row, f'{value:.2f}{stars}', ha='center', va='center', fontsize=9,
                color='white' if abs(value) >= 0.65 else '#172B3A')
    ax.spines[:].set_visible(False)
    cax = fig.add_axes((0.925, 0.24, 0.012, 0.32))
    fig.colorbar(image, cax=cax, ticks=[-1, -0.5, 0, 0.5, 1], label='Spearman ρ')
    fig.text(0.135, 0.065, '* P < 0.05,  ** P < 0.01  |  Two-sided permutation P, before BH correction',
             fontsize=9, color='#4A565F')
    fig.savefig(figure_path, dpi=240)
    plt.close(fig)


def run_analysis(output):
    ranking = pd.read_csv(UNION_OUTPUT / 'tables/feature_ranking.tsv', sep='\t', index_col='基因ID')
    genes = ranking.nsmallest(8, '排名')['基因符号']
    metadata = pd.read_csv(STAGE_ONE / 'tables/sample_inclusion.csv', index_col='GSM编号')
    days = metadata.loc[metadata['主分析纳入'].eq('是'), '时间点（天）']
    vst = pd.read_csv(STAGE_ONE / 'matrices/timecourse_vst.tsv.gz', sep='\t', index_col='基因ID')
    if days.value_counts().to_dict() != MAIN_LIBRARY_COUNTS or set(days.index) != set(vst.columns):
        raise ValueError('VST 矩阵与阶段 01_1 的 22 份主分析文库不一致')
    vst = vst.loc[:, days.index]

    membership = pd.read_csv(HERE / 'ssgsea_pathway_membership.tsv', sep='\t')
    pathways = membership[['PathwayID', 'Pathway', 'Group', 'Database']].drop_duplicates().set_index('PathwayID')
    mapping = pd.read_csv(STAGE_ONE / 'tables/enrichment_gene_mapping.tsv', sep='\t',
                          usecols=['基因ID', 'EntrezID']).rename(columns={'基因ID': 'EnsemblID'})
    members = membership.merge(mapping, on='EntrezID', how='left', validate='many_to_one')
    members['In_VST'] = members['EnsemblID'].isin(vst.index)
    members['Core_gene'] = members['EnsemblID'].map(genes)
    gene_sets = {pathway: sorted(members.loc[members['PathwayID'].eq(pathway) & members['In_VST'], 'EnsemblID'])
                 for pathway in pathways.index}
    coverage = members.groupby('PathwayID')['In_VST'].agg(['size', 'sum']).rename(
        columns={'size': 'Source_genes', 'sum': 'Genes_in_VST'}
    )
    scores = score_pathways(vst, gene_sets)
    expression = vst.loc[genes.index].T
    rho, p = spearman_permutation(expression, scores)
    q = false_discovery_control(p.ravel()).reshape(p.shape)

    # Re-score after removing the eight queried genes from the ranking universe.
    without_core = {pathway: sorted(set(identifiers) - set(genes.index)) for pathway, identifiers in gene_sets.items()}
    sensitivity = score_pathways(vst.drop(index=genes.index), without_core)
    sensitivity_rho = standardized_ranks(expression.to_numpy()).T @ standardized_ranks(sensitivity.to_numpy())
    rows = []
    for i, (identifier, symbol) in enumerate(genes.items()):
        for j, (pathway, annotation) in enumerate(pathways.iterrows()):
            rows.append((identifier, symbol, pathway, annotation['Pathway'], annotation['Group'],
                         rho[i, j], p[i, j], q[i, j], identifier in gene_sets[pathway], sensitivity_rho[i, j]))
    correlations = pd.DataFrame(rows, columns=['EnsemblID', 'Gene', 'PathwayID', 'Pathway', 'Group',
                                              'Spearman_rho', 'P', 'q_BH', 'Gene_in_pathway', 'Rho_without_core'])
    tables, figures, reports = (output / name for name in ('tables', 'figures', 'reports'))
    for directory in (tables, figures, reports):
        directory.mkdir(parents=True, exist_ok=True)
    matrices = {
        'ssgsea_scores': scores.assign(Day=days)[['Day', *pathways.index]].rename_axis('Library'),
        'ssgsea_spearman_rho': pd.DataFrame(rho, index=pd.Index(genes, name='Gene'), columns=pathways.index),
        'ssgsea_pathway_coverage': pathways.join(coverage),
    }
    for name, matrix in matrices.items():
        matrix.to_csv(tables / f'{name}.tsv', sep='\t', encoding='utf-8-sig', float_format='%.10g')
    correlations.to_csv(tables / 'ssgsea_correlations.tsv', sep='\t', index=False,
                        encoding='utf-8-sig', float_format='%.10g')
    summary = {'GSEApy': gseapy.__version__, 'VST_shape': list(vst.shape), 'Core_genes': genes.to_dict(),
               'library_counts': days.value_counts().sort_index().to_dict(),
               'score': 'ssGSEA NES, rank/rank, weight 0.25, permutation_num 0',
               'gene_set_sizes': [MIN_SIZE, MAX_SIZE], 'permutations': N_PERMUTATIONS, 'seed': SEED,
               'P_value': 'Two-sided Spearman pairing permutation, all 22 libraries pooled',
               'BH_tests': int(p.size), 'pathway_sources': json.loads((HERE / 'ssgsea_pathway_sources.json').read_text())}
    (reports / 'ssgsea_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    plot_correlations(genes, pathways, rho, p, figures / '04_ssgsea_correlation.png', len(days))
    print(f'{len(days)} libraries × {len(pathways)} pathways, {p.size} gene–pathway correlations: {output}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=UNION_OUTPUT)
    run_analysis(parser.parse_args().output_dir)
