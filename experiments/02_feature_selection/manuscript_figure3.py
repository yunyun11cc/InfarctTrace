"""Generate the nine-candidate ssGSEA and marker coexpression figure."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

from manuscript_style import BLUE, CORE_GENES, EXPORT_DIR, INK, WARM, WIDTH_MM
from manuscript_style import configure_style, panel_heading, save_figure
from pathway_activity import score_pathways, spearman_permutation
from select_features import MAIN_LIBRARY_COUNTS, STAGE_ONE


HERE = Path(__file__).resolve().parent
MARKERS = {'Nlrp3': 'Inflammation', 'Casp1': 'Inflammation', 'Il1b': 'Inflammation',
           'Tgfb1': 'Fibrosis', 'Smad3': 'Fibrosis', 'Col1a1': 'Fibrosis', 'Acta2': 'Fibrosis'}
PATHWAY_LABELS = {
    'mmu04064': 'NF-κB', 'R-MMU-844456': 'NLRP3\ninflammasome', 'mmu04668': 'TNF',
    'mmu04060': 'Cytokine–\nreceptor', 'mmu04350': 'TGF-β',
    'R-MMU-1474244': 'ECM\norganization', 'R-MMU-1474290': 'Collagen\nformation',
}
NORM = TwoSlopeNorm(0, vmin=-1, vmax=1)
CORRELATION_CMAP = LinearSegmentedColormap.from_list('blue_white_rust', [BLUE, '#FAF9F5', WARM])
CORRELATION_CMAP.set_bad('#E6E6E6')


def get_expression(vst, annotation, symbols):
    """Resolve requested symbols against uniquely quantified Ensembl IDs."""
    relevant = annotation.loc[annotation['基因ID'].isin(vst.index)].set_index('基因符号')['基因ID']
    ids = relevant.loc[relevant.index.isin(symbols)].reindex(symbols).rename('gene_id')
    values = vst.reindex(index=ids).set_axis(symbols).T
    return ids, values


def correlate(expression, features):
    """Retain missing features in the output and adjust only the estimable tests."""
    rows, cols = len(expression.columns), len(features.columns)
    rho, p, q = (np.full((rows, cols), np.nan) for _ in range(3))
    valid = features.notna().all() & features.nunique().gt(1)
    if valid.any():
        rho[:, valid], p[:, valid] = spearman_permutation(expression, features.loc[:, valid])
        q[:, valid] = false_discovery_control(p[:, valid].ravel()).reshape(rows, int(valid.sum()))
    return rho, p, q


def as_rows(genes, features, statistics):
    columns = pd.MultiIndex.from_product([genes, features], names=['gene', 'feature'])
    return pd.DataFrame({name: values.ravel() for name, values in zip(('rho', 'p', 'q'), statistics)},
                        index=columns).reset_index()


def pathway_data(vst, expression, ids, days, tables):
    membership = pd.read_csv(HERE / 'ssgsea_pathway_membership.tsv', sep='\t')
    pathways = membership[['PathwayID', 'Pathway', 'Group']].drop_duplicates().set_index('PathwayID')
    pathways = pathways.loc[list(PATHWAY_LABELS)]
    mapping = pd.read_csv(STAGE_ONE / 'tables/enrichment_gene_mapping.tsv', sep='\t',
                          usecols=['基因ID', 'EntrezID'])
    members = membership.merge(mapping, on='EntrezID', how='left', validate='many_to_one')
    matched = members.loc[members['基因ID'].isin(vst.index)]
    gene_sets = {pathway: sorted(set(matched.loc[matched['PathwayID'].eq(pathway), '基因ID']))
                 for pathway in pathways.index}
    reduced_sets = {pathway: sorted(set(genes) - set(ids)) for pathway, genes in gene_sets.items()}
    scores = score_pathways(vst, gene_sets)
    reduced_scores = score_pathways(vst.drop(index=ids), reduced_sets)
    names = list(expression.columns)
    comparisons = as_rows(names, scores.columns, correlate(expression, scores))
    sensitivity = as_rows(names, scores.columns, correlate(expression, reduced_scores))
    result = comparisons.merge(sensitivity, on=['gene', 'feature'], suffixes=('', '_without_core'))
    result['gene_id'] = result['gene'].map(ids)
    result = result.assign(pathway=result['feature'].map(pathways['Pathway']),
                           group=result['feature'].map(pathways['Group']),
                           genes_in_pathway=result['feature'].map({key: len(value) for key, value in gene_sets.items()}),
                           genes_in_pathway_without_core=result['feature'].map(
                               {key: len(value) for key, value in reduced_sets.items()}),
                           candidate_in_pathway=[gene_id in gene_sets[pathway]
                                                 for gene_id, pathway in zip(result['gene_id'], result['feature'])])
    result = result.rename(columns={'feature': 'pathway_id'})
    columns = ['gene', 'gene_id', 'pathway_id', 'pathway', 'group', 'rho', 'p', 'q',
               'rho_without_core', 'p_without_core', 'q_without_core', 'genes_in_pathway',
               'genes_in_pathway_without_core', 'candidate_in_pathway']
    result[columns].to_csv(tables / 'figure_3a_pathway_correlations.tsv', sep='\t', index=False, float_format='%.12g')
    scores.insert(0, 'day', days)
    scores.join(reduced_scores.add_suffix('_without_core')).to_csv(
        tables / 'figure_3a_pathway_scores.tsv', sep='\t', index_label='library', float_format='%.12g',
    )
    return result


def marker_data(expression, ids, markers, marker_ids, days, tables):
    results = as_rows(expression.columns, markers.columns, correlate(expression, markers))
    marker_status = pd.Series(np.select([markers.isna().all(), markers.nunique().le(1)],
                                        ['not_quantified', 'constant'], default='measured'), index=markers.columns)
    results = results.assign(gene_id=results['gene'].map(ids), marker_id=results['feature'].map(marker_ids),
                             group=results['feature'].map(MARKERS), status=results['feature'].map(marker_status))
    results = results.rename(columns={'feature': 'marker'})
    results[['gene', 'gene_id', 'marker', 'marker_id', 'group', 'rho', 'p', 'q', 'status']].to_csv(
        tables / 'figure_3b_marker_correlations.tsv', sep='\t', index=False, float_format='%.12g', na_rep='NA',
    )
    expression.join(markers).assign(day=days)[['day', *expression.columns, *markers.columns]].to_csv(
        tables / 'figure_3_expression.tsv', sep='\t', index_label='library', float_format='%.12g', na_rep='NA',
    )
    return results


def draw_heatmap(ax, results, feature, labels, inflammation_columns, italic_columns=False, show_genes=True):
    rho = results.pivot(index='gene', columns=feature, values='rho')
    q = results.pivot(index='gene', columns=feature, values='q')
    ordered_columns = results[feature].unique()
    rho = rho.loc[list(CORE_GENES), ordered_columns].to_numpy()
    q = q.loc[list(CORE_GENES), ordered_columns].to_numpy()
    image = ax.imshow(rho, cmap=CORRELATION_CMAP, norm=NORM, interpolation='nearest', aspect='auto')
    ax.set_xticks(range(len(labels)), labels, fontstyle='italic' if italic_columns else 'normal')
    ax.set_yticks(range(len(CORE_GENES)), CORE_GENES if show_genes else [''] * len(CORE_GENES), fontstyle='italic')
    ax.tick_params(axis='both', length=0, pad=3, labelsize=6)
    ax.xaxis.tick_top()
    for tick in ax.get_yticklabels():
        if tick.get_text() == 'Mrc1':
            tick.set_color(WARM)
            tick.set_fontweight('bold')
    ax.spines[:].set_visible(False)
    for x in np.arange(0.5, rho.shape[1] - 0.5):
        ax.axvline(x, color='white', lw=1.5 if x == inflammation_columns - 0.5 else 0.45)
    for y in np.arange(0.5, rho.shape[0] - 0.5):
        ax.axhline(y, color='white', lw=0.45)
    for (row, col), value in np.ndenumerate(rho):
        stars = '**' if q[row, col] < 0.01 else '*' if q[row, col] < 0.05 else ''
        color = 'white' if abs(value) >= 0.78 else INK
        label = f'{value:.2f}{stars}'.replace('-', '−') if np.isfinite(value) else 'NA'
        ax.text(col, row, label, ha='center', va='center', color=color, fontsize=5.5)
    ax.add_patch(Rectangle((-0.48, len(CORE_GENES) - 1.48), len(labels) - 0.04, 0.96,
                           fill=False, edgecolor=WARM, lw=0.85, clip_on=False))
    for left, right, label, color in [(-0.45, inflammation_columns - 0.55, 'Inflammation', BLUE),
                                      (inflammation_columns - 0.45, len(labels) - 0.55, 'Fibrosis', WARM)]:
        transform = ax.get_xaxis_transform()
        ax.plot([left, right], [1.18, 1.18], transform=transform, color=color, lw=1, clip_on=False)
        ax.text((left + right) / 2, 1.205, label, transform=transform, ha='center', va='bottom',
                color=color, fontsize=6.5, fontweight='bold')
    return image


def plot_figure(pathways, markers, output, n_libraries):
    configure_style()
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, 91 / 25.4))
    left = fig.add_axes((0.11, 0.245, 0.42, 0.50))
    right = fig.add_axes((0.55, 0.245, 0.42, 0.50))
    image = draw_heatmap(left, pathways, 'pathway_id', tuple(PATHWAY_LABELS.values()), 4)
    draw_heatmap(right, markers, 'marker', list(MARKERS), 3, italic_columns=True, show_genes=False)
    panel_heading(fig, 0.07, 0.927, 'a', 'ssGSEA pathway scores')
    panel_heading(fig, 0.54, 0.927, 'b', 'Marker coexpression')
    cax = fig.add_axes((0.35, 0.166, 0.34, 0.024))
    cb = fig.colorbar(image, cax=cax, orientation='horizontal', ticks=[-1, -0.5, 0, 0.5, 1])
    cb.ax.tick_params(labelsize=6, length=2, width=0.5, pad=2)
    cb.outline.set_linewidth(0.45)
    fig.text(0.32, 0.176, 'Spearman ρ', ha='right', va='center', fontsize=6.5)
    fig.text(0.54, 0.063,
             f'{n_libraries} libraries pooled across six time points  ·  BH q < 0.05 (*) or < 0.01 (**)',
             ha='center', va='center', fontsize=6.5, color='#465865')
    save_figure(fig, output, 'figure_3')
    plt.close(fig)


def run(output):
    metadata = pd.read_csv(STAGE_ONE / 'tables/sample_inclusion.csv', index_col='GSM编号')
    days = metadata.loc[metadata['主分析纳入'].eq('是'), '时间点（天）']
    vst = pd.read_csv(STAGE_ONE / 'matrices/timecourse_vst.tsv.gz', sep='\t', index_col='基因ID')
    if days.value_counts().to_dict() != MAIN_LIBRARY_COUNTS or set(days.index) != set(vst.columns):
        raise ValueError('Retained library metadata do not match the Stage 01 VST matrix')
    vst = vst.loc[:, days.index]
    annotation = pd.read_csv(STAGE_ONE / 'matrices/gene_annotation.tsv.gz', sep='\t')
    ids, expression = get_expression(vst, annotation, list(CORE_GENES))
    if ids.isna().any() or expression.nunique().le(1).any():
        raise ValueError('A shortlisted gene is missing or invariant across the retained libraries')
    marker_ids, markers = get_expression(vst, annotation, list(MARKERS))
    tables = output / 'source_data'
    tables.mkdir(parents=True, exist_ok=True)
    pathways = pathway_data(vst, expression, ids, days, tables)
    marker_results = marker_data(expression, ids, markers, marker_ids, days, tables)
    plot_figure(pathways, marker_results, output, len(days))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=EXPORT_DIR)
    run(parser.parse_args().output_dir)
