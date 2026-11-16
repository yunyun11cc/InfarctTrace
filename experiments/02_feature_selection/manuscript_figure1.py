"""Draw the GSE252027 library trajectory, candidate screen and expression patterns."""

import argparse
from collections import Counter
from itertools import product
from pathlib import Path

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import StrMethodFormatter
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cut_tree, leaves_list, linkage
from sklearn.decomposition import PCA

from manuscript_style import BLUE, EXPORT_DIR, WARM, WIDTH_MM, configure_style, panel_heading, save_figure
from select_features import MAIN_LIBRARY_COUNTS, STAGE_ONE


DAYS = list(MAIN_LIBRARY_COUNTS)
DAY_LABELS = [f'D{day}' for day in DAYS]
DAY_COLORS = ['#475B6D', '#6D97AD', '#31668C', '#B85B3F', '#99724D', '#656888']
DAY_MARKERS = ['o', 's', '^', 'D', 'P', 'v']
CLUSTER_COLORS = ['#317593', '#6397AB', '#B2543B', '#956639', '#635D84']
GRAY, LIGHT = '#536471', '#E4EBEE'
SETS = ['D3 vs D0', 'D7 vs D0', 'Time LRT']


def read_analysis():
    metadata = pd.read_csv(STAGE_ONE / 'tables/sample_inclusion.csv', index_col='GSM编号')
    vst = pd.read_csv(STAGE_ONE / 'matrices/timecourse_vst.tsv.gz', sep='\t', index_col='基因ID')
    libraries = metadata.loc[metadata['主分析纳入'].eq('是')].rename(columns={
        '时间点（天）': 'Day', '组内重复编号': 'Replicate',
    })
    if libraries['Day'].value_counts().to_dict() != MAIN_LIBRARY_COUNTS:
        raise ValueError('Stage 01 library counts differ from those of the primary analysis')
    if set(libraries.index) != set(vst.columns):
        raise ValueError('Stage 01 VST libraries and sample inclusion table disagree')
    return vst.loc[:, libraries.index], libraries[['Day', 'Replicate']], len(metadata)


def library_pca(vst, libraries):
    model = PCA(n_components=2, svd_solver='full')
    pcs = model.fit_transform(vst.to_numpy().T)
    points = libraries.assign(PC1=pcs[:, 0], PC2=pcs[:, 1]).rename_axis('Library')
    centroids = points.groupby('Day')[['PC1', 'PC2']].mean().assign(
        Libraries=points.groupby('Day').size(),
    )
    return points, centroids, model.explained_variance_ratio_


def screen_membership(gene_ids):
    wald = pd.read_csv(STAGE_ONE / 'tables/timecourse_wald.tsv.gz', sep='\t', dtype={'基因ID': str})
    lrt = pd.read_csv(STAGE_ONE / 'tables/timecourse_lrt.tsv.gz', sep='\t', index_col='基因ID')
    comparisons = [wald.loc[wald['比较'].eq(f'D{day}_vs_D0')].set_index('基因ID').reindex(gene_ids)
                   for day in (3, 7)]
    flags = [frame['校正后P值'].lt(0.05) & frame['log2表达倍数变化'].abs().ge(1)
             for frame in comparisons]
    membership = pd.DataFrame(dict(zip(SETS, [*flags, lrt['校正后P值'].reindex(gene_ids).lt(0.05)])))
    candidates = membership[SETS[0]] & membership[SETS[1]]
    membership['D3 and D7 candidate'] = candidates
    membership.insert(0, 'Gene_symbol', lrt['基因符号'].reindex(gene_ids))
    counts = Counter(membership[SETS].itertuples(index=False, name=None))
    combinations = list(product([True, False], repeat=3))
    overlap = pd.DataFrame(combinations, columns=SETS).assign(
        Genes=[counts[combination] for combination in combinations],
    )
    return membership.rename_axis('Gene_ID'), overlap


def cluster_candidates(vst, days, membership):
    keep = membership.index[membership['D3 and D7 candidate']].sort_values()
    means = vst.loc[keep].T.groupby(days).mean().T.reindex(columns=DAYS)
    centered = means.sub(means.mean(axis=1), axis=0)
    # Ward operates on equally weighted, gene-wise standardized time-point means.
    scaled = centered.div(means.std(axis=1, ddof=0), axis=0)
    if not np.isfinite(scaled.to_numpy()).all():
        raise ValueError('A candidate has a non-finite or constant time course')
    tree = linkage(scaled.to_numpy(), method='ward', optimal_ordering=True)
    groups = pd.Series(cut_tree(tree, n_clusters=5).ravel(), index=keep)
    profiles = scaled.groupby(groups).mean()
    ordered = sorted(profiles.index, key=lambda group: (
        profiles.loc[group].idxmax() == DAYS[0], profiles.loc[group].idxmax(), -profiles.loc[group].max(),
    ))
    names = {group: f'C{i}' for i, group in enumerate(ordered, start=1)}
    groups = groups.map(names)
    leaves = keep[leaves_list(tree)]
    ordered_genes = [gene for cluster in names.values() for gene in leaves if groups[gene] == cluster]
    metadata = membership.loc[ordered_genes, ['Gene_symbol', 'Time LRT']].assign(Cluster=groups)
    means = means.loc[ordered_genes].set_axis(DAY_LABELS, axis=1)
    scaled = scaled.loc[ordered_genes].set_axis(DAY_LABELS, axis=1)
    profiles = scaled.groupby(groups).mean().assign(Genes=groups.value_counts())
    return metadata, means, scaled, profiles


def panel_title(ax, letter, title):
    position = ax.get_position()
    panel_heading(ax.figure, position.x0 - 0.022, position.y1 + 0.027, letter, title)


def draw_pca(ax, points, centroids, variance, submitted, genes):
    panel_title(ax, 'a', 'Library trajectories')
    ax.text(0.0, 1.005, f'{len(points)} of {submitted} libraries · {genes:,} VST genes',
            transform=ax.transAxes, fontsize=6.5, color=GRAY)
    handles = []
    for day, color, marker in zip(DAYS, DAY_COLORS, DAY_MARKERS):
        group = points.loc[points['Day'].eq(day)]
        ax.scatter(group['PC1'], group['PC2'], s=19, marker=marker, facecolors='white',
                   edgecolors=color, linewidths=1, zorder=3)
        centroid = centroids.loc[day]
        ax.scatter(centroid['PC1'], centroid['PC2'], s=26, marker=marker, facecolors=color,
                   edgecolors='white', linewidths=0.55, zorder=5)
        handles.append(Line2D([], [], marker=marker, markersize=4, linestyle='none', markerfacecolor=color,
                              markeredgecolor=color, label=f'D{day} ({len(group)})'))
    for position in range(len(centroids) - 1):
        start, end = centroids[['PC1', 'PC2']].iloc[[position, position + 1]].to_numpy()
        ax.annotate('', xy=end, xytext=start, arrowprops={'arrowstyle': '->', 'color': '#8797A1',
                    'lw': 0.9, 'shrinkA': 2.5, 'shrinkB': 2.5, 'mutation_scale': 6}, zorder=2)
    ax.legend(handles=handles, ncols=3, loc='upper right', fontsize=5.8, labelspacing=0.38,
              handlelength=0.9, handletextpad=0.35, columnspacing=0.65, borderpad=0.1, frameon=True,
              facecolor='white', edgecolor='none', framealpha=0.88)
    ax.set(xlabel=f'PC1 ({variance[0]:.1%})', ylabel=f'PC2 ({variance[1]:.1%})')
    ax.margins(x=0.10, y=0.23)
    ax.axhline(0, color=LIGHT, lw=0.5, zorder=0)
    ax.axvline(0, color=LIGHT, lw=0.5, zorder=0)


def draw_overlap(parent, fig, membership, overlap):
    spec = parent.subgridspec(2, 1, height_ratios=[3.6, 1], hspace=0.06)
    bars = fig.add_subplot(spec[0])
    dots = fig.add_subplot(spec[1], sharex=bars)
    panel_title(bars, 'b', 'Wald and LRT screening overlap')
    pool_n = membership['D3 and D7 candidate'].sum()
    bars.text(0, 1.008, f'D3 ∩ D7  ·  {pool_n:,} candidate genes', transform=bars.transAxes,
              color=WARM, fontsize=6.5, fontweight='bold')
    shown = overlap.loc[overlap[SETS].any(axis=1)].sort_values('Genes', ascending=False, kind='stable')
    x = np.arange(len(shown))
    highlighted = shown[SETS[:2]].all(axis=1)
    colors = [WARM if selected else BLUE for selected in highlighted]
    bars.bar(x, shown['Genes'], color=colors, width=0.54, zorder=3)
    for position, count, color in zip(x, shown['Genes'], colors):
        bars.annotate(f'{count:,}', (position, count), xytext=(0, 2.8), textcoords='offset points',
                      ha='center', color=color, fontsize=6, fontweight='bold' if color == WARM else 'normal')
    bars.set(ylabel='Genes', ylim=(0, 9300), xlim=(-0.52, len(x) - 0.48))
    bars.set_yticks([0, 4000, 8000])
    bars.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    bars.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    bars.grid(axis='y', color=LIGHT, linewidth=0.5)
    for position, (_, row), color in zip(x, shown.iterrows(), colors):
        members = np.flatnonzero(row[SETS].to_numpy(dtype=bool))
        dots.scatter([position] * len(SETS), range(len(SETS)), s=11, color=LIGHT, edgecolors='none', zorder=1)
        dots.plot([position, position], [members.min(), members.max()], color=color, lw=1.1, zorder=2)
        dots.scatter([position] * len(members), members, s=12, color=color, edgecolors='none', zorder=3)
    dots.set(ylim=(2.42, -0.42), xticks=[])
    dots.set_yticks(range(len(SETS)), [f'{label} ({membership[label].sum():,})' for label in SETS], fontsize=6)
    dots.tick_params(axis='y', length=0, pad=4)
    dots.spines[:].set_visible(False)


def draw_heatmap(parent, fig, scaled, profiles):
    spec = parent.subgridspec(1, 3, width_ratios=[5.5, 0.11, 2.45], wspace=0.17)
    heat = fig.add_subplot(spec[0])
    panel_title(heat, 'c', 'Candidate expression patterns')
    heat.text(0.0, 1.008, f'{len(scaled):,} genes grouped by time-course similarity', transform=heat.transAxes,
              color=GRAY, fontsize=6.5)
    image = heat.imshow(scaled, cmap='RdBu_r', vmin=-2.25, vmax=2.25, aspect='auto', interpolation='nearest')
    counts = profiles['Genes'].to_numpy()
    boundaries = np.r_[0, counts.cumsum()]
    for i, (first, last) in enumerate(zip(boundaries[:-1], boundaries[1:])):
        heat.axhline(first - 0.5, lw=0.85, color='white')
        heat.plot([-0.61, -0.61], [first - 0.25, last - 0.74], color=CLUSTER_COLORS[i], lw=2.0, clip_on=False)
    heat.set_yticks([(first + last - 1) / 2 for first, last in zip(boundaries[:-1], boundaries[1:])], profiles.index)
    heat.set_xticks(range(len(DAYS)), DAY_LABELS)
    heat.tick_params(axis='both', length=0, pad=4)
    heat.spines[:].set_visible(False)
    colorbar_ax = fig.add_subplot(spec[1].subgridspec(3, 1, height_ratios=[1, 1.3, 1])[1])
    bar = fig.colorbar(image, cax=colorbar_ax, ticks=[-2, 0, 2])
    bar.ax.tick_params(labelsize=6, length=1.4, pad=2)
    bar.outline.set_visible(False)
    bar.ax.set_title('z', fontsize=6.5, pad=3)
    inset = spec[2].subgridspec(len(profiles), 1, hspace=0.11)
    for i, (cluster, row) in enumerate(profiles.iterrows()):
        ax = fig.add_subplot(inset[i])
        ax.axhline(0, lw=0.6, color=LIGHT)
        ax.plot(range(len(DAYS)), row[DAY_LABELS], color=CLUSTER_COLORS[i], lw=1.05, marker='o', ms=2.2)
        ax.text(0.02, 0.91, f'{cluster}  ·  {int(row.Genes)} genes', color=CLUSTER_COLORS[i],
                transform=ax.transAxes, va='top', fontsize=6, fontweight='bold')
        ax.set(xlim=(-0.15, 5.15), ylim=(-2, 3.1), yticks=[])
        ax.set_xticks(range(len(DAYS)), DAY_LABELS if i == len(profiles) - 1 else [])
        ax.tick_params(axis='x', labelsize=5.9, length=0, pad=3)
        ax.spines[:].set_visible(False)
        if i == 0:
            ax.set_title('Cluster mean z', fontsize=7, loc='left', fontweight='bold', pad=5)


def render(output):
    vst, libraries, submitted = read_analysis()
    points, centroids, variance = library_pca(vst, libraries)
    membership, overlap = screen_membership(vst.index)
    metadata, means, scaled, profiles = cluster_candidates(vst, libraries['Day'], membership)
    sources = output / 'source_data'
    sources.mkdir(parents=True, exist_ok=True)
    points.to_csv(sources / 'figure_1_library_pca.tsv', sep='\t')
    centroids.to_csv(sources / 'figure_1_pca_centroids.tsv', sep='\t')
    pd.Series(variance, index=['PC1', 'PC2'], name='Explained_variance_fraction').to_csv(
        sources / 'figure_1_pca_variance.tsv', sep='\t', index_label='Component',
    )
    membership.to_csv(sources / 'figure_1_screen_membership.tsv', sep='\t')
    overlap.to_csv(sources / 'figure_1_screen_overlap.tsv', sep='\t', index=False)
    metadata.join(means.add_prefix('Mean_VST_')).join(scaled.add_prefix('Z_')).to_csv(
        sources / 'figure_1_candidate_profiles.tsv', sep='\t',
    )
    profiles.to_csv(sources / 'figure_1_cluster_centroids.tsv', sep='\t', index_label='Cluster')

    configure_style()
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, 151 / 25.4))
    grid = fig.add_gridspec(2, 2, left=0.103, right=0.974, top=0.94, bottom=0.064,
                           width_ratios=[1, 1.09], height_ratios=[1, 1.32], wspace=0.37, hspace=0.41)
    draw_pca(fig.add_subplot(grid[0, 0]), points, centroids, variance, submitted, len(vst))
    draw_overlap(grid[0, 1], fig, membership, overlap)
    draw_heatmap(grid[1, :], fig, scaled, profiles)
    save_figure(fig, output, 'figure_1')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=EXPORT_DIR)
    render(parser.parse_args().output_dir)
