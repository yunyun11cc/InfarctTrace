"""Plot Stage 02 selection stability, within-cohort CV, and candidate trajectories."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

from manuscript_style import (BLUE, CORE_GENES, EXPORT_DIR, INK, WARM, WIDTH_MM,
                              configure_style, panel_heading, save_figure)
from select_features import MAIN_LIBRARY_COUNTS, STAGE_ONE, UNION_OUTPUT


SELECTORS = {
    'VST mean difference': ('均值差', BLUE, 'o'),
    'LASSO': ('LASSO', '#6E5D85', 's'),
    'Elastic Net': ('ElasticNet', '#39776E', '^'),
    'Random forest': ('随机森林', '#687582', 'D'),
}
ENSEMBLES = {'EFS-RA': ('#8F3D2A', 'o'), 'EFS-SA': ('#C67D39', 's'), 'EFS-IA': ('#432B60', '^')}
DAYS = tuple(MAIN_LIBRARY_COUNTS)


def candidate_data():
    ranking = pd.read_csv(UNION_OUTPUT / 'tables/feature_ranking.tsv', sep='\t')
    core = ranking.loc[ranking['基因符号'].isin(CORE_GENES)].set_index('基因符号')
    if not core.index.is_unique:
        raise ValueError('Candidate symbols must map uniquely to Ensembl IDs')
    core = core.loc[list(CORE_GENES)].rename(columns={'基因ID': 'gene_id'}).rename_axis('gene')
    bootstrap = pd.read_csv(UNION_OUTPUT / 'resampling/tables/bootstrap_selection_probability.tsv',
                            sep='\t', index_col='基因ID')
    columns = [method[0] for method in SELECTORS.values()] + ['Mean']
    core = core.join(bootstrap[columns], on='gene_id', validate='one_to_one')
    if not core[columns].apply(lambda series: series.between(0, 1).all()).all():
        raise ValueError('Missing or out-of-range bootstrap frequencies for the nine candidates')
    metadata = pd.read_csv(STAGE_ONE / 'tables/sample_inclusion.csv', index_col='GSM编号')
    days = metadata.loc[metadata['主分析纳入'].eq('是'), '时间点（天）'].sort_index()
    if days.value_counts().to_dict() != MAIN_LIBRARY_COUNTS:
        raise ValueError('Library counts do not match the Stage 01 main analysis')
    vst = pd.read_csv(STAGE_ONE / 'matrices/timecourse_vst.tsv.gz', sep='\t', index_col='基因ID')
    values = vst.loc[core['gene_id'], days.index].copy()
    values.index = core.index
    return core, values, days


def plot_bootstrap(ax, source):
    details = pd.read_csv(UNION_OUTPUT / 'resampling/tables/bootstrap_jaccard.tsv', sep='\t')
    groups = [details.loc[details['Method'].eq(method), 'Jaccard'].dropna() for method, _, _ in SELECTORS.values()]
    if any(len(group) != 1000 for group in groups):
        raise ValueError('Expected 1,000 defined Jaccard values per selector')
    details.assign(Method=details['Method'].replace({value[0]: key for key, value in SELECTORS.items()})).to_csv(
        source / 'figure_2a_bootstrap_jaccard.tsv', sep='\t', index=False,
    )
    boxes = ax.boxplot(groups, widths=0.53, patch_artist=True, showfliers=True, zorder=2,
                       medianprops={'color': INK, 'linewidth': 1.05},
                       boxprops={'linewidth': 0.65}, whiskerprops={'linewidth': 0.65},
                       capprops={'linewidth': 0.65}, flierprops={'markersize': 1.6, 'alpha': 0.48})
    for box, (_, color, _) in zip(boxes['boxes'], SELECTORS.values()):
        box.set_facecolor(color)
        box.set_alpha(0.33)
    ax.set_xticks(range(1, 5), ['VST mean\ndifference', 'LASSO', 'Elastic\nNet', 'Random\nforest'])
    ax.set(ylim=(-0.025, 1.025), yticks=np.linspace(0, 1, 5), ylabel='Jaccard overlap')


def plot_cv(ax, source):
    summary = pd.read_csv(source / 'figure_2b_cv_summary.tsv', sep='\t')
    methods = {name: (color, marker) for name, (_, color, marker) in SELECTORS.items()} | ENSEMBLES
    measured = summary['auroc_mean'].to_numpy()
    if set(summary['method']) != set(methods) or not np.isfinite(measured).all():
        raise ValueError('Incomplete selector comparison')
    for method, (color, marker) in methods.items():
        values = summary.loc[summary['method'].eq(method)].sort_values('k')
        if values['k'].tolist() != list(range(5, 31)) or not values['n_folds'].eq(200).all():
            raise ValueError(f'Expected 200 folds per feature count from 5 to 30 for {method}')
        ensemble = method in ENSEMBLES
        ax.plot(values['k'], values['auroc_mean'], color=color, label=method,
                linestyle='-' if ensemble else '--', linewidth=1.20 if ensemble else 0.82,
                marker=marker, markevery=5, markersize=2.5, markeredgewidth=0.30)
    legacy = pd.read_csv(UNION_OUTPUT / 'resampling/tables/cv_summary.tsv', sep='\t')
    legacy_en = legacy.set_index('Method').at['ElasticNet', 'AUROC_mean']
    pd.DataFrame([{'method': 'Legacy Elastic Net selection and refit', 'k_target': 20,
                   'auroc_mean': legacy_en}]).to_csv(source / 'figure_2b_legacy_reference.tsv', sep='\t', index=False)
    ax.scatter([20], [legacy_en], marker='*', s=34, facecolor='white', edgecolor=INK,
               linewidth=0.8, zorder=6, label=f'Legacy EN ({legacy_en:.3f})')
    low = max(0, min(np.floor(measured.min() * 20) / 20 - 0.015, legacy_en - 0.025))
    ax.set(xlim=(4.4, 30.6), ylim=(low, 1.005), xticks=range(5, 31, 5),
           xlabel='Features selected (k)', ylabel='Mean validation AUROC')
    ax.legend(loc='lower left', bbox_to_anchor=(-0.025, 1.014), ncol=4, handlelength=1.55,
              columnspacing=0.80, handletextpad=0.46, borderaxespad=0, fontsize=5.7)


def plot_candidate_selection(ax, core, source):
    ranks = pd.read_csv(source / 'figure_2c_efs_rank.tsv', sep='\t', index_col='gene_id')
    plotted = core.join(ranks[['efs_ra_rank']], on='gene_id', validate='one_to_one')
    columns = [item[0] for item in SELECTORS.values()] + ['Mean', 'efs_ra_rank']
    plotted[['gene_id', *columns]].to_csv(source / 'figure_2c_candidate_frequencies.tsv', sep='\t')
    plotted = plotted.sort_values('efs_ra_rank')
    rows = np.arange(len(plotted))
    mrc1_row = plotted.index.get_loc('Mrc1')
    ax.axhspan(mrc1_row - 0.44, mrc1_row + 0.44, color='#F5E6DE', zorder=0, linewidth=0)
    ax.hlines(rows, 0, plotted['Mean'], color='#A6B0B6', linewidth=0.8, zorder=1)
    ax.scatter(plotted['Mean'], rows, color=INK, s=12, linewidths=0, zorder=2)
    offsets = (-0.24, -0.08, 0.08, 0.24)
    for offset, (column, color, marker) in zip(offsets, SELECTORS.values()):
        ax.scatter(plotted[column], rows + offset, color=color, marker=marker, s=10,
                   linewidths=0.28, edgecolors='white', zorder=3, clip_on=False)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticks(rows, plotted.index, fontstyle='italic')
    ax.yaxis.set_tick_params(length=0, pad=4)
    ax.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    ax.set(xlim=(0, 1), ylim=(len(rows) - 0.48, -0.5), xlabel='Bootstrap selection frequency')
    ax.spines['left'].set_visible(False)
    for row, rank in zip(rows, plotted['efs_ra_rank']):
        ax.text(1.14, row, f'{rank:,}', transform=ax.get_yaxis_transform(),
                ha='right', va='center', fontsize=6.0, fontweight='bold' if row == mrc1_row else 'normal')
    ax.text(1.14, 1.015, 'RA rank', transform=ax.transAxes, ha='right', va='bottom', fontsize=6.0)
    handle = Line2D([], [], color='#A6B0B6', marker='o', markerfacecolor=INK, markeredgewidth=0,
                    markersize=3.2, linewidth=0.8, label='Mean of four selectors')
    ax.legend(handles=[handle], loc='lower left', bbox_to_anchor=(-0.04, 1.011), borderaxespad=0,
              handlelength=1.6, fontsize=6)


def plot_trajectories(fig, grid, core, values, days, source):
    table = values.T.assign(day=days).rename_axis('library').reset_index().melt(
        id_vars=['library', 'day'], var_name='gene', value_name='vst',
    )
    table['gene_id'] = table['gene'].map(core['gene_id'])
    table.to_csv(source / 'figure_2d_candidate_vst.tsv', sep='\t', index=False)
    positions = np.arange(len(DAYS))
    limits = (np.floor(values.min().min()), np.ceil(values.max().max()))
    ticks = range(int(limits[0]) + 1, int(limits[1]), 2)
    for index, gene in enumerate(CORE_GENES):
        row, column = divmod(index, 3)
        ax = fig.add_subplot(grid[row, column])
        ax.axvspan(2.73, 3.27, color='#F2EBE4', linewidth=0, zorder=0)
        means = [values.loc[gene, days.index[days.eq(day)]].mean() for day in DAYS]
        color = WARM if gene == 'Mrc1' else BLUE
        ax.plot(positions, means, linewidth=1.05, color=color, zorder=2)
        for pos, day in zip(positions, DAYS):
            library_values = values.loc[gene, days.index[days.eq(day)]].to_numpy()
            jitter = np.linspace(-0.16, 0.16, len(library_values))
            ax.scatter(pos + jitter, library_values, s=8, color=color,
                       edgecolor='white', linewidth=0.27, zorder=3)
        ax.set_title(gene, fontstyle='italic', fontweight='bold' if gene == 'Mrc1' else 'normal', fontsize=7, pad=2)
        ax.set(xlim=(-0.35, 5.35), ylim=limits, yticks=ticks)
        ax.set_xticks(positions, [str(day) for day in DAYS])
        ax.tick_params(labelbottom=row == 2, labelleft=column == 0, labelsize=5.7, length=1.8, pad=1.5)


def render(output):
    configure_style()
    source = output / 'source_data'
    if not (source / 'figure_2b_cv_summary.tsv').is_file():
        raise FileNotFoundError('Run manuscript_efs.py to create the CV comparison before drawing Figure 2')
    core, values, days = candidate_data()
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, 161 / 25.4))
    top = fig.add_gridspec(1, 2, left=0.095, right=0.977, bottom=0.645, top=0.905,
                           width_ratios=(0.87, 1.3), wspace=0.30)
    bottom = fig.add_gridspec(1, 2, left=0.105, right=0.977, bottom=0.087, top=0.515,
                              width_ratios=(0.72, 1.56), wspace=0.41)
    a, b = [fig.add_subplot(cell) for cell in top]
    c = fig.add_subplot(bottom[0, 0])
    d = bottom[0, 1].subgridspec(3, 3, hspace=0.49, wspace=0.21)
    for x, y, label, title in [
        (0.033, 0.97, 'a', 'Stability across 1,000 bootstrap draws'),
        (0.458, 0.97, 'b', 'Selector comparison across feature counts'),
        (0.033, 0.573, 'c', 'Selection frequency and rank'),
        (0.443, 0.573, 'd', 'Candidate expression trajectories'),
    ]:
        panel_heading(fig, x, y, label, title)
    plot_bootstrap(a, source)
    plot_cv(b, source)
    plot_candidate_selection(c, core, source)
    plot_trajectories(fig, d, core, values, days, source)
    fig.text(0.712, 0.046, 'Days after injury', ha='center', va='center', fontsize=7)
    fig.text(0.438, 0.299, 'VST expression', rotation=90, ha='center', va='center', fontsize=7)
    save_figure(fig, output, 'figure_2')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=EXPORT_DIR)
    render(parser.parse_args().output_dir)
