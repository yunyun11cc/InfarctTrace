"""Plot the Stage 02 intersection, candidate expression, and selection stability."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
from matplotlib import font_manager
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.patheffects import withStroke
from matplotlib.ticker import PercentFormatter, StrMethodFormatter
import numpy as np
import pandas as pd

from figure_data import CORE_GENES, SETS, efs_ra_ranks, screen_membership
from select_features import OUTPUT, STAGE_ONE, UNION_OUTPUT


BLUE = '#235F84'
ORANGE = '#B4552A'
GRAY = '#62717D'
LIGHT = '#DFE5EA'
INK = '#192D3B'
CHINESE_FONTS = ('Noto Sans CJK SC', 'PingFang SC', 'Microsoft YaHei', 'Source Han Sans SC', 'SimHei')
UNION_BLUE, UNION_RED, UNION_GRAY = '#1f5c7e', '#a63a2a', '#8c9ba5'
UNION_METHODS = {
    '均值差': ('VST mean difference', UNION_BLUE, 'o'),
    'LASSO': ('LASSO', '#596f80', 's'),
    'ElasticNet': ('Elastic Net', UNION_GRAY, '^'),
    '随机森林': ('Random forest', UNION_RED, 'D'),
}


def configure_style():
    matplotlib.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7,
        'text.color': INK,
        'axes.labelcolor': INK,
        'axes.edgecolor': INK,
        'axes.linewidth': 0.65,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.color': INK,
        'ytick.color': INK,
        'xtick.labelsize': 6,
        'ytick.labelsize': 6,
        'xtick.major.size': 2.5,
        'ytick.major.size': 2.5,
        'legend.frameon': False,
        'legend.fontsize': 6,
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
    })


def union_figure(height, title, subtitle):
    fig = plt.figure(figsize=(10.6, height))
    fig.text(0.135, 0.95, title, fontsize=19, fontweight='bold', color=INK)
    fig.text(0.135, 0.901, subtitle, fontsize=11, color=GRAY)
    return fig


def union_axes(ax, grid_axis):
    ax.set_axisbelow(True)
    ax.grid(axis=grid_axis, color=LIGHT, linewidth=0.7)
    ax.spines[['left', 'bottom']].set_color(UNION_GRAY)
    ax.tick_params(labelsize=9)


def plot_union(ranking, overlap):
    d3, d7 = (ranking[f'D{day}通过差异筛选'] for day in (3, 7))
    shared = d3 & d7
    fig = union_figure(6.0, 'Differential gene overlap and effect trajectories',
                       f'22 libraries  |  {len(ranking):,} genes passed the D3 or D7 Wald screen')
    grid = fig.add_gridspec(2, 2, left=0.135, right=0.963, top=0.79, bottom=0.22,
                            height_ratios=(3.8, 1), width_ratios=(0.92, 1.17), hspace=0.05, wspace=0.23)
    bars = fig.add_subplot(grid[0, 0])
    matrix = fig.add_subplot(grid[1, 0], sharex=bars)
    scatter = fig.add_subplot(grid[:, 1])
    bars.set_title('A   Wald and time-course LRT intersections', loc='left', fontweight='bold', pad=16)
    scatter.set_title('B   D3 and D7 effects', loc='left', fontweight='bold', pad=16)

    exact = overlap.loc[overlap[SETS].any(axis=1)]
    colors = np.where(exact[SETS[:2]].all(axis=1), UNION_RED, UNION_GRAY)
    x = np.arange(len(exact))
    bars.bar(x, exact['Genes'], width=0.60, color=colors)
    for column, count, color in zip(x, exact['Genes'], colors, strict=True):
        bars.annotate(f'{count:,}', (column, count), xytext=(0, 4), textcoords='offset points',
                      ha='center', va='bottom', fontsize=9, color=UNION_RED if color == UNION_RED else INK)
    bars.text(0.04, 0.92, f'D3 ∩ D7   {shared.sum():,} genes  ·  {shared.mean():.1%} of union',
              transform=bars.transAxes, va='top', color=UNION_RED, fontsize=10, fontweight='bold')
    bars.set(xlim=(-0.55, len(x) - 0.45), ylim=(0, exact['Genes'].max() * 1.18),
             ylabel='Genes in exact intersection', yticks=[0, 2000, 4000, 6000, 8000])
    bars.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    bars.tick_params(axis='x', bottom=False, labelbottom=False)
    union_axes(bars, 'y')
    for row in range(len(SETS)):
        if row % 2 == 0:
            matrix.axhspan(row - 0.47, row + 0.47, color='#f2f5f7', zorder=0)
        matrix.scatter(x, np.full(len(x), row), s=22, color='#d5dfe4', linewidths=0, zorder=1)
    for column, membership, color in zip(x, exact[SETS].to_numpy(), colors, strict=True):
        rows = np.flatnonzero(membership)
        matrix.vlines(column, rows.min(), rows.max(), color=color, lw=1.7, zorder=2)
        matrix.scatter(np.full(len(rows), column), rows, s=34, color=color, linewidths=0, zorder=3)
    matrix.set(ylim=(2.5, -0.5), yticks=range(len(SETS)), yticklabels=SETS, xticks=[])
    matrix.tick_params(axis='both', length=0, labelsize=9)
    matrix.spines[:].set_visible(False)

    x_field, y_field = 'D3对D0_log2FC', 'D7对D0_log2FC'
    x_min, x_max = np.floor(ranking[x_field].min()) - 0.5, np.ceil(ranking[x_field].max()) + 0.5
    y_min, y_max = np.floor(ranking[y_field].min()), np.ceil(ranking[y_field].max()) + 0.6
    shade = np.array([x_min, y_min, y_max, x_max])
    scatter.fill_between(shade, np.clip(shade, y_min, y_max), y_max, color='#f0f5f8', zorder=0)
    scatter.axline((0, 0), slope=1, color=UNION_BLUE, ls='--', lw=0.9, alpha=0.7, zorder=1)
    scatter.axhline(0, color=UNION_GRAY, lw=0.7, zorder=1)
    scatter.axvline(0, color=UNION_GRAY, lw=0.7, zorder=1)
    scatter.text(0.035, 0.95, 'D7 effect > D3 effect', transform=scatter.transAxes,
                 va='top', color=UNION_BLUE, fontsize=9)
    direction = ranking['相对D0方向']
    same = shared & direction.isin(['均上调', '均下调'])
    for mask, label, color, marker, size, alpha in [
        (~shared, 'One contrast', UNION_GRAY, 'o', 7, 0.28),
        (same, 'Both, same sign', UNION_BLUE, 'o', 10, 0.48),
        (shared & ~same, 'Both, sign reversal', UNION_RED, 'D', 29, 1),
    ]:
        points = ranking.loc[mask]
        scatter.scatter(points[x_field], points[y_field], color=color, marker=marker, s=size,
                        alpha=alpha, linewidths=0, label=label, zorder=3)
    spotlight = ranking.set_index('基因符号')
    for gene, offset in [('C1qtnf3', (12, 20)), ('Postn', (17, 9)), ('Cthrc1', (9, 21))]:
        point = tuple(spotlight.loc[gene, [x_field, y_field]])
        scatter.plot(*point, marker='o', ms=5.5, mfc=UNION_BLUE, mec='white', mew=0.8, zorder=4)
        scatter.annotate(gene, point, xytext=offset, textcoords='offset points', color=UNION_BLUE,
                         fontsize=10, fontweight='bold', fontstyle='italic', va='center', zorder=5,
                         path_effects=[withStroke(linewidth=2, foreground='white')],
                         arrowprops={'arrowstyle': '-', 'color': UNION_BLUE, 'lw': 0.8})
    scatter.set(xlim=(x_min, x_max), ylim=(y_min, y_max), xlabel='D3 vs D0 (log$_2$ fold change)',
                ylabel='D7 vs D0 (log$_2$ fold change)')
    scatter.spines[['left', 'bottom']].set_color(UNION_GRAY)
    scatter.tick_params(labelsize=9)
    scatter.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), ncols=3, fontsize=8.5,
                   handletextpad=0.25, columnspacing=1.0, scatterpoints=1, markerscale=1.4)
    fig.text(0.135, 0.073, 'Wald screens: BH q < 0.05 and |log$_2$FC| ≥ 1.   LRT: BH q < 0.05.   Dashed line: y = x.',
             fontsize=9, color=GRAY)
    return fig


def plot_union_expression(core, vst, days):
    grouped = vst.loc[core.index, days.index].T.groupby(days, sort=True)
    mean, sd = grouped.mean(), grouped.std(ddof=1)
    counts = grouped.size()
    fig = union_figure(6.6, 'Expression trajectories of nine core genes',
                       f'{len(days)} libraries across {len(mean)} sampling days  |  Lines show means, ribbons show ±1 SD')
    grid = fig.add_gridspec(3, 3, left=0.135, right=0.963, top=0.80, bottom=0.16, wspace=0.23, hspace=0.66)
    positions = np.arange(len(mean))
    peak = mean.index.get_loc(7)
    low, high = np.floor((mean - sd).min().min()), np.ceil((mean + sd).max().max())
    for index, (gene_id, gene) in enumerate(core['基因符号'].items()):
        row, column = divmod(index, 3)
        ax = fig.add_subplot(grid[row, column])
        accent = gene == 'Mrc1'
        color = UNION_RED if accent else UNION_BLUE
        ax.set_facecolor('#fcf7f4' if accent else 'white')
        ax.axvspan(peak - 0.31, peak + 0.31, color='#f1e7df', lw=0, zorder=0)
        ax.fill_between(positions, mean[gene_id] - sd[gene_id], mean[gene_id] + sd[gene_id],
                        color=color, alpha=0.16, lw=0, zorder=2)
        ax.plot(positions, mean[gene_id], color=color, lw=2.05, marker='o', ms=4.1,
                markeredgecolor='white', markeredgewidth=0.55, zorder=3)
        ax.set_title(gene, loc='left', fontstyle='italic', fontweight='bold', color=color, pad=7)
        ax.set(xlim=(-0.30, len(mean) - 0.70), ylim=(low, high),
               xticks=positions, xticklabels=[f'D{day}' for day in mean.index],
               yticks=np.arange(np.ceil(low / 2) * 2, high + 0.1, 4))
        union_axes(ax, 'y')
        ax.tick_params(axis='y', labelleft=column == 0, left=column == 0)
        ax.tick_params(axis='x', length=0, pad=4)
        peak_label = ax.get_xticklabels()[peak]
        peak_label.set_fontweight('bold')
        if accent:
            ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, fill=False, ls='--',
                                   lw=0.8, edgecolor=UNION_RED, clip_on=False, zorder=4))
            ax.text(1, 1.065, 'Multivariable addition', transform=ax.transAxes, ha='right',
                    va='bottom', fontsize=9, color=UNION_RED)
    fig.text(0.074, 0.48, 'VST expression', rotation=90, ha='center', va='center', fontsize=10)
    labels = '   ·   '.join(f'D{day}: {count}' for day, count in counts.items())
    fig.text(0.135, 0.066, f'Libraries per day   {labels}', color=GRAY, fontsize=9)
    return fig


def plot_union_stability(detail, evidence, universe_size):
    groups = [detail.loc[detail['Method'].eq(method), 'Jaccard'].to_numpy() for method in UNION_METHODS]
    draws = detail['Draw'].nunique()
    if any(len(group) != draws or not np.isfinite(group).all() for group in groups):
        raise ValueError('Expected one defined Jaccard value per bootstrap draw and method')
    fig = union_figure(6.15, 'Selection stability and core gene support',
                       f'{draws:,} stratified bootstrap draws  |  D7 vs all other days')
    distribution = fig.add_subplot(fig.add_gridspec(1, 1, left=0.135, right=0.424, bottom=0.199, top=0.775)[0])
    grid = fig.add_gridspec(1, 2, left=0.588, right=0.963, bottom=0.199, top=0.775,
                            width_ratios=(5, 0.86), wspace=0.19)
    frequency = fig.add_subplot(grid[0, 0])
    ranks = fig.add_subplot(grid[0, 1], sharey=frequency)
    distribution.set_title('A   Overlap with full-cohort selection', loc='left', fontweight='bold', pad=39)
    frequency.set_title('B   Core gene selection and rank', loc='left', fontweight='bold', pad=39)

    positions = np.arange(len(groups))
    violins = distribution.violinplot(groups, positions=positions, orientation='horizontal',
                                     side='low', widths=0.68, showextrema=False)
    for density, (_, color, _) in zip(violins['bodies'], UNION_METHODS.values(), strict=True):
        density.set_facecolor(color)
        density.set_edgecolor(color)
        density.set_alpha(0.30)
        density.set_linewidth(0.7)
    boxes = distribution.boxplot(groups, positions=positions + 0.11, orientation='horizontal', widths=0.19,
                                 patch_artist=True, manage_ticks=False, zorder=3,
                                 medianprops={'color': INK, 'linewidth': 1.4},
                                 flierprops={'marker': '.', 'markersize': 2.3, 'alpha': 0.35,
                                             'markeredgecolor': UNION_GRAY})
    for position, values, box, (_, color, _) in zip(positions, groups, boxes['boxes'], UNION_METHODS.values(), strict=True):
        box.set(facecolor=color, edgecolor=INK, alpha=0.62, linewidth=0.8)
        distribution.text(1.05, position, f'Median = {np.median(values):.2f}',
                          transform=distribution.get_yaxis_transform(), color=color,
                          fontsize=9, va='center', fontweight='bold')
    distribution.set(xlim=(0, 1), ylim=(len(groups) - 0.52, -0.52), yticks=positions,
                     yticklabels=['VST mean\ndifference', 'LASSO', 'Elastic Net', 'Random forest'],
                     xticks=np.linspace(0, 1, 6), xlabel='Jaccard overlap')
    union_axes(distribution, 'x')
    distribution.tick_params(axis='y', length=0, pad=6)

    evidence = evidence.sort_values('efs_ra_rank')
    rows = np.arange(len(evidence))
    mrc1 = evidence.index.get_loc('Mrc1')
    frequency.axhspan(mrc1 - 0.45, mrc1 + 0.45, color='#f8eeea', lw=0, zorder=0)
    frequency.hlines(rows, 0, evidence['Mean'], color=UNION_GRAY, lw=0.9, zorder=2)
    frequency.scatter(evidence['Mean'], rows, s=19, c=INK, linewidths=0.4, edgecolors='white', zorder=3)
    handles = []
    for offset, (key, (label, color, marker)) in zip((-0.29, -0.10, 0.10, 0.29), UNION_METHODS.items(), strict=True):
        frequency.scatter(evidence[key], rows + offset, color=color, marker=marker, s=24,
                          linewidths=0.45, edgecolors='white', clip_on=False, zorder=4)
        handles.append(Line2D([], [], marker=marker, color=color, linestyle='none', markersize=5.1, label=label))
    handles.append(Line2D([], [], marker='o', color=UNION_GRAY, markerfacecolor=INK, markeredgecolor=INK,
                          markersize=3.6, lw=0.8, label='Four-method mean'))
    frequency.set(xlim=(0, 1), ylim=(len(rows) - 0.5, -0.5), xticks=np.linspace(0, 1, 5),
                  xlabel='Bootstrap selection frequency', yticks=rows, yticklabels=evidence.index)
    frequency.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    union_axes(frequency, 'x')
    frequency.spines['left'].set_visible(False)
    frequency.tick_params(axis='y', length=0, pad=8)
    for row, label in enumerate(frequency.get_yticklabels()):
        label.set_fontstyle('italic')
        if row == mrc1:
            label.set(color=UNION_RED, fontweight='bold')
    frequency.legend(handles=handles, loc='lower left', bbox_to_anchor=(-0.03, 1.008), ncols=3,
                     fontsize=8, handlelength=1.2, columnspacing=1.05, handletextpad=0.36, borderpad=0.15)

    ranks.set_axis_off()
    ranks.set_title('EFS-RA\nrank', fontsize=9, fontweight='bold', pad=6)
    for row, rank in zip(rows, evidence['efs_ra_rank'], strict=True):
        ranks.text(0.5, row, f'{rank:,}', transform=ranks.get_yaxis_transform(), ha='center', va='center',
                   fontsize=9, fontweight='bold' if row == mrc1 else 'normal',
                   color=UNION_RED if row == mrc1 else INK,
                   bbox={'boxstyle': 'round,pad=0.27', 'facecolor': '#edf1f3', 'edgecolor': 'none'})
    fig.text(0.135, 0.083, 'Each draw resampled 4 D7 and 18 non-D7 libraries with replacement. Top 20 retains boundary ties.',
             color=GRAY, fontsize=9)
    fig.text(0.135, 0.052, f'Bootstrap pool: 3,907 Wald candidates. EFS-RA ranks: {universe_size:,} VST genes, full cohort.',
             color=GRAY, fontsize=9)
    return fig


def build_union_figures(output, ranking, vst, libraries):
    resampling = output / 'resampling/tables'
    detail = pd.read_csv(resampling / 'bootstrap_jaccard.tsv', sep='\t')
    frequency = pd.read_csv(resampling / 'bootstrap_selection_probability.tsv', sep='\t', index_col='基因ID')
    if len(frequency) != len(ranking) or not frequency.index.difference(ranking.index).empty:
        raise ValueError('Bootstrap probabilities and Wald candidate IDs differ')

    days = pd.Series({library: int(day) for day, group in libraries.items() for library in group}).sort_index()
    core = ranking.loc[ranking['基因符号'].isin(CORE_GENES)].sort_values('排名')
    if len(core) != len(CORE_GENES) or set(core['基因符号']) != set(CORE_GENES):
        raise ValueError('Expected one Ensembl ID per core gene')
    membership, overlap = screen_membership(vst.index)
    for screen, day in zip(SETS[:2], (3, 7), strict=True):
        passed = ranking.index[ranking[f'D{day}通过差异筛选']]
        if set(passed) != set(membership.index[membership[screen]]):
            raise ValueError(f'D{day} Wald results differ from the recorded candidate screen')

    # Bootstrap frequencies use the fixed Wald pool. Refit the same rankers on all VST genes for EFS-RA.
    rank = efs_ra_ranks(vst.sort_index().loc[:, days.index].T, days)
    evidence = core[['基因符号']].join(frequency[[*UNION_METHODS, 'Mean']]).join(rank).set_index('基因符号')
    if evidence.isna().any().any():
        raise ValueError('Missing bootstrap frequency or EFS-RA rank for a core gene')
    return {
        '01_candidate_union': plot_union(ranking, overlap),
        '02_candidate_expression': plot_union_expression(core, vst, days),
        '03_selection_stability': plot_union_stability(detail, evidence, len(vst)),
    }


def plot_candidate_pool(ranking, wald, use_union=False):
    passed = wald[wald['通过差异筛选'].eq('是')]['比较'].value_counts()
    shared = ranking.loc[ranking['D3通过差异筛选'] & ranking['D7通过差异筛选']] if use_union else ranking
    overlap = len(shared)
    counts = [passed['D3_vs_D0'] - overlap, overlap, passed['D7_vs_D0'] - overlap]
    fig, (bars, scatter) = plt.subplots(1, 2, figsize=(10.4, 5.1), layout='constrained',
                                       gridspec_kw={'width_ratios': [1, 1.65]})
    bars.barh(['仅 D3 对 D0', '两组交集', '仅 D7 对 D0'], counts, color=[LIGHT, BLUE, LIGHT], height=0.58)
    for i, count in enumerate(counts):
        bars.text(count + 55, i, f'{count:,}', va='center')
    bars.set(xlim=(0, max(counts) * 1.25), xlabel='基因数', title='A  通过差异筛选的基因')
    bars.invert_yaxis()
    bars.tick_params(axis='y', length=0)

    directions = [('均上调', '均上调', BLUE), ('均下调', '均下调', GRAY),
                  ('D3下调_D7上调', 'D3 下调、D7 上调', ORANGE),
                  ('D3上调_D7下调', 'D3 上调、D7 下调', ORANGE)]
    if use_union:
        single = ranking.loc[~ranking.index.isin(shared.index)]
        scatter.scatter(single['D3对D0_log2FC'], single['D7对D0_log2FC'], color=GRAY,
                        s=8, alpha=0.18, edgecolor='none', label=f'仅一个时点达标（{len(single):,}）')
    for key, label, color in directions:
        genes = shared[shared['相对D0方向'].eq(key)]
        reverse = key.startswith('D3')
        scatter.scatter(genes['D3对D0_log2FC'], genes['D7对D0_log2FC'], color=color,
                        s=32 if reverse else 10, marker='D' if reverse else 'o',
                        alpha=0.95 if reverse else 0.48, edgecolor='none',
                        label=f'{label}（{len(genes):,}）')
    scatter.axhline(0, color=LIGHT, zorder=0)
    scatter.axvline(0, color=LIGHT, zorder=0)
    pool = '并集' if use_union else '交集'
    scatter.set(xlabel='D3 对 D0 的 log₂FC', ylabel='D7 对 D0 的 log₂FC', title=f'B  {pool}内的效应方向')
    scatter.legend(loc='lower right', fontsize=8, frameon=True, facecolor='white', framealpha=0.93)
    fig.suptitle(f'D3 与 D7 相对 D0 的差异基因{pool}共有 {len(ranking):,} 个基因', fontweight='bold')
    return fig


def plot_expression(ranking, vst, libraries, repeats):
    genes = ranking.nsmallest(8, '排名')
    fig, panels = plt.subplots(2, 4, figsize=(11, 6), layout='constrained', sharey=True)
    days = sorted(libraries, key=int)
    reference_label = 'D3' if days == ['3', '7'] else '非D7'
    for ax, (gene_id, gene) in zip(panels.flat, genes.iterrows()):
        for x, day in enumerate(days):
            color = BLUE if day == '7' else GRAY
            values = vst.loc[gene_id, libraries[day]].to_numpy()
            offsets = np.linspace(-0.13, 0.13, len(values))
            ax.scatter(x + offsets, values, color=color, s=27, zorder=3)
            ax.hlines(values.mean(), x - 0.24, x + 0.24, color='#182833', linewidth=1.5, zorder=4)
        effect = gene[f'D7减{reference_label}_VST均值差']
        count = round(gene['均值差入选频率'] * repeats)
        label = gene['基因符号'] if pd.notna(gene['基因符号']) else gene_id
        ax.set_title(f'{label}\nD7−{reference_label} {effect:+.2f}，入选 {count}/{repeats}', fontsize=10)
        ax.set(xticks=range(len(days)), xticklabels=[f'D{day}' for day in days], xlim=(-0.45, len(days) - 0.55))
        ax.grid(axis='y', color=LIGHT, linewidth=0.7)
        ax.tick_params(axis='x', length=0)
    fig.suptitle('主排序前八个基因的逐文库表达', fontweight='bold')
    fig.supylabel('VST 表达值')
    counts = '、'.join(f'D{day} {len(libraries[day])} 份' for day in days)
    fig.supxlabel(f'每点为一份文库，短线为各时点均值。{counts}。', fontsize=10)
    return fig


def plot_stability(summary):
    methods = list(summary['全样本入选基因数'])
    repeats = summary['重抽样详情']
    offsets = np.linspace(-0.20, 0.20, len(repeats))
    fig, (overlap, sizes) = plt.subplots(1, 2, figsize=(10.4, 4.6), layout='constrained', sharey=True,
                                        gridspec_kw={'width_ratios': [1.3, 1]})
    for i, method in enumerate(methods):
        color = ORANGE if method == '随机森林' else BLUE if method == '均值差' else GRAY
        metrics = [row['方法'][method] for row in repeats]
        for ax, field in [(overlap, '与全样本名单Jaccard'), (sizes, '入选基因数')]:
            values = [metric[field] for metric in metrics]
            ax.scatter(values, i + offsets, color=color, s=27, zorder=3)
    labels = {'均值差': 'VST 均值差', 'ElasticNet': 'Elastic Net'}
    overlap.set(yticks=range(len(methods)), yticklabels=[labels.get(method, method) for method in methods],
                xlim=(-0.02, 1.02), xlabel='与全样本入选名单的 Jaccard', title='A  名单重合度')
    overlap.invert_yaxis()
    target = summary['每次目标入选数']
    sizes.axvline(target, color=GRAY, linestyle='--', linewidth=1)
    sizes.set(xlim=(0, None), xlabel='实际入选基因数', title='B  名单大小')
    for ax in (overlap, sizes):
        ax.grid(axis='x', color=LIGHT, linewidth=0.7)
        ax.tick_params(axis='y', length=0)
    fig.suptitle('八种文库扰动下的特征筛选稳定性', fontweight='bold')
    fig.supxlabel(f'每点对应一组文库扰动。虚线为目标 {target} 个，边界同分的基因一并入选。', fontsize=10)
    return fig


def visualize(output):
    ranking = pd.read_csv(output / 'tables/feature_ranking.tsv', sep='\t', index_col='基因ID')
    summary = json.loads((output / 'reports/selection_summary.json').read_text(encoding='utf-8'))
    vst = pd.read_csv(STAGE_ONE / 'matrices/timecourse_vst.tsv.gz', sep='\t', index_col='基因ID')
    figures = output / 'figures'
    figures.mkdir(parents=True, exist_ok=True)
    if '并集基因数' in summary:
        configure_style()
        plt.rcParams.update({'axes.titlesize': 12, 'axes.labelsize': 10})
        for name, fig in build_union_figures(output, ranking, vst, summary['文库']).items():
            fig.savefig(figures / f'{name}.png', dpi=300)
            plt.close(fig)
        return

    fonts = {font.name for font in font_manager.fontManager.ttflist}
    family = next((name for name in CHINESE_FONTS if name in fonts), None)
    if family is None:
        raise RuntimeError('未找到中文字体。请安装 Noto Sans CJK SC 后重新运行。')
    plt.rcParams.update({'font.family': [family, 'DejaVu Sans'], 'font.size': 10,
                         'axes.unicode_minus': False, 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.linewidth': 0.8})
    wald = pd.read_csv(STAGE_ONE / 'tables/timecourse_wald.tsv.gz', sep='\t')
    for name, fig in [('01_candidate_intersection', plot_candidate_pool(ranking, wald)),
                      ('02_candidate_expression', plot_expression(ranking, vst, summary['文库'], summary['重抽样次数'])),
                      ('03_selection_stability', plot_stability(summary))]:
        fig.savefig(figures / f'{name}.png', dpi=180)
        plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--union', action='store_true', help='绘制全部 22 份文库的并集筛选结果')
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    visualize(args.output_dir or (UNION_OUTPUT if args.union else OUTPUT))
