"""组合阶段二的筛选、重采样和生物学互证图。

读取已导出的统计表和表达矩阵，写入相应目录下的 PNG。
"""

import json

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patheffects import withStroke
from matplotlib.ticker import MultipleLocator, PercentFormatter, StrMethodFormatter
import numpy as np
import pandas as pd

from .data_loader import SETS, read_ranking, read_vst, read_wald, window
from .plot_style import (
    BLUE,
    GRAY,
    INK,
    LIGHT,
    ORANGE,
    UNION_BLUE,
    UNION_GRAY,
    UNION_METHODS,
    UNION_RED,
    WIDTH_MM,
    configure_intersection_style,
    configure_panel_style,
    configure_style,
    panel_heading,
    save_panel_figure,
    union_axes,
    union_figure,
)
from ..bio_validation.marker_coexpression import draw_markers
from ..bio_validation.ssgsea_pathway import draw_pathways
from ..config import CORE_GENES, METHODS
from ..methods.m05_efs_rank_agg import efs_ra_ranks


def plot_union(ranking, overlap):
    """用 Wald 重叠统计和候选表生成并集效应图。"""
    d3, d7 = (ranking[f'D{day}通过差异筛选'] for day in (3, 7))
    shared = d3 & d7
    fig = union_figure(6.0, 'Differential gene overlap and effect trajectories')
    grid = fig.add_gridspec(
        2,
        2,
        left=0.135,
        right=0.963,
        top=0.84,
        bottom=0.22,
        height_ratios=(3.8, 1),
        width_ratios=(0.92, 1.17),
        hspace=0.05,
        wspace=0.23,
    )
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
        bars.annotate(
            f'{count:,}',
            (column, count),
            xytext=(0, 4),
            textcoords='offset points',
            ha='center',
            va='bottom',
            fontsize=9,
            color=UNION_RED if color == UNION_RED else INK,
        )
    bars.text(
        0.04,
        0.92,
        f'D3 ∩ D7   {shared.sum():,} genes  ·  {shared.mean():.1%} of union',
        transform=bars.transAxes,
        va='top',
        color=UNION_RED,
        fontsize=10,
        fontweight='bold',
    )
    bars.set(
        xlim=(-0.55, len(x) - 0.45),
        ylim=(0, exact['Genes'].max() * 1.18),
        ylabel='Genes in exact intersection',
        yticks=[0, 2000, 4000, 6000, 8000],
    )
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
    scatter.text(
        0.035,
        0.95,
        'D7 effect > D3 effect',
        transform=scatter.transAxes,
        va='top',
        color=UNION_BLUE,
        fontsize=9,
    )
    direction = ranking['相对D0方向']
    same = shared & direction.isin(['均上调', '均下调'])
    for mask, label, color, marker, size, alpha in [
        (~shared, 'One contrast', UNION_GRAY, 'o', 7, 0.28),
        (same, 'Both, same sign', UNION_BLUE, 'o', 10, 0.48),
        (shared & ~same, 'Both, sign reversal', UNION_RED, 'D', 29, 1),
    ]:
        points = ranking.loc[mask]
        scatter.scatter(
            points[x_field],
            points[y_field],
            color=color,
            marker=marker,
            s=size,
            alpha=alpha,
            linewidths=0,
            label=label,
            zorder=3,
        )
    spotlight = ranking.set_index('基因符号')
    for gene, offset in [('C1qtnf3', (12, 20)), ('Postn', (17, 9)), ('Cthrc1', (9, 21))]:
        point = tuple(spotlight.loc[gene, [x_field, y_field]])
        scatter.plot(*point, marker='o', ms=5.5, mfc=UNION_BLUE, mec='white', mew=0.8, zorder=4)
        scatter.annotate(
            gene,
            point,
            xytext=offset,
            textcoords='offset points',
            color=UNION_BLUE,
            fontsize=10,
            fontweight='bold',
            fontstyle='italic',
            va='center',
            zorder=5,
            path_effects=[withStroke(linewidth=2, foreground='white')],
            arrowprops={'arrowstyle': '-', 'color': UNION_BLUE, 'lw': 0.8},
        )
    scatter.set(
        xlim=(x_min, x_max),
        ylim=(y_min, y_max),
        xlabel='D3 vs D0 (log$_2$ fold change)',
        ylabel='D7 vs D0 (log$_2$ fold change)',
    )
    scatter.spines[['left', 'bottom']].set_color(UNION_GRAY)
    scatter.tick_params(labelsize=9)
    scatter.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.13),
        ncols=3,
        fontsize=8.5,
        handletextpad=0.25,
        columnspacing=1.0,
        scatterpoints=1,
        markerscale=1.4,
    )
    fig.text(
        0.135,
        0.073,
        'Wald screens: BH q < 0.05 and |log$_2$FC| ≥ 1.   LRT: BH q < 0.05.   Dashed line: y = x.',
        fontsize=9,
        color=GRAY,
    )
    return fig


def plot_expression_trajectories(genes, vst, days):
    """共享 VST 纵轴，显示逐文库表达、各采样日均值和 D7 窗口。"""
    time_points = sorted(days.unique())
    groups = [days.index[days.eq(day)] for day in time_points]
    positions = np.arange(len(time_points))
    d7 = time_points.index(7)
    expression = vst.loc[genes.index, days.index]
    low = np.floor(expression.min().min() - 0.15)
    high = np.ceil(expression.max().max() + 0.15)

    fig, panels = plt.subplots(3, 3, figsize=(11, 8.8), sharex=True, sharey=True, facecolor='white')
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.13, top=0.865, wspace=0.21, hspace=0.49)
    fig.suptitle(
        'Candidate expression trajectories',
        x=0.10,
        y=0.955,
        ha='left',
        fontsize=19,
        fontweight='bold',
        color=INK,
    )
    for index, (gene_id, gene) in enumerate(genes['基因符号'].items()):
        row, column = divmod(index, 3)
        ax = panels[row, column]
        gene_name = gene if pd.notna(gene) else gene_id
        accent = gene_name == 'Mrc1'
        color = UNION_RED if accent else UNION_BLUE

        # 横轴按采样顺序等距排列。同一组内轻微错开散点，避免样本相互遮挡。
        values_by_day = [expression.loc[gene_id, libraries].to_numpy() for libraries in groups]
        ax.axvspan(d7 - 0.27, d7 + 0.27, color='#efe9e2', lw=0, zorder=0)
        ax.plot(positions, [values.mean() for values in values_by_day], color=color, lw=2.1)
        for position, values in zip(positions, values_by_day, strict=True):
            offsets = np.linspace(-0.16, 0.16, len(values)) if len(values) > 1 else 0.0
            ax.scatter(
                position + offsets,
                values,
                s=48,
                color=color,
                edgecolors='white',
                linewidths=0.7,
                zorder=3,
            )

        ax.set_title(
            gene_name,
            fontstyle='italic',
            fontweight='bold' if accent else 'normal',
            color=INK,
            fontsize=17,
            pad=9,
        )
        ax.grid(False)
        ax.spines[['top', 'right']].set_visible(False)
        for side in ('left', 'bottom'):
            ax.spines[side].set(color=INK, linewidth=1.45)
        ax.tick_params(
            axis='both',
            direction='out',
            length=5,
            width=1.45,
            colors=INK,
            labelsize=12,
            labelbottom=row == 2,
            labelleft=column == 0,
        )

    for ax in panels.flat[len(genes):]:
        ax.set_visible(False)
    panels[0, 0].set(
        xlim=(-0.35, len(time_points) - 0.65),
        ylim=(low, high),
        xticks=positions,
        xticklabels=[str(day) for day in time_points],
    )
    panels[0, 0].yaxis.set_major_locator(MultipleLocator(2))
    fig.supxlabel('Days after injury', x=0.54, y=0.04, fontsize=16, color=INK)
    fig.supylabel('VST expression', x=0.03, y=0.50, fontsize=16, color=INK)
    return fig


def plot_union_stability(detail, evidence, universe_size):
    """绘制四方法的 Bootstrap Jaccard、候选入选率与全 VST 排名。"""
    groups = [detail.loc[detail['Method'].eq(method), 'Jaccard'].to_numpy() for method in UNION_METHODS]
    draws = detail['Draw'].nunique()
    if any(len(group) != draws or not np.isfinite(group).all() for group in groups):
        raise ValueError('Expected one defined Jaccard value per bootstrap draw and method')
    fig = union_figure(6.15, 'Selection stability and core gene support')
    distribution = fig.add_subplot(
        fig.add_gridspec(1, 1, left=0.135, right=0.424, bottom=0.199, top=0.815)[0]
    )
    grid = fig.add_gridspec(
        1, 2, left=0.588, right=0.963, bottom=0.199, top=0.815, width_ratios=(5, 0.86), wspace=0.19
    )
    frequency = fig.add_subplot(grid[0, 0])
    ranks = fig.add_subplot(grid[0, 1], sharey=frequency)
    distribution.set_title(
        'A   Overlap with full-cohort selection', loc='left', fontweight='bold', pad=39
    )
    frequency.set_title('B   Core gene selection and rank', loc='left', fontweight='bold', pad=39)

    positions = np.arange(len(groups))
    violins = distribution.violinplot(
        groups, positions=positions, orientation='horizontal', side='low', widths=0.68, showextrema=False
    )
    for density, (_, color, _) in zip(violins['bodies'], UNION_METHODS.values(), strict=True):
        density.set_facecolor(color)
        density.set_edgecolor(color)
        density.set_alpha(0.30)
        density.set_linewidth(0.7)
    boxes = distribution.boxplot(
        groups,
        positions=positions + 0.11,
        orientation='horizontal',
        widths=0.19,
        patch_artist=True,
        manage_ticks=False,
        zorder=3,
        medianprops={'color': INK, 'linewidth': 1.4},
        flierprops={'marker': '.', 'markersize': 2.3, 'alpha': 0.35, 'markeredgecolor': UNION_GRAY},
    )
    for position, values, box, (_, color, _) in zip(
        positions, groups, boxes['boxes'], UNION_METHODS.values(), strict=True
    ):
        box.set(facecolor=color, edgecolor=INK, alpha=0.62, linewidth=0.8)
        distribution.text(
            1.05,
            position,
            f'Median = {np.median(values):.2f}',
            transform=distribution.get_yaxis_transform(),
            color=color,
            fontsize=9,
            va='center',
            fontweight='bold',
        )
    distribution.set(
        xlim=(0, 1),
        ylim=(len(groups) - 0.52, -0.52),
        yticks=positions,
        yticklabels=['VST mean\ndifference', 'LASSO', 'Elastic Net', 'Random forest'],
        xticks=np.linspace(0, 1, 6),
        xlabel='Jaccard overlap',
    )
    union_axes(distribution, 'x')
    distribution.tick_params(axis='y', length=0, pad=6)

    evidence = evidence.sort_values('efs_ra_rank')
    rows = np.arange(len(evidence))
    mrc1 = evidence.index.get_loc('Mrc1')
    frequency.axhspan(mrc1 - 0.45, mrc1 + 0.45, color='#f8eeea', lw=0, zorder=0)
    frequency.hlines(rows, 0, evidence['Mean'], color=UNION_GRAY, lw=0.9, zorder=2)
    frequency.scatter(evidence['Mean'], rows, s=19, c=INK, linewidths=0.4, edgecolors='white', zorder=3)
    handles = []
    for offset, (key, (label, color, marker)) in zip(
        (-0.29, -0.10, 0.10, 0.29), UNION_METHODS.items(), strict=True
    ):
        frequency.scatter(
            evidence[key],
            rows + offset,
            color=color,
            marker=marker,
            s=24,
            linewidths=0.45,
            edgecolors='white',
            clip_on=False,
            zorder=4,
        )
        handles.append(
            Line2D([], [], marker=marker, color=color, linestyle='none', markersize=5.1, label=label)
        )
    handles.append(
        Line2D(
            [],
            [],
            marker='o',
            color=UNION_GRAY,
            markerfacecolor=INK,
            markeredgecolor=INK,
            markersize=3.6,
            lw=0.8,
            label='Four-method mean',
        )
    )
    frequency.set(
        xlim=(0, 1),
        ylim=(len(rows) - 0.5, -0.5),
        xticks=np.linspace(0, 1, 5),
        xlabel='Bootstrap selection frequency',
        yticks=rows,
        yticklabels=evidence.index,
    )
    frequency.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    union_axes(frequency, 'x')
    frequency.spines['left'].set_visible(False)
    frequency.tick_params(axis='y', length=0, pad=8)
    for row, label in enumerate(frequency.get_yticklabels()):
        label.set_fontstyle('italic')
        if row == mrc1:
            label.set(color=UNION_RED, fontweight='bold')
    frequency.legend(
        handles=handles,
        loc='lower left',
        bbox_to_anchor=(-0.03, 1.008),
        ncols=3,
        fontsize=8,
        handlelength=1.2,
        columnspacing=1.05,
        handletextpad=0.36,
        borderpad=0.15,
    )

    ranks.set_axis_off()
    ranks.set_title('EFS-RA\nrank', fontsize=9, fontweight='bold', pad=6)
    for row, rank in zip(rows, evidence['efs_ra_rank'], strict=True):
        ranks.text(
            0.5,
            row,
            f'{rank:,}',
            transform=ranks.get_yaxis_transform(),
            ha='center',
            va='center',
            fontsize=9,
            fontweight='bold' if row == mrc1 else 'normal',
            color=UNION_RED if row == mrc1 else INK,
            bbox={'boxstyle': 'round,pad=0.27', 'facecolor': '#edf1f3', 'edgecolor': 'none'},
        )
    fig.text(
        0.135,
        0.083,
        'Each draw resampled 4 D7 and 18 non-D7 libraries with replacement. Top 20 retains boundary ties.',
        color=GRAY,
        fontsize=9,
    )
    fig.text(
        0.135,
        0.052,
        f'Bootstrap pool: 3,907 Wald candidates. EFS-RA ranks: {universe_size:,} VST genes, full cohort.',
        color=GRAY,
        fontsize=9,
    )
    return fig


def build_union_figures(output, ranking, vst, days):
    """核对并集、VST 和重采样表，返回三张待保存的并集图。"""
    resampling = output / 'resampling/tables'
    detail = pd.read_csv(resampling / 'bootstrap_jaccard.tsv', sep='\t')
    frequency = pd.read_csv(
        resampling / 'bootstrap_selection_probability.tsv', sep='\t', index_col='基因ID'
    )
    if len(frequency) != len(ranking) or not frequency.index.difference(ranking.index).empty:
        raise ValueError('Bootstrap probabilities and Wald candidate IDs differ')

    gene_order = {gene: index for index, gene in enumerate(CORE_GENES)}
    core = ranking.loc[ranking['基因符号'].isin(CORE_GENES)].sort_values(
        '基因符号', key=lambda symbols: symbols.map(gene_order)
    )
    if len(core) != len(CORE_GENES) or set(core['基因符号']) != set(CORE_GENES):
        raise ValueError('Expected one Ensembl ID per core gene')
    membership, overlap = window.screen_membership(vst.index)
    for screen, day in zip(SETS[:2], (3, 7), strict=True):
        passed = ranking.index[ranking[f'D{day}通过差异筛选']]
        if set(passed) != set(membership.index[membership[screen]]):
            raise ValueError(f'D{day} Wald results differ from the recorded candidate screen')

    # Bootstrap frequencies use the fixed Wald pool. Refit the same rankers on all VST genes for EFS-RA.
    rank = efs_ra_ranks(vst.sort_index().loc[:, days.index].T, days)
    evidence = (
        core[['基因符号']].join(frequency[[*UNION_METHODS, 'Mean']]).join(rank).set_index('基因符号')
    )
    if evidence.isna().any().any():
        raise ValueError('Missing bootstrap frequency or EFS-RA rank for a core gene')
    return {
        '01_candidate_union': plot_union(ranking, overlap),
        '02_candidate_expression': plot_expression_trajectories(core, vst, days),
        '03_selection_stability': plot_union_stability(detail, evidence, len(vst)),
    }


def plot_candidate_pool(ranking, wald, use_union=False):
    """绘制 D3/D7 交集规模和两时点效应方向。"""
    passed = wald[wald['通过差异筛选'].eq('是')]['比较'].value_counts()
    shared = ranking.loc[ranking['D3通过差异筛选'] & ranking['D7通过差异筛选']] if use_union else ranking
    overlap = len(shared)
    counts = [passed['D3_vs_D0'] - overlap, overlap, passed['D7_vs_D0'] - overlap]
    fig, (bars, scatter) = plt.subplots(
        1, 2, figsize=(10.4, 5.1), layout='constrained', gridspec_kw={'width_ratios': [1, 1.65]}
    )
    bars.barh(
        ['仅 D3 对 D0', '两组交集', '仅 D7 对 D0'], counts, color=[LIGHT, BLUE, LIGHT], height=0.58
    )
    for i, count in enumerate(counts):
        bars.text(count + 55, i, f'{count:,}', va='center')
    bars.set(xlim=(0, max(counts) * 1.25), xlabel='基因数', title='A  通过差异筛选的基因')
    bars.invert_yaxis()
    bars.tick_params(axis='y', length=0)

    directions = [
        ('均上调', '均上调', BLUE),
        ('均下调', '均下调', GRAY),
        ('D3下调_D7上调', 'D3 下调、D7 上调', ORANGE),
        ('D3上调_D7下调', 'D3 上调、D7 下调', ORANGE),
    ]
    if use_union:
        single = ranking.loc[~ranking.index.isin(shared.index)]
        scatter.scatter(
            single['D3对D0_log2FC'],
            single['D7对D0_log2FC'],
            color=GRAY,
            s=8,
            alpha=0.18,
            edgecolor='none',
            label=f'仅一个时点达标（{len(single):,}）',
        )
    for key, label, color in directions:
        genes = shared[shared['相对D0方向'].eq(key)]
        reverse = key.startswith('D3')
        scatter.scatter(
            genes['D3对D0_log2FC'],
            genes['D7对D0_log2FC'],
            color=color,
            s=32 if reverse else 10,
            marker='D' if reverse else 'o',
            alpha=0.95 if reverse else 0.48,
            edgecolor='none',
            label=f'{label}（{len(genes):,}）',
        )
    scatter.axhline(0, color=LIGHT, zorder=0)
    scatter.axvline(0, color=LIGHT, zorder=0)
    pool = '并集' if use_union else '交集'
    scatter.set(xlabel='D3 对 D0 的 log₂FC', ylabel='D7 对 D0 的 log₂FC', title=f'B  {pool}内的效应方向')
    scatter.legend(loc='lower right', fontsize=8, frameon=True, facecolor='white', framealpha=0.93)
    fig.suptitle(f'D3 与 D7 相对 D0 的差异基因{pool}共有 {len(ranking):,} 个基因', fontweight='bold')
    return fig


def plot_stability(summary):
    """绘制八次文库排除下的名单 Jaccard 和入选数量。"""
    methods = list(summary['全样本入选基因数'])
    repeats = summary['重抽样详情']
    offsets = np.linspace(-0.20, 0.20, len(repeats))
    fig, (overlap, sizes) = plt.subplots(
        1,
        2,
        figsize=(10.4, 4.6),
        layout='constrained',
        sharey=True,
        gridspec_kw={'width_ratios': [1.3, 1]},
    )
    for i, method in enumerate(methods):
        color = ORANGE if method == '随机森林' else BLUE if method == '均值差' else GRAY
        metrics = [row['方法'][method] for row in repeats]
        for ax, field in [(overlap, '与全样本名单Jaccard'), (sizes, '入选基因数')]:
            values = [metric[field] for metric in metrics]
            ax.scatter(values, i + offsets, color=color, s=27, zorder=3)
    labels = {'均值差': 'VST 均值差', 'ElasticNet': 'Elastic Net'}
    overlap.set(
        yticks=range(len(methods)),
        yticklabels=[labels.get(method, method) for method in methods],
        xlim=(-0.02, 1.02),
        xlabel='与全样本入选名单的 Jaccard',
        title='A  名单重合度',
    )
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


def plot_selection(output):
    """根据现有筛选结果生成对应交集或并集版本的图。"""
    ranking = read_ranking(output)
    summary = json.loads((output / 'reports/selection_summary.json').read_text(encoding='utf-8'))
    vst = read_vst()
    days = pd.Series(
        {library: int(day) for day, group in summary['文库'].items() for library in group}
    ).sort_index()
    figures = output / 'figures'
    figures.mkdir(parents=True, exist_ok=True)
    if '并集基因数' in summary:
        configure_style()
        plt.rcParams.update({'axes.titlesize': 12, 'axes.labelsize': 10})
        for name, fig in build_union_figures(output, ranking, vst, days).items():
            fig.savefig(figures / f'{name}.png', dpi=300)
            plt.close(fig)
        return

    configure_intersection_style()
    wald = read_wald()
    for name, fig in [
        ('01_candidate_intersection', plot_candidate_pool(ranking, wald)),
        (
            '02_candidate_expression',
            plot_expression_trajectories(ranking.nsmallest(8, '排名'), vst, days),
        ),
        ('03_selection_stability', plot_stability(summary)),
    ]:
        fig.savefig(figures / f'{name}.png', dpi=180)
        plt.close(fig)


def plot_resampling(output):
    """读取 Bootstrap 结果并绘制分布及平均入选概率最高的候选。"""
    detail = pd.read_csv(output / 'tables/bootstrap_jaccard.tsv', sep='\t')
    probabilities = pd.read_csv(
        output / 'tables/bootstrap_selection_probability.tsv', sep='\t', skipinitialspace=True
    )
    summary = json.loads((output / 'reports/resampling_summary.json').read_text())
    top = probabilities.head(20)
    labels = top['基因符号'].fillna(top['基因ID'])
    values = top['Mean'].to_numpy()
    groups = [detail.loc[detail['Method'].eq(method), 'Jaccard'].dropna() for method in METHODS]
    colors = ['#286485', '#8B6456', '#786894', '#547B71']
    fig, (box, lollipop) = plt.subplots(
        1, 2, figsize=(12, 7.3), layout='constrained', gridspec_kw={'width_ratios': [0.95, 1]}
    )
    boxes = box.boxplot(
        groups,
        tick_labels=['VST mean\ndifference', 'LASSO', 'Elastic Net', 'Random\nforest'],
        patch_artist=True,
        medianprops={'color': '#243744', 'linewidth': 1.6},
        flierprops={'markersize': 2, 'alpha': 0.48},
    )
    for patch, color in zip(boxes['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
    box.set(ylim=(-0.025, 1.025), ylabel='Jaccard with full-cohort selection')
    box.set_title(
        'a   Selection overlap across bootstrap draws', loc='left', fontsize=12, fontweight='bold'
    )
    box.grid(axis='y', color='#E1E8ED', linewidth=0.75)

    rows = np.arange(len(top))
    lollipop.hlines(rows, 0, values, color='#A7C0CE', linewidth=1.5)
    lollipop.scatter(values, rows, color='#286485', s=33, zorder=3)
    for row, value in zip(rows, values):
        lollipop.annotate(
            f'{value:.1%}',
            (value, row),
            xytext=(5, 0),
            textcoords='offset points',
            va='center',
            fontsize=8,
            color='#243744',
        )
    lollipop.set(
        yticks=rows, yticklabels=labels, xlim=(0, 1.12), xlabel='Mean Top 20 selection probability'
    )
    lollipop.set_xticks(np.linspace(0, 1, 6))
    lollipop.xaxis.set_major_formatter(PercentFormatter(1))
    lollipop.invert_yaxis()
    lollipop.set_title(
        'b   Top genes by selection probability', loc='left', fontsize=12, fontweight='bold'
    )
    lollipop.grid(axis='x', color='#E1E8ED', linewidth=0.75)
    for ax in (box, lollipop):
        ax.set_axisbelow(True)
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(labelsize=9)
    count = summary['bootstrap_repeats']
    fig.suptitle(
        f'Selection stability across {count:,} stratified bootstrap draws',
        fontsize=16,
        fontweight='bold',
    )
    fig.supxlabel(
        'Each draw samples 4 D7 and 18 non-D7 libraries with replacement.\n'
        'Selection includes boundary ties. Probabilities are averaged across four methods.',
        fontsize=9,
    )
    figures = output / 'figures'
    figures.mkdir(parents=True, exist_ok=True)
    fig.savefig(figures / '03_selection_stability.png', dpi=220)
    plt.close(fig)


def plot_biology(output):
    """读入通路和标志物结果，绘制两类相关系数的对照热图。"""
    tables = output / 'tables'
    pathways = pd.read_csv(tables / 'candidate_pathway_correlations.tsv', sep='\t')
    markers = pd.read_csv(tables / 'candidate_marker_correlations.tsv', sep='\t')
    n_libraries = len(pd.read_csv(tables / 'candidate_marker_expression.tsv', sep='\t'))
    configure_panel_style()
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, 91 / 25.4))
    left = fig.add_axes((0.11, 0.245, 0.42, 0.50))
    right = fig.add_axes((0.55, 0.245, 0.42, 0.50))
    image = draw_pathways(left, pathways)
    draw_markers(right, markers)
    panel_heading(fig, 0.07, 0.927, 'a', 'ssGSEA pathway scores')
    panel_heading(fig, 0.54, 0.927, 'b', 'Marker coexpression')
    cax = fig.add_axes((0.35, 0.166, 0.34, 0.024))
    cb = fig.colorbar(image, cax=cax, orientation='horizontal', ticks=[-1, -0.5, 0, 0.5, 1])
    cb.ax.tick_params(labelsize=6, length=2, width=0.5, pad=2)
    cb.outline.set_linewidth(0.45)
    fig.text(0.32, 0.176, 'Spearman ρ', ha='right', va='center', fontsize=6.5)
    fig.text(
        0.54,
        0.063,
        f'{n_libraries} libraries pooled across six time points  ·  BH q < 0.05 (*) or < 0.01 (**)',
        ha='center',
        va='center',
        fontsize=6.5,
        color='#465865',
    )
    save_panel_figure(fig, output / 'figures', 'candidate_biology_correlations')
    plt.close(fig)
