"""从已有外部评估结果绘制表达分布、AUC 和相关性图。

读取步骤 01 的样本纳入表及步骤 02 的 TSV。
保存两队列的 3×3 表达图、AUC 对照图及通路和标志物相关点图。
"""

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator
import numpy as np

from .config import COHORTS, COHORT_SCALES, GENES
from .tables import read


def boxplot(cohort, x, meta, scale, tests, output):
    """绘制所给队列的九候选分组箱线图，使用各时间组表达及主要比较的 Welch P。"""
    fig, axes = plt.subplots(3, 3, figsize=(11, 8), layout='constrained')
    order = list(meta.Group.unique())
    for ax, gene in zip(axes.flat, GENES):
        groups = [x.loc[gene, meta.Group.eq(group)].to_numpy() for group in order]
        ax.boxplot(groups, tick_labels=[group.replace('MI_', '') for group in order], showfliers=False)
        for j, values in enumerate(groups, 1):
            ax.scatter(
                j + np.linspace(-0.11, 0.11, len(values)),
                values,
                s=21,
                color='#a65439' if j > 1 else '#276687',
                zorder=3,
            )
        p = tests.loc[tests.Gene.eq(gene) & tests.Primary, 'P'].iat[0]
        ax.set_title(f'{gene}   D7 P={p:.3g}', fontsize=10)
        ax.tick_params(axis='x', labelsize=8)
    fig.suptitle(f'{cohort} | candidate expression by time')
    fig.supylabel(scale)
    fig.savefig(output / 'figures' / f'{cohort}_expression.png', dpi=200)
    plt.close(fig)


def plot_expression_auc(scores_by_cohort, metrics, output):
    """用已保存的单基因预测分数画表达分布及对应 AUC。"""
    available = {f.name for f in font_manager.fontManager.ttflist}
    font = next(
        (
            name
            for name in ('Microsoft YaHei', 'Noto Sans CJK SC', 'PingFang SC', 'SimHei')
            if name in available
        ),
        None,
    )
    ink, blue, orange = '#192d3b', '#256080', '#b5543a'
    auc = metrics.set_index(['Cohort', 'Score']).AUROC
    style = {
        'font.family': ['DejaVu Sans', font] if font else ['DejaVu Sans'],
        'text.color': ink,
        'axes.labelcolor': ink,
        'xtick.color': ink,
        'ytick.color': ink,
    }
    with plt.rc_context(style):
        fig = plt.figure(figsize=(14, 10), facecolor='white')
        grid = fig.add_gridspec(
            3, 3, left=0.055, right=0.98, top=0.90, bottom=0.085, wspace=0.31, hspace=0.32
        )
        fig.suptitle(
            '九个候选基因的外部表达验证'
            if font
            else 'External expression validation of nine candidates',
            y=0.971,
            fontsize=24,
            fontweight='bold',
        )
        for cell, gene in zip(grid, GENES):
            inner = cell.subgridspec(2, 2, height_ratios=[0.19, 1], hspace=0.20, wspace=0.42)
            title = fig.add_subplot(inner[0, :])
            title.set_axis_off()
            title.text(
                0.5,
                0.47,
                gene,
                ha='center',
                va='center',
                fontsize=15,
                fontstyle='italic',
                fontweight='bold',
                color=orange if gene == 'Cxadr' else ink,
            )
            for j, (cohort, control) in enumerate(COHORTS):
                ax = fig.add_subplot(inner[1, j])
                score = scores_by_cohort[cohort]
                groups = [score.loc[score.Label.eq(label), gene].to_numpy() for label in (0, 1)]
                boxes = ax.boxplot(
                    groups,
                    tick_labels=['未损伤' if control == 'uninjured' else 'Sham', 'D7']
                    if font
                    else [control, 'D7'],
                    showfliers=False,
                    showmeans=True,
                    meanline=True,
                    patch_artist=True,
                    widths=0.47,
                    medianprops={'alpha': 0},
                    meanprops={'linestyle': '-', 'linewidth': 2, 'color': ink},
                )
                for pos, values, color, box in zip((1, 2), groups, (blue, orange), boxes['boxes']):
                    box.set(facecolor=color, edgecolor=color, alpha=0.20)
                    ax.scatter(
                        pos + np.linspace(-0.14, 0.14, len(values)),
                        values,
                        s=23,
                        color=color,
                        edgecolor='white',
                        linewidth=0.5,
                        zorder=3,
                    )
                values = np.concatenate(groups)
                span = np.ptp(values) or 1
                ax.set_ylim(values.min() - 0.12 * span, values.max() + 0.38 * span)
                ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
                ax.ticklabel_format(axis='y', style='plain', useOffset=False)
                ax.set_title(cohort, fontsize=9.5, pad=5)
                ax.tick_params(labelsize=8, length=3, pad=2)
                ax.set_axisbelow(True)
                ax.grid(axis='y', color='#e4e9ef', linewidth=0.65)
                ax.spines[['top', 'right']].set_visible(False)
                ax.text(
                    0.97,
                    0.95,
                    f'AUC = {auc.loc[(cohort, gene)]:.2f}',
                    transform=ax.transAxes,
                    ha='right',
                    va='top',
                    fontsize=9,
                    fontweight='bold',
                    color=orange if gene == 'Cxadr' else ink,
                    bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': 1.4},
                )
        note = (
            '散点为单份样本，深色横线为组均值。各队列纵轴独立（左：SST-RMA log₂ 信号，右：log₂(CPM+1)）。'
            if font
            else 'Dots show samples; dark lines show means. Scales: left SST-RMA log2, right log2(CPM+1).'
        )
        fig.text(0.5, 0.025, note, ha='center', fontsize=10, color='#455b6b')
        fig.savefig(output / 'figures/bulk_roc.png', dpi=225)
        plt.close(fig)


def plot_correlation_dots(tables, name, output):
    """并排呈现两个外部队列的 Spearman 相关系数和 BH 校正显著性。"""
    available = {f.name for f in font_manager.fontManager.ttflist}
    font = next(
        (f for f in ('Microsoft YaHei', 'Noto Sans CJK SC', 'PingFang SC', 'SimHei') if f in available),
        None,
    )
    ink, muted = '#172b3a', '#4a565f'
    cmap = LinearSegmentedColormap.from_list('correlation', ['#256080', '#f3f3f3', '#a6492e'])
    is_pathway = name == 'pathway_correlations'
    targets = tables[0].Target.drop_duplicates().tolist()
    pathway_labels = {
        'NLRP3 inflammasome': 'NLRP3 炎性小体' if font else 'NLRP3 inflammasome',
        'Extracellular matrix organization': '细胞外基质组织'
        if font
        else 'Extracellular matrix\norganization',
        'Collagen formation': '胶原形成' if font else 'Collagen formation',
    }
    labels = [pathway_labels.get(target, target) for target in targets]
    title = '候选基因与 ssGSEA 通路活性的相关性' if is_pathway else '候选基因与经典标志物的共表达'
    if not font:
        title = (
            'Candidate genes and ssGSEA pathway activity'
            if is_pathway
            else 'Candidate gene and marker co-expression'
        )
    family = ['DejaVu Sans', font] if font else ['DejaVu Sans']
    with plt.rc_context(
        {'font.family': family, 'text.color': ink, 'xtick.color': ink, 'ytick.color': ink}
    ):
        fig, axes = plt.subplots(1, 2, figsize=(12, 6.6), sharey=True, facecolor='white')
        fig.subplots_adjust(left=0.125, right=0.96, bottom=0.19, top=0.80, wspace=0.18)
        fig.text(0.125, 0.94, title, fontsize=20, fontweight='bold')
        for ax, table in zip(axes, tables):
            rho = (
                table.pivot(index='Gene', columns='Target', values='Spearman_rho')
                .reindex(index=GENES, columns=targets)
                .to_numpy()
            )
            q = (
                table.pivot(index='Gene', columns='Target', values='BH_q')
                .reindex(index=GENES, columns=targets)
                .to_numpy()
            )
            rows, columns = np.where(np.isfinite(rho))
            values = rho[rows, columns]
            # 用同一色标和面积比例比较两个队列，保留零相关的小圆点。
            dots = ax.scatter(
                columns,
                rows,
                s=16 + 132 * np.abs(values),
                c=values,
                cmap=cmap,
                vmin=-1,
                vmax=1,
                edgecolor='#7b8790',
                linewidth=0.45,
            )
            for i, j in np.ndindex(rho.shape):
                value = rho[i, j]
                stars = '**' if q[i, j] < 0.01 else '*' if q[i, j] < 0.05 else ''
                if np.isfinite(value):
                    ax.annotate(
                        f'{value:.2f}{stars}',
                        (j, i),
                        xytext=(0, -7),
                        textcoords='offset points',
                        ha='center',
                        va='top',
                        fontsize=9,
                    )
                else:
                    ax.text(j, i, 'NA', ha='center', va='center', fontsize=9, color=muted)
            ax.set(xlim=(-0.5, len(targets) - 0.5), ylim=(len(GENES) - 0.33, -0.58))
            ax.set_xticks(range(len(targets)), labels)
            ax.set_yticks(range(len(GENES)), GENES, fontstyle='italic')
            ax.xaxis.tick_top()
            ax.tick_params(axis='both', length=0, labelsize=10, pad=8)
            ax.set_title(table.Cohort.iat[0], y=1.095, fontsize=13, fontweight='bold')
            ax.spines[:].set_visible(False)
            ax.spines['top'].set_visible(True)
            ax.spines['top'].set_color('#909da6')
        cax = fig.add_axes((0.125, 0.125, 0.36, 0.026))
        fig.colorbar(dots, cax=cax, orientation='horizontal', ticks=[-1, -0.5, 0, 0.5, 1])
        cax.tick_params(labelsize=9)
        fig.text(
            0.51,
            0.138,
            'Spearman ρ  |  圆点越大，相关越强'
            if font
            else 'Spearman ρ  |  Dot area increases with |ρ|',
            fontsize=10,
            va='center',
        )
        fig.text(
            0.125,
            0.038,
            '* BH q < 0.05，** BH q < 0.01' if font else '* BH q < 0.05, ** BH q < 0.01',
            fontsize=10,
            color=muted,
        )
        fig.savefig(output / 'figures' / f'{name}.png', dpi=240)
        plt.close(fig)


def run(preprocessing, output):
    """读取外部验证 TSV，使用原分组、取值尺度和版式生成图。"""
    (output / 'figures').mkdir(parents=True, exist_ok=True)
    tables = output / 'tables'
    # 可往返解析保留步骤 02 写出的浮点数，避免绘图重新计算统计量。
    tests = read(tables / 'bulk_expression_tests.tsv', float_precision='round_trip')
    metrics = read(tables / 'bulk_roc_metrics.tsv', float_precision='round_trip')
    scores = {}
    for cohort, _ in COHORTS:
        expression = read(
            tables / f'{cohort}_frozen_expression.tsv', index_col='Gene', float_precision='round_trip'
        )
        meta = read(preprocessing / cohort / 'tables/sample_inclusion.tsv', index_col='GSM')
        cohort_tests = tests[tests.Cohort.eq(cohort)]
        boxplot(cohort, expression, meta, COHORT_SCALES[cohort], cohort_tests, output)
        scores[cohort] = read(tables / f'{cohort}_roc_scores.tsv', float_precision='round_trip')
    plot_expression_auc(scores, metrics, output)

    for name in ('marker_correlations', 'pathway_correlations'):
        rows = read(tables / f'{name}.tsv', float_precision='round_trip')
        cohorts = [rows[rows.Cohort.eq(cohort)] for cohort, _ in COHORTS]
        plot_correlation_dots(cohorts, name, output)
