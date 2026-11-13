"""Plot the Stage 02 intersection, candidate expression, and selection stability."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
from matplotlib import font_manager
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from select_features import OUTPUT, STAGE_ONE


BLUE = '#235F84'
ORANGE = '#B4552A'
GRAY = '#62717D'
LIGHT = '#DFE5EA'
CHINESE_FONTS = ('Noto Sans CJK SC', 'PingFang SC', 'Microsoft YaHei', 'Source Han Sans SC', 'SimHei')


def plot_intersection(ranking, wald):
    passed = wald[wald['通过差异筛选'].eq('是')]['比较'].value_counts()
    overlap = len(ranking)
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
    for key, label, color in directions:
        genes = ranking[ranking['相对D0方向'].eq(key)]
        reverse = key.startswith('D3')
        scatter.scatter(genes['D3对D0_log2FC'], genes['D7对D0_log2FC'], color=color,
                        s=32 if reverse else 10, marker='D' if reverse else 'o',
                        alpha=0.95 if reverse else 0.48, edgecolor='none',
                        label=f'{label}（{len(genes):,}）')
    scatter.axhline(0, color=LIGHT, zorder=0)
    scatter.axvline(0, color=LIGHT, zorder=0)
    scatter.set(xlabel='D3 对 D0 的 log₂FC', ylabel='D7 对 D0 的 log₂FC', title='B  交集内的效应方向')
    scatter.legend(loc='lower right', fontsize=8, frameon=True, facecolor='white', framealpha=0.93)
    fig.suptitle(f'D3 与 D7 相对 D0 的差异基因交集共有 {overlap:,} 个基因', fontweight='bold')
    return fig


def plot_expression(ranking, vst, libraries, repeats):
    genes = ranking.nsmallest(8, '排名')
    fig, panels = plt.subplots(2, 4, figsize=(11, 6), layout='constrained', sharey=True)
    days = ('3', '7')
    for ax, (gene_id, gene) in zip(panels.flat, genes.iterrows()):
        for x, day, color in zip(range(2), days, (GRAY, BLUE)):
            values = vst.loc[gene_id, libraries[day]].to_numpy()
            offsets = np.linspace(-0.13, 0.13, len(values))
            ax.scatter(x + offsets, values, color=color, s=27, zorder=3)
            ax.hlines(values.mean(), x - 0.24, x + 0.24, color='#182833', linewidth=1.5, zorder=4)
        effect = gene['D7减D3_VST均值差']
        count = round(gene['均值差入选频率'] * repeats)
        label = gene['基因符号'] if pd.notna(gene['基因符号']) else gene_id
        ax.set_title(f'{label}\nD7−D3 {effect:+.2f}，入选 {count}/{repeats}', fontsize=10)
        ax.set(xticks=[0, 1], xticklabels=['D3', 'D7'], xlim=(-0.45, 1.45))
        ax.grid(axis='y', color=LIGHT, linewidth=0.7)
        ax.tick_params(axis='x', length=0)
    fig.suptitle('主排序前八个基因的逐文库表达', fontweight='bold')
    fig.supylabel('VST 表达值')
    fig.supxlabel('每点为一份文库，短线为组均值。D3 有 2 份，D7 有 4 份。', fontsize=10)
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
    fonts = {font.name for font in font_manager.fontManager.ttflist}
    family = next((name for name in CHINESE_FONTS if name in fonts), None)
    if family is None:
        raise RuntimeError('未找到中文字体。请安装 Noto Sans CJK SC 后重新运行。')
    plt.rcParams.update({'font.family': [family, 'DejaVu Sans'], 'font.size': 10,
                         'axes.unicode_minus': False, 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.linewidth': 0.8})
    ranking = pd.read_csv(output / 'tables/feature_ranking.tsv', sep='\t', index_col='基因ID')
    summary = json.loads((output / 'reports/selection_summary.json').read_text(encoding='utf-8'))
    wald = pd.read_csv(STAGE_ONE / 'tables/timecourse_wald.tsv.gz', sep='\t')
    vst = pd.read_csv(STAGE_ONE / 'matrices/timecourse_vst.tsv.gz', sep='\t', index_col='基因ID')
    figures = output / 'figures'
    figures.mkdir(parents=True, exist_ok=True)
    for name, fig in [('01_candidate_intersection', plot_intersection(ranking, wald)),
                      ('02_candidate_expression', plot_expression(ranking, vst, summary['文库'], summary['重抽样次数'])),
                      ('03_selection_stability', plot_stability(summary))]:
        fig.savefig(figures / f'{name}.png', dpi=180)
        plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    visualize(parser.parse_args().output_dir)
