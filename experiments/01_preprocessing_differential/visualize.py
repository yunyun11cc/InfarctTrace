"""绘制阶段 01 的中文质控、差异与富集图。"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from paths import ROOT, FIGURES, MATRICES, TABLES

TERM_ZH = {
    'GO:0001501': '骨骼系统发育', 'GO:0007094': '有丝分裂纺锤体组装检查点信号',
    'GO:0031577': '纺锤体检查点信号', 'GO:0033046': '姐妹染色单体分离的负调控',
    'GO:0033048': '有丝分裂姐妹染色单体分离的负调控', 'GO:0045841': '有丝分裂中期至后期转变的负调控',
    'mmu04512': '细胞外基质与受体相互作用', 'mmu05322': '系统性红斑狼疮',
    'mmu04060': '细胞因子与受体相互作用', 'mmu04514': '细胞黏附分子相互作用',
    'mmu04974': '蛋白质消化与吸收', 'mmu05150': '金黄色葡萄球菌感染',
    'mmu00190': '氧化磷酸化', 'mmu05014': '肌萎缩侧索硬化', 'mmu05208': '化学致癌中的活性氧',
    'mmu04714': '产热', 'mmu05012': '帕金森病', 'mmu05020': '朊病毒病',
}


def configure_plot_font():
    """在常见系统中选择已安装的中文字体。"""
    available = {font.name for font in font_manager.fontManager.ttflist}
    for name in ('Microsoft YaHei', 'Noto Sans CJK SC', 'Source Han Sans SC',
                 'PingFang SC', 'SimHei', 'WenQuanYi Zen Hei'):
        if name in available:
            plt.rcParams.update({'font.family': name, 'font.size': 10, 'axes.unicode_minus': False})
            return
    raise RuntimeError('未找到中文字体。请安装 Noto Sans CJK SC 后重新运行。')


def submitted_library_pca(meta, tpm):
    # 图 1 使用全部提交文库，标出需要结合作者记录核对的样本。
    if not tpm.columns.equals(meta.index):
        raise ValueError('提交文库的 TPM 列和样本表未对齐。')
    logtpm = np.log2(tpm + 1)
    variable = logtpm.var(axis=1).nlargest(1000).index
    pca = PCA(n_components=2)
    coordinates = pca.fit_transform(logtpm.loc[variable].T)
    fig, ax = plt.subplots(figsize=(8.8, 5.4), layout='constrained')
    for day in sorted(meta.day.unique()):
        mask = meta.day.to_numpy() == day
        ax.scatter(coordinates[mask, 0], coordinates[mask, 1], label=f'D{day}', s=42)
    for i, (_, row) in enumerate(meta.iterrows()):
        if (row.day, row.replicate) in {(3, 2), (3, 3), (28, 4)}:
            ax.annotate(f'D{row.day}r{row.replicate}', coordinates[i], fontsize=8,
                        xytext=(3, 5), textcoords='offset points')
    ax.margins(x=0.12)
    ax.set(xlabel=f'第一主成分（{pca.explained_variance_ratio_[0]:.1%}）',
           ylabel=f'第二主成分（{pca.explained_variance_ratio_[1]:.1%}）',
           title='图 1  成纤维细胞提交文库的主成分分析')
    ax.legend(ncol=2, fontsize=9)
    fig.supxlabel('每个点代表一个提交文库。基于 log2（TPM + 1）方差最大的 1,000 个基因进行 PCA。', fontsize=9)
    fig.savefig(FIGURES / '01_submitted_library_pca.png', dpi=180)
    plt.close(fig)


def read_result(path):
    return pd.read_csv(path, sep='\t', index_col=0)


def timecourse_qc(meta, counts, normalized, vst):
    # 两个箱线图使用相同的基因和文库，便于比较标准化前后的分布。
    days = meta['时间点（天）'].to_numpy()
    labels = [f'D{day}r{rep}' for day, rep in zip(days, meta['组内重复编号'])]
    colors = dict(zip(sorted(set(days)), plt.get_cmap('tab10').colors))
    fig, axes = plt.subplots(2, 1, figsize=(11.8, 7), layout='constrained')
    for ax, matrix, title in zip(axes, (counts, normalized),
                                 ('A  长度缩放计数', 'B  文库大小归一化计数')):
        b = ax.boxplot(np.log2(matrix.to_numpy() + 1), tick_labels=labels, showfliers=False,
                       patch_artist=True, medianprops={'color': '#263238'})
        for box, day in zip(b['boxes'], days):
            box.set_facecolor(colors[day])
            box.set_alpha(0.45)
        ax.set_title(title, loc='left')
        ax.set_ylabel('log2（计数 + 1）')
        ax.tick_params(axis='x', labelrotation=45)
    fig.suptitle('图 2  主分析文库标准化前后的计数分布', fontsize=14)
    fig.supxlabel('每个箱体对应一个主分析文库。两图采用相同基因，大小因子由 PyDESeq2 估计。')
    fig.savefig(FIGURES / '02_timecourse_normalization.png', dpi=180)
    plt.close(fig)

    variable = vst.var(axis=1).nlargest(500).index
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(vst.loc[variable].T)
    fig, ax = plt.subplots(figsize=(8.2, 5.3), layout='constrained')
    for day, color in colors.items():
        mask = days == day
        ax.scatter(pcs[mask, 0], pcs[mask, 1], s=52, label=f'D{day}', color=color)
    ax.set(xlabel=f'第一主成分（{pca.explained_variance_ratio_[0]:.1%}）',
           ylabel=f'第二主成分（{pca.explained_variance_ratio_[1]:.1%}）',
           title='图 3  标准化后成纤维细胞文库的主成分分析')
    ax.legend(ncol=3)
    fig.supxlabel('22 个主分析文库。按 PyDESeq2 VST 方差选取前 500 个基因，颜色表示采样组。')
    fig.savefig(FIGURES / '03_timecourse_pca.png', dpi=180)
    plt.close(fig)


def timecourse_heatmap(meta, vst, lrt):
    # 仅展示总体检验排名靠前的基因，按行标准化便于读时间模式。
    top = lrt[lrt['总体时间效应显著'].eq('是')].nsmallest(25, '校正后P值')
    matrix = vst.loc[top.index]
    z = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1), axis=0)
    days = meta['时间点（天）'].to_numpy()
    labels = [f'D{day}r{rep}' for day, rep in zip(days, meta['组内重复编号'])]
    genes = top['基因符号'].fillna(top.index.to_series())
    group = pd.factorize(days, sort=True)[0]
    fig, axes = plt.subplots(2, 1, figsize=(11, 9.6), height_ratios=(0.18, 8), layout='constrained')
    axes[0].imshow(group[None, :], aspect='auto', cmap=ListedColormap(plt.get_cmap('tab10').colors[:6]))
    axes[0].axis('off')
    image = axes[1].imshow(z, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5)
    axes[1].set(xticks=range(len(labels)), xticklabels=labels,
                yticks=range(len(genes)), yticklabels=genes, ylabel='基因符号')
    axes[1].tick_params(axis='x', labelrotation=55)
    fig.colorbar(image, ax=axes, shrink=0.84, label='基因行内标准分数')
    fig.suptitle('图 4  总体时间效应排名前 25 个基因的表达模式', fontsize=14)
    fig.supxlabel('按 Python 总体 LRT 校正后 P 值排序。颜色比较同一基因在不同文库中的 VST 表达。')
    fig.savefig(FIGURES / '04_timecourse_heatmap.png', dpi=180)
    plt.close(fig)


def timecourse_volcano(wald):
    # 火山图直接以校正后 P 值作纵轴，水平线对应预设统计门槛。
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), layout='constrained')
    colors = {'上调': '#bc4539', '下调': '#286b9e', '未达到门槛': '#b6bec6'}
    for ax, (contrast, frame) in zip(axes.flat, wald.groupby('比较', sort=False)):
        frame = frame.dropna(subset=['校正后P值', 'log2表达倍数变化'])
        x, q = frame['log2表达倍数变化'], frame['校正后P值']
        selected = frame['通过差异筛选'].eq('是')
        states = np.where(selected, np.where(x > 0, '上调', '下调'), '未达到门槛')
        y = -np.log10(np.maximum(q, 1e-300))
        for state in ('未达到门槛', '下调', '上调'):
            mask = states == state
            ax.scatter(x[mask], y[mask], s=4.5, alpha=0.52, color=colors[state], rasterized=True)
        ax.axhline(-np.log10(0.05), color='#717171', linestyle='--', linewidth=0.7)
        for cutoff in (-1, 1):
            ax.axvline(cutoff, color='#717171', linestyle='--', linewidth=0.7)
        ax.set(title=f'{contrast.split("_")[0]} 对 D0', xlabel='log2 表达倍数变化', ylabel='−log10（校正后 P 值）')
    axes.flat[-1].axis('off')
    axes.flat[-1].legend(handles=[Line2D([], [], marker='o', linestyle='', color=color, label=label)
                                 for label, color in colors.items()], loc='center', frameon=False)
    fig.suptitle('图 5  心梗后各时间点相对未损伤对照的差异表达', fontsize=14)
    fig.supxlabel('每点为一个 PyDESeq2 Wald 检验基因。筛选门槛为 |log2FC| ≥ 1 且校正后 P 值 < 0.05。')
    fig.savefig(FIGURES / '05_timecourse_volcano.png', dpi=180)
    plt.close(fig)


def enrichment_overview():
    # 代表图展示 D7，所有时间点及上调、下调方向的结果均保存在表中。
    go, kegg, gsea = [pd.read_csv(TABLES / name, sep='\t') for name in (
        'enrichment_go_ora.tsv.gz', 'enrichment_kegg_ora.tsv.gz', 'enrichment_kegg_gsea.tsv.gz')]
    go = go[(go['注释集'] == 'BP') & (go['方向'] == '上调')]
    kegg = kegg[kegg['方向'] == '上调']
    panels = [(go, 'A  GO 生物过程（上调基因）', '富集倍数'),
              (kegg, 'B  KEGG 通路（上调基因）', '富集倍数'),
              (gsea, 'C  KEGG 基因集富集分析', '标准化富集分数NES')]
    fig, axes = plt.subplots(3, 1, figsize=(10.7, 10.2), layout='constrained')
    for ax, (frame, title, value) in zip(axes, panels):
        significant = frame[(frame['比较'] == 'D7_vs_D0') & (frame['校正后P值'] < 0.05)]
        if value == '标准化富集分数NES':
            # GSEA 按方向分别取前三项，避免单一方向占满展示位置。
            top = pd.concat([significant[significant[value] > 0].nsmallest(3, '校正后P值'),
                             significant[significant[value] < 0].nsmallest(3, '校正后P值')]).iloc[::-1]
            ax.axvline(0, color='#687076', linewidth=0.8)
        else:
            top = significant.nsmallest(5, '校正后P值').iloc[::-1]
        ax.set_title(title, loc='left')
        if top.empty:
            ax.text(0.5, 0.5, '该图所选范围内无达到显著性门槛的条目', ha='center', va='center', transform=ax.transAxes)
            ax.set_axis_off()
            continue
        labels = [f'{TERM_ZH.get(term, name[:58])}（{term}）' for term, name in zip(top['条目ID'], top['条目原名'])]
        ax.barh(range(len(top)), top[value], color=np.where(top[value] > 0, '#b24b3f', '#32789f'))
        ax.set(yticks=range(len(top)), yticklabels=labels, xlabel=value)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
    fig.suptitle('图 6  心梗后第 7 天的功能与通路富集', fontsize=14)
    fig.supxlabel('基于 Python 差异结果。校正后 P 值 < 0.05 的条目中，A、B 各取前 5 项，C 正负方向各取前 3 项。')
    fig.savefig(FIGURES / '06_enrichment_D7.png', dpi=180)
    plt.close(fig)


def main():
    configure_plot_font()
    FIGURES.mkdir(parents=True, exist_ok=True)
    submitted = pd.read_csv(ROOT / 'datasets/GSE252027/metadata/samples_FB.csv').set_index('gsm')
    submitted_library_pca(submitted, read_result(MATRICES / 'FB_TPM.tsv.gz'))
    counts, normalized, vst = [read_result(MATRICES / name) for name in (
        'timecourse_length_scaled_counts.tsv.gz', 'timecourse_normalized_counts.tsv.gz', 'timecourse_vst.tsv.gz')]
    meta = pd.read_csv(TABLES / 'sample_inclusion.csv').set_index('GSM编号').loc[counts.columns]
    lrt = read_result(TABLES / 'timecourse_lrt.tsv.gz')
    wald = pd.read_csv(TABLES / 'timecourse_wald.tsv.gz', sep='\t')
    timecourse_qc(meta, counts, normalized, vst)
    timecourse_heatmap(meta, vst, lrt)
    timecourse_volcano(wald)
    enrichment_overview()
    print('已生成全时段中文质控、差异和富集图。')


if __name__ == '__main__':
    main()


