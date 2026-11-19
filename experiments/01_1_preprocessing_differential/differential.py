"""用 Python 完成 GSE252027 的全时段差异分析。

流程衔接 ``prepare.py``：后者产出三个定量矩阵、基因注释和样本纳入表，本模块
在 22 个（或 21 个）文库上拟合负二项模型，回答两个不同的问题：

1. **总体时间效应**：六个时间组的表达均值是否不全相等（LRT，5 个自由度）。
2. **逐时间点比较**：每个损伤后时间点相对 D0 的方向和幅度（Wald 检验）。

主分析使用 22 个文库；另用排除 D28r4 的 21 个文库复算总体效应与 D28 对比，
量化单个样本取舍对结论的影响。
"""

import json

import numpy as np
import pandas as pd
import pydeseq2
import scipy
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from scipy.special import xlogy
from scipy.stats import chi2, false_discovery_control

from paths import MATRICES, REPORTS, TABLES
from prepare import length_scaled_counts, save_gene_table


# 六个采样时间点，顺序同时决定分类因子的水平顺序和 Wald 比较的先后
DAYS = [0, 1, 3, 7, 14, 28]

# PyDESeq2 结果列到中文表头的映射；五个比较共用同一套命名
WALD_LABELS = {'baseMean': '平均标准化计数', 'log2FoldChange': 'log2表达倍数变化',
               'lfcSE': '效应标准误', 'stat': 'Wald统计量', 'pvalue': 'P值', 'padj': '校正后P值'}


# ---------------------------------------------------------------------------
# 1. 拟合全时段模型
# ---------------------------------------------------------------------------

def fit_time_model(quant, samples, inference):
    """在同一批计数、大小因子和离散度上拟合 ``~day`` 负二项模型。

    参数
    ----
    quant : list[DataFrame]
        ``[TPM, EffectiveLength, NumReads]`` 三个矩阵，列名是 GSM。
    samples : DataFrame
        本次纳入的样本子表，索引为 GSM，含 ``day``。
    inference : DefaultInference
        PyDESeq2 的并行估计器。

    返回
    ----
    DeseqDataSet
        已运行 ``deseq2()``、含拟合结果的对象；主分析与敏感性分析各调用一次。
    """
    # 长度缩放只使用本次纳入的文库，平均有效长度与列总量都随纳入集合变化
    counts = length_scaled_counts(*(matrix.loc[:, samples.index] for matrix in quant))

    # 低表达过滤：至少存在一个时间组，组内有 2 个文库的计数达到 10。
    # 这样只在单个时间窗口表达的基因也有机会进入模型。
    keep = np.logical_or.reduce([(counts.loc[:, samples.index[samples.day.eq(day)]] >= 10).sum(axis=1) >= 2
                                 for day in DAYS])

    # 时间点作为分类因子，水平顺序固定为 D0 → D28，D0 是后续对比的参照
    design = pd.DataFrame({'day': pd.Categorical([f'D{d}' for d in samples.day],
                                                categories=[f'D{d}' for d in DAYS])}, index=samples.index)

    # 每组至少两个文库才可估计组内离散度；D3 纳入 2 个文库，正好卡在下限
    if design['day'].value_counts().min() < 2:
        raise ValueError('每个时间组至少需要两个纳入文库。')

    # refit_cooks=False：关闭 Cook 异常计数的自动替换（替换会改变计数尺度），
    # 异常影响改由总体检验里把对应基因的 P 值置空来表达
    dds = DeseqDataSet(counts=counts.loc[keep].T, metadata=design, design='~day',
                       refit_cooks=False, inference=inference, quiet=False)
    dds.deseq2()

    return dds


# ---------------------------------------------------------------------------
# 2. 总体时间效应：固定离散度的负二项似然比检验
# ---------------------------------------------------------------------------

def profile_loglik(counts, sizes, dispersion):
    """在固定离散度和文库大小因子时，求一组文库的最大负二项对数似然。

    记文库 j 的期望为 μ_j = s_j·q，q 为该组共同均值。固定离散度 α 后，对数
    似然中随 q 变化的部分为

        ℓ(q) = Σ_j [ C_j·ln μ_j − (C_j + 1/α)·ln(1 + α·μ_j) ]。

    由于得分方程单调，代码在“归一化计数”的取值区间内用二分法求 q 的极大值
    点，因此不需要外部优化器。返回每个基因各自的最大对数似然。
    """
    # 组内期望是 size_factor × group_mean。得分方程单调，解位于归一化计数的范围内。
    scaled = counts / sizes[:, None]
    lower, upper = scaled.min(axis=0), scaled.max(axis=0)

    # 48 次二分后区间宽度约为初始范围的 2^-48，对似然比精度足够
    for _ in range(48):
        middle = (lower + upper) / 2
        mu = sizes[:, None] * middle

        # 得分函数符号 > 0 表示最优点在当前中点右侧
        positive = ((counts - mu) / (1 + dispersion * mu)).sum(axis=0) > 0
        lower, upper = np.where(positive, middle, lower), np.where(positive, upper, middle)

    # 用区间中点作为组均值，回代得到该组最大对数似然
    mu = sizes[:, None] * ((lower + upper) / 2)

    # 与模型无关的常数在似然比中抵消。xlogy 让全零组按极限正确贡献 0。
    return (xlogy(counts, mu) - (counts + 1 / dispersion) * np.log1p(dispersion * mu)).sum(axis=0)


def overall_lrt(dds, samples, annotation):
    """比较“六组各自均值”与“六组共用均值”的负二项似然比检验。

    全模型每组一个均值，简化模型共用一个均值。两者使用同一份由 PyDESeq2 估计
    的离散度与文库大小因子，因此统计量只反映时间组带来的拟合改善：

        D = 2(ℓ_full − ℓ_reduced)，p = P(χ²_5 ≥ D)。

    Cook 异常标记的基因不参与 BH 家族，其校正后 P 值留空。
    """
    # 全模型每组一个均值，简化模型共用一个均值。两者使用同一份离散度估计。
    counts = dds.X.astype(float)
    sizes = dds.obs['size_factors'].to_numpy()
    dispersion = dds.var['dispersions'].to_numpy()

    # 六组分别求最大似然，求和即完整模型的最大对数似然
    full = sum(profile_loglik(counts[mask], sizes[mask], dispersion)
               for mask in (samples.day.eq(day).to_numpy() for day in DAYS))

    # 简化模型忽略分组，所有文库共用一个均值
    reduced = profile_loglik(counts, sizes, dispersion)

    statistic = 2 * (full - reduced)

    # 数值误差可能让本应为 0 的差变成极小负数，超过容差说明二分或输入有问题
    if np.any(statistic < -1e-5) or not np.isfinite(statistic).all():
        raise ValueError('总体时间检验得到异常似然差。')

    statistic = np.maximum(statistic, 0)

    # 自由度 = 时间组数 − 1 = 5
    p = chi2.sf(statistic, df=len(DAYS) - 1)

    # Cook 异常基因在总体检验中不参与多重校正
    p[np.asarray(dds.cooks_outlier(), dtype=bool)] = np.nan

    # 只对有限 P 值做 BH；其余保持 NaN，避免把“未检验”当成“不显著”
    tested = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    q[tested] = false_discovery_control(p[tested], method='bh')

    return pd.DataFrame({'基因符号': annotation['基因符号'].reindex(dds.var_names),
                         '平均标准化计数': dds.layers['normed_counts'].mean(axis=0),
                         'LRT统计量': statistic, 'P值': p, '校正后P值': q,
                         '总体时间效应显著': np.where(q < 0.05, '是', '否')}, index=dds.var_names)


# ---------------------------------------------------------------------------
# 3. 逐时间点比较：Wald 检验
# ---------------------------------------------------------------------------

def wald_contrast(dds, day, annotation, inference):
    """用同一个六组模型做 ``D{day} vs D0`` 的 Wald 对比。

    与 LRT 不同，Wald 给出**有方向**的效应量：log2FC 为正表示相对 D0 上调。
    PyDESeq2 的独立过滤会按 baseMean 分布调整可检验基因集，因此同一基因在
    LRT 与 Wald 表中的校正后 P 值可能一个有限、一个为空。
    """
    # 同一个六组模型提供各时点对 D0 的方向和效应，保留 PyDESeq2 的独立过滤。
    stats = DeseqStats(dds, contrast=['day', f'D{day}', 'D0'], alpha=0.05,
                       inference=inference, quiet=True)
    stats.summary()

    table = stats.results_df.rename(columns=WALD_LABELS).copy()

    # 两个标识列放在最左侧：比较名与基因符号
    table.insert(0, '比较', f'D{day}_vs_D0')
    table.insert(0, '基因符号', annotation['基因符号'].reindex(table.index))

    # 预设门槛：校正后 P 值 < 0.05 且 |log2FC| ≥ 1
    table['通过差异筛选'] = np.where(table['校正后P值'].lt(0.05) & table['log2表达倍数变化'].abs().ge(1), '是', '否')

    return table


# ---------------------------------------------------------------------------
# 4. D28r4 敏感性对照
# ---------------------------------------------------------------------------

def sensitivity_tables(primary_lrt, primary_d28, other_lrt, other_d28):
    """写出主分析与敏感性分析的两类对照表。

    1. 逐基因表：同一基因在两套模型中的总体校正后 P 值、D28 效应与校正后 P 值。
    2. 汇总表：两份显著名单的入选数、交集数与 Jaccard 指数。
    """
    # 两次拟合的基因集可能因低表达过滤不同而略有差异，只对照共同基因
    shared = primary_d28.index.intersection(other_d28.index)

    joined = pd.DataFrame({'主分析总体校正后P值': primary_lrt.loc[shared, '校正后P值'],
                           '敏感性总体校正后P值': other_lrt.loc[shared, '校正后P值'],
                           '主分析D28效应': primary_d28.loc[shared, 'log2表达倍数变化'],
                           '敏感性D28效应': other_d28.loc[shared, 'log2表达倍数变化'],
                           '主分析D28校正后P值': primary_d28.loc[shared, '校正后P值'],
                           '敏感性D28校正后P值': other_d28.loc[shared, '校正后P值']})
    save_gene_table(joined, TABLES / 'timecourse_d28_sensitivity.tsv.gz')

    rows = []
    for name, primary, other, key in [
        ('总体时间效应', primary_lrt, other_lrt, '总体时间效应显著'),
        ('D28_vs_D0', primary_d28, other_d28, '通过差异筛选'),
    ]:
        # 两份“显著”名单的交集与并集之比即 Jaccard 指数
        a = set(primary.index[primary[key].eq('是')])
        b = set(other.index[other[key].eq('是')])
        rows.append({'检验': name, '主分析入选基因数': len(a), '敏感性分析入选基因数': len(b),
                     '共同入选基因数': len(a & b), 'Jaccard指数': len(a & b) / len(a | b) if a | b else np.nan})

    pd.DataFrame(rows).to_csv(TABLES / 'timecourse_sensitivity_summary.tsv', sep='\t', index=False, encoding='utf-8-sig')


# ---------------------------------------------------------------------------
# 5. 入口
# ---------------------------------------------------------------------------

def main():
    for directory in (MATRICES, TABLES, REPORTS):
        directory.mkdir(parents=True, exist_ok=True)

    # 样本纳入表是 prepare.py 的产物，这里只读取“是否纳入”的结论
    meta = pd.read_csv(TABLES / 'sample_inclusion.csv').set_index('GSM编号').rename(columns={'时间点（天）': 'day'})

    # 三个定量矩阵的列顺序应与样本表完全一致
    quant = [pd.read_csv(MATRICES / f'FB_{field}.tsv.gz', sep='\t', index_col=0)
             for field in ('TPM', 'EffectiveLength', 'NumReads')]
    if any(not matrix.columns.equals(meta.index) or not matrix.index.equals(quant[0].index) for matrix in quant):
        raise ValueError('样本或基因索引不一致。')

    annotation = pd.read_csv(MATRICES / 'gene_annotation.tsv.gz', sep='\t', index_col=0)
    inference = DefaultInference(n_cpus=2)

    # 主分析：22 个文库（D3 只保留 r1、r4）
    samples = meta[meta['主分析纳入'].eq('是')]
    print(f'拟合 {len(samples)} 个文库的 Python 主模型。', flush=True)

    dds = fit_time_model(quant, samples, inference)
    lrt = overall_lrt(dds, samples, annotation)

    # 五个对比共用同一个 dds，避免重复拟合
    comparisons = [wald_contrast(dds, day, annotation, inference) for day in DAYS[1:]]

    # 总体表按校正后 P 值升序，最显著的基因在最上面；Wald 表把五个对比纵向堆叠
    save_gene_table(lrt.sort_values('校正后P值'), TABLES / 'timecourse_lrt.tsv.gz')
    save_gene_table(pd.concat(comparisons), TABLES / 'timecourse_wald.tsv.gz')

    # 每个对比的上调/下调计数，便于快速核对表 2 的结论
    rows = []
    for day, table in zip(DAYS[1:], comparisons):
        passed = table['通过差异筛选'].eq('是')
        rows.append({'比较': f'D{day}_vs_D0', '文库数': int(samples.day.eq(day).sum()),
                     '有限校正P值基因数': int(table['校正后P值'].notna().sum()),
                     '上调基因数': int((passed & table['log2表达倍数变化'].gt(0)).sum()),
                     '下调基因数': int((passed & table['log2表达倍数变化'].lt(0)).sum())})

    summary = pd.DataFrame(rows)
    summary.to_csv(TABLES / 'timecourse_contrast_summary.tsv', sep='\t', index=False, encoding='utf-8-sig')

    # VST 让高表达基因的方差趋于稳定，用于 PCA 与热图；计数矩阵仍按原始尺度保存
    dds.vst(use_design=True)
    for name, values in [('length_scaled_counts', dds.X), ('normalized_counts', dds.layers['normed_counts']),
                         ('vst', dds.layers['vst_counts'])]:
        matrix = pd.DataFrame(values.T, index=dds.var_names, columns=dds.obs_names)
        save_gene_table(matrix, MATRICES / f'timecourse_{name}.tsv.gz')

    # 大小因子反映各文库总体测序量的相对高低
    pd.DataFrame({'GSM编号': dds.obs_names, '时间组': dds.obs['day'].to_numpy(),
                  '文库大小因子': dds.obs['size_factors'].to_numpy()}).to_csv(
        TABLES / 'timecourse_size_factors.tsv', sep='\t', index=False, encoding='utf-8-sig')

    # 敏感性分析：排除 D28r4，只用 21 个文库重跑同一套模型
    sensitivity = meta[meta['敏感性分析纳入'].eq('是')]
    print(f'拟合 {len(sensitivity)} 个文库的 D28 敏感性模型。', flush=True)

    other = fit_time_model(quant, sensitivity, inference)
    sensitivity_tables(lrt, comparisons[-1], overall_lrt(other, sensitivity, annotation),
                       wald_contrast(other, 28, annotation, inference))

    # 运行摘要记录版本与口径，方便复现和跨版本对照
    report = {'pydeseq2': pydeseq2.__version__, 'scipy': scipy.__version__,
              '纳入文库': len(samples), '敏感性纳入文库': len(sensitivity), '进入模型基因数': dds.n_vars,
              '总体显著基因数': int(lrt['总体时间效应显著'].eq('是').sum()),
              '总体检验': '固定离散度负二项 LRT（六组与仅截距对比，5 自由度），Cook 过滤后 BH 校正',
              '逐组检验': 'PyDESeq2 Wald，含独立过滤'}
    (REPORTS / 'differential_python_session.json').write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

    print(summary.to_string(index=False))
    print(f"总体时间效应显著基因数 {report['总体显著基因数']:,}")


if __name__ == '__main__':
    main()
