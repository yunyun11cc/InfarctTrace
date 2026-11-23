"""对共享的发现队列拟合模型进行总体 LRT、逐时间点 Wald 与 D28 敏感性分析。

读取步骤 01 的样本纳入决定和步骤 02 的拟合对象。
输出主检验统计表、D28r4 排除敏感性表和运行报告。
"""

import json

import numpy as np
import pandas as pd
import pydeseq2
import scipy
from pydeseq2.ds import DeseqStats
from scipy.special import xlogy
from scipy.stats import chi2, false_discovery_control

from .config import ALPHA, DAYS, LRT_BISECTION_STEPS, MIN_ABS_LOG2FC, REPORTS, TABLES
from .step01_load_and_qc import save_gene_table
from .step02_normalization_vst import fit_time_model


WALD_LABELS = {
    'baseMean': '平均标准化计数',
    'log2FoldChange': 'log2表达倍数变化',
    'lfcSE': '效应标准误',
    'stat': 'Wald统计量',
    'pvalue': 'P值',
    'padj': '校正后P值',
}


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
    for _ in range(LRT_BISECTION_STEPS):
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
    full = sum(
        profile_loglik(counts[mask], sizes[mask], dispersion)
        for mask in (samples.day.eq(day).to_numpy() for day in DAYS)
    )

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

    return pd.DataFrame(
        {
            '基因符号': annotation['基因符号'].reindex(dds.var_names),
            '平均标准化计数': dds.layers['normed_counts'].mean(axis=0),
            'LRT统计量': statistic,
            'P值': p,
            '校正后P值': q,
            '总体时间效应显著': np.where(q < ALPHA, '是', '否'),
        },
        index=dds.var_names,
    )


def wald_contrast(dds, day, annotation, inference):
    """用同一个六组模型做 ``D{day} vs D0`` 的 Wald 对比。

    与 LRT 不同，Wald 给出**有方向**的效应量：log2FC 为正表示相对 D0 上调。
    PyDESeq2 的独立过滤会按 baseMean 分布调整可检验基因集，因此同一基因在
    LRT 与 Wald 表中的校正后 P 值可能一个有限、一个为空。
    """
    # 同一个六组模型提供各时点对 D0 的方向和效应，保留 PyDESeq2 的独立过滤。
    stats = DeseqStats(
        dds, contrast=['day', f'D{day}', 'D0'], alpha=ALPHA, inference=inference, quiet=True
    )
    stats.summary()

    table = stats.results_df.rename(columns=WALD_LABELS).copy()

    # 两个标识列放在最左侧：比较名与基因符号
    table.insert(0, '比较', f'D{day}_vs_D0')
    table.insert(0, '基因符号', annotation['基因符号'].reindex(table.index))

    # 预设门槛：校正后 P 值 < 0.05 且 |log2FC| ≥ 1
    significant = table['校正后P值'].lt(ALPHA)
    large_effect = table['log2表达倍数变化'].abs().ge(MIN_ABS_LOG2FC)
    table['通过差异筛选'] = np.where(significant & large_effect, '是', '否')

    return table


def sensitivity_tables(primary_lrt, primary_d28, other_lrt, other_d28):
    """写出主分析与敏感性分析的两类对照表。

    1. 逐基因表：同一基因在两套模型中的总体校正后 P 值、D28 效应与校正后 P 值。
    2. 汇总表：两份显著名单的入选数、交集数与 Jaccard 指数。
    """
    # 两次拟合的基因集可能因低表达过滤不同而略有差异，只对照共同基因
    shared = primary_d28.index.intersection(other_d28.index)

    joined = pd.DataFrame(
        {
            '主分析总体校正后P值': primary_lrt.loc[shared, '校正后P值'],
            '敏感性总体校正后P值': other_lrt.loc[shared, '校正后P值'],
            '主分析D28效应': primary_d28.loc[shared, 'log2表达倍数变化'],
            '敏感性D28效应': other_d28.loc[shared, 'log2表达倍数变化'],
            '主分析D28校正后P值': primary_d28.loc[shared, '校正后P值'],
            '敏感性D28校正后P值': other_d28.loc[shared, '校正后P值'],
        }
    )
    save_gene_table(joined, TABLES / 'timecourse_d28_sensitivity.tsv.gz')

    rows = []
    for name, primary, other, key in [
        ('总体时间效应', primary_lrt, other_lrt, '总体时间效应显著'),
        ('D28_vs_D0', primary_d28, other_d28, '通过差异筛选'),
    ]:
        # 两份“显著”名单的交集与并集之比即 Jaccard 指数
        a = set(primary.index[primary[key].eq('是')])
        b = set(other.index[other[key].eq('是')])
        rows.append(
            {
                '检验': name,
                '主分析入选基因数': len(a),
                '敏感性分析入选基因数': len(b),
                '共同入选基因数': len(a & b),
                'Jaccard指数': len(a & b) / len(a | b) if a | b else np.nan,
            }
        )

    pd.DataFrame(rows).to_csv(
        TABLES / 'timecourse_sensitivity_summary.tsv', sep='\t', index=False, encoding='utf-8-sig'
    )


def write_primary(dds, samples, annotation, inference):
    """写出主模型结果，并返回 LRT、逐时点对比和汇总表供敏感性分析使用。"""
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
        rows.append(
            {
                '比较': f'D{day}_vs_D0',
                '文库数': int(samples.day.eq(day).sum()),
                '有限校正P值基因数': int(table['校正后P值'].notna().sum()),
                '上调基因数': int((passed & table['log2表达倍数变化'].gt(0)).sum()),
                '下调基因数': int((passed & table['log2表达倍数变化'].lt(0)).sum()),
            }
        )

    summary = pd.DataFrame(rows)
    summary.to_csv(
        TABLES / 'timecourse_contrast_summary.tsv', sep='\t', index=False, encoding='utf-8-sig'
    )
    return lrt, comparisons, summary


def write_sensitivity(dds, samples, sensitivity, quant, annotation, inference, primary):
    """重拟合 21 文库模型，对照主分析并写出敏感性结果和软件版本。"""
    lrt, comparisons, summary = primary
    print(f'拟合 {len(sensitivity)} 个文库的 D28 敏感性模型。', flush=True)
    other = fit_time_model(quant, sensitivity, inference)
    sensitivity_tables(
        lrt,
        comparisons[-1],
        overall_lrt(other, sensitivity, annotation),
        wald_contrast(other, 28, annotation, inference),
    )

    # 运行摘要记录版本与口径，方便复现和跨版本对照
    report = {
        'pydeseq2': pydeseq2.__version__,
        'scipy': scipy.__version__,
        '纳入文库': len(samples),
        '敏感性纳入文库': len(sensitivity),
        '进入模型基因数': dds.n_vars,
        '总体显著基因数': int(lrt['总体时间效应显著'].eq('是').sum()),
        '总体检验': '固定离散度负二项 LRT（六组与仅截距对比，5 自由度），Cook 过滤后 BH 校正',
        '逐组检验': 'PyDESeq2 Wald，含独立过滤',
    }
    (REPORTS / 'differential_python_session.json').write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8'
    )

    print(summary.to_string(index=False))
    print(f"总体时间效应显著基因数 {report['总体显著基因数']:,}")
