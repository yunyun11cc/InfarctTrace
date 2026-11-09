"""用 Python 完成 GSE252027 的全时段差异分析。"""
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

DAYS = [0, 1, 3, 7, 14, 28]
WALD_LABELS = {'baseMean': '平均标准化计数', 'log2FoldChange': 'log2表达倍数变化',
               'lfcSE': '效应标准误', 'stat': 'Wald统计量', 'pvalue': 'P值', 'padj': '校正后P值'}


def fit_time_model(quant, samples, inference):
    # 保留在至少一个时间组内、有两份文库计数达到 10 的基因。
    counts = length_scaled_counts(*(matrix.loc[:, samples.index] for matrix in quant))
    keep = np.logical_or.reduce([(counts.loc[:, samples.index[samples.day.eq(day)]] >= 10).sum(axis=1) >= 2
                                 for day in DAYS])
    design = pd.DataFrame({'day': pd.Categorical([f'D{d}' for d in samples.day],
                                                categories=[f'D{d}' for d in DAYS])}, index=samples.index)
    if design['day'].value_counts().min() < 2:
        raise ValueError('每个时间组至少需要两个纳入文库。')
    dds = DeseqDataSet(counts=counts.loc[keep].T, metadata=design, design='~day',
                      refit_cooks=False, inference=inference, quiet=False)
    dds.deseq2()
    return dds


def profile_loglik(counts, sizes, dispersion):
    """在固定离散度和文库大小因子时，求一组文库的最大负二项对数似然。"""
    # 组内期望是 size_factor × group_mean。得分方程单调，解位于归一化计数的范围内。
    scaled = counts / sizes[:, None]
    lower, upper = scaled.min(axis=0), scaled.max(axis=0)
    for _ in range(48):
        middle = (lower + upper) / 2
        mu = sizes[:, None] * middle
        positive = ((counts - mu) / (1 + dispersion * mu)).sum(axis=0) > 0
        lower, upper = np.where(positive, middle, lower), np.where(positive, upper, middle)
    mu = sizes[:, None] * ((lower + upper) / 2)
    # 与模型无关的常数在似然比中抵消。xlogy 让全零组按极限正确贡献 0。
    return (xlogy(counts, mu) - (counts + 1 / dispersion) * np.log1p(dispersion * mu)).sum(axis=0)


def overall_lrt(dds, samples, annotation):
    # 全模型每组一个均值，简化模型共用一个均值。两者使用同一份离散度估计。
    counts = dds.X.astype(float)
    sizes = dds.obs['size_factors'].to_numpy()
    dispersion = dds.var['dispersions'].to_numpy()
    full = sum(profile_loglik(counts[mask], sizes[mask], dispersion)
               for mask in (samples.day.eq(day).to_numpy() for day in DAYS))
    reduced = profile_loglik(counts, sizes, dispersion)
    statistic = 2 * (full - reduced)
    if np.any(statistic < -1e-5) or not np.isfinite(statistic).all():
        raise ValueError('总体时间检验得到异常似然差。')
    statistic = np.maximum(statistic, 0)
    p = chi2.sf(statistic, df=len(DAYS) - 1)
    p[np.asarray(dds.cooks_outlier(), dtype=bool)] = np.nan
    tested = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    q[tested] = false_discovery_control(p[tested], method='bh')
    return pd.DataFrame({'基因符号': annotation['基因符号'].reindex(dds.var_names),
                         '平均标准化计数': dds.layers['normed_counts'].mean(axis=0),
                         'LRT统计量': statistic, 'P值': p, '校正后P值': q,
                         '总体时间效应显著': np.where(q < 0.05, '是', '否')}, index=dds.var_names)


def wald_contrast(dds, day, annotation, inference):
    # 同一个六组模型提供各时点对 D0 的方向和效应，保留 PyDESeq2 的独立过滤。
    stats = DeseqStats(dds, contrast=['day', f'D{day}', 'D0'], alpha=0.05,
                       inference=inference, quiet=True)
    stats.summary()
    table = stats.results_df.rename(columns=WALD_LABELS).copy()
    table.insert(0, '比较', f'D{day}_vs_D0')
    table.insert(0, '基因符号', annotation['基因符号'].reindex(table.index))
    table['通过差异筛选'] = np.where(table['校正后P值'].lt(0.05) & table['log2表达倍数变化'].abs().ge(1), '是', '否')
    return table


def sensitivity_tables(primary_lrt, primary_d28, other_lrt, other_d28):
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
        a = set(primary.index[primary[key].eq('是')])
        b = set(other.index[other[key].eq('是')])
        rows.append({'检验': name, '主分析入选基因数': len(a), '敏感性分析入选基因数': len(b),
                     '共同入选基因数': len(a & b), 'Jaccard指数': len(a & b) / len(a | b) if a | b else np.nan})
    pd.DataFrame(rows).to_csv(TABLES / 'timecourse_sensitivity_summary.tsv', sep='\t', index=False, encoding='utf-8-sig')


def main():
    for directory in (MATRICES, TABLES, REPORTS):
        directory.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(TABLES / 'sample_inclusion.csv').set_index('GSM编号').rename(columns={'时间点（天）': 'day'})
    quant = [pd.read_csv(MATRICES / f'FB_{field}.tsv.gz', sep='\t', index_col=0)
             for field in ('TPM', 'EffectiveLength', 'NumReads')]
    if any(not matrix.columns.equals(meta.index) or not matrix.index.equals(quant[0].index) for matrix in quant):
        raise ValueError('样本或基因索引不一致。')
    annotation = pd.read_csv(MATRICES / 'gene_annotation.tsv.gz', sep='\t', index_col=0)
    inference = DefaultInference(n_cpus=2)
    samples = meta[meta['主分析纳入'].eq('是')]
    print(f'拟合 {len(samples)} 个文库的 Python 主模型。', flush=True)
    dds = fit_time_model(quant, samples, inference)
    lrt = overall_lrt(dds, samples, annotation)
    comparisons = [wald_contrast(dds, day, annotation, inference) for day in DAYS[1:]]
    save_gene_table(lrt.sort_values('校正后P值'), TABLES / 'timecourse_lrt.tsv.gz')
    save_gene_table(pd.concat(comparisons), TABLES / 'timecourse_wald.tsv.gz')
    rows = []
    for day, table in zip(DAYS[1:], comparisons):
        passed = table['通过差异筛选'].eq('是')
        rows.append({'比较': f'D{day}_vs_D0', '文库数': int(samples.day.eq(day).sum()),
                     '有限校正P值基因数': int(table['校正后P值'].notna().sum()),
                     '上调基因数': int((passed & table['log2表达倍数变化'].gt(0)).sum()),
                     '下调基因数': int((passed & table['log2表达倍数变化'].lt(0)).sum())})
    summary = pd.DataFrame(rows)
    summary.to_csv(TABLES / 'timecourse_contrast_summary.tsv', sep='\t', index=False, encoding='utf-8-sig')
    dds.vst(use_design=True)
    for name, values in [('length_scaled_counts', dds.X), ('normalized_counts', dds.layers['normed_counts']),
                         ('vst', dds.layers['vst_counts'])]:
        matrix = pd.DataFrame(values.T, index=dds.var_names, columns=dds.obs_names)
        save_gene_table(matrix, MATRICES / f'timecourse_{name}.tsv.gz')
    pd.DataFrame({'GSM编号': dds.obs_names, '时间组': dds.obs['day'].to_numpy(),
                  '文库大小因子': dds.obs['size_factors'].to_numpy()}).to_csv(
        TABLES / 'timecourse_size_factors.tsv', sep='\t', index=False, encoding='utf-8-sig')

    sensitivity = meta[meta['敏感性分析纳入'].eq('是')]
    print(f'拟合 {len(sensitivity)} 个文库的 D28 敏感性模型。', flush=True)
    other = fit_time_model(quant, sensitivity, inference)
    sensitivity_tables(lrt, comparisons[-1], overall_lrt(other, sensitivity, annotation),
                       wald_contrast(other, 28, annotation, inference))
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
