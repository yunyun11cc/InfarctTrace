"""拟合发现队列时间模型并导出计数、大小因子及 VST。

读取步骤 01 的 GSM 纳入表、基因注释及 TPM、有效长度、估计计数矩阵。
输出主分析的 timecourse_* 计数和 VST 矩阵及文库大小因子。
"""

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference

from .config import DAYS, INFERENCE_CPUS, MATRICES, MIN_LIBRARIES, REPORTS, TABLES
from .step01_load_and_qc import expressed_genes, save_gene_table


def length_scaled_counts(tpm, effective_length, estimated_reads):
    """按 tximport 的 lengthScaledTPM 算法生成计数尺度输入。

    记基因为 i、文库为 j。先用本次纳入文库的平均有效长度 L̄_i 把 TPM 换回
    长度校正后的丰度 U_ij = TPM_ij × L̄_i，再按每列“估计片段总数 / U 的列和”
    缩放，使每列总量与保留基因行的原始 NumReads 总量一致，最后四舍五入。

    这样得到的整数计数既保留了 TPM 的相对关系，又具备负二项模型所需的
    文库总量尺度。
    """
    # 每个基因用所纳入文库的平均有效长度，再恢复每列原估计片段总量
    unscaled = tpm.mul(effective_length.mean(axis=1), axis=0)
    counts = unscaled.mul(estimated_reads.sum(axis=0) / unscaled.sum(axis=0), axis=1)

    # 取整前做一次守恒检查：各列缩放后的总量应等于原始估计片段总量
    if not np.allclose(counts.sum(), estimated_reads.sum(), rtol=1e-10):
        raise ValueError('Length scaling changed library totals')

    return counts.round().astype('int64')


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
    keep = expressed_genes(counts, samples)

    # 时间点作为分类因子，水平顺序固定为 D0 → D28，D0 是后续对比的参照
    design = pd.DataFrame(
        {'day': pd.Categorical([f'D{d}' for d in samples.day], categories=[f'D{d}' for d in DAYS])},
        index=samples.index,
    )

    # 每组至少两个文库才可估计组内离散度；D3 纳入 2 个文库，正好卡在下限
    if design['day'].value_counts().min() < MIN_LIBRARIES:
        raise ValueError('每个时间组至少需要两个纳入文库。')

    # refit_cooks=False：关闭 Cook 异常计数的自动替换（替换会改变计数尺度），
    # 异常影响改由总体检验里把对应基因的 P 值置空来表达
    # 定量矩阵从基因 × 文库转置为 PyDESeq2 所需的文库 × 基因。
    dds = DeseqDataSet(
        counts=counts.loc[keep].T,
        metadata=design,
        design='~day',
        refit_cooks=False,
        inference=inference,
        quiet=False,
    )
    dds.deseq2()

    return dds


def read_inputs():
    """读取并核对阶段 01 矩阵，返回元数据、定量矩阵、注释和估计器。"""
    for directory in (MATRICES, TABLES, REPORTS):
        directory.mkdir(parents=True, exist_ok=True)

    # 样本纳入表是 prepare.py 的产物，这里只读取“是否纳入”的结论
    meta = (
        pd.read_csv(TABLES / 'sample_inclusion.csv')
        .set_index('GSM编号')
        .rename(columns={'时间点（天）': 'day'})
    )

    # 三个定量矩阵的列顺序应与样本表完全一致
    quant = [
        pd.read_csv(MATRICES / f'FB_{field}.tsv.gz', sep='\t', index_col=0)
        for field in ('TPM', 'EffectiveLength', 'NumReads')
    ]
    if any(
        not matrix.columns.equals(meta.index) or not matrix.index.equals(quant[0].index)
        for matrix in quant
    ):
        raise ValueError('样本或基因索引不一致。')

    annotation = pd.read_csv(MATRICES / 'gene_annotation.tsv.gz', sep='\t', index_col=0)
    inference = DefaultInference(n_cpus=INFERENCE_CPUS)
    return meta, quant, annotation, inference


def export_normalization(dds):
    """对已拟合的主模型执行 VST，将三个矩阵和大小因子写入磁盘。"""
    # VST 让高表达基因的方差趋于稳定，用于 PCA 与热图；计数矩阵仍按原始尺度保存
    dds.vst(use_design=True)
    for name, values in [
        ('length_scaled_counts', dds.X),
        ('normalized_counts', dds.layers['normed_counts']),
        ('vst', dds.layers['vst_counts']),
    ]:
        matrix = pd.DataFrame(values.T, index=dds.var_names, columns=dds.obs_names)
        save_gene_table(matrix, MATRICES / f'timecourse_{name}.tsv.gz')

    # 大小因子反映各文库总体测序量的相对高低
    pd.DataFrame(
        {
            'GSM编号': dds.obs_names,
            '时间组': dds.obs['day'].to_numpy(),
            '文库大小因子': dds.obs['size_factors'].to_numpy(),
        }
    ).to_csv(TABLES / 'timecourse_size_factors.tsv', sep='\t', index=False, encoding='utf-8-sig')
