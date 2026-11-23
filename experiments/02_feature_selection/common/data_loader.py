"""读取阶段一的 Wald、VST 和样本纳入表，以及并集排序表。

输出为按基因 ID 和 GSM 精确对齐的 pandas 数据，供各选择器及相关分析使用。
"""

from importlib import import_module

import pandas as pd

from ..config import STAGE_ONE


window = import_module('experiments.01_preprocessing_differential.step04_window_selection')
SETS = window.SETS


def read_wald():
    """读取带基因 ID 与通过标记的 Wald 长表。"""
    return pd.read_csv(STAGE_ONE / 'tables/timecourse_wald.tsv.gz', sep='\t', encoding='utf-8-sig')


def read_vst():
    """读取基因 × 文库表达矩阵，保留原行列顺序。"""
    path = STAGE_ONE / 'matrices/timecourse_vst.tsv.gz'
    return pd.read_csv(path, sep='\t', encoding='utf-8-sig', index_col='基因ID')


def read_main_days():
    """返回主分析文库的采样日，索引为原样本表中的 GSM。"""
    path = STAGE_ONE / 'tables/sample_inclusion.csv'
    metadata = pd.read_csv(path, encoding='utf-8-sig', index_col='GSM编号')
    return metadata.loc[metadata['主分析纳入'].eq('是'), '时间点（天）']


def read_ranking(output):
    """从给定筛选目录读取特征排序表并以基因 ID 索引。"""
    return pd.read_csv(output / 'tables/feature_ranking.tsv', sep='\t', index_col='基因ID')


def read_union_inputs(output):
    """返回按 ID 排序的并集、GSM 采样日与文库 × 基因 VST。"""
    candidates = read_ranking(output).sort_index()
    days = read_main_days().sort_index()
    expression = read_vst().loc[candidates.index, days.index].T
    return candidates, days, expression
