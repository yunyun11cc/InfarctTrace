"""根据已有 Wald 和 LRT 表确定 D3/D7 候选窗口。

读取 timecourse_wald.tsv.gz 和 timecourse_lrt.tsv.gz。
写出交集、并集、筛选成员和交叠计数表，特征筛选也可从旧 Wald 表直接构建候选池。
"""

from collections import Counter
from itertools import product

import numpy as np
import pandas as pd

from .config import ALPHA, MIN_ABS_LOG2FC, TABLES


SETS = ['D3 vs D0', 'D7 vs D0', 'Time LRT']


def select_candidates(wald, use_union=False):
    """按 Wald 达标标记构建交集或并集，返回按基因 ID 排序的效应和标记表。"""
    wald = wald.set_index('基因ID')
    evidence = ['log2表达倍数变化', '校正后P值', '通过差异筛选']
    d3, d7 = (
        wald.loc[wald['比较'].eq(f'D{day}_vs_D0'), ['基因符号', *evidence]].rename(
            columns={
                'log2表达倍数变化': f'D{day}对D0_log2FC',
                '校正后P值': f'D{day}对D0_padj',
                '通过差异筛选': f'D{day}通过差异筛选',
            }
        )
        for day in (3, 7)
    )
    merged = d3.join(d7.rename(columns={'基因符号': '基因符号_D7'}), how='outer', validate='one_to_one')
    pass3 = merged['D3通过差异筛选'].eq('是')
    pass7 = merged['D7通过差异筛选'].eq('是')
    keep = pass3 | pass7 if use_union else pass3 & pass7
    candidates = merged.loc[keep].sort_index()
    up3 = candidates['D3对D0_log2FC'].gt(0)
    up7 = candidates['D7对D0_log2FC'].gt(0)
    candidates['相对D0方向'] = np.select(
        [up3 & up7, ~up3 & ~up7, ~up3 & up7],
        ['均上调', '均下调', 'D3下调_D7上调'],
        default='D3上调_D7下调',
    )
    if use_union:
        candidates['基因符号'] = candidates['基因符号'].fillna(candidates['基因符号_D7'])
        for day, up, other_day in [(3, up3, 7), (7, up7, 3)]:
            column = f'D{day}通过差异筛选'
            only = candidates[column].eq('是') & ~candidates[f'D{other_day}通过差异筛选'].eq('是')
            candidates.loc[only, '相对D0方向'] = np.where(up[only], f'仅D{day}上调', f'仅D{day}下调')
        return candidates.drop(columns='基因符号_D7').assign(
            D3通过差异筛选=pass3.loc[candidates.index],
            D7通过差异筛选=pass7.loc[candidates.index],
        )
    return candidates.drop(columns=['基因符号_D7', 'D3通过差异筛选', 'D7通过差异筛选'])


def screen_membership(gene_ids):
    """按指定基因 ID 对齐 Wald 与 LRT，返回逐基因成员掩码和各交集计数。"""
    wald = pd.read_csv(TABLES / 'timecourse_wald.tsv.gz', sep='\t', dtype={'基因ID': str})
    lrt = pd.read_csv(TABLES / 'timecourse_lrt.tsv.gz', sep='\t', index_col='基因ID')
    comparisons = [
        wald.loc[wald['比较'].eq(f'D{day}_vs_D0')].set_index('基因ID').reindex(gene_ids)
        for day in (3, 7)
    ]
    flags = [
        frame['校正后P值'].lt(ALPHA) & frame['log2表达倍数变化'].abs().ge(MIN_ABS_LOG2FC)
        for frame in comparisons
    ]
    membership = pd.DataFrame(dict(zip(SETS, [*flags, lrt['校正后P值'].reindex(gene_ids).lt(ALPHA)])))
    candidates = membership[SETS[0]] & membership[SETS[1]]
    membership['D3 and D7 candidate'] = candidates
    membership.insert(0, 'Gene_symbol', lrt['基因符号'].reindex(gene_ids))
    counts = Counter(membership[SETS].itertuples(index=False, name=None))
    combinations = list(product([True, False], repeat=3))
    overlap = pd.DataFrame(combinations, columns=SETS).assign(
        Genes=[counts[combination] for combination in combinations],
    )
    return membership.rename_axis('Gene_ID'), overlap


def run():
    """写出两个候选池和绘图用的 Wald/LRT 交叠统计。"""
    wald = pd.read_csv(TABLES / 'timecourse_wald.tsv.gz', sep='\t', encoding='utf-8-sig')
    for name, use_union in (('intersection', False), ('union', True)):
        selection = select_candidates(wald, use_union=use_union)
        selection.to_csv(TABLES / f'd3_d7_{name}.tsv', sep='\t', encoding='utf-8-sig')

    ids = pd.Index(wald['基因ID'].unique()).sort_values()
    membership, overlap = screen_membership(ids)
    membership.to_csv(TABLES / 'd3_d7_screen_membership.tsv', sep='\t')
    overlap.to_csv(TABLES / 'd3_d7_screen_overlap.tsv', sep='\t', index=False)
