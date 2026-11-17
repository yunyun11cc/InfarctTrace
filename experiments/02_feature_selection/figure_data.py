"""Prepare overlap counts and gene ranks for the union figures."""

from collections import Counter
from itertools import product

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from select_features import STAGE_ONE, fit_rankers


SETS = ['D3 vs D0', 'D7 vs D0', 'Time LRT']
CORE_GENES = ('C1qtnf3', 'Cthrc1', 'Postn', 'Wisp1', 'Cemip', 'Cxadr', 'Kif26b', 'Col24a1', 'Mrc1')


def efs_ra_ranks(expression, days):
    """Average descending ranks from the four selectors; break ties by gene ID."""
    scores = np.abs(np.column_stack(list(fit_rankers(expression, days).values())))
    if not np.isfinite(scores).all():
        raise ValueError('Feature scores must be finite')
    ids = expression.columns.to_numpy()
    mean_ranks = rankdata(-scores, axis=0, method='average').mean(axis=1)
    order = np.lexsort((ids, mean_ranks))
    return pd.Series(np.arange(1, len(ids) + 1), index=ids[order], name='efs_ra_rank')


def screen_membership(gene_ids):
    wald = pd.read_csv(STAGE_ONE / 'tables/timecourse_wald.tsv.gz', sep='\t', dtype={'基因ID': str})
    lrt = pd.read_csv(STAGE_ONE / 'tables/timecourse_lrt.tsv.gz', sep='\t', index_col='基因ID')
    comparisons = [wald.loc[wald['比较'].eq(f'D{day}_vs_D0')].set_index('基因ID').reindex(gene_ids)
                   for day in (3, 7)]
    flags = [frame['校正后P值'].lt(0.05) & frame['log2表达倍数变化'].abs().ge(1)
             for frame in comparisons]
    membership = pd.DataFrame(dict(zip(SETS, [*flags, lrt['校正后P值'].reindex(gene_ids).lt(0.05)])))
    candidates = membership[SETS[0]] & membership[SETS[1]]
    membership['D3 and D7 candidate'] = candidates
    membership.insert(0, 'Gene_symbol', lrt['基因符号'].reindex(gene_ids))
    counts = Counter(membership[SETS].itertuples(index=False, name=None))
    combinations = list(product([True, False], repeat=3))
    overlap = pd.DataFrame(combinations, columns=SETS).assign(
        Genes=[counts[combination] for combination in combinations],
    )
    return membership.rename_axis('Gene_ID'), overlap
