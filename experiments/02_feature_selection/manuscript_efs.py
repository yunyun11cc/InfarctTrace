"""Compare manuscript feature selectors under one within-cohort CV protocol.

RA, SA and IA are local adaptations of the three aggregation concepts in
Ding et al., Pattern Recognition (2026), doi:10.1016/j.patcog.2025.112809.
"""

import argparse
from pathlib import Path

from joblib import Parallel, delayed, parallel_config
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

from resampling import CV_REPEATS, FOLDS, cv_splits
from select_features import STAGE_ONE, STAGE_TWO, UNION_OUTPUT, fit_rankers


SOURCE_DATA = STAGE_TWO / 'manuscript/source_data'
SELECTORS = {
    '均值差': 'VST mean difference', 'LASSO': 'LASSO',
    'ElasticNet': 'Elastic Net', '随机森林': 'Random forest',
}
METHODS = (*SELECTORS.values(), 'EFS-RA', 'EFS-SA', 'EFS-IA')
FEATURE_COUNTS = range(5, 31)


def score_and_order(expression, days):
    """Rank four base scores and three aggregates, using gene ID as the final tie break."""
    raw_scores = fit_rankers(expression, days)
    scores = np.abs(np.column_stack([raw_scores[method] for method in SELECTORS]).astype(float))
    if not np.isfinite(scores).all():
        raise ValueError('Feature scores must be finite')
    ids = expression.columns.to_numpy()
    ranks = rankdata(-scores, axis=0, method='average')
    mean_rank = ranks.mean(axis=1)
    span = scores.max(axis=0) - scores.min(axis=0)
    normalized_scores = np.divide(
        scores - scores.min(axis=0), span, out=np.zeros_like(scores), where=span > 0,
    )
    mean_score = normalized_scores.mean(axis=1)
    # A prefix that retains a tied block in full must reach the block's last rank.
    common_prefix = rankdata(-scores, axis=0, method='max').max(axis=1)
    order = {name: np.lexsort((ids, -scores[:, column])) for column, name in enumerate(SELECTORS.values())}
    order['EFS-RA'] = np.lexsort((ids, mean_rank))
    order['EFS-SA'] = np.lexsort((ids, -mean_score))
    order['EFS-IA'] = np.lexsort((ids, mean_rank, common_prefix))

    evidence = pd.DataFrame({
        'gene_id': ids, 'mean_base_rank': mean_rank,
        'mean_minmax_score': mean_score, 'common_prefix_size': common_prefix,
    })
    return order, evidence


def score_fold(number, expression, days, train, validation):
    """Fit all selectors and the common L2 model using the training libraries only."""
    train_expression = expression.iloc[train]
    train_days = days.iloc[train]
    train_is_d7 = train_days.eq(7).to_numpy()
    valid_is_d7 = days.iloc[validation].eq(7).to_numpy()
    order, _ = score_and_order(train_expression, train_days)
    scaler = StandardScaler().fit(train_expression)
    train_values = scaler.transform(train_expression)
    valid_values = scaler.transform(expression.iloc[validation])
    rows = []
    for method in METHODS:
        for k in FEATURE_COUNTS:
            keep = order[method][:k]
            fit = LogisticRegression(C=1, solver='lbfgs', l1_ratio=0, max_iter=1000,
                                     class_weight='balanced').fit(train_values[:, keep], train_is_d7)
            predicted = fit.predict_proba(valid_values[:, keep])[:, 1]
            rows.append((number, method, k, len(train), len(validation), int(valid_is_d7.sum()),
                         roc_auc_score(valid_is_d7, predicted)))
    return rows


def run_analysis(output, repeats, jobs):
    candidates = pd.read_csv(UNION_OUTPUT / 'tables/feature_ranking.tsv', sep='\t', index_col='基因ID').sort_index()
    metadata = pd.read_csv(STAGE_ONE / 'tables/sample_inclusion.csv', index_col='GSM编号')
    days = metadata.loc[metadata['主分析纳入'].eq('是'), '时间点（天）'].sort_index()
    vst = pd.read_csv(STAGE_ONE / 'matrices/timecourse_vst.tsv.gz', sep='\t', index_col='基因ID')
    expression = vst.loc[candidates.index, days.index].T
    is_d7 = days.eq(7).to_numpy()
    if expression.shape != (22, 3907) or is_d7.sum() != 4:
        raise ValueError('Expected 22 libraries, four D7 positives, and 3,907 union genes')

    output.mkdir(parents=True, exist_ok=True)
    order, evidence = score_and_order(expression, days)
    evidence = evidence.iloc[order['EFS-RA']].copy()
    evidence.insert(1, 'efs_ra_rank', np.arange(1, len(evidence) + 1))
    evidence.insert(1, 'gene_symbol', candidates.loc[evidence['gene_id'], '基因符号'].to_numpy())
    evidence.to_csv(output / 'figure_2c_efs_rank.tsv', sep='\t', index=False)

    splits = list(cv_splits(is_d7, repeats))
    tasks = (delayed(score_fold)(number, expression, days, train, valid)
             for number, (train, valid) in enumerate(splits, 1))
    rows = []
    with parallel_config(backend='loky', inner_max_num_threads=1):
        for done, fold in enumerate(Parallel(n_jobs=jobs, return_as='generator_unordered')(tasks), 1):
            rows.extend(fold)
            if done % 25 == 0 or done == len(splits):
                print(f'Aggregation CV: {done}/{len(splits)} folds', flush=True)
    detail = pd.DataFrame(rows, columns=['split', 'method', 'k', 'n_train', 'n_validation',
                                         'n_d7_validation', 'auroc']).sort_values(['split', 'method', 'k'])
    detail['repeat'] = (detail['split'] - 1) // FOLDS + 1
    detail['fold'] = (detail['split'] - 1) % FOLDS + 1
    summary = detail.groupby(['method', 'k'])['auroc'].agg(
        auroc_mean='mean', auroc_sd='std', n_folds='size',
    ).reset_index()
    detail.to_csv(output / 'figure_2b_cv_folds.tsv', sep='\t', index=False)
    summary.to_csv(output / 'figure_2b_cv_summary.tsv', sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=SOURCE_DATA)
    parser.add_argument('--cv-repeats', type=int, default=CV_REPEATS)
    parser.add_argument('--jobs', type=int, default=6)
    args = parser.parse_args()
    if min(args.cv_repeats, args.jobs) < 1:
        parser.error('repeat and job counts must be positive')
    run_analysis(args.output_dir, args.cv_repeats, args.jobs)
