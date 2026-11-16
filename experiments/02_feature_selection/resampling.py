"""Repeated stratified CV and bootstrap on the Stage 02 union candidate pool."""

import argparse
import json
from pathlib import Path

from joblib import Parallel, delayed, parallel_config
import numpy as np
import pandas as pd
import sklearn
from sklearn.metrics import accuracy_score, average_precision_score, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler

from select_features import (CLASSIFIERS, SEED, STAGE_ONE, TOP_K, UNION_OUTPUT,
                             classifier_scores, fit_classifier, mean_difference, top_membership)


FOLDS = 4
CV_REPEATS = 50
BOOTSTRAPS = 1000
METHODS = ('均值差', *CLASSIFIERS)
OUTPUT = UNION_OUTPUT / 'resampling'


def cv_splits(labels, repeats):
    splitter = RepeatedStratifiedKFold(n_splits=FOLDS, n_repeats=repeats, random_state=SEED)
    return splitter.split(np.empty(len(labels)), labels)


def bootstrap_draws(labels, repeats):
    rng = np.random.default_rng(SEED)
    classes = [np.flatnonzero(labels == group) for group in (False, True)]
    draws = np.column_stack([rng.choice(indices, size=(repeats, len(indices)), replace=True)
                             for indices in classes])
    return np.sort(draws, axis=1)


def resampled_fit(number, expression, labels, train, validation=None):
    vst, y_train = expression[train], labels[train]
    scaler = StandardScaler().fit(vst)
    scaled = scaler.transform(vst)
    selected = {}
    rows = []
    if validation is None:
        selected['均值差'] = np.flatnonzero(top_membership(np.abs(mean_difference(vst, y_train)), TOP_K))
    for method in CLASSIFIERS:
        train_values = vst if method == '随机森林' else scaled
        fit = fit_classifier(method, train_values, y_train)
        keep = top_membership(np.abs(classifier_scores(method, fit)), TOP_K)
        selected[method] = np.flatnonzero(keep)
        if validation is not None:
            y_valid = labels[validation]
            valid_values = expression[validation]
            if method != '随机森林':
                valid_values = scaler.transform(valid_values)
            if keep.any():
                refit = fit_classifier(method, train_values[:, keep], y_train)
                probability = refit.predict_proba(valid_values[:, keep])[:, 1]
                prediction = refit.predict(valid_values[:, keep])
            else:
                probability = np.full(len(validation), y_train.mean())
                prediction = probability > 0.5
            rows.append((number, method, len(train), len(validation), int(y_train.sum()), int(y_valid.sum()),
                         int(keep.sum()), roc_auc_score(y_valid, probability),
                         average_precision_score(y_valid, probability), accuracy_score(y_valid, prediction)))
    return number, selected, rows


def cv_results(parallel, expression, labels, candidates, repeats, tables):
    splits = list(cv_splits(labels, repeats))
    counts = {method: np.zeros(len(candidates), dtype=int) for method in CLASSIFIERS}
    rows = []
    jobs = (delayed(resampled_fit)(i, expression, labels, train, valid)
            for i, (train, valid) in enumerate(splits, 1))
    for done, (_, selected, metrics) in enumerate(parallel(jobs), 1):
        for method, indices in selected.items():
            counts[method][indices] += 1
        rows.extend(metrics)
        if done % 25 == 0:
            print(f'CV: {done}/{len(splits)} folds', flush=True)
    detail = pd.DataFrame(rows, columns=['Split', 'Method', 'N_train', 'N_valid', 'D7_train', 'D7_valid',
                                         'Selected', 'AUROC', 'AUPRC_AP', 'Accuracy']).sort_values(['Split', 'Method'])
    detail['Repeat'] = (detail['Split'] - 1) // FOLDS + 1
    detail['Fold'] = (detail['Split'] - 1) % FOLDS + 1
    detail.to_csv(tables / 'cv_fold_metrics.tsv', sep='\t', index=False)
    summary = detail.groupby('Method', sort=False)[['AUROC', 'AUPRC_AP', 'Accuracy']].agg(['mean', 'std'])
    summary.columns = ['_'.join(column) for column in summary.columns]
    summary.loc[list(CLASSIFIERS)].to_csv(tables / 'cv_summary.tsv', sep='\t')
    frequencies = candidates[['基因符号']].assign(**{name: count / len(splits) for name, count in counts.items()})
    frequencies.to_csv(tables / 'cv_selection_frequency.tsv', sep='\t')


def bootstrap_results(parallel, expression, labels, candidates, repeats, tables):
    draws = bootstrap_draws(labels, repeats)
    counts = {method: np.zeros(len(candidates), dtype=int) for method in METHODS}
    reference = {method: candidates[f'{method}全样本入选'].to_numpy() for method in METHODS}
    rows = []
    jobs = (delayed(resampled_fit)(i, expression, labels, train) for i, train in enumerate(draws, 1))
    for done, (number, selected, _) in enumerate(parallel(jobs), 1):
        unique = len(np.unique(draws[number - 1]))
        for method, indices in selected.items():
            counts[method][indices] += 1
            overlap = reference[method][indices].sum()
            union = reference[method].sum() + len(indices) - overlap
            rows.append((number, method, len(indices), unique, float(overlap / union) if union else np.nan))
        if done % 50 == 0:
            print(f'Bootstrap: {done}/{repeats} draws', flush=True)
    detail = pd.DataFrame(rows, columns=['Draw', 'Method', 'Selected', 'Unique_libraries', 'Jaccard'])
    detail.sort_values(['Draw', 'Method']).to_csv(tables / 'bootstrap_jaccard.tsv', sep='\t', index=False)
    frequencies = candidates[['基因符号']].assign(**{name: count / repeats for name, count in counts.items()})
    frequencies['Mean'] = frequencies[list(METHODS)].mean(axis=1)
    frequencies.sort_values(['Mean', '基因ID'], ascending=[False, True]).to_csv(
        tables / 'bootstrap_selection_probability.tsv', sep='\t',
    )


def run_analysis(output, cv_repeats, bootstrap_repeats, jobs):
    candidates = pd.read_csv(UNION_OUTPUT / 'tables/feature_ranking.tsv', sep='\t', index_col='基因ID').sort_index()
    meta = pd.read_csv(STAGE_ONE / 'tables/sample_inclusion.csv', index_col='GSM编号')
    days = meta.loc[meta['主分析纳入'].eq('是'), '时间点（天）'].sort_index()
    labels = days.eq(7).to_numpy()
    vst = pd.read_csv(STAGE_ONE / 'matrices/timecourse_vst.tsv.gz', sep='\t', index_col='基因ID')
    expression = np.ascontiguousarray(vst.loc[candidates.index, days.index].to_numpy().T)
    if expression.shape != (22, 3907) or labels.sum() != 4:
        raise ValueError('需要 22 份文库、4 个 D7 正类与 3,907 个并集基因')
    tables = output / 'tables'
    reports = output / 'reports'
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    with parallel_config(backend='loky', inner_max_num_threads=1):
        with Parallel(n_jobs=jobs, return_as='generator_unordered') as parallel:
            cv_results(parallel, expression, labels, candidates, cv_repeats, tables)
            bootstrap_results(parallel, expression, labels, candidates, bootstrap_repeats, tables)
    summary = {'scikit-learn': sklearn.__version__, 'GSM': days.index.tolist(), 'D7': days.index[labels].tolist(),
               'CV_folds': FOLDS, 'CV_repeats': cv_repeats, 'CV_splits': FOLDS * cv_repeats,
               'bootstrap_repeats': bootstrap_repeats, 'bootstrap_class_sizes': {'D7': 4, 'non_D7': 18},
               'seed': SEED, 'top_k': TOP_K, 'top_k_ties': 'keep boundary ties, exclude zero scores',
               'CV_prediction': 'refit same estimator on training-fold Top 20, then score validation fold',
               'AUPRC': 'average precision', 'sd': 'sample SD across validation folds',
               'candidate_pool_and_VST': 'fixed full-cohort preprocessing, including validation libraries',
               'bootstrap_Jaccard_reference': 'union_22 full-library selected list, separately per method'}
    (reports / 'resampling_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n')
    print(f'Completed CV and bootstrap: {output}', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    parser.add_argument('--cv-repeats', type=int, default=CV_REPEATS)
    parser.add_argument('--bootstrap', type=int, default=BOOTSTRAPS)
    parser.add_argument('--jobs', type=int, default=8)
    args = parser.parse_args()
    if min(args.cv_repeats, args.bootstrap, args.jobs) < 1:
        parser.error('repeat counts and jobs must be positive')
    run_analysis(args.output_dir, args.cv_repeats, args.bootstrap, args.jobs)
