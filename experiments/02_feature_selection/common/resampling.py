"""执行并集版的原分类器重复 CV 和 Bootstrap。

从并集排序、VST 和纳入样本表读取固定输入。
写出逐折指标、入选频率、Bootstrap Jaccard 和重采样报告。
"""

import json

from joblib import Parallel, delayed, parallel_config
import numpy as np
import pandas as pd
import sklearn
from sklearn.metrics import accuracy_score, average_precision_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

from .cv_splitter import bootstrap_draws, cv_splits
from .data_loader import read_union_inputs
from .selection import CLASSIFIER_MODULES, top_membership
from ..config import CLASSIFIERS, D7_POSITIVES, FOLDS, METHODS, SEED, TOP_K, UNION_SHAPE
from ..methods import m01_vst_mean_diff


def resampled_fit(number, expression, labels, train, validation=None):
    """在训练下标上筛选；若给定验证下标，用原分类器在入选基因上重拟合并打分。"""
    vst, y_train = expression[train], labels[train]
    scaler = StandardScaler().fit(vst)
    scaled = scaler.transform(vst)
    selected = {}
    rows = []
    if validation is None:
        selected['均值差'] = np.flatnonzero(
            top_membership(np.abs(m01_vst_mean_diff.score(vst, y_train)), TOP_K)
        )
    for method in CLASSIFIER_MODULES:
        train_values = scaled if method.USES_SCALING else vst
        fit = method.fit(train_values, y_train)
        keep = top_membership(np.abs(method.scores(fit)), TOP_K)
        selected[method.METHOD] = np.flatnonzero(keep)
        if validation is not None:
            y_valid = labels[validation]
            valid_values = expression[validation]
            if method.USES_SCALING:
                valid_values = scaler.transform(valid_values)
            if keep.any():
                refit = method.fit(train_values[:, keep], y_train)
                probability = refit.predict_proba(valid_values[:, keep])[:, 1]
                prediction = refit.predict(valid_values[:, keep])
            else:
                probability = np.full(len(validation), y_train.mean())
                prediction = probability > 0.5
            rows.append(
                (
                    number,
                    method.METHOD,
                    len(train),
                    len(validation),
                    int(y_train.sum()),
                    int(y_valid.sum()),
                    int(keep.sum()),
                    roc_auc_score(y_valid, probability),
                    average_precision_score(y_valid, probability),
                    accuracy_score(y_valid, prediction),
                )
            )
    return number, selected, rows


def cv_results(parallel, expression, labels, candidates, repeats, tables):
    """执行给定重复次数的 CV，写出每折指标、聚合指标与入选频率。"""
    splits = list(cv_splits(labels, repeats))
    counts = {method: np.zeros(len(candidates), dtype=int) for method in CLASSIFIERS}
    rows = []
    jobs = (
        delayed(resampled_fit)(i, expression, labels, train, valid)
        for i, (train, valid) in enumerate(splits, 1)
    )
    for done, (_, selected, metrics) in enumerate(parallel(jobs), 1):
        for method, indices in selected.items():
            counts[method][indices] += 1
        rows.extend(metrics)
        if done % 25 == 0:
            print(f'CV: {done}/{len(splits)} folds', flush=True)
    detail = pd.DataFrame(
        rows,
        columns=[
            'Split',
            'Method',
            'N_train',
            'N_valid',
            'D7_train',
            'D7_valid',
            'Selected',
            'AUROC',
            'AUPRC_AP',
            'Accuracy',
        ],
    ).sort_values(['Split', 'Method'])
    detail['Repeat'] = (detail['Split'] - 1) // FOLDS + 1
    detail['Fold'] = (detail['Split'] - 1) % FOLDS + 1
    detail.to_csv(tables / 'cv_fold_metrics.tsv', sep='\t', index=False)
    summary = detail.groupby('Method', sort=False)[['AUROC', 'AUPRC_AP', 'Accuracy']].agg(
        ['mean', 'std']
    )
    summary.columns = ['_'.join(column) for column in summary.columns]
    summary.loc[list(CLASSIFIERS)].to_csv(tables / 'cv_summary.tsv', sep='\t')
    frequencies = candidates[['基因符号']].assign(
        **{name: count / len(splits) for name, count in counts.items()}
    )
    frequencies.to_csv(tables / 'cv_selection_frequency.tsv', sep='\t')


def bootstrap_results(parallel, expression, labels, candidates, repeats, tables):
    """执行配对标签分层 Bootstrap，写出名单重合度和入选概率。"""
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
            rows.append(
                (number, method, len(indices), unique, float(overlap / union) if union else np.nan)
            )
        if done % 50 == 0:
            print(f'Bootstrap: {done}/{repeats} draws', flush=True)
    detail = pd.DataFrame(rows, columns=['Draw', 'Method', 'Selected', 'Unique_libraries', 'Jaccard'])
    detail.sort_values(['Draw', 'Method']).to_csv(
        tables / 'bootstrap_jaccard.tsv', sep='\t', index=False
    )
    frequencies = candidates[['基因符号']].assign(
        **{name: count / repeats for name, count in counts.items()}
    )
    frequencies['Mean'] = frequencies[list(METHODS)].mean(axis=1)
    frequencies.sort_values(['Mean', '基因ID'], ascending=[False, True]).to_csv(
        tables / 'bootstrap_selection_probability.tsv',
        sep='\t',
    )


def run_analysis(selection_dir, output, cv_repeats, bootstrap_repeats, jobs):
    """加载固定并集并运行 CV 和 Bootstrap，将结果写入所选目录。"""
    candidates, days, matrix = read_union_inputs(selection_dir)
    labels = days.eq(7).to_numpy()
    # 保持与原重采样一致的 C 连续布局，行按 GSM、列按候选基因 ID。
    expression = np.ascontiguousarray(matrix.to_numpy())
    if expression.shape != UNION_SHAPE or labels.sum() != D7_POSITIVES:
        raise ValueError('需要 22 份文库、4 个 D7 正类与 3,907 个并集基因')
    tables = output / 'tables'
    reports = output / 'reports'
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    with parallel_config(backend='loky', inner_max_num_threads=1):
        with Parallel(n_jobs=jobs, return_as='generator_unordered') as parallel:
            cv_results(parallel, expression, labels, candidates, cv_repeats, tables)
            bootstrap_results(parallel, expression, labels, candidates, bootstrap_repeats, tables)
    summary = {
        'scikit-learn': sklearn.__version__,
        'GSM': days.index.tolist(),
        'D7': days.index[labels].tolist(),
        'CV_folds': FOLDS,
        'CV_repeats': cv_repeats,
        'CV_splits': FOLDS * cv_repeats,
        'bootstrap_repeats': bootstrap_repeats,
        'bootstrap_class_sizes': {'D7': 4, 'non_D7': 18},
        'seed': SEED,
        'top_k': TOP_K,
        'top_k_ties': 'keep boundary ties, exclude zero scores',
        'CV_prediction': 'refit same estimator on training-fold Top 20, then score validation fold',
        'AUPRC': 'average precision',
        'sd': 'sample SD across validation folds',
        'candidate_pool_and_VST': 'fixed full-cohort preprocessing, including validation libraries',
        'bootstrap_Jaccard_reference': 'union_22 full-library selected list, separately per method',
    }
    (reports / 'resampling_summary.json').write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + '\n'
    )
    print(f'Completed CV and bootstrap: {output}', flush=True)
