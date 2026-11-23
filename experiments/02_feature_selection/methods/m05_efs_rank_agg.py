"""EFS-RA、EFS-SA、EFS-IA 及其与四种基础选择器的比较。

输入为阶段一 VST、并集排序、Bootstrap 结果和九候选通路数据。
输出为 EFS 排名、逐折和汇总 AUROC、分析图及冻结候选表。
"""

import json

from joblib import Parallel, delayed, parallel_config
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

from ..common.cv_splitter import cv_splits
from ..common.data_loader import read_union_inputs
from ..common.plot_style import (
    PANEL_BLUE as BLUE,
    INK,
    configure_panel_style,
    save_panel_figure,
)
from ..common.selection import fit_rankers
from ..config import (
    CORE_GENES,
    D7_POSITIVES,
    EFS_FEATURE_COUNTS,
    EFS_MAX_ITER,
    FOLDS,
    LOGISTIC_C,
    FROZEN_CANDIDATES_FILE,
    STAGE_ONE,
    UNION_SHAPE,
)


SELECTORS = {
    '均值差': 'VST mean difference',
    'LASSO': 'LASSO',
    'ElasticNet': 'Elastic Net',
    '随机森林': 'Random forest',
}
METHODS = (*SELECTORS.values(), 'EFS-RA', 'EFS-SA', 'EFS-IA')
SELECTOR_STYLES = {
    'VST mean difference': ('均值差', BLUE, 'o'),
    'LASSO': ('LASSO', '#6E5D85', 's'),
    'Elastic Net': ('ElasticNet', '#39776E', '^'),
    'Random forest': ('随机森林', '#687582', 'D'),
}
ENSEMBLES = {'EFS-RA': ('#8F3D2A', 'o'), 'EFS-SA': ('#C67D39', 's'), 'EFS-IA': ('#432B60', '^')}


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
        scores - scores.min(axis=0),
        span,
        out=np.zeros_like(scores),
        where=span > 0,
    )
    mean_score = normalized_scores.mean(axis=1)
    # A prefix that retains a tied block in full must reach the block's last rank.
    common_prefix = rankdata(-scores, axis=0, method='max').max(axis=1)
    order = {
        name: np.lexsort((ids, -scores[:, column])) for column, name in enumerate(SELECTORS.values())
    }
    order['EFS-RA'] = np.lexsort((ids, mean_rank))
    order['EFS-SA'] = np.lexsort((ids, -mean_score))
    order['EFS-IA'] = np.lexsort((ids, mean_rank, common_prefix))

    evidence = pd.DataFrame(
        {
            'gene_id': ids,
            'mean_base_rank': mean_rank,
            'mean_minmax_score': mean_score,
            'common_prefix_size': common_prefix,
        }
    )
    return order, evidence


def efs_ra_ranks(expression, days):
    """Average descending ranks from the four selectors; break ties by gene ID."""
    scores = np.abs(np.column_stack(list(fit_rankers(expression, days).values())))
    if not np.isfinite(scores).all():
        raise ValueError('Feature scores must be finite')
    ids = expression.columns.to_numpy()
    mean_ranks = rankdata(-scores, axis=0, method='average').mean(axis=1)
    order = np.lexsort((ids, mean_ranks))
    return pd.Series(np.arange(1, len(ids) + 1), index=ids[order], name='efs_ra_rank')


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
        for k in EFS_FEATURE_COUNTS:
            keep = order[method][:k]
            fit = LogisticRegression(
                C=LOGISTIC_C, solver='lbfgs', l1_ratio=0, max_iter=EFS_MAX_ITER, class_weight='balanced'
            ).fit(train_values[:, keep], train_is_d7)
            predicted = fit.predict_proba(valid_values[:, keep])[:, 1]
            rows.append(
                (
                    number,
                    method,
                    k,
                    len(train),
                    len(validation),
                    int(valid_is_d7.sum()),
                    roc_auc_score(valid_is_d7, predicted),
                )
            )
    return rows


def run_analysis(selection_dir, output, repeats, jobs):
    """在固定并集上运行全队列聚合排序和重复 CV，并导出逐折与汇总矩阵。"""
    candidates, days, expression = read_union_inputs(selection_dir)
    is_d7 = days.eq(7).to_numpy()
    if expression.shape != UNION_SHAPE or is_d7.sum() != D7_POSITIVES:
        raise ValueError('Expected 22 libraries, four D7 positives, and 3,907 union genes')

    output.mkdir(parents=True, exist_ok=True)
    order, evidence = score_and_order(expression, days)
    evidence = evidence.iloc[order['EFS-RA']].copy()
    evidence.insert(1, 'efs_ra_rank', np.arange(1, len(evidence) + 1))
    evidence.insert(1, 'gene_symbol', candidates.loc[evidence['gene_id'], '基因符号'].to_numpy())
    evidence.to_csv(output / 'efs_ra_ranking.tsv', sep='\t', index=False)

    splits = list(cv_splits(is_d7, repeats))
    tasks = (
        delayed(score_fold)(number, expression, days, train, valid)
        for number, (train, valid) in enumerate(splits, 1)
    )
    rows = []
    with parallel_config(backend='loky', inner_max_num_threads=1):
        for done, fold in enumerate(Parallel(n_jobs=jobs, return_as='generator_unordered')(tasks), 1):
            rows.extend(fold)
            if done % 25 == 0 or done == len(splits):
                print(f'Aggregation CV: {done}/{len(splits)} folds', flush=True)
    detail = pd.DataFrame(
        rows, columns=['split', 'method', 'k', 'n_train', 'n_validation', 'n_d7_validation', 'auroc']
    ).sort_values(['split', 'method', 'k'])
    detail['repeat'] = (detail['split'] - 1) // FOLDS + 1
    detail['fold'] = (detail['split'] - 1) % FOLDS + 1
    summary = (
        detail.groupby(['method', 'k'])['auroc']
        .agg(
            auroc_mean='mean',
            auroc_sd='std',
            n_folds='size',
        )
        .reset_index()
    )
    detail.to_csv(output / 'efs_cv_folds.tsv', sep='\t', index=False)
    summary.to_csv(output / 'efs_cv_summary.tsv', sep='\t', index=False)


def write_candidate_frequencies(selection_dir, tables):
    """对齐九个冻结候选的 Bootstrap 入选率与并集 EFS-RA 名次。"""
    ranking = pd.read_csv(selection_dir / 'tables/feature_ranking.tsv', sep='\t')
    core = ranking.loc[ranking['基因符号'].isin(CORE_GENES)].set_index('基因符号')
    if not core.index.is_unique:
        raise ValueError('Candidate symbols must map uniquely to Ensembl IDs')
    core = core.loc[list(CORE_GENES)].rename(columns={'基因ID': 'gene_id'}).rename_axis('gene')
    bootstrap = pd.read_csv(
        selection_dir / 'resampling/tables/bootstrap_selection_probability.tsv',
        sep='\t',
        index_col='基因ID',
    )
    frequency_columns = [method[0] for method in SELECTOR_STYLES.values()] + ['Mean']
    core = core.join(bootstrap[frequency_columns], on='gene_id', validate='one_to_one')
    if not core[frequency_columns].apply(lambda series: series.between(0, 1).all()).all():
        raise ValueError('Missing or out-of-range bootstrap frequencies for the nine candidates')

    ranks = pd.read_csv(tables / 'efs_ra_ranking.tsv', sep='\t', index_col='gene_id')
    frequencies = core.join(ranks[['efs_ra_rank']], on='gene_id', validate='one_to_one')
    columns = ['gene_id', *frequency_columns, 'efs_ra_rank']
    frequencies[columns].to_csv(tables / 'candidate_frequencies.tsv', sep='\t')


def plot_cv(ax, source, selection_dir, folds):
    """绘制基础选择器和聚合方法的 AUROC，再标注原分类器重拟合参考值。"""
    summary = pd.read_csv(source / 'efs_cv_summary.tsv', sep='\t')
    methods = {name: (color, marker) for name, (_, color, marker) in SELECTOR_STYLES.items()} | ENSEMBLES
    measured = summary['auroc_mean'].to_numpy()
    if set(summary['method']) != set(methods) or not np.isfinite(measured).all():
        raise ValueError('Incomplete selector comparison')
    for method, (color, marker) in methods.items():
        values = summary.loc[summary['method'].eq(method)].sort_values('k')
        if values['k'].tolist() != list(EFS_FEATURE_COUNTS) or not values['n_folds'].eq(folds).all():
            raise ValueError(f'Expected {folds} folds per feature count for {method}')
        ensemble = method in ENSEMBLES
        ax.plot(
            values['k'],
            values['auroc_mean'],
            color=color,
            label=method,
            linestyle='-' if ensemble else '--',
            linewidth=1.20 if ensemble else 0.82,
            marker=marker,
            markevery=5,
            markersize=2.5,
            markeredgewidth=0.30,
        )
    legacy = pd.read_csv(selection_dir / 'resampling/tables/cv_summary.tsv', sep='\t')
    legacy_en = legacy.set_index('Method').at['ElasticNet', 'AUROC_mean']
    ax.scatter(
        [20],
        [legacy_en],
        marker='*',
        s=34,
        facecolor='white',
        edgecolor=INK,
        linewidth=0.8,
        zorder=6,
        label=f'Elastic Net refit ({legacy_en:.3f})',
    )
    low = max(0, min(np.floor(measured.min() * 20) / 20 - 0.015, legacy_en - 0.025))
    ax.set(
        xlim=(4.4, 30.6),
        ylim=(low, 1.005),
        xticks=range(5, 31, 5),
        xlabel='Features selected (k)',
        ylabel='Mean validation AUROC',
    )
    ax.legend(
        loc='lower left',
        bbox_to_anchor=(-0.025, 1.014),
        ncol=4,
        handlelength=1.55,
        columnspacing=0.80,
        handletextpad=0.46,
        borderaxespad=0,
        fontsize=5.7,
    )


def plot_candidate_selection(ax, frequencies):
    """绘制九候选入选频率和 EFS-RA 排名。"""
    plotted = frequencies.sort_values('efs_ra_rank')
    rows = np.arange(len(plotted))
    mrc1_row = plotted.index.get_loc('Mrc1')
    ax.axhspan(mrc1_row - 0.44, mrc1_row + 0.44, color='#F5E6DE', zorder=0, linewidth=0)
    ax.hlines(rows, 0, plotted['Mean'], color='#A6B0B6', linewidth=0.8, zorder=1)
    ax.scatter(plotted['Mean'], rows, color=INK, s=12, linewidths=0, zorder=2)
    offsets = (-0.24, -0.08, 0.08, 0.24)
    method_handles = []
    for offset, (name, (column, color, marker)) in zip(offsets, SELECTOR_STYLES.items(), strict=True):
        markers = ax.scatter(
            plotted[column],
            rows + offset,
            color=color,
            marker=marker,
            s=10,
            linewidths=0.28,
            edgecolors='white',
            zorder=3,
            clip_on=False,
            label=name,
        )
        method_handles.append(markers)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticks(rows, plotted.index, fontstyle='italic')
    ax.yaxis.set_tick_params(length=0, pad=4)
    ax.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    ax.set(xlim=(0, 1), ylim=(len(rows) - 0.48, -0.5), xlabel='Bootstrap selection frequency')
    ax.spines['left'].set_visible(False)
    for row, rank in zip(rows, plotted['efs_ra_rank']):
        ax.text(
            1.14,
            row,
            f'{rank:,}',
            transform=ax.get_yaxis_transform(),
            ha='right',
            va='center',
            fontsize=6.0,
            fontweight='bold' if row == mrc1_row else 'normal',
        )
    ax.text(1.14, 1.015, 'RA rank', transform=ax.transAxes, ha='right', va='bottom', fontsize=6.0)
    handle = Line2D(
        [],
        [],
        color='#A6B0B6',
        marker='o',
        markerfacecolor=INK,
        markeredgewidth=0,
        markersize=3.2,
        linewidth=0.8,
        label='Mean of four selectors',
    )
    ax.legend(
        handles=[handle, *method_handles],
        loc='lower left',
        bbox_to_anchor=(-0.04, 1.011),
        borderaxespad=0,
        ncol=3,
        handlelength=1.6,
        fontsize=6,
    )


def plot_results(selection_dir, output):
    """用已有统计表生成 EFS 交叉验证比较图和冻结候选入选率图。"""
    configure_panel_style()
    tables = output / 'tables'
    figures = output / 'figures'
    cv_table = tables / 'efs_cv_summary.tsv'
    if not cv_table.is_file():
        raise FileNotFoundError(f'{cv_table} is missing. Run Stage 02 with --steps efs.')
    report = selection_dir / 'resampling/reports/resampling_summary.json'
    folds = json.loads(report.read_text(encoding='utf-8'))['CV_splits']

    fig, ax = plt.subplots(figsize=(8.4, 4.8), layout='constrained')
    fig.suptitle('Cross-validation by feature selector', fontsize=12, fontweight='bold')
    plot_cv(ax, tables, selection_dir, folds)
    save_panel_figure(fig, figures, 'efs_cv_auroc')
    plt.close(fig)

    frequencies = pd.read_csv(tables / 'candidate_frequencies.tsv', sep='\t', index_col='gene')
    fig, ax = plt.subplots(figsize=(7.2, 4.8), layout='constrained')
    fig.suptitle('Bootstrap selection frequency and EFS-RA rank', fontsize=12, fontweight='bold')
    plot_candidate_selection(ax, frequencies)
    save_panel_figure(fig, figures, 'candidate_selection_frequency')
    plt.close(fig)


def write_table(selection_dir, aggregation_tables, biology_tables):
    """对齐入选频率、LRT 和通路表，写出外部验证使用的冻结九基因表。"""
    frequencies = pd.read_csv(
        aggregation_tables / 'candidate_frequencies.tsv',
        sep='\t',
        index_col='gene',
        dtype={'gene_id': str},
    ).loc[list(CORE_GENES)]
    ranking = pd.read_csv(selection_dir / 'tables/feature_ranking.tsv', sep='\t', index_col='基因ID')
    lrt = pd.read_csv(STAGE_ONE / 'tables/timecourse_lrt.tsv.gz', sep='\t', index_col='基因ID')
    pathways = pd.read_csv(biology_tables / 'candidate_pathway_correlations.tsv', sep='\t')
    if len(pathways) != 63 or pathways.groupby('gene_id')['pathway'].nunique().ne(7).any():
        raise ValueError('Expected seven pathways per candidate for all nine candidates')
    associations = pathways.loc[pathways['rho'].gt(0) & pathways['q'].lt(0.05)]
    groups = associations.groupby('gene_id')['group'].agg(frozenset)
    pathway_labels = {
        frozenset(): '未达阈值',
        frozenset(['Inflammation']): '炎症',
        frozenset(['Fibrosis']): '纤维化',
        frozenset(['Inflammation', 'Fibrosis']): '双通路',
    }
    ids = frequencies['gene_id']
    result = pd.DataFrame(
        {
            '基因符号': frequencies.index,
            '基因ID': ids.to_numpy(),
            'D7 - 非D7 均值差': ranking.loc[ids, 'D7减非D7_VST均值差'].to_numpy(),
            'LRT 校正后 P 值': lrt.loc[ids, '校正后P值'].to_numpy(),
            'VST 均值差 Bootstrap 选中率': frequencies['均值差'].to_numpy(),
            'LASSO Bootstrap 选中率': frequencies['LASSO'].to_numpy(),
            'EN Bootstrap 选中率': frequencies['ElasticNet'].to_numpy(),
            'RF Bootstrap 选中率': frequencies['随机森林'].to_numpy(),
            'EFS-RA 综合秩聚合排名': frequencies['efs_ra_rank'].to_numpy(),
            '主要关联通路': [pathway_labels[groups.get(gene_id, frozenset())] for gene_id in ids],
        }
    )
    if result.isna().any().any() or not ids.is_unique or set(ids) != set(pathways['gene_id']):
        raise ValueError('Candidate table requires nine unique and fully annotated Ensembl IDs')
    result.to_csv(
        aggregation_tables / FROZEN_CANDIDATES_FILE, sep='\t', index=False, float_format='%.10g'
    )
