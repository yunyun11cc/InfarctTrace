"""对两外部队列的冻结候选执行表达、ROC/AUC 与共表达评估。

输入为步骤 01 矩阵和覆盖表、阶段二冻结基因与固定通路成员、发现队列注释。
输出表达差异、逐样本预测分数、ROC 指标、通路得分及相关性 TSV。
"""

from contextlib import closing
import sqlite3

import gseapy
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, spearmanr, ttest_ind
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from .config import (
    AUC_BOOTSTRAPS,
    AUC_QUANTILES,
    COHORTS,
    FROZEN,
    GENES,
    LOOCV_C,
    MARKERS,
    MOUSE_ANNOTATION,
    PATHWAY_MEMBERSHIP,
    PATHWAYS,
    PRIMARY_GROUP,
    SEED,
    SSGSEA_MAX_SIZE,
    SSGSEA_MIN_SIZE,
    SSGSEA_WEIGHT,
    WELCH_GATE_ALPHA,
)
from .tables import read, save


def aliases(index, symbol):
    """按冻结名称匹配，Wisp1 另查当前符号 Ccn4。"""
    options = [symbol] if symbol != 'Wisp1' else ['Wisp1', 'Ccn4']
    for name in options:
        hits = [x for x in index if x.casefold() == name.casefold()]
        if len(hits) > 1:
            raise ValueError(f'Ambiguous symbol {symbol}: {hits}')
        if hits:
            return hits[0]
    return None


def load_cohort(cohort, coverage, mapping, preprocessing, output):
    """从清洗矩阵对齐外部队列文库，返回全矩阵、冻结基因与标志物切片、样本表和映射。"""
    directory = preprocessing / cohort
    meta = read(directory / 'tables/sample_inclusion.tsv').set_index('GSM')
    if cohort == 'GSE111059':
        full = read(directory / 'matrices/gene_sst_rma.tsv.gz', index_col=0)
        full.index = full.index.astype(str)
        by_symbol = dict(zip(coverage.Gene_symbol, coverage.EntrezID.astype(str)))
        for marker in MARKERS:
            found = mapping.loc[mapping.symbol.eq(marker), 'entrez'].unique()
            if len(found) == 1:
                by_symbol[marker] = found[0]
    else:
        raw = read(directory / 'matrices/raw_counts.tsv.gz', index_col=0)
        full = np.log2(raw.div(raw.sum()) * 1e6 + 1)
        by_symbol = {g: aliases(full.index, g) for g in GENES + MARKERS}
    assert full.index.is_unique and set(full.columns) == set(meta.index)
    # 基因 × 文库矩阵按 GEO 样本表的 GSM 顺序对齐。
    full = full.loc[:, meta.index]
    labels = []
    sliced = pd.DataFrame(index=GENES + MARKERS, columns=full.columns, dtype=float)
    for symbol in sliced.index:
        key = by_symbol.get(symbol)
        present = key in full.index
        if present:
            sliced.loc[symbol] = full.loc[key]
        labels.append(dict(Cohort=cohort, Gene=symbol, MatrixID=key, Present=present))
    if sliced.loc[GENES].isna().any().any():
        raise ValueError(f'{cohort}: frozen panel is incomplete')
    sliced.loc[GENES].rename_axis('Gene').reset_index().pipe(save, f'{cohort}_frozen_expression', output)
    return full, sliced, meta, labels


def expression_tests(cohort, control, x, meta):
    """比较各时间组与该队列对照，返回 Welch 效应、P 值、BH q 及晋级标记。"""
    rows = []
    cases = {group: meta.Group.eq(group) for group in meta.Group.unique() if group != control}
    cases['MI_all_days'] = meta.Group.ne(control)
    for group, case in cases.items():
        for gene in GENES:
            a = x.loc[gene, case].to_numpy()
            b = x.loc[gene, meta.Group.eq(control)].to_numpy()
            effect = a.mean() - b.mean()
            rows.append(
                dict(
                    Cohort=cohort,
                    Gene=gene,
                    Contrast=f'{group}_vs_{control}',
                    Primary=group == PRIMARY_GROUP,
                    N_MI=len(a),
                    N_control=len(b),
                    Mean_MI=a.mean(),
                    Mean_control=b.mean(),
                    Log2_scale_difference=effect,
                    Direction='up' if effect > 0 else 'down' if effect < 0 else 'flat',
                    P=ttest_ind(a, b, equal_var=False).pvalue,
                )
            )
    table = pd.DataFrame(rows)
    for _, block in table[table.P.notna()].groupby('Contrast'):
        table.loc[block.index, 'BH_q'] = false_discovery_control(block.P)
    table['Pass'] = table.P.lt(WELCH_GATE_ALPHA) & table.Log2_scale_difference.gt(0)
    return table


def roc_scores(cohort, control, sliced, meta, passed, output):
    """以表达值评估单基因，并拟合原有 LOOCV 面板，返回逐文库分数和指标。"""
    groups = meta.Group[meta.Group.isin([control, PRIMARY_GROUP])]
    y = groups.eq(PRIMARY_GROUP).astype(int).to_numpy()
    x = sliced.loc[GENES, groups.index].T
    rng = np.random.default_rng(SEED)
    scores = {}
    rows = []
    for gene in GENES:
        score = x[gene].to_numpy()
        # 分层重采样保留正、负样本数，使用 ROC 的成对比较定义。
        a = rng.choice(score[y == 1], (AUC_BOOTSTRAPS, sum(y == 1)))
        b = rng.choice(score[y == 0], (AUC_BOOTSTRAPS, sum(y == 0)))
        bootstrap = ((a[:, :, None] > b[:, None, :]) + 0.5 * (a[:, :, None] == b[:, None, :])).mean(
            axis=(1, 2)
        )
        low, high = np.quantile(bootstrap, AUC_QUANTILES)
        scores[gene] = score
        rows.append(
            (
                cohort,
                gene,
                'fixed expression',
                roc_auc_score(y, score),
                low,
                high,
                average_precision_score(y, score),
            )
        )
    panels = [('Frozen 9 (LOOCV)', GENES, 'fixed panel, within-cohort LOOCV')]
    if passed:
        panels.append(('Gate subset (LOOCV)', passed, 'post-selection, within-cohort LOOCV'))
    for name, features, design in panels:
        pipe = make_pipeline(StandardScaler(), LogisticRegression(C=LOOCV_C, solver='lbfgs'))
        score = cross_val_predict(pipe, x[features], y, cv=LeaveOneOut(), method='predict_proba')[:, 1]
        scores[name] = score
        rows.append(
            (
                cohort,
                name,
                design,
                roc_auc_score(y, score),
                np.nan,
                np.nan,
                average_precision_score(y, score),
            )
        )
    scores = pd.DataFrame(scores, index=groups.index).assign(Cohort=cohort, Label=y, Group=groups)
    save(scores.rename_axis('GSM').reset_index(), f'{cohort}_roc_scores', output)
    return scores, pd.DataFrame(
        rows, columns=['Cohort', 'Score', 'Design', 'AUROC', 'CI_low', 'CI_high', 'AUPRC_AP']
    )


def correlations(cohort, x, other):
    """对每个冻结基因计算成对 Spearman 相关并在该队列内做 BH 校正。"""
    rows = []
    for gene in GENES:
        for name, values in other.iterrows():
            a = x.loc[gene]
            valid = a.notna() & values.notna()
            if valid.sum() >= 3 and a[valid].nunique() > 1 and values[valid].nunique() > 1:
                rho, p = spearmanr(a[valid], values[valid])
            else:
                rho = p = np.nan
            rows.append((cohort, gene, name, int(valid.sum()), rho, p))
    table = pd.DataFrame(rows, columns=['Cohort', 'Gene', 'Target', 'N', 'Spearman_rho', 'P'])
    valid = table.P.notna()
    table.loc[valid, 'BH_q'] = false_discovery_control(table.loc[valid, 'P'])
    return table


def pathway_scores(cohort, full, mapping, membership, output):
    """映射外部全矩阵的通路成员，写出 ssGSEA NES 并返回覆盖表。"""
    sets = {}
    for name, frame in membership.groupby('Pathway', sort=False):
        ids = set(frame.EntrezID.astype(str))
        members = ids if cohort == 'GSE111059' else set(mapping.loc[mapping.entrez.isin(ids), 'symbol'])
        sets[name] = sorted(members & set(full.index))
    coverage = pd.DataFrame(
        {'Cohort': cohort, 'Pathway': list(sets), 'Matched_features': [len(x) for x in sets.values()]}
    )
    print(coverage.to_string(index=False), flush=True)
    usable = {name: ids for name, ids in sets.items() if SSGSEA_MIN_SIZE <= len(ids) <= SSGSEA_MAX_SIZE}
    # GSEApy 将纯数字索引视为行序号。前缀同时应用于矩阵和基因集。
    score_input = full.rename(index=lambda gene: f'gene:{gene}')
    usable = {name: [f'gene:{gene}' for gene in ids] for name, ids in usable.items()}
    result = gseapy.ssgsea(
        data=score_input,
        gene_sets=usable,
        organism='mouse',
        outdir=None,
        sample_norm_method='rank',
        correl_norm_type='rank',
        min_size=SSGSEA_MIN_SIZE,
        max_size=SSGSEA_MAX_SIZE,
        weight=SSGSEA_WEIGHT,
        permutation_num=0,
        threads=1,
        seed=SEED,
        no_plot=True,
    )
    score = result.res2d.pivot(index='Term', columns='Name', values='NES').astype(float)
    score = score.reindex(index=sets, columns=full.columns)
    save(score.rename_axis('Pathway').reset_index(), f'{cohort}_ssgsea_scores', output)
    return score, coverage


def run(preprocessing, output):
    """分别评估两队列；按主比较双队列晋级标记保存全部表格。"""
    (output / 'tables').mkdir(parents=True, exist_ok=True)
    frozen = read(FROZEN)
    assert frozen['基因符号'].tolist() == GENES
    coverage = read(preprocessing / 'tables/candidate_coverage.tsv')
    assert frozen['基因ID'].tolist() == coverage.EnsemblID.tolist()
    with closing(sqlite3.connect(f"{MOUSE_ANNOTATION.as_uri()}?mode=ro", uri=True)) as db:
        mapping = pd.read_sql_query(
            'SELECT gene_id AS entrez, symbol FROM genes JOIN gene_info USING(_id)', db
        )
    membership = read(PATHWAY_MEMBERSHIP)
    membership = membership[membership.Pathway.isin(PATHWAYS)]
    assert membership.Pathway.nunique() == len(PATHWAYS)
    data, tests, mapping_rows = {}, [], []
    for cohort, control in COHORTS:
        full, sliced, meta, labels = load_cohort(cohort, coverage, mapping, preprocessing, output)
        data[cohort] = (full, sliced, meta)
        mapping_rows.extend(labels)
        tests.append(expression_tests(cohort, control, sliced, meta))
    combined = pd.concat(tests, ignore_index=True)
    save(combined, 'bulk_expression_tests', output)
    save(pd.DataFrame(mapping_rows), 'bulk_gene_mapping', output)
    passed = combined[combined.Primary].groupby('Gene').Pass.all()
    selected = [gene for gene in GENES if passed.loc[gene]]
    print('D7 gate:', selected, flush=True)
    roc, marker, pathway, align = [], [], [], []
    for cohort, control in COHORTS:
        full, x, meta = data[cohort]
        _, metrics = roc_scores(cohort, control, x, meta, selected, output)
        roc.append(metrics)
        marker.append(correlations(cohort, x, x.loc[MARKERS]))
        score, availability = pathway_scores(cohort, full, mapping, membership, output)
        align.append(availability)
        pathway.append(correlations(cohort, x, score))
    roc = pd.concat(roc, ignore_index=True)
    save(roc, 'bulk_roc_metrics', output)
    save(pd.concat(marker, ignore_index=True), 'marker_correlations', output)
    save(pd.concat(pathway, ignore_index=True), 'pathway_correlations', output)
    save(pd.concat(align, ignore_index=True), 'pathway_coverage', output)
