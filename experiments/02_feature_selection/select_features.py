"""Rank a fixed Stage 01 D3/D7 candidate pool and audit sample sensitivity."""

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[2]
STAGE_ONE = ROOT / 'outputs/01_preprocessing_differential'
STAGE_TWO = ROOT / 'outputs/02_feature_selection'
D3_D7_OUTPUT = STAGE_TWO / 'intersection_6'
OUTPUT = STAGE_TWO / 'intersection_22'
UNION_OUTPUT = STAGE_TWO / 'union_22'
MAIN_LIBRARY_COUNTS = {0: 4, 1: 4, 3: 2, 7: 4, 14: 4, 28: 4}
TOP_K = 20
SEED = 252027
MAX_ITER = 20_000
CLASSIFIERS = ('LASSO', 'ElasticNet', '随机森林')


def select_candidates(wald, use_union=False):
    wald = wald.set_index('基因ID')
    evidence = ['log2表达倍数变化', '校正后P值', '通过差异筛选']
    d3, d7 = (wald.loc[wald['比较'].eq(f'D{day}_vs_D0'), ['基因符号', *evidence]].rename(columns={
        'log2表达倍数变化': f'D{day}对D0_log2FC', '校正后P值': f'D{day}对D0_padj',
        '通过差异筛选': f'D{day}通过差异筛选',
    }) for day in (3, 7))
    merged = d3.join(d7.rename(columns={'基因符号': '基因符号_D7'}), how='outer', validate='one_to_one')
    pass3 = merged['D3通过差异筛选'].eq('是')
    pass7 = merged['D7通过差异筛选'].eq('是')
    keep = pass3 | pass7 if use_union else pass3 & pass7
    candidates = merged.loc[keep].sort_index()
    up3 = candidates['D3对D0_log2FC'].gt(0)
    up7 = candidates['D7对D0_log2FC'].gt(0)
    candidates['相对D0方向'] = np.select(
        [up3 & up7, ~up3 & ~up7, ~up3 & up7],
        ['均上调', '均下调', 'D3下调_D7上调'], default='D3上调_D7下调',
    )
    if use_union:
        candidates['基因符号'] = candidates['基因符号'].fillna(candidates['基因符号_D7'])
        for day, up, other_day in [(3, up3, 7), (7, up7, 3)]:
            column = f'D{day}通过差异筛选'
            only = candidates[column].eq('是') & ~candidates[f'D{other_day}通过差异筛选'].eq('是')
            candidates.loc[only, '相对D0方向'] = np.where(up[only], f'仅D{day}上调', f'仅D{day}下调')
        return candidates.drop(columns='基因符号_D7').assign(
            D3通过差异筛选=pass3.loc[candidates.index], D7通过差异筛选=pass7.loc[candidates.index],
        )
    return candidates.drop(columns=['基因符号_D7', 'D3通过差异筛选', 'D7通过差异筛选'])


def paired_omissions(days):
    return product(days.index[days.eq(3)], days.index[days.eq(7)])


def top_membership(scores, k):
    """Retain positive scores at or above the kth score, including boundary ties."""
    boundary = max(0, len(scores) - k)
    cutoff = np.partition(scores, boundary)[boundary]
    return (scores > 0) & (scores >= cutoff)


def fit_classifier(method, values, labels):
    if method == '随机森林':
        return RandomForestClassifier(n_estimators=512, max_features='sqrt',
                                      class_weight='balanced', random_state=SEED).fit(values, labels)
    l1_ratio = {'LASSO': 1.0, 'ElasticNet': 0.5}[method]
    fit = LogisticRegression(solver='saga', l1_ratio=l1_ratio, C=1, class_weight='balanced',
                             random_state=SEED, max_iter=MAX_ITER).fit(values, labels)
    if fit.n_iter_.max() >= MAX_ITER:
        raise RuntimeError(f'{method} 逻辑回归未收敛')
    return fit


def classifier_scores(method, classifier):
    return classifier.feature_importances_ if method == '随机森林' else classifier.coef_[0]


def mean_difference(vst, is_d7):
    return vst[is_d7].mean(axis=0) - vst[~is_d7].mean(axis=0)


def fit_rankers(expression, days):
    """Fit fixed-parameter rankers for D7 versus all other input libraries."""
    is_d7 = days.eq(7).to_numpy()
    vst = expression.to_numpy()
    scaled = StandardScaler().fit_transform(vst)
    scores = {'均值差': mean_difference(vst, is_d7)}
    for method in CLASSIFIERS:
        values = vst if method == '随机森林' else scaled
        scores[method] = classifier_scores(method, fit_classifier(method, values, is_d7))
    return scores


def run_selection(output, d3_d7_only=False, use_union=False):
    wald = pd.read_csv(STAGE_ONE / 'tables/timecourse_wald.tsv.gz', sep='\t', encoding='utf-8-sig')
    candidates = select_candidates(wald, use_union=use_union)
    metadata = pd.read_csv(STAGE_ONE / 'tables/sample_inclusion.csv', encoding='utf-8-sig', index_col='GSM编号')
    days = metadata.loc[metadata['主分析纳入'].eq('是'), '时间点（天）']
    if d3_d7_only:
        days = days.loc[days.isin([3, 7])]
    days = days.sort_index()
    expected_counts = {3: 2, 7: 4} if d3_d7_only else MAIN_LIBRARY_COUNTS
    if days.value_counts().to_dict() != expected_counts:
        raise ValueError(f'文库数量与阶段 01 主分析不符，预期各时点文库数为 {expected_counts}')

    vst = pd.read_csv(STAGE_ONE / 'matrices/timecourse_vst.tsv.gz', sep='\t', encoding='utf-8-sig', index_col='基因ID')
    expression = vst.loc[candidates.index, days.index].T
    if len(candidates) < TOP_K or not np.isfinite(expression.to_numpy()).all():
        raise ValueError('候选基因数量不足或 VST 矩阵中存在非有限数值')

    full = fit_rankers(expression, days)
    full_selected = {method: top_membership(np.abs(scores), TOP_K) for method, scores in full.items()}
    omissions = []
    detail = []
    deltas = []
    for fit_id, (d3, d7) in enumerate(paired_omissions(days), start=1):
        subset = days.drop([d3, d7])
        scores = fit_rankers(expression.loc[subset.index], subset)
        deltas.append(scores['均值差'])
        fit_summary = {'重抽样编号': fit_id, '排除D3': d3, '排除D7': d7, '方法': {}}
        for method, values in scores.items():
            selected = top_membership(np.abs(values), TOP_K)
            detail.extend((fit_id, method, gene, float(score))
                          for gene, score in zip(candidates.index[selected], values[selected]))
            reference = full_selected[method]
            union = np.count_nonzero(selected | reference)
            fit_summary['方法'][method] = {
                '非零分数基因数': int(np.count_nonzero(values)),
                '入选基因数': int(selected.sum()),
                '与全样本名单Jaccard': float(np.count_nonzero(selected & reference) / union) if union else None,
            }
        omissions.append(fit_summary)

    selected = pd.DataFrame(detail, columns=['重抽样编号', '方法', '基因ID', '分数'])
    reference_label = 'D3' if d3_d7_only else '非D7'
    ranking = candidates.assign(
        **{f'{reference_label}平均VST': expression.loc[days.ne(7)].mean(axis=0)},
        D7平均VST=expression.loc[days.eq(7)].mean(axis=0),
        **{f'D7减{reference_label}_VST均值差': full['均值差']},
        均值差绝对值=np.abs(full['均值差']),
        扰动方向一致率=((np.sign(deltas) == np.sign(full['均值差'])) & (full['均值差'] != 0)).mean(axis=0),
    )
    for method, scores in full.items():
        counts = selected.loc[selected['方法'].eq(method), '基因ID'].value_counts()
        ranking[f'{method}全样本分数'] = scores
        ranking[f'{method}全样本入选'] = full_selected[method]
        ranking[f'{method}入选频率'] = counts.reindex(ranking.index, fill_value=0) / len(omissions)
    ranking = ranking.sort_values(['均值差入选频率', '均值差绝对值', '基因ID'],
                                  ascending=[False, False, True])
    ranking.insert(0, '排名', np.arange(1, len(ranking) + 1))

    pool = '并集' if use_union else '交集'
    summary = {
        'scikit-learn': sklearn.__version__, '候选池': f'全队列阶段01的D3对D0和D7对D0差异基因{pool}',
        '表达矩阵': '全队列阶段01 VST', '评估范围': '固定候选池内文库扰动的描述性稳定性',
        '文库': {str(day): days.index[days.eq(day)].tolist() for day in sorted(days.unique())},
        f'{pool}基因数': len(ranking), '方向计数': ranking['相对D0方向'].value_counts().to_dict(),
        '每次目标入选数': TOP_K, '并列规则': '保留边界并列的正分数（用绝对值排序）',
        '随机种子': SEED, '逻辑回归参数': {'C': 1, 'solver': 'saga', 'max_iter': MAX_ITER,
                                           'class_weight': 'balanced', 'l1_ratio': [1.0, 0.5]},
        '随机森林参数': {'n_estimators': 512, 'max_features': 'sqrt', 'class_weight': 'balanced'},
        '全样本非零分数基因数': {method: int(np.count_nonzero(scores)) for method, scores in full.items()},
        '全样本入选基因数': {method: int(mask.sum()) for method, mask in full_selected.items()},
        '重抽样次数': len(omissions), '重抽样详情': omissions,
    }
    if not d3_d7_only:
        summary.update({'模型对比': 'D7 对其余五个时点（按文库）',
                        '输入矩阵形状': list(expression.shape), '扰动后文库数': len(days) - 2})
    tables = output / 'tables'
    reports = output / 'reports'
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    ranking.to_csv(tables / 'feature_ranking.tsv', sep='\t', encoding='utf-8-sig', float_format='%.10g')
    selected.to_csv(tables / 'perturbation_selections.tsv', sep='\t', encoding='utf-8-sig', index=False,
                    float_format='%.10g')
    (reports / 'selection_summary.json').write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8',
    )
    print(f'{len(ranking)} 个候选基因，{len(omissions)} 次文库扰动，输出位于 {output}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--d3-d7-only', action='store_true', help='在 D3 和 D7 的六份文库上筛选交集基因')
    modes.add_argument('--union', action='store_true', help='在全部 22 份文库上筛选 D3/D7 并集基因')
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    default_output = UNION_OUTPUT if args.union else D3_D7_OUTPUT if args.d3_d7_only else OUTPUT
    run_selection(args.output_dir or default_output, d3_d7_only=args.d3_d7_only, use_union=args.union)
