"""Rank the fixed Stage 01 D3/D7 intersection and audit sample sensitivity."""

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
OUTPUT = ROOT / 'outputs/02_feature_selection'
TOP_K = 20
SEED = 252027
MAX_ITER = 20_000


def intersect_candidates(wald):
    passed = wald.loc[wald['通过差异筛选'].eq('是')].set_index('基因ID')
    evidence = ['log2表达倍数变化', '校正后P值']
    d3 = passed.loc[passed['比较'].eq('D3_vs_D0'), ['基因符号', *evidence]]
    d7 = passed.loc[passed['比较'].eq('D7_vs_D0'), evidence]
    candidates = d3.rename(columns=dict(zip(evidence, ['D3对D0_log2FC', 'D3对D0_padj']))).join(
        d7.rename(columns=dict(zip(evidence, ['D7对D0_log2FC', 'D7对D0_padj']))),
        how='inner', validate='one_to_one',
    ).sort_index()
    up3 = candidates['D3对D0_log2FC'].gt(0)
    up7 = candidates['D7对D0_log2FC'].gt(0)
    return candidates.assign(相对D0方向=np.select(
        [up3 & up7, ~up3 & ~up7, ~up3 & up7],
        ['均上调', '均下调', 'D3下调_D7上调'], default='D3上调_D7下调',
    ))


def paired_omissions(days):
    return product(days.index[days.eq(3)], days.index[days.eq(7)])


def top_membership(scores, k):
    """Retain positive scores at or above the kth score, including boundary ties."""
    boundary = max(0, len(scores) - k)
    cutoff = np.partition(scores, boundary)[boundary]
    return (scores > 0) & (scores >= cutoff)


def fit_rankers(expression, days):
    """Fit fixed-parameter rankers on one subset of the six libraries."""
    is_d7 = days.eq(7).to_numpy()
    vst = expression.to_numpy()
    scaled = StandardScaler().fit_transform(vst)
    mean_difference = vst[is_d7].mean(axis=0) - vst[~is_d7].mean(axis=0)

    def logistic_coefficients(l1_ratio):
        fit = LogisticRegression(solver='saga', l1_ratio=l1_ratio, C=1,
                                 class_weight='balanced', random_state=SEED,
                                 max_iter=MAX_ITER).fit(scaled, is_d7)
        if fit.n_iter_.max() >= MAX_ITER:
            raise RuntimeError('稀疏逻辑回归未收敛')
        return fit.coef_[0]

    forest = RandomForestClassifier(n_estimators=512, max_features='sqrt',
                                    class_weight='balanced', random_state=SEED).fit(vst, is_d7)
    return {'均值差': mean_difference, 'LASSO': logistic_coefficients(1.0),
            'ElasticNet': logistic_coefficients(0.5), '随机森林': forest.feature_importances_}


def run_selection(output):
    wald = pd.read_csv(STAGE_ONE / 'tables/timecourse_wald.tsv.gz', sep='\t', encoding='utf-8-sig')
    candidates = intersect_candidates(wald)
    metadata = pd.read_csv(STAGE_ONE / 'tables/sample_inclusion.csv', encoding='utf-8-sig', index_col='GSM编号')
    days = metadata.loc[metadata['主分析纳入'].eq('是'), '时间点（天）']
    days = days.loc[days.isin([3, 7])].sort_index()
    if days.value_counts().to_dict() != {3: 2, 7: 4}:
        raise ValueError('当前实验需要阶段 01 主分析中的 2 份 D3 和 4 份 D7 文库')

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
    ranking = candidates.assign(
        D3平均VST=expression.loc[days.eq(3)].mean(axis=0),
        D7平均VST=expression.loc[days.eq(7)].mean(axis=0),
        D7减D3_VST均值差=full['均值差'],
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

    summary = {
        'scikit-learn': sklearn.__version__, '候选池': '全队列阶段01的D3对D0和D7对D0差异基因交集',
        '表达矩阵': '全队列阶段01 VST', '评估范围': '固定候选池内文库扰动的描述性稳定性',
        '文库': {str(day): days.index[days.eq(day)].tolist() for day in (3, 7)},
        '交集基因数': len(ranking), '方向计数': ranking['相对D0方向'].value_counts().to_dict(),
        '每次目标入选数': TOP_K, '并列规则': '保留边界并列的正分数（用绝对值排序）',
        '随机种子': SEED, '逻辑回归参数': {'C': 1, 'solver': 'saga', 'max_iter': MAX_ITER,
                                           'class_weight': 'balanced', 'l1_ratio': [1.0, 0.5]},
        '随机森林参数': {'n_estimators': 512, 'max_features': 'sqrt', 'class_weight': 'balanced'},
        '全样本非零分数基因数': {method: int(np.count_nonzero(scores)) for method, scores in full.items()},
        '全样本入选基因数': {method: int(mask.sum()) for method, mask in full_selected.items()},
        '重抽样次数': len(omissions), '重抽样详情': omissions,
    }
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
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    run_selection(parser.parse_args().output_dir)
