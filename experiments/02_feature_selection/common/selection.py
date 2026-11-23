"""组合四种基础选择器，计算并写出文库配对排除的稳定性。

读取阶段一的候选窗口、GSM 纳入表和 VST。
输出 feature_ranking.tsv、perturbation_selections.tsv 和 selection_summary.json。
"""

from itertools import product
import json

import numpy as np
import pandas as pd
import sklearn
from sklearn.preprocessing import StandardScaler

from .data_loader import read_main_days, read_vst, read_wald, window
from ..config import (
    D3_D7_LIBRARY_COUNTS,
    ELASTIC_NET_L1_RATIO,
    FOREST_MAX_FEATURES,
    FOREST_TREES,
    LASSO_L1_RATIO,
    LOGISTIC_C,
    MAIN_LIBRARY_COUNTS,
    MAX_ITER,
    SEED,
    TOP_K,
)
from ..methods import m01_vst_mean_diff, m02_lasso, m03_elastic_net, m04_random_forest


RANKERS = (m01_vst_mean_diff, m02_lasso, m03_elastic_net, m04_random_forest)
CLASSIFIER_MODULES = RANKERS[1:]


def paired_omissions(days):
    """逐一生成各排除一份 D3 和 D7 文库的 GSM 对。"""
    return product(days.index[days.eq(3)], days.index[days.eq(7)])


def top_membership(scores, k):
    """Retain positive scores at or above the kth score, including boundary ties."""
    boundary = max(0, len(scores) - k)
    cutoff = np.partition(scores, boundary)[boundary]
    return (scores > 0) & (scores >= cutoff)


def fit_rankers(expression, days):
    """返回按方法顺序组织的全样本或子集得分，行标签决定 D7 阳性文库。"""
    is_d7 = days.eq(7).to_numpy()
    vst = expression.to_numpy()
    scaled = StandardScaler().fit_transform(vst)
    scores = {m01_vst_mean_diff.METHOD: m01_vst_mean_diff.score(vst, is_d7)}

    # 逻辑回归使用同一批训练文库的标准化矩阵；随机森林直接使用原 VST。
    for method in CLASSIFIER_MODULES:
        values = scaled if method.USES_SCALING else vst
        model = method.fit(values, is_d7)
        scores[method.METHOD] = method.scores(model)
    return scores


def run_selection(output, d3_d7_only=False, use_union=False):
    """构建所选模式的候选排序和八次扰动明细，写表并返回结果矩阵。"""
    candidates = window.select_candidates(read_wald(), use_union=use_union)
    days = read_main_days()
    if d3_d7_only:
        days = days.loc[days.isin([3, 7])]
    days = days.sort_index()
    expected_counts = D3_D7_LIBRARY_COUNTS if d3_d7_only else MAIN_LIBRARY_COUNTS
    if days.value_counts().to_dict() != expected_counts:
        raise ValueError(f'文库数量与阶段 01 主分析不符，预期各时点文库数为 {expected_counts}')

    # 阶段一为基因 × 文库，选择器以文库为行、候选基因为列。
    expression = read_vst().loc[candidates.index, days.index].T
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
            detail.extend(
                (fit_id, method, gene, float(score))
                for gene, score in zip(candidates.index[selected], values[selected])
            )
            reference = full_selected[method]
            union = np.count_nonzero(selected | reference)
            fit_summary['方法'][method] = {
                '非零分数基因数': int(np.count_nonzero(values)),
                '入选基因数': int(selected.sum()),
                '与全样本名单Jaccard': float(np.count_nonzero(selected & reference) / union)
                if union
                else None,
            }
        omissions.append(fit_summary)

    selected = pd.DataFrame(detail, columns=['重抽样编号', '方法', '基因ID', '分数'])
    reference_label = 'D3' if d3_d7_only else '非D7'
    ranking = candidates.assign(
        **{f'{reference_label}平均VST': expression.loc[days.ne(7)].mean(axis=0)},
        D7平均VST=expression.loc[days.eq(7)].mean(axis=0),
        **{f'D7减{reference_label}_VST均值差': full['均值差']},
        均值差绝对值=np.abs(full['均值差']),
        扰动方向一致率=((np.sign(deltas) == np.sign(full['均值差'])) & (full['均值差'] != 0)).mean(
            axis=0
        ),
    )
    for method, scores in full.items():
        counts = selected.loc[selected['方法'].eq(method), '基因ID'].value_counts()
        ranking[f'{method}全样本分数'] = scores
        ranking[f'{method}全样本入选'] = full_selected[method]
        ranking[f'{method}入选频率'] = counts.reindex(ranking.index, fill_value=0) / len(omissions)
    ranking = ranking.sort_values(
        ['均值差入选频率', '均值差绝对值', '基因ID'], ascending=[False, False, True]
    )
    ranking.insert(0, '排名', np.arange(1, len(ranking) + 1))

    pool = '并集' if use_union else '交集'
    summary = {
        'scikit-learn': sklearn.__version__,
        '候选池': f'全队列阶段01的D3对D0和D7对D0差异基因{pool}',
        '表达矩阵': '全队列阶段01 VST',
        '评估范围': '固定候选池内文库扰动的描述性稳定性',
        '文库': {str(day): days.index[days.eq(day)].tolist() for day in sorted(days.unique())},
        f'{pool}基因数': len(ranking),
        '方向计数': ranking['相对D0方向'].value_counts().to_dict(),
        '每次目标入选数': TOP_K,
        '并列规则': '保留边界并列的正分数（用绝对值排序）',
        '随机种子': SEED,
        '逻辑回归参数': {
            'C': LOGISTIC_C,
            'solver': 'saga',
            'max_iter': MAX_ITER,
            'class_weight': 'balanced',
            'l1_ratio': [LASSO_L1_RATIO, ELASTIC_NET_L1_RATIO],
        },
        '随机森林参数': {
            'n_estimators': FOREST_TREES,
            'max_features': FOREST_MAX_FEATURES,
            'class_weight': 'balanced',
        },
        '全样本非零分数基因数': {
            method: int(np.count_nonzero(scores)) for method, scores in full.items()
        },
        '全样本入选基因数': {method: int(mask.sum()) for method, mask in full_selected.items()},
        '重抽样次数': len(omissions),
        '重抽样详情': omissions,
    }
    if not d3_d7_only:
        summary.update(
            {
                '模型对比': 'D7 对其余五个时点（按文库）',
                '输入矩阵形状': list(expression.shape),
                '扰动后文库数': len(days) - 2,
            }
        )
    tables = output / 'tables'
    reports = output / 'reports'
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    ranking.to_csv(tables / 'feature_ranking.tsv', sep='\t', encoding='utf-8-sig', float_format='%.10g')
    selected.to_csv(
        tables / 'perturbation_selections.tsv',
        sep='\t',
        encoding='utf-8-sig',
        index=False,
        float_format='%.10g',
    )
    (reports / 'selection_summary.json').write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + '\n',
        encoding='utf-8',
    )
    print(f'{len(ranking)} 个候选基因，{len(omissions)} 次文库扰动，输出位于 {output}')
    return ranking
