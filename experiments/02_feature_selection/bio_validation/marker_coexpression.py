"""用发现队列 VST 评估冻结九候选与七个标志物的共表达。

读取主分析 VST 与注释，输出候选与标志物相关性及表达源表和互证图面板。
"""

import numpy as np
import pandas as pd

from .ssgsea_pathway import as_rows, correlate, get_expression, load_candidate_expression
from ..common.plot_style import draw_heatmap
from ..config import MARKERS


def marker_data(expression, ids, markers, marker_ids, days, tables):
    """计算逐候选–标志物的相关和 BH 校正，导出统计及逐文库表达。"""
    results = as_rows(expression.columns, markers.columns, correlate(expression, markers))
    marker_status = pd.Series(
        np.select(
            [markers.isna().all(), markers.nunique().le(1)],
            ['not_quantified', 'constant'],
            default='measured',
        ),
        index=markers.columns,
    )
    results = results.assign(
        gene_id=results['gene'].map(ids),
        marker_id=results['feature'].map(marker_ids),
        group=results['feature'].map(MARKERS),
        status=results['feature'].map(marker_status),
    )
    results = results.rename(columns={'feature': 'marker'})
    results[['gene', 'gene_id', 'marker', 'marker_id', 'group', 'rho', 'p', 'q', 'status']].to_csv(
        tables / 'candidate_marker_correlations.tsv',
        sep='\t',
        index=False,
        float_format='%.12g',
        na_rep='NA',
    )
    expression.join(markers).assign(day=days)[['day', *expression.columns, *markers.columns]].to_csv(
        tables / 'candidate_marker_expression.tsv',
        sep='\t',
        index_label='library',
        float_format='%.12g',
        na_rep='NA',
    )
    return results


def run(tables):
    """解析七个标志物对应的基因 ID，写出共表达源表。"""
    vst, expression, ids, days, annotation = load_candidate_expression()
    marker_ids, markers = get_expression(vst, annotation, list(MARKERS))
    tables.mkdir(parents=True, exist_ok=True)
    return marker_data(expression, ids, markers, marker_ids, days, tables)


def draw_markers(ax, results):
    """绘制共表达分析图的标志物面板，共用通路面板的轴标签与色标。"""
    return draw_heatmap(ax, results, 'marker', list(MARKERS), 3, italic_columns=True, show_genes=False)
