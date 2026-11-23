"""将四个选择器在联合排序和文库扰动中的结果单独导出。

读取内存中的 feature_ranking 表，写入各方法的 TSV 和频率图。
"""

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

from .plot_style import UNION_METHODS, configure_style
from ..config import TOP_K


def write_method_outputs(ranking, output, method, stem):
    """导出一个方法的全队列得分、入选标记和配对排除入选频率。"""
    label, color, _ = UNION_METHODS[method]
    columns = ['基因符号', f'{method}全样本分数', f'{method}全样本入选', f'{method}入选频率']
    tables = output / 'tables'
    figures = output / 'figures'
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    table = ranking[columns].copy()
    table.to_csv(tables / f'{stem}.tsv', sep='\t', encoding='utf-8-sig', float_format='%.10g')

    # 同频率时按绝对分数和基因 ID 决定展示顺序，不修改原联合排序。
    table['绝对分数'] = table[columns[1]].abs()
    ranked = table.sort_values([columns[3], '绝对分数', '基因ID'], ascending=[False, False, True])
    top = ranked.head(TOP_K).iloc[::-1]
    symbols = top[columns[0]].fillna(top.index.to_series())
    with plt.rc_context():
        configure_style()
        fig, ax = plt.subplots(figsize=(6.6, 5.6), layout='constrained')
        order = range(len(top))
        ax.barh(order, top[columns[3]], color=color)
        ax.set(
            yticks=order,
            yticklabels=symbols,
            xlim=(0, 1),
            xlabel='Selection frequency across paired omissions',
            title=label,
        )
        ax.xaxis.set_major_formatter(PercentFormatter(1))
        fig.savefig(figures / f'{stem}.png', dpi=200)
        plt.close(fig)
