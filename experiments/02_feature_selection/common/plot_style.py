"""阶段二现有单图、分析组图和相关热图的绘图配置。

函数设置 Matplotlib 字体、颜色和 DPI，绘图会写入调用方指定的目录。
"""

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle
import numpy as np

from ..config import CORE_GENES


BLUE = '#235F84'
ORANGE = '#B4552A'
GRAY = '#62717D'
LIGHT = '#DFE5EA'
INK = '#192D3B'
CHINESE_FONTS = ('Noto Sans CJK SC', 'PingFang SC', 'Microsoft YaHei', 'Source Han Sans SC', 'SimHei')
UNION_BLUE, UNION_RED, UNION_GRAY = '#1f5c7e', '#a63a2a', '#8c9ba5'
UNION_METHODS = {
    '均值差': ('VST mean difference', UNION_BLUE, 'o'),
    'LASSO': ('LASSO', '#596f80', 's'),
    'ElasticNet': ('Elastic Net', UNION_GRAY, '^'),
    '随机森林': ('Random forest', UNION_RED, 'D'),
}

PANEL_BLUE = '#2A6385'
PANEL_WARM = '#A64B32'
WIDTH_MM = 183
PATHWAY_LABELS = {
    'mmu04064': 'NF-κB',
    'R-MMU-844456': 'NLRP3\ninflammasome',
    'mmu04668': 'TNF',
    'mmu04060': 'Cytokine–\nreceptor',
    'mmu04350': 'TGF-β',
    'R-MMU-1474244': 'ECM\norganization',
    'R-MMU-1474290': 'Collagen\nformation',
}
NORM = TwoSlopeNorm(0, vmin=-1, vmax=1)
CORRELATION_CMAP = LinearSegmentedColormap.from_list(
    'blue_white_rust', [PANEL_BLUE, '#FAF9F5', PANEL_WARM]
)
CORRELATION_CMAP.set_bad('#E6E6E6')


def configure_style():
    """应用并集图的字体、边框和配色。"""
    matplotlib.rcParams.update(
        {
            'font.family': 'sans-serif',
            'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
            'font.size': 7,
            'text.color': INK,
            'axes.labelcolor': INK,
            'axes.edgecolor': INK,
            'axes.linewidth': 0.65,
            'axes.spines.top': False,
            'axes.spines.right': False,
            'xtick.color': INK,
            'ytick.color': INK,
            'xtick.labelsize': 6,
            'ytick.labelsize': 6,
            'xtick.major.size': 2.5,
            'ytick.major.size': 2.5,
            'legend.frameon': False,
            'legend.fontsize': 6,
            'figure.facecolor': 'white',
            'savefig.facecolor': 'white',
        }
    )


def configure_intersection_style():
    """选取已安装中文字体并设置交集图的轴样式。"""
    fonts = {font.name for font in font_manager.fontManager.ttflist}
    family = next((name for name in CHINESE_FONTS if name in fonts), None)
    if family is None:
        raise RuntimeError('未找到中文字体。请安装 Noto Sans CJK SC 后重新运行。')
    plt.rcParams.update(
        {
            'font.family': [family, 'DejaVu Sans'],
            'font.size': 10,
            'axes.unicode_minus': False,
            'axes.spines.top': False,
            'axes.spines.right': False,
            'axes.linewidth': 0.8,
        }
    )


def union_figure(height, title):
    """建立带主标题的并集图画布。"""
    fig = plt.figure(figsize=(10.6, height))
    fig.text(0.135, 0.95, title, fontsize=19, fontweight='bold', color=INK)
    return fig


def union_axes(ax, grid_axis):
    """为并集图坐标轴设置网格方向、边框和刻度。"""
    ax.set_axisbelow(True)
    ax.grid(axis=grid_axis, color=LIGHT, linewidth=0.7)
    ax.spines[['left', 'bottom']].set_color(UNION_GRAY)
    ax.tick_params(labelsize=9)


def configure_panel_style():
    """应用分析面板的排版和字体设置。"""
    matplotlib.rcParams.update(
        {
            'font.family': 'sans-serif',
            'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
            'font.size': 7,
            'text.color': INK,
            'axes.labelcolor': INK,
            'axes.edgecolor': INK,
            'axes.linewidth': 0.65,
            'axes.spines.top': False,
            'axes.spines.right': False,
            'xtick.color': INK,
            'ytick.color': INK,
            'xtick.labelsize': 6,
            'ytick.labelsize': 6,
            'xtick.major.size': 2.5,
            'ytick.major.size': 2.5,
            'legend.frameon': False,
            'legend.fontsize': 6,
            'figure.facecolor': 'white',
            'savefig.facecolor': 'white',
        }
    )


def panel_heading(fig, x, y, label, title):
    """在图的指定坐标添加面板序号和标题。"""
    fig.text(x, y, label, fontsize=8, fontweight='bold', va='baseline')
    fig.text(x + 0.022, y, title, fontsize=7, fontweight='bold', va='baseline')


def save_panel_figure(fig, output, stem):
    """将分析图保存为 240 DPI 的 PNG。"""
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / f'{stem}.png', dpi=240)


def draw_heatmap(
    ax, results, feature, labels, inflammation_columns, italic_columns=False, show_genes=True
):
    """用共同色标绘制 Spearman ρ，按传入的 q 值标注显著性。"""
    rho = results.pivot(index='gene', columns=feature, values='rho')
    q = results.pivot(index='gene', columns=feature, values='q')
    ordered_columns = results[feature].unique()
    rho = rho.loc[list(CORE_GENES), ordered_columns].to_numpy()
    q = q.loc[list(CORE_GENES), ordered_columns].to_numpy()
    image = ax.imshow(rho, cmap=CORRELATION_CMAP, norm=NORM, interpolation='nearest', aspect='auto')
    ax.set_xticks(range(len(labels)), labels, fontstyle='italic' if italic_columns else 'normal')
    ax.set_yticks(
        range(len(CORE_GENES)), CORE_GENES if show_genes else [''] * len(CORE_GENES), fontstyle='italic'
    )
    ax.tick_params(axis='both', length=0, pad=3, labelsize=6)
    ax.xaxis.tick_top()
    for tick in ax.get_yticklabels():
        if tick.get_text() == 'Mrc1':
            tick.set_color(PANEL_WARM)
            tick.set_fontweight('bold')
    ax.spines[:].set_visible(False)
    for x in np.arange(0.5, rho.shape[1] - 0.5):
        ax.axvline(x, color='white', lw=1.5 if x == inflammation_columns - 0.5 else 0.45)
    for y in np.arange(0.5, rho.shape[0] - 0.5):
        ax.axhline(y, color='white', lw=0.45)
    for (row, col), value in np.ndenumerate(rho):
        stars = '**' if q[row, col] < 0.01 else '*' if q[row, col] < 0.05 else ''
        color = 'white' if abs(value) >= 0.78 else INK
        label = f'{value:.2f}{stars}'.replace('-', '−') if np.isfinite(value) else 'NA'
        ax.text(col, row, label, ha='center', va='center', color=color, fontsize=5.5)
    ax.add_patch(
        Rectangle(
            (-0.48, len(CORE_GENES) - 1.48),
            len(labels) - 0.04,
            0.96,
            fill=False,
            edgecolor=PANEL_WARM,
            lw=0.85,
            clip_on=False,
        )
    )
    for left, right, label, color in [
        (-0.45, inflammation_columns - 0.55, 'Inflammation', BLUE),
        (inflammation_columns - 0.45, len(labels) - 0.55, 'Fibrosis', PANEL_WARM),
    ]:
        transform = ax.get_xaxis_transform()
        ax.plot([left, right], [1.18, 1.18], transform=transform, color=color, lw=1, clip_on=False)
        ax.text(
            (left + right) / 2,
            1.205,
            label,
            transform=transform,
            ha='center',
            va='bottom',
            color=color,
            fontsize=6.5,
            fontweight='bold',
        )
    return image
