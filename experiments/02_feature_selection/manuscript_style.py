"""Shared formatting and export settings for the Stage 01–02 manuscript figures."""

import matplotlib as mpl

from select_features import STAGE_TWO


EXPORT_DIR = STAGE_TWO / 'manuscript'
CORE_GENES = ('C1qtnf3', 'Cthrc1', 'Postn', 'Wisp1', 'Cemip', 'Cxadr', 'Kif26b', 'Col24a1', 'Mrc1')
INK = '#192D3B'
BLUE = '#2A6385'
WARM = '#A64B32'
WIDTH_MM = 183


def configure_style():
    mpl.rcParams.update({
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
        'svg.fonttype': 'none',
        'svg.hashsalt': 'infarcttrace-manuscript',
        'pdf.fonttype': 42,
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
    })


def panel_heading(fig, x, y, label, title):
    fig.text(x, y, label, fontsize=8, fontweight='bold', va='baseline')
    fig.text(x + 0.022, y, title, fontsize=7, fontweight='bold', va='baseline')


def save_figure(fig, output, stem):
    output.mkdir(parents=True, exist_ok=True)
    for extension in ('pdf', 'svg', 'png', 'tiff'):
        kwargs = {'dpi': 600}
        if extension == 'png':
            kwargs['dpi'] = 240
        if extension == 'tiff':
            kwargs['pil_kwargs'] = {'compression': 'tiff_lzw'}
        fig.savefig(output / f'{stem}.{extension}', **kwargs)
