"""Draw the bootstrap Jaccard and mean Top 20 selection probability plots."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

from resampling import METHODS, OUTPUT


def visualize(output):
    detail = pd.read_csv(output / 'tables/bootstrap_jaccard.tsv', sep='\t')
    probabilities = pd.read_csv(output / 'tables/bootstrap_selection_probability.tsv', sep='\t')
    summary = json.loads((output / 'reports/resampling_summary.json').read_text())
    top = probabilities.head(20)
    labels = top['基因符号'].fillna(top['基因ID'])
    values = top['Mean'].to_numpy()
    groups = [detail.loc[detail['Method'].eq(method), 'Jaccard'].dropna() for method in METHODS]
    colors = ['#286485', '#8B6456', '#786894', '#547B71']
    fig, (box, lollipop) = plt.subplots(1, 2, figsize=(12, 7.3), layout='constrained',
                                       gridspec_kw={'width_ratios': [0.95, 1]})
    boxes = box.boxplot(groups, tick_labels=['VST mean\ndifference', 'LASSO', 'Elastic Net', 'Random\nforest'],
                        patch_artist=True, medianprops={'color': '#243744', 'linewidth': 1.6},
                        flierprops={'markersize': 2, 'alpha': 0.48})
    for patch, color in zip(boxes['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
    box.set(ylim=(-0.025, 1.025), ylabel='Jaccard with full-cohort selection')
    box.set_title('a   Selection overlap across bootstrap draws', loc='left', fontsize=12, fontweight='bold')
    box.grid(axis='y', color='#E1E8ED', linewidth=0.75)

    rows = np.arange(len(top))
    lollipop.hlines(rows, 0, values, color='#A7C0CE', linewidth=1.5)
    lollipop.scatter(values, rows, color='#286485', s=33, zorder=3)
    for row, value in zip(rows, values):
        lollipop.annotate(f'{value:.1%}', (value, row), xytext=(5, 0), textcoords='offset points',
                          va='center', fontsize=8, color='#243744')
    lollipop.set(yticks=rows, yticklabels=labels, xlim=(0, 1.12), xlabel='Mean Top 20 selection probability')
    lollipop.set_xticks(np.linspace(0, 1, 6))
    lollipop.xaxis.set_major_formatter(PercentFormatter(1))
    lollipop.invert_yaxis()
    lollipop.set_title('b   Top genes by selection probability', loc='left', fontsize=12, fontweight='bold')
    lollipop.grid(axis='x', color='#E1E8ED', linewidth=0.75)
    for ax in (box, lollipop):
        ax.set_axisbelow(True)
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(labelsize=9)
    count = summary['bootstrap_repeats']
    fig.suptitle(f'Selection stability across {count:,} stratified bootstrap draws',
                 fontsize=16, fontweight='bold')
    fig.supxlabel('Each draw samples 4 D7 and 18 non-D7 libraries with replacement.\n'
                  'Selection includes boundary ties. Probabilities are averaged across four methods.', fontsize=9)
    figures = output / 'figures'
    figures.mkdir(parents=True, exist_ok=True)
    fig.savefig(figures / '03_selection_stability.png', dpi=220)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    visualize(parser.parse_args().output_dir)
