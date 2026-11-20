"""阶段 03 的冻结候选外部验证，以共同采样日 D7 为主要比较。"""
from pathlib import Path
import sqlite3
from contextlib import closing

import gseapy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, spearmanr, ttest_ind
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/03_validation'
PREP = ROOT / 'outputs/01_2_validation_preprocessing'
FROZEN = ROOT / 'outputs/02_feature_selection/manuscript/table_02_candidate_genes_summary.tsv'
GENES = ['C1qtnf3', 'Cthrc1', 'Postn', 'Wisp1', 'Cemip', 'Cxadr', 'Kif26b', 'Col24a1', 'Mrc1']
MARKERS = ['Col1a1', 'Col3a1', 'Acta2', 'Tgfb1', 'Nlrp3', 'Il1b']
COHORTS = [('GSE111059', 'uninjured'), ('GSE236374', 'sham')]


def read(path, **kwargs):
    return pd.read_csv(path, sep='\t', **kwargs)


def save(frame, name, **kwargs):
    frame.to_csv(OUT / 'tables' / f'{name}.tsv', sep='\t', index=False, na_rep='NA', **kwargs)


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


def load_cohort(cohort, coverage, mapping):
    directory = PREP / cohort
    meta = read(directory / 'tables/sample_inclusion.tsv').set_index('GSM')
    if cohort == 'GSE111059':
        full = read(directory / 'matrices/gene_sst_rma.tsv.gz', index_col=0)
        full.index = full.index.astype(str)
        by_symbol = dict(zip(coverage.Gene_symbol, coverage.EntrezID.astype(str)))
        for marker in MARKERS:
            found = mapping.loc[mapping.symbol.eq(marker), 'entrez'].unique()
            if len(found) == 1:
                by_symbol[marker] = found[0]
        scale = 'SST-RMA log2 signal'
    else:
        raw = read(directory / 'matrices/raw_counts.tsv.gz', index_col=0)
        full = np.log2(raw.div(raw.sum()) * 1e6 + 1)
        by_symbol = {g: aliases(full.index, g) for g in GENES + MARKERS}
        scale = 'log2(CPM+1)'
    assert full.index.is_unique and set(full.columns) == set(meta.index)
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
    sliced.loc[GENES].rename_axis('Gene').reset_index().pipe(save, f'{cohort}_frozen_expression')
    return full, sliced, meta, scale, labels


def expression_tests(cohort, control, x, meta):
    rows = []
    cases = {group: meta.Group.eq(group) for group in meta.Group.unique() if group != control}
    cases['MI_all_days'] = meta.Group.ne(control)
    for group, case in cases.items():
        for gene in GENES:
            a = x.loc[gene, case].to_numpy()
            b = x.loc[gene, meta.Group.eq(control)].to_numpy()
            effect = a.mean() - b.mean()
            rows.append(dict(Cohort=cohort, Gene=gene, Contrast=f'{group}_vs_{control}',
                             Primary=group == 'MI_D7', N_MI=len(a), N_control=len(b),
                             Mean_MI=a.mean(), Mean_control=b.mean(), Log2_scale_difference=effect,
                             Direction='up' if effect > 0 else 'down' if effect < 0 else 'flat',
                             P=ttest_ind(a, b, equal_var=False).pvalue))
    table = pd.DataFrame(rows)
    for _, block in table[table.P.notna()].groupby('Contrast'):
        table.loc[block.index, 'BH_q'] = false_discovery_control(block.P)
    table['Pass'] = table.P.lt(.05) & table.Log2_scale_difference.gt(0)
    return table


def roc_scores(cohort, control, sliced, meta, passed):
    groups = meta.Group[meta.Group.isin([control, 'MI_D7'])]
    y = groups.eq('MI_D7').astype(int).to_numpy()
    x = sliced.loc[GENES, groups.index].T
    rng = np.random.default_rng(252027)
    scores = {}
    rows = []
    for gene in GENES:
        score = x[gene].to_numpy()
        # 分层重采样保留正、负样本数，使用 ROC 的成对比较定义。
        a = rng.choice(score[y == 1], (10000, sum(y == 1)))
        b = rng.choice(score[y == 0], (10000, sum(y == 0)))
        bootstrap = ((a[:, :, None] > b[:, None, :]) + .5 *
                     (a[:, :, None] == b[:, None, :])).mean(axis=(1, 2))
        low, high = np.quantile(bootstrap, [.025, .975])
        scores[gene] = score
        rows.append((cohort, gene, 'fixed expression', roc_auc_score(y, score),
                     low, high, average_precision_score(y, score)))
    panels = [('Frozen 9 (LOOCV)', GENES, 'fixed panel, within-cohort LOOCV')]
    if passed:
        panels.append(('Gate subset (LOOCV)', passed, 'post-selection, within-cohort LOOCV'))
    for name, features, design in panels:
        pipe = make_pipeline(StandardScaler(), LogisticRegression(C=1, solver='lbfgs'))
        score = cross_val_predict(pipe, x[features], y, cv=LeaveOneOut(), method='predict_proba')[:, 1]
        scores[name] = score
        rows.append((cohort, name, design, roc_auc_score(y, score), np.nan, np.nan,
                     average_precision_score(y, score)))
    scores = pd.DataFrame(scores, index=groups.index).assign(Cohort=cohort, Label=y, Group=groups)
    save(scores.rename_axis('GSM').reset_index(), f'{cohort}_roc_scores')
    return scores, pd.DataFrame(rows, columns=['Cohort', 'Score', 'Design', 'AUROC', 'CI_low', 'CI_high', 'AUPRC_AP'])


def correlations(cohort, x, other):
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


def pathway_scores(cohort, full, mapping, membership):
    sets = {}
    for name, frame in membership.groupby('Pathway', sort=False):
        ids = set(frame.EntrezID.astype(str))
        members = ids if cohort == 'GSE111059' else set(mapping.loc[mapping.entrez.isin(ids), 'symbol'])
        sets[name] = sorted(members & set(full.index))
    coverage = pd.DataFrame({'Cohort': cohort, 'Pathway': list(sets),
                             'Matched_features': [len(x) for x in sets.values()]})
    print(coverage.to_string(index=False), flush=True)
    usable = {name: ids for name, ids in sets.items() if 5 <= len(ids) <= 500}
    # GSEApy 将纯数字索引视为行序号。前缀同时应用于矩阵和基因集。
    score_input = full.rename(index=lambda gene: f'gene:{gene}')
    usable = {name: [f'gene:{gene}' for gene in ids] for name, ids in usable.items()}
    result = gseapy.ssgsea(data=score_input, gene_sets=usable, organism='mouse', outdir=None,
                          sample_norm_method='rank', correl_norm_type='rank', min_size=5, max_size=500,
                          weight=.25, permutation_num=0, threads=1, seed=252027, no_plot=True)
    score = result.res2d.pivot(index='Term', columns='Name', values='NES').astype(float)
    score = score.reindex(index=sets, columns=full.columns)
    save(score.rename_axis('Pathway').reset_index(), f'{cohort}_ssgsea_scores')
    return score, coverage


def boxplot(cohort, x, meta, scale, tests, control):
    fig, axes = plt.subplots(3, 3, figsize=(11, 8), layout='constrained')
    order = list(meta.Group.unique())
    for ax, gene in zip(axes.flat, GENES):
        groups = [x.loc[gene, meta.Group.eq(group)].to_numpy() for group in order]
        ax.boxplot(groups, tick_labels=[group.replace('MI_', '') for group in order], showfliers=False)
        for j, values in enumerate(groups, 1):
            ax.scatter(j + np.linspace(-.11, .11, len(values)), values, s=21, color='#a65439' if j > 1 else '#276687', zorder=3)
        p = tests.loc[tests.Gene.eq(gene) & tests.Primary, 'P'].iat[0]
        ax.set_title(f'{gene}   D7 P={p:.3g}', fontsize=10)
        ax.tick_params(axis='x', labelsize=8)
    fig.suptitle(f'{cohort} | candidate expression by time\nWelch P: D7 vs {control}')
    fig.supylabel(scale)
    fig.savefig(OUT / 'figures' / f'{cohort}_expression.png', dpi=200)
    plt.close(fig)


def plot_expression_auc(scores_by_cohort, metrics):
    """用已保存的单基因预测分数画表达分布及对应 AUC。"""
    available = {f.name for f in font_manager.fontManager.ttflist}
    font = next((name for name in ('Microsoft YaHei', 'Noto Sans CJK SC', 'PingFang SC', 'SimHei')
                 if name in available), None)
    ink, blue, orange = '#192d3b', '#256080', '#b5543a'
    auc = metrics.set_index(['Cohort', 'Score']).AUROC
    style = {'font.family': ['DejaVu Sans', font] if font else ['DejaVu Sans'],
             'text.color': ink, 'axes.labelcolor': ink, 'xtick.color': ink, 'ytick.color': ink}
    with plt.rc_context(style):
        fig = plt.figure(figsize=(14, 10), facecolor='white')
        grid = fig.add_gridspec(3, 3, left=.055, right=.98, top=.86, bottom=.085, wspace=.31, hspace=.32)
        fig.suptitle('九个候选基因的外部表达验证' if font else 'External expression validation of nine candidates',
                     y=.971, fontsize=24, fontweight='bold')
        fig.text(.5, .918, '对照 vs 心梗后 D7  |  每队列每组 3 份样本' if font
                 else 'Control vs day 7 after MI  |  3 samples per group and cohort',
                 ha='center', fontsize=11, color='#455b6b')
        for cell, gene in zip(grid, GENES):
            inner = cell.subgridspec(2, 2, height_ratios=[.19, 1], hspace=.20, wspace=.42)
            title = fig.add_subplot(inner[0, :])
            title.set_axis_off()
            title.text(.5, .47, gene, ha='center', va='center', fontsize=15, fontstyle='italic',
                       fontweight='bold', color=orange if gene == 'Cxadr' else ink)
            for j, (cohort, control) in enumerate(COHORTS):
                ax = fig.add_subplot(inner[1, j])
                score = scores_by_cohort[cohort]
                groups = [score.loc[score.Label.eq(label), gene].to_numpy() for label in (0, 1)]
                boxes = ax.boxplot(groups, tick_labels=['未损伤' if control == 'uninjured' else 'Sham', 'D7']
                                   if font else [control, 'D7'], showfliers=False, showmeans=True,
                                   meanline=True, patch_artist=True, widths=.47,
                                   medianprops={'alpha': 0},
                                   meanprops={'linestyle': '-', 'linewidth': 2, 'color': ink})
                for pos, values, color, box in zip((1, 2), groups, (blue, orange), boxes['boxes']):
                    box.set(facecolor=color, edgecolor=color, alpha=.20)
                    ax.scatter(pos + np.linspace(-.14, .14, len(values)), values, s=23, color=color,
                               edgecolor='white', linewidth=.5, zorder=3)
                values = np.concatenate(groups)
                span = np.ptp(values) or 1
                ax.set_ylim(values.min() - .12 * span, values.max() + .38 * span)
                ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
                ax.ticklabel_format(axis='y', style='plain', useOffset=False)
                ax.set_title(cohort, fontsize=9.5, pad=5)
                ax.tick_params(labelsize=8, length=3, pad=2)
                ax.set_axisbelow(True)
                ax.grid(axis='y', color='#e4e9ef', linewidth=.65)
                ax.spines[['top', 'right']].set_visible(False)
                ax.text(.97, .95, f'AUC = {auc.loc[(cohort, gene)]:.2f}', transform=ax.transAxes,
                        ha='right', va='top', fontsize=9, fontweight='bold',
                        color=orange if gene == 'Cxadr' else ink,
                        bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': 1.4})
        note = ('散点为单份样本，深色横线为组均值。各队列纵轴独立（左：SST-RMA log₂ 信号，右：log₂(CPM+1)）。'
                if font else 'Dots show samples; dark lines show means. Scales: left SST-RMA log2, right log2(CPM+1).')
        fig.text(.5, .025, note, ha='center', fontsize=10, color='#455b6b')
        fig.savefig(OUT / 'figures/bulk_roc.png', dpi=225)
        plt.close(fig)


def plot_correlation_dots(tables, name):
    """并排呈现两个外部队列的 Spearman 相关系数和 BH 校正显著性。"""
    available = {f.name for f in font_manager.fontManager.ttflist}
    font = next((f for f in ('Microsoft YaHei', 'Noto Sans CJK SC', 'PingFang SC', 'SimHei')
                 if f in available), None)
    ink, muted = '#172b3a', '#4a565f'
    cmap = LinearSegmentedColormap.from_list('correlation', ['#256080', '#f3f3f3', '#a6492e'])
    is_pathway = name == 'pathway_correlations'
    targets = tables[0].Target.drop_duplicates().tolist()
    pathway_labels = {
        'NLRP3 inflammasome': 'NLRP3 炎性小体' if font else 'NLRP3 inflammasome',
        'Extracellular matrix organization': '细胞外基质组织' if font else 'Extracellular matrix\norganization',
        'Collagen formation': '胶原形成' if font else 'Collagen formation',
    }
    labels = [pathway_labels.get(target, target) for target in targets]
    title = ('候选基因与 ssGSEA 通路活性的相关性' if is_pathway else '候选基因与经典标志物的共表达')
    if not font:
        title = 'Candidate genes and ssGSEA pathway activity' if is_pathway else 'Candidate gene and marker co-expression'
    family = ['DejaVu Sans', font] if font else ['DejaVu Sans']
    with plt.rc_context({'font.family': family, 'text.color': ink, 'xtick.color': ink, 'ytick.color': ink}):
        fig, axes = plt.subplots(1, 2, figsize=(12, 6.6), sharey=True, facecolor='white')
        fig.subplots_adjust(left=.125, right=.96, bottom=.19, top=.76, wspace=.18)
        fig.text(.125, .94, title, fontsize=20, fontweight='bold')
        fig.text(.125, .885, '外部队列内合并各采样日计算 Spearman 相关系数' if font
                 else 'Spearman correlations across all sampling days within each external cohort',
                 fontsize=11, color=muted)
        for ax, table in zip(axes, tables):
            rho = table.pivot(index='Gene', columns='Target', values='Spearman_rho').reindex(
                index=GENES, columns=targets).to_numpy()
            q = table.pivot(index='Gene', columns='Target', values='BH_q').reindex(
                index=GENES, columns=targets).to_numpy()
            rows, columns = np.where(np.isfinite(rho))
            values = rho[rows, columns]
            # 用同一色标和面积比例比较两个队列，保留零相关的小圆点。
            dots = ax.scatter(columns, rows, s=16 + 132 * np.abs(values), c=values, cmap=cmap,
                              vmin=-1, vmax=1, edgecolor='#7b8790', linewidth=.45)
            for i, j in np.ndindex(rho.shape):
                value = rho[i, j]
                stars = '**' if q[i, j] < .01 else '*' if q[i, j] < .05 else ''
                if np.isfinite(value):
                    ax.annotate(f'{value:.2f}{stars}', (j, i), xytext=(0, -7), textcoords='offset points',
                                ha='center', va='top', fontsize=9)
                else:
                    ax.text(j, i, 'NA', ha='center', va='center', fontsize=9, color=muted)
            ax.set(xlim=(-.5, len(targets) - .5), ylim=(len(GENES) - .33, -.58))
            ax.set_xticks(range(len(targets)), labels)
            ax.set_yticks(range(len(GENES)), GENES, fontstyle='italic')
            ax.xaxis.tick_top()
            ax.tick_params(axis='both', length=0, labelsize=10, pad=8)
            ax.set_title(table.Cohort.iat[0], y=1.095, fontsize=13, fontweight='bold')
            ax.spines[:].set_visible(False)
            ax.spines['top'].set_visible(True)
            ax.spines['top'].set_color('#909da6')
        cax = fig.add_axes((.125, .125, .36, .026))
        fig.colorbar(dots, cax=cax, orientation='horizontal', ticks=[-1, -.5, 0, .5, 1])
        cax.tick_params(labelsize=9)
        fig.text(.51, .138, 'Spearman ρ  |  圆点越大，相关越强' if font
                 else 'Spearman ρ  |  Dot area increases with |ρ|', fontsize=10, va='center')
        fig.text(.125, .038, '* BH q < 0.05，** BH q < 0.01' if font
                 else '* BH q < 0.05, ** BH q < 0.01', fontsize=10, color=muted)
        fig.savefig(OUT / 'figures' / f'{name}.png', dpi=240)
        plt.close(fig)


def main():
    for directory in ('tables', 'figures'):
        (OUT / directory).mkdir(parents=True, exist_ok=True)
    frozen = read(FROZEN)
    assert frozen['基因符号'].tolist() == GENES
    coverage = read(PREP / 'tables/candidate_coverage.tsv')
    assert frozen['基因ID'].tolist() == coverage.EnsemblID.tolist()
    with closing(sqlite3.connect(f"{(ROOT / '.tools/org.Mm.eg.sqlite').as_uri()}?mode=ro", uri=True)) as db:
        mapping = pd.read_sql_query(
            'SELECT gene_id AS entrez, symbol FROM genes JOIN gene_info USING(_id)', db)
    membership = read(ROOT / 'experiments/02_feature_selection/ssgsea_pathway_membership.tsv')
    membership = membership[membership.Pathway.isin(['NLRP3 inflammasome', 'Extracellular matrix organization', 'Collagen formation'])]
    assert membership.Pathway.nunique() == 3
    data, tests, mapping_rows = {}, [], []
    for cohort, control in COHORTS:
        full, sliced, meta, scale, labels = load_cohort(cohort, coverage, mapping)
        data[cohort] = (full, sliced, meta, scale)
        mapping_rows.extend(labels)
        tests.append(expression_tests(cohort, control, sliced, meta))
    combined = pd.concat(tests, ignore_index=True)
    save(combined, 'bulk_expression_tests')
    save(pd.DataFrame(mapping_rows), 'bulk_gene_mapping')
    passed = combined[combined.Primary].groupby('Gene').Pass.all()
    selected = [gene for gene in GENES if passed.loc[gene]]
    print('D7 gate:', selected, flush=True)
    roc, results, marker, pathway, align = [], {}, [], [], []
    for cohort, control in COHORTS:
        full, x, meta, scale = data[cohort]
        boxplot(cohort, x, meta, scale, combined[combined.Cohort.eq(cohort)], control)
        predictions, metrics = roc_scores(cohort, control, x, meta, selected)
        results[cohort] = predictions
        roc.append(metrics)
        marker.append(correlations(cohort, x, x.loc[MARKERS]))
        score, availability = pathway_scores(cohort, full, mapping, membership)
        align.append(availability)
        pathway.append(correlations(cohort, x, score))
    roc = pd.concat(roc, ignore_index=True)
    save(roc, 'bulk_roc_metrics')
    save(pd.concat(marker, ignore_index=True), 'marker_correlations')
    save(pd.concat(pathway, ignore_index=True), 'pathway_correlations')
    save(pd.concat(align, ignore_index=True), 'pathway_coverage')
    plot_expression_auc(results, roc)
    plot_correlation_dots(marker, 'marker_correlations')
    plot_correlation_dots(pathway, 'pathway_correlations')


if __name__ == '__main__':
    main()
