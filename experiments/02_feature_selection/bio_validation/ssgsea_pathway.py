"""发现队列前八基因的ssGSEA 分析及冻结九候选的通路互证。

读取阶段一 VST、GSM 纳入表、基因 ID 映射、并集排序和固定通路成员。
导出 NES、Spearman 相关、BH 校正和敏感性源表，以及相应热图。
"""

import json
from textwrap import wrap

import gseapy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, rankdata

from ..common.data_loader import read_main_days, read_ranking, read_vst
from ..common.plot_style import PATHWAY_LABELS, draw_heatmap
from ..config import (
    CORE_GENES,
    CORRELATION_PERMUTATIONS,
    TOP_PATHWAY_GENES,
    MAIN_LIBRARY_COUNTS,
    PATHWAY_MEMBERSHIP,
    PATHWAY_SOURCES,
    SEED,
    SSGSEA_MIN_SIZE,
    SSGSEA_MAX_SIZE,
    SSGSEA_WEIGHT,
    STAGE_ONE,
)


def score_pathways(vst, gene_sets):
    """从基因 × 文库 VST 和基因集计算 NES，返回文库 × 通路矩阵。"""
    result = gseapy.ssgsea(
        data=vst,
        gene_sets=gene_sets,
        organism='mouse',
        outdir=None,
        sample_norm_method='rank',
        correl_norm_type='rank',
        min_size=SSGSEA_MIN_SIZE,
        max_size=SSGSEA_MAX_SIZE,
        permutation_num=0,
        weight=SSGSEA_WEIGHT,
        threads=1,
        seed=SEED,
        no_plot=True,
    )
    scores = result.res2d.pivot(index='Name', columns='Term', values='NES')
    return scores.loc[vst.columns, list(gene_sets)].astype(float)


def standardized_ranks(values):
    """按列转换为中心化、单位范数的秩，返回同形状数组。"""
    ranks = rankdata(values, axis=0)
    ranks = ranks - ranks.mean(axis=0)
    norms = np.linalg.norm(ranks, axis=0)
    if (norms == 0).any():
        raise ValueError('核心基因表达或通路得分恒定，无法计算 Spearman 相关')
    return ranks / norms


def spearman_permutation(expression, scores, permutations=CORRELATION_PERMUTATIONS):
    """Two-sided marginal tests with shared random permutations of the library pairing."""
    if not expression.index.equals(scores.index):
        raise ValueError('核心基因表达与 ssGSEA 得分的文库顺序不一致')
    gene = standardized_ranks(expression.to_numpy())
    pathway = standardized_ranks(scores.to_numpy())
    observed = gene.T @ pathway
    rng = np.random.default_rng(SEED)
    order = rng.permuted(np.tile(np.arange(len(expression)), (permutations, 1)), axis=1)
    null = np.einsum('pna,nb->pab', gene[order], pathway)
    p = (1 + (np.abs(null) >= np.abs(observed) - 1e-12).sum(axis=0)) / (permutations + 1)
    return observed, p


def get_expression(vst, annotation, symbols):
    """Resolve requested symbols against uniquely quantified Ensembl IDs."""
    relevant = annotation.loc[annotation['基因ID'].isin(vst.index)].set_index('基因符号')['基因ID']
    ids = relevant.loc[relevant.index.isin(symbols)].reindex(symbols).rename('gene_id')
    values = vst.reindex(index=ids).set_axis(symbols).T
    return ids, values


def correlate(expression, features):
    """Retain missing features in the output and adjust only the estimable tests."""
    rows, cols = len(expression.columns), len(features.columns)
    rho, p, q = (np.full((rows, cols), np.nan) for _ in range(3))
    valid = features.notna().all() & features.nunique().gt(1)
    if valid.any():
        rho[:, valid], p[:, valid] = spearman_permutation(expression, features.loc[:, valid])
        q[:, valid] = false_discovery_control(p[:, valid].ravel()).reshape(rows, int(valid.sum()))
    return rho, p, q


def as_rows(genes, features, statistics):
    """将 ρ、P 和 q 数组展开为逐基因–特征长表。"""
    columns = pd.MultiIndex.from_product([genes, features], names=['gene', 'feature'])
    return pd.DataFrame(
        {name: values.ravel() for name, values in zip(('rho', 'p', 'q'), statistics)}, index=columns
    ).reset_index()


def plot_correlations(genes, pathways, rho, p, figure_path):
    """绘制前八基因通路热图，标注双侧置换原始 P 值。"""
    fig, ax = plt.subplots(figsize=(10.6, 6.0))
    fig.subplots_adjust(left=0.135, right=0.90, bottom=0.13, top=0.735)
    fig.text(0.135, 0.95, 'Core gene expression and pathway enrichment', fontsize=17, fontweight='bold')
    image = ax.imshow(rho, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    labels = ['\n'.join(wrap(name, width=15, break_long_words=False)) for name in pathways['Pathway']]
    labels = [
        label.replace('Cytokine–receptor\ninteraction', 'Cytokine–\nreceptor interaction')
        for label in labels
    ]
    ax.set_xticks(range(len(pathways)), labels)
    ax.set_yticks(range(len(genes)), genes)
    ax.xaxis.tick_top()
    ax.tick_params(axis='both', length=0, labelsize=9, pad=8)
    for tick in ax.get_yticklabels():
        tick.set_fontstyle('italic')
    for group, color in [('Inflammation', '#35566B'), ('Fibrosis', '#855743')]:
        columns = np.flatnonzero(pathways['Group'].eq(group))
        left, right = columns[0] - 0.44, columns[-1] + 0.44
        ax.plot(
            [left, right],
            [1.18, 1.18],
            transform=ax.get_xaxis_transform(),
            color=color,
            linewidth=2,
            clip_on=False,
        )
        ax.text(
            (left + right) / 2,
            1.205,
            group,
            transform=ax.get_xaxis_transform(),
            ha='center',
            va='bottom',
            color=color,
            fontsize=10,
            fontweight='bold',
        )
    ax.axvline(pathways['Group'].eq('Inflammation').sum() - 0.5, color='white', linewidth=2)
    for row, column in np.ndindex(rho.shape):
        stars = '**' if p[row, column] < 0.01 else '*' if p[row, column] < 0.05 else ''
        value = rho[row, column]
        ax.text(
            column,
            row,
            f'{value:.2f}{stars}',
            ha='center',
            va='center',
            fontsize=9,
            color='white' if abs(value) >= 0.65 else '#172B3A',
        )
    ax.spines[:].set_visible(False)
    cax = fig.add_axes((0.925, 0.24, 0.012, 0.32))
    fig.colorbar(image, cax=cax, ticks=[-1, -0.5, 0, 0.5, 1], label='Spearman ρ')
    fig.text(
        0.135,
        0.065,
        '* P < 0.05,  ** P < 0.01  |  Two-sided permutation P, before BH correction',
        fontsize=9,
        color='#4A565F',
    )
    fig.savefig(figure_path, dpi=240)
    plt.close(fig)


def run_top_candidates(selection_dir, output):
    """为并集排名前八的候选生成ssGSEA 表、敏感性结果与热图。"""
    ranking = read_ranking(selection_dir)
    genes = ranking.nsmallest(TOP_PATHWAY_GENES, '排名')['基因符号']
    days = read_main_days()
    vst = read_vst()
    if days.value_counts().to_dict() != MAIN_LIBRARY_COUNTS or set(days.index) != set(vst.columns):
        raise ValueError('VST 矩阵与阶段 01 的 22 份主分析文库不一致')
    vst = vst.loc[:, days.index]

    membership = pd.read_csv(PATHWAY_MEMBERSHIP, sep='\t')
    pathways = (
        membership[['PathwayID', 'Pathway', 'Group', 'Database']]
        .drop_duplicates()
        .set_index('PathwayID')
    )
    mapping = pd.read_csv(
        STAGE_ONE / 'tables/enrichment_gene_mapping.tsv', sep='\t', usecols=['基因ID', 'EntrezID']
    ).rename(columns={'基因ID': 'EnsemblID'})
    members = membership.merge(mapping, on='EntrezID', how='left', validate='many_to_one')
    members['In_VST'] = members['EnsemblID'].isin(vst.index)
    members['Core_gene'] = members['EnsemblID'].map(genes)
    gene_sets = {
        pathway: sorted(members.loc[members['PathwayID'].eq(pathway) & members['In_VST'], 'EnsemblID'])
        for pathway in pathways.index
    }
    coverage = (
        members.groupby('PathwayID')['In_VST']
        .agg(['size', 'sum'])
        .rename(columns={'size': 'Source_genes', 'sum': 'Genes_in_VST'})
    )
    scores = score_pathways(vst, gene_sets)
    expression = vst.loc[genes.index].T
    rho, p = spearman_permutation(expression, scores)
    q = false_discovery_control(p.ravel()).reshape(p.shape)

    # Re-score after removing the eight queried genes from the ranking universe.
    without_core = {
        pathway: sorted(set(identifiers) - set(genes.index))
        for pathway, identifiers in gene_sets.items()
    }
    sensitivity = score_pathways(vst.drop(index=genes.index), without_core)
    sensitivity_rho = standardized_ranks(expression.to_numpy()).T @ standardized_ranks(
        sensitivity.to_numpy()
    )
    rows = []
    for i, (identifier, symbol) in enumerate(genes.items()):
        for j, (pathway, annotation) in enumerate(pathways.iterrows()):
            rows.append(
                (
                    identifier,
                    symbol,
                    pathway,
                    annotation['Pathway'],
                    annotation['Group'],
                    rho[i, j],
                    p[i, j],
                    q[i, j],
                    identifier in gene_sets[pathway],
                    sensitivity_rho[i, j],
                )
            )
    correlations = pd.DataFrame(
        rows,
        columns=[
            'EnsemblID',
            'Gene',
            'PathwayID',
            'Pathway',
            'Group',
            'Spearman_rho',
            'P',
            'q_BH',
            'Gene_in_pathway',
            'Rho_without_core',
        ],
    )
    tables, figures, reports = (output / name for name in ('tables', 'figures', 'reports'))
    for directory in (tables, figures, reports):
        directory.mkdir(parents=True, exist_ok=True)
    matrices = {
        'ssgsea_scores': scores.assign(Day=days)[['Day', *pathways.index]].rename_axis('Library'),
        'ssgsea_spearman_rho': pd.DataFrame(
            rho, index=pd.Index(genes, name='Gene'), columns=pathways.index
        ),
        'ssgsea_pathway_coverage': pathways.join(coverage),
    }
    for name, matrix in matrices.items():
        matrix.to_csv(tables / f'{name}.tsv', sep='\t', encoding='utf-8-sig', float_format='%.10g')
    correlations.to_csv(
        tables / 'ssgsea_correlations.tsv',
        sep='\t',
        index=False,
        encoding='utf-8-sig',
        float_format='%.10g',
    )
    summary = {
        'GSEApy': gseapy.__version__,
        'VST_shape': list(vst.shape),
        'Core_genes': genes.to_dict(),
        'library_counts': days.value_counts().sort_index().to_dict(),
        'score': 'ssGSEA NES, rank/rank, weight 0.25, permutation_num 0',
        'gene_set_sizes': [SSGSEA_MIN_SIZE, SSGSEA_MAX_SIZE],
        'permutations': CORRELATION_PERMUTATIONS,
        'seed': SEED,
        'P_value': 'Two-sided Spearman pairing permutation, all 22 libraries pooled',
        'BH_tests': int(p.size),
        'pathway_sources': json.loads(PATHWAY_SOURCES.read_text(encoding='utf-8')),
    }
    (reports / 'ssgsea_summary.json').write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8'
    )
    plot_correlations(genes, pathways, rho, p, figures / 'ssgsea_top_candidates.png')
    print(
        f'{len(days)} libraries × {len(pathways)} pathways, {p.size} gene–pathway correlations: {output}'
    )


def pathway_data(vst, expression, ids, days, tables):
    """计算冻结九候选的七通路 NES 和关联，导出全量与去候选敏感性结果。"""
    membership = pd.read_csv(PATHWAY_MEMBERSHIP, sep='\t')
    pathways = membership[['PathwayID', 'Pathway', 'Group']].drop_duplicates().set_index('PathwayID')
    pathways = pathways.loc[list(PATHWAY_LABELS)]
    mapping = pd.read_csv(
        STAGE_ONE / 'tables/enrichment_gene_mapping.tsv', sep='\t', usecols=['基因ID', 'EntrezID']
    )
    members = membership.merge(mapping, on='EntrezID', how='left', validate='many_to_one')
    matched = members.loc[members['基因ID'].isin(vst.index)]
    gene_sets = {
        pathway: sorted(set(matched.loc[matched['PathwayID'].eq(pathway), '基因ID']))
        for pathway in pathways.index
    }
    reduced_sets = {pathway: sorted(set(genes) - set(ids)) for pathway, genes in gene_sets.items()}
    scores = score_pathways(vst, gene_sets)
    reduced_scores = score_pathways(vst.drop(index=ids), reduced_sets)
    names = list(expression.columns)
    comparisons = as_rows(names, scores.columns, correlate(expression, scores))
    sensitivity = as_rows(names, scores.columns, correlate(expression, reduced_scores))
    result = comparisons.merge(sensitivity, on=['gene', 'feature'], suffixes=('', '_without_core'))
    result['gene_id'] = result['gene'].map(ids)
    result = result.assign(
        pathway=result['feature'].map(pathways['Pathway']),
        group=result['feature'].map(pathways['Group']),
        genes_in_pathway=result['feature'].map({key: len(value) for key, value in gene_sets.items()}),
        genes_in_pathway_without_core=result['feature'].map(
            {key: len(value) for key, value in reduced_sets.items()}
        ),
        candidate_in_pathway=[
            gene_id in gene_sets[pathway]
            for gene_id, pathway in zip(result['gene_id'], result['feature'])
        ],
    )
    result = result.rename(columns={'feature': 'pathway_id'})
    columns = [
        'gene',
        'gene_id',
        'pathway_id',
        'pathway',
        'group',
        'rho',
        'p',
        'q',
        'rho_without_core',
        'p_without_core',
        'q_without_core',
        'genes_in_pathway',
        'genes_in_pathway_without_core',
        'candidate_in_pathway',
    ]
    result[columns].to_csv(
        tables / 'candidate_pathway_correlations.tsv', sep='\t', index=False, float_format='%.12g'
    )
    scores.insert(0, 'day', days)
    scores.join(reduced_scores.add_suffix('_without_core')).to_csv(
        tables / 'candidate_pathway_scores.tsv',
        sep='\t',
        index_label='library',
        float_format='%.12g',
    )
    return result


def load_candidate_expression():
    """返回按样本表对齐的 VST、九候选表达、候选 ID、采样日和注释。"""
    days = read_main_days()
    vst = read_vst()
    if days.value_counts().to_dict() != MAIN_LIBRARY_COUNTS or set(days.index) != set(vst.columns):
        raise ValueError('Retained library metadata do not match the Stage 01 VST matrix')

    # VST 以基因为行；get_expression 将候选矩阵转置成文库 × 候选。
    vst = vst.loc[:, days.index]
    annotation = pd.read_csv(STAGE_ONE / 'matrices/gene_annotation.tsv.gz', sep='\t')
    ids, expression = get_expression(vst, annotation, list(CORE_GENES))
    if ids.isna().any() or expression.nunique().le(1).any():
        raise ValueError('A shortlisted gene is missing or invariant across the retained libraries')
    return vst, expression, ids, days, annotation


def run_frozen_candidates(tables):
    """为冻结九候选运行通路分析，返回热图的长表。"""
    vst, expression, ids, days, _ = load_candidate_expression()
    tables.mkdir(parents=True, exist_ok=True)
    return pathway_data(vst, expression, ids, days, tables)


def draw_pathways(ax, results):
    """为共表达分析图的通路面板绘图，并返回可供共享色标使用的图像。"""
    return draw_heatmap(ax, results, 'pathway_id', tuple(PATHWAY_LABELS.values()), 4)
