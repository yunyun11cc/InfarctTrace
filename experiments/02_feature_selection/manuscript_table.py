"""Assemble the nine nominated genes and their discovery-cohort evidence."""

import argparse
from pathlib import Path

import pandas as pd

from manuscript_style import CORE_GENES, EXPORT_DIR
from select_features import STAGE_ONE, UNION_OUTPUT


def write_table(output):
    source = output / 'source_data'
    frequencies = pd.read_csv(source / 'figure_2c_candidate_frequencies.tsv', sep='\t',
                             index_col='gene', dtype={'gene_id': str}).loc[list(CORE_GENES)]
    ranking = pd.read_csv(UNION_OUTPUT / 'tables/feature_ranking.tsv', sep='\t', index_col='基因ID')
    lrt = pd.read_csv(STAGE_ONE / 'tables/timecourse_lrt.tsv.gz', sep='\t', index_col='基因ID')
    pathways = pd.read_csv(source / 'figure_3a_pathway_correlations.tsv', sep='\t')
    if len(pathways) != 63 or pathways.groupby('gene_id')['pathway'].nunique().ne(7).any():
        raise ValueError('Expected seven pathways per candidate for all nine candidates')
    associations = pathways.loc[pathways['rho'].gt(0) & pathways['q'].lt(0.05)]
    groups = associations.groupby('gene_id')['group'].agg(frozenset)
    pathway_labels = {
        frozenset(): '未达阈值',
        frozenset(['Inflammation']): '炎症',
        frozenset(['Fibrosis']): '纤维化',
        frozenset(['Inflammation', 'Fibrosis']): '双通路',
    }
    ids = frequencies['gene_id']
    result = pd.DataFrame({
        '基因符号': frequencies.index,
        '基因ID': ids.to_numpy(),
        'D7 - 非D7 均值差': ranking.loc[ids, 'D7减非D7_VST均值差'].to_numpy(),
        'LRT 校正后 P 值': lrt.loc[ids, '校正后P值'].to_numpy(),
        'VST 均值差 Bootstrap 选中率': frequencies['均值差'].to_numpy(),
        'LASSO Bootstrap 选中率': frequencies['LASSO'].to_numpy(),
        'EN Bootstrap 选中率': frequencies['ElasticNet'].to_numpy(),
        'RF Bootstrap 选中率': frequencies['随机森林'].to_numpy(),
        'EFS-RA 综合秩聚合排名': frequencies['efs_ra_rank'].to_numpy(),
        '主要关联通路': [pathway_labels[groups.get(gene_id, frozenset())] for gene_id in ids],
    })
    if result.isna().any().any() or not ids.is_unique or set(ids) != set(pathways['gene_id']):
        raise ValueError('Candidate table requires nine unique and fully annotated Ensembl IDs')
    result.to_csv(output / 'table_02_candidate_genes_summary.tsv', sep='\t', index=False, float_format='%.10g')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=EXPORT_DIR)
    write_table(parser.parse_args().output_dir)
