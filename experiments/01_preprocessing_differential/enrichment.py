"""对五个时间对比进行 GO、KEGG 过度富集检验和 KEGG GSEA。"""
from collections import defaultdict
from contextlib import closing
import json
import sqlite3

import gseapy
import numpy as np
import pandas as pd
import scipy
from scipy.stats import false_discovery_control, hypergeom

from paths import ROOT, REPORTS, TABLES


def annotation_sets(genes):
    # 将主队列稳定 Ensembl ID 对应到小鼠 Entrez ID，双向歧义均排除。
    with closing(sqlite3.connect((ROOT / '.tools/org.Mm.eg.sqlite').as_uri() + '?mode=ro', uri=True)) as db:
        mapping = pd.read_sql_query('SELECT DISTINCT ensembl_id, gene_id FROM ensembl JOIN genes USING (_id)', db)
        mapping = mapping[mapping.ensembl_id.isin(genes)]
        unique = ~mapping.ensembl_id.duplicated(keep=False) & ~mapping.gene_id.duplicated(keep=False)
        mapping = mapping.loc[unique].copy()
        db.execute('CREATE TEMP TABLE included (gene_id TEXT PRIMARY KEY)')
        db.executemany('INSERT INTO included VALUES (?)', ((gene,) for gene in mapping.gene_id))
        go_sets = {}
        for ontology in ('BP', 'CC', 'MF'):
            sets = defaultdict(set)
            query = f'SELECT g.gene_id, a.go_id FROM go_{ontology.lower()}_all a JOIN genes g USING (_id) JOIN included i ON g.gene_id=i.gene_id'
            for gene, term in db.execute(query):
                sets[term].add(gene)
            go_sets[ontology] = dict(sets)
        metadata = dict(db.execute('SELECT name, value FROM metadata'))
    with closing(sqlite3.connect((ROOT / '.tools/GO.sqlite').as_uri() + '?mode=ro', uri=True)) as db:
        go_titles = dict(db.execute('SELECT go_id, term FROM go_term'))
        go_metadata = dict(db.execute('SELECT name, value FROM metadata'))
    for ontology, sets in go_sets.items():
        go_sets[ontology] = {term: members for term, members in sets.items() if term in go_titles}
    return mapping, go_sets, go_titles, {'小鼠注释': metadata, 'GO注释': go_metadata}


def kegg_sets():
    source = ROOT / 'datasets/annotation/kegg_mmu'
    edge = pd.read_csv(source / 'pathway_genes.tsv', sep='\t', names=['gene', 'term'], dtype=str)
    title = pd.read_csv(source / 'pathway_names.tsv', sep='\t', names=['term', 'name'], dtype=str)
    edge['gene'] = edge.gene.str.removeprefix('mmu:')
    edge['term'] = edge.term.str.removeprefix('path:')
    title['term'] = title.term.str.removeprefix('path:')
    title['name'] = title.name.str.split(' - Mus musculus').str[0]
    sets, names = edge.groupby('term').gene.agg(set).to_dict(), title.set_index('term').name.to_dict()
    if not len(sets) > 100 or not sets.keys() <= names.keys():
        raise ValueError('KEGG 小鼠通路名称与基因映射不一致。')
    return sets, names


def overrepresentation(sets, names, background, selected, contrast, direction, collection):
    # 各注释集分别定义背景，保留零命中的条目参与多重检验校正。
    universe = background & set().union(*sets.values())
    selected = selected & universe
    rows = []
    if not selected:
        return pd.DataFrame()
    for term, members in sets.items():
        size = len(members & universe)
        if not 10 <= size <= 500:
            continue
        overlap = members & selected
        rows.append({'比较': contrast, '方向': direction, '注释集': collection, '条目ID': term,
                     '条目原名': names[term], '命中数': len(overlap), '候选数': len(selected),
                     '背景条目大小': size, '背景基因数': len(universe),
                     '富集倍数': (len(overlap) / len(selected)) / (size / len(universe)),
                     '命中基因EntrezID': ','.join(sorted(overlap))})
    result = pd.DataFrame(rows)
    if not result.empty:
        p = hypergeom.sf(result['命中数'] - 1, result['背景基因数'], result['背景条目大小'], result['候选数'])
        result['P值'] = p
        result['校正后P值'] = false_discovery_control(p, method='bh')
    return result


def prerank(part, sets, names, contrast, seed):
    # 对所有具备有限 Wald 统计量与 P 值的映射基因排序。
    rank = part.dropna(subset=['EntrezID', 'Wald统计量', 'P值']).copy()
    rank = rank[np.isfinite(rank['Wald统计量']) & np.isfinite(rank['P值'])]
    rank = rank.sort_values(['Wald统计量', 'EntrezID'], ascending=[False, True])
    rank = rank.set_index('EntrezID')['Wald统计量']
    if not rank.index.is_unique or len(rank) < 1000:
        raise ValueError(f'{contrast} 的排序基因不足或存在重复。')
    analysis = gseapy.prerank(rnk=rank, gene_sets={k: sorted(v) for k, v in sets.items()},
                              min_size=15, max_size=500, permutation_num=1000, weight=1,
                              threads=2, ascending=None, method='multilevel', sample_size=101,
                              eps=1e-10, seed=seed, outdir=None, no_plot=True)
    rows = analysis.res2d.rename(columns={'Term': '条目ID', 'ES': '富集分数ES', 'NES': '标准化富集分数NES',
                                           'NOM p-val': 'P值', 'Lead_genes': '核心基因EntrezID'}).copy()
    rows['比较'] = contrast
    rows['条目原名'] = rows['条目ID'].map(names)
    rows['排序基因数'] = len(rank)
    rank_ids = set(rank.index)
    rows['集合大小'] = rows['条目ID'].map(lambda term: len(sets[term] & rank_ids))
    rows['P值'] = pd.to_numeric(rows['P值'], errors='coerce')
    valid = rows['P值'].notna()
    rows['校正后P值'] = np.nan
    rows.loc[valid, '校正后P值'] = false_discovery_control(rows.loc[valid, 'P值'], method='bh')
    return rows[['比较', '条目ID', '条目原名', '集合大小', '排序基因数', '富集分数ES', '标准化富集分数NES',
                 'P值', '校正后P值', '核心基因EntrezID']]


def main():
    wald = pd.read_csv(TABLES / 'timecourse_wald.tsv.gz', sep='\t')
    mapping, go, go_names, metadata = annotation_sets(set(wald['基因ID']))
    kegg, kegg_names = kegg_sets()
    ids = mapping.set_index('ensembl_id').gene_id
    mapped = mapping.rename(columns={'ensembl_id': '基因ID', 'gene_id': 'EntrezID'})
    symbols = wald.drop_duplicates('基因ID').set_index('基因ID')['基因符号']
    mapped['基因符号'] = mapped['基因ID'].map(symbols)
    mapped.to_csv(TABLES / 'enrichment_gene_mapping.tsv', sep='\t', index=False, encoding='utf-8-sig')
    wald['EntrezID'] = wald['基因ID'].map(ids)
    outputs = {'GO_ORA': [], 'KEGG_ORA': [], 'KEGG_GSEA': []}
    for i, (contrast, part) in enumerate(wald.groupby('比较', sort=False)):
        print(f'富集分析 {contrast}', flush=True)
        background = set(part.loc[part['校正后P值'].notna(), 'EntrezID'].dropna())
        for direction, sign in [('上调', 1), ('下调', -1)]:
            mask = part['通过差异筛选'].eq('是') & (sign * part['log2表达倍数变化'] > 0)
            selected = set(part.loc[mask, 'EntrezID'].dropna())
            for ontology, sets in go.items():
                outputs['GO_ORA'].append(overrepresentation(sets, go_names, background, selected, contrast, direction, ontology))
            outputs['KEGG_ORA'].append(overrepresentation(kegg, kegg_names, background, selected, contrast, direction, 'KEGG'))
        outputs['KEGG_GSEA'].append(prerank(part, kegg, kegg_names, contrast, 252028 + i))
    summary = []
    for kind, frames in outputs.items():
        table = pd.concat(frames, ignore_index=True).sort_values(['比较', '校正后P值', '条目ID'])
        table.to_csv(TABLES / f'enrichment_{kind.lower()}.tsv.gz', sep='\t', index=False, encoding='utf-8-sig')
        for contrast, rows in table.groupby('比较', sort=False):
            summary.append({'分析': kind, '比较': contrast, '检验条目数': len(rows),
                            '显著条目数': int(rows['校正后P值'].lt(0.05).sum())})
    summary = pd.DataFrame(summary)
    summary.to_csv(TABLES / 'enrichment_summary.tsv', sep='\t', index=False, encoding='utf-8-sig')
    report = {'差异模型基因数': len(symbols), '一对一映射基因数': len(mapped), 'gseapy': gseapy.__version__,
              'scipy': scipy.__version__, '注释来源': metadata,
              'GSEA方法': 'multilevel', 'GSEA集合大小范围': [15, 500], 'GSEA基因集': 'KEGG mmu',
              'ORA集合大小范围': [10, 500], '校正方法': '各检验家族内 Benjamini–Hochberg'}
    (REPORTS / 'enrichment_python_session.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
