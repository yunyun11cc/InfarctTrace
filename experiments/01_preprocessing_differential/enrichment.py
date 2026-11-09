"""对五个时间对比进行 GO、KEGG 过度富集检验和 KEGG GSEA。

两种富集回答不同问题，使用不同的基因名单：

- **过度富集（ORA）**：先按 ``|log2FC| ≥ 1 且 padj < 0.05`` 筛出差异基因，再看
  这些基因里某个 GO/KEGG 条目的命中比例是否高于背景（超几何检验）。
- **基因集富集（GSEA）**：不预先筛基因，而是把**所有**有有效 Wald 统计量的
  基因按带符号的统计量排序，看某个通路成员是否集中在排序的一端。

注释来源固定为本次快照：Bioconductor 3.22 的 ``org.Mm.eg.db``、``GO.db``
（由 ``go_reference.py`` 下载并校验）和 KEGG ``mmu`` 通路文件（由
``kegg_reference.py`` 下载并校验）。
"""

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


# ---------------------------------------------------------------------------
# 1. 注释基因集
# ---------------------------------------------------------------------------

def annotation_sets(genes):
    """构建 GO 基因集并把 Ensembl ID 一对一映射到 Entrez ID。

    ``org.Mm.eg.db`` 与 ``GO.db`` 都是 SQLite 数据库，这里用只读 URI 打开，
    避免分析过程意外改动注释文件。

    映射规则：任一 Ensembl ID 或 Entrez ID 在两侧出现多次都整体排除，只保留
    双向唯一的对应关系，防止一个基因的计数被重复摊到多个条目上。

    返回
    ----
    mapping : DataFrame
        列 ``ensembl_id``、``gene_id``（Entrez），只含本次检验基因。
    go_sets : dict
        ``{'BP': {term: {entrez, ...}}, 'CC': ..., 'MF': ...}``。
    go_titles : dict
        GO 条目 ID → 名称。
    metadata : dict
        两个数据库的版本/构建信息，写入运行报告。
    """
    # 将主队列稳定 Ensembl ID 对应到小鼠 Entrez ID，双向歧义均排除。
    with closing(sqlite3.connect((ROOT / '.tools/org.Mm.eg.sqlite').as_uri() + '?mode=ro', uri=True)) as db:
        mapping = pd.read_sql_query('SELECT DISTINCT ensembl_id, gene_id FROM ensembl JOIN genes USING (_id)', db)

        # 只保留本次进入差异模型的基因
        mapping = mapping[mapping.ensembl_id.isin(genes)]

        # 双向唯一：没有任何 ID 在映射表中出现两次
        unique = ~mapping.ensembl_id.duplicated(keep=False) & ~mapping.gene_id.duplicated(keep=False)
        mapping = mapping.loc[unique].copy()

        # 用临时表把允许的 Entrez ID 传入后续查询，避免拼接超长 IN 列表
        db.execute('CREATE TEMP TABLE included (gene_id TEXT PRIMARY KEY)')
        db.executemany('INSERT INTO included VALUES (?)', ((gene,) for gene in mapping.gene_id))

        # GO 按三个本体分开收集基因集
        go_sets = {}
        for ontology in ('BP', 'CC', 'MF'):
            sets = defaultdict(set)
            query = f'SELECT g.gene_id, a.go_id FROM go_{ontology.lower()}_all a JOIN genes g USING (_id) JOIN included i ON g.gene_id=i.gene_id'
            for gene, term in db.execute(query):
                sets[term].add(gene)
            go_sets[ontology] = dict(sets)

        metadata = dict(db.execute('SELECT name, value FROM metadata'))

    # GO 条目名称存放在另一个数据库
    with closing(sqlite3.connect((ROOT / '.tools/GO.sqlite').as_uri() + '?mode=ro', uri=True)) as db:
        go_titles = dict(db.execute('SELECT go_id, term FROM go_term'))
        go_metadata = dict(db.execute('SELECT name, value FROM metadata'))

    # 丢弃没有名称的条目（通常是已废弃编号），保证写表时都有可读标签
    for ontology, sets in go_sets.items():
        go_sets[ontology] = {term: members for term, members in sets.items() if term in go_titles}

    return mapping, go_sets, go_titles, {'小鼠注释': metadata, 'GO注释': go_metadata}


def kegg_sets():
    """读取 KEGG 小鼠通路快照，返回 ``(通路基因集, 通路名称)``。

    ``pathway_genes.tsv`` 每行是 ``mmu:基因ID`` 与 ``path:mmuXXXXX``；
    ``pathway_names.tsv`` 每行是通路 ID 与名称。KEGG 的名称带
    `` - Mus musculus`` 后缀，这里去掉以便中文图注显示。
    """
    source = ROOT / 'datasets/annotation/kegg_mmu'

    edge = pd.read_csv(source / 'pathway_genes.tsv', sep='\t', names=['gene', 'term'], dtype=str)
    title = pd.read_csv(source / 'pathway_names.tsv', sep='\t', names=['term', 'name'], dtype=str)

    # 去掉 KEGG 的 mmu: 与 path: 前缀，只留编号
    edge['gene'] = edge.gene.str.removeprefix('mmu:')
    edge['term'] = edge.term.str.removeprefix('path:')
    title['term'] = title.term.str.removeprefix('path:')
    title['name'] = title.name.str.split(' - Mus musculus').str[0]

    # 按通路聚合成基因集合，并建立 ID → 名称的字典
    sets, names = edge.groupby('term').gene.agg(set).to_dict(), title.set_index('term').name.to_dict()

    # 完整性检查：通路数量应达到常规规模，且每条通路都有名称
    if not len(sets) > 100 or not sets.keys() <= names.keys():
        raise ValueError('KEGG 小鼠通路名称与基因映射不一致。')

    return sets, names


# ---------------------------------------------------------------------------
# 2. 过度富集检验（ORA）
# ---------------------------------------------------------------------------

def overrepresentation(sets, names, background, selected, contrast, direction, collection):
    """对一个基因集集合做超几何右尾检验。

    记号：N 为背景基因数，K 为背景中属于该条目的基因数，n 为候选（差异）基因数，
    k 为候选中命中的基因数。零假设下 X ~ 超几何(N, K, n)，P 值为 P(X ≥ k)。

    ``background`` 是该对比里所有可检验的 Entrez ID；``selected`` 是其中达到
    差异门槛且方向一致的部分。仅检验背景内大小在 10–500 的条目，但**零命中
    的条目仍保留**，以便它们进入同一 BH 家族。
    """
    # 各注释集分别定义背景，保留零命中的条目参与多重检验校正。
    universe = background & set().union(*sets.values())
    selected = selected & universe

    rows = []
    if not selected:
        return pd.DataFrame()

    for term, members in sets.items():
        size = len(members & universe)

        # 过小的条目统计不稳定，过大的条目接近全基因组，均不检验
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
        # surf 的右尾从 k-1 开始累计，等价于 P(X ≥ k)
        p = hypergeom.sf(result['命中数'] - 1, result['背景基因数'], result['背景条目大小'], result['候选数'])
        result['P值'] = p
        result['校正后P值'] = false_discovery_control(p, method='bh')

    return result


# ---------------------------------------------------------------------------
# 3. 基因集富集分析（GSEA）
# ---------------------------------------------------------------------------

def prerank(part, sets, names, contrast, seed):
    """用带符号的 Wald 统计量排序，跑 GSEApy 的 prerank 分析。

    排序只用“能映射到 Entrez 且有有限统计量”的基因，因此每个对比的排序
    基因数略少于进入模型的基因数。固定种子保证同一条命令可复现。
    """
    # 对所有具备有限 Wald 统计量与 P 值的映射基因排序。
    rank = part.dropna(subset=['EntrezID', 'Wald统计量', 'P值']).copy()
    rank = rank[np.isfinite(rank['Wald统计量']) & np.isfinite(rank['P值'])]

    # 主键是统计量降序；并列时用 Entrez ID 升序，避免排序不稳定
    rank = rank.sort_values(['Wald统计量', 'EntrezID'], ascending=[False, True])
    rank = rank.set_index('EntrezID')['Wald统计量']

    # 至少 1000 个基因才做 GSEA，且 ID 必须唯一
    if not rank.index.is_unique or len(rank) < 1000:
        raise ValueError(f'{contrast} 的排序基因不足或存在重复。')

    # multilevel 用少量置换估计极小的名义 P 值；weight=1 对应 |r| 权重
    analysis = gseapy.prerank(rnk=rank, gene_sets={k: sorted(v) for k, v in sets.items()},
                              min_size=15, max_size=500, permutation_num=1000, weight=1,
                              threads=2, ascending=None, method='multilevel', sample_size=101,
                              eps=1e-10, seed=seed, outdir=None, no_plot=True)

    # GSEApy 的结果列改成中文，并补上比较名、通路名和集合大小
    rows = analysis.res2d.rename(columns={'Term': '条目ID', 'ES': '富集分数ES', 'NES': '标准化富集分数NES',
                                           'NOM p-val': 'P值', 'Lead_genes': '核心基因EntrezID'}).copy()
    rows['比较'] = contrast
    rows['条目原名'] = rows['条目ID'].map(names)
    rows['排序基因数'] = len(rank)

    # 集合大小按“排序列表中实际出现的通路成员”重新计数
    rank_ids = set(rank.index)
    rows['集合大小'] = rows['条目ID'].map(lambda term: len(sets[term] & rank_ids))

    # 名义 P 值可能为空（置换不足），只对有限值做 BH
    rows['P值'] = pd.to_numeric(rows['P值'], errors='coerce')
    valid = rows['P值'].notna()
    rows['校正后P值'] = np.nan
    rows.loc[valid, '校正后P值'] = false_discovery_control(rows.loc[valid, 'P值'], method='bh')

    return rows[['比较', '条目ID', '条目原名', '集合大小', '排序基因数', '富集分数ES', '标准化富集分数NES',
                 'P值', '校正后P值', '核心基因EntrezID']]


# ---------------------------------------------------------------------------
# 4. 入口
# ---------------------------------------------------------------------------

def main():
    # Wald 表每个基因在每个对比各一行，是富集分析的输入
    wald = pd.read_csv(TABLES / 'timecourse_wald.tsv.gz', sep='\t')

    mapping, go, go_names, metadata = annotation_sets(set(wald['基因ID']))
    kegg, kegg_names = kegg_sets()

    # 生成 Ensembl ↔ Entrez ↔ 基因符号的对照表，供人工核对
    ids = mapping.set_index('ensembl_id').gene_id
    mapped = mapping.rename(columns={'ensembl_id': '基因ID', 'gene_id': 'EntrezID'})
    symbols = wald.drop_duplicates('基因ID').set_index('基因ID')['基因符号']
    mapped['基因符号'] = mapped['基因ID'].map(symbols)
    mapped.to_csv(TABLES / 'enrichment_gene_mapping.tsv', sep='\t', index=False, encoding='utf-8-sig')

    # 把 Entrez ID 作为新列拼回 Wald 表，供后续按对比分组
    wald['EntrezID'] = wald['基因ID'].map(ids)

    outputs = {'GO_ORA': [], 'KEGG_ORA': [], 'KEGG_GSEA': []}

    for i, (contrast, part) in enumerate(wald.groupby('比较', sort=False)):
        print(f'富集分析 {contrast}', flush=True)

        # 背景 = 该对比中校正后 P 值有限且能映射到 Entrez 的基因
        background = set(part.loc[part['校正后P值'].notna(), 'EntrezID'].dropna())

        # 上调、下调分别检验，方向由 log2FC 的符号决定
        for direction, sign in [('上调', 1), ('下调', -1)]:
            mask = part['通过差异筛选'].eq('是') & (sign * part['log2表达倍数变化'] > 0)
            selected = set(part.loc[mask, 'EntrezID'].dropna())

            # GO 三个本体与 KEGG 分别做 ORA，各自独立校正
            for ontology, sets in go.items():
                outputs['GO_ORA'].append(overrepresentation(sets, go_names, background, selected, contrast, direction, ontology))
            outputs['KEGG_ORA'].append(overrepresentation(kegg, kegg_names, background, selected, contrast, direction, 'KEGG'))

        # 每个对比一条 GSEA；种子随对比序号递增，方便复现单个结果
        outputs['KEGG_GSEA'].append(prerank(part, kegg, kegg_names, contrast, 252028 + i))

    # 汇总每个分析类型的全部对比，并在各对比内按校正后 P 值排序
    summary = []
    for kind, frames in outputs.items():
        table = pd.concat(frames, ignore_index=True).sort_values(['比较', '校正后P值', '条目ID'])
        table.to_csv(TABLES / f'enrichment_{kind.lower()}.tsv.gz', sep='\t', index=False, encoding='utf-8-sig')

        for contrast, rows in table.groupby('比较', sort=False):
            summary.append({'分析': kind, '比较': contrast, '检验条目数': len(rows),
                            '显著条目数': int(rows['校正后P值'].lt(0.05).sum())})

    summary = pd.DataFrame(summary)
    summary.to_csv(TABLES / 'enrichment_summary.tsv', sep='\t', index=False, encoding='utf-8-sig')

    # 记录版本、注释来源与统计口径，便于复核
    report = {'差异模型基因数': len(symbols), '一对一映射基因数': len(mapped), 'gseapy': gseapy.__version__,
              'scipy': scipy.__version__, '注释来源': metadata,
              'GSEA方法': 'multilevel', 'GSEA集合大小范围': [15, 500], 'GSEA基因集': 'KEGG mmu',
              'ORA集合大小范围': [10, 500], '校正方法': '各检验家族内 Benjamini–Hochberg'}
    (REPORTS / 'enrichment_python_session.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
