"""保留原 D7 对 D0 的 Python 教学基线。"""
import json

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

from paths import MATRICES, TABLES, FIGURES, REPORTS
from prepare import annotation, length_scaled_counts, load_fb, save_gene_table
from visualize import configure_plot_font, submitted_library_pca


def discovery(meta, matrices, gene_info):
    """保持既有 D7 对 D0 的对比、过滤门槛与模型参数。"""
    selected = meta.index[meta.day.isin([0, 7])]
    counts = length_scaled_counts(*(matrices[k][selected] for k in
                                   ('TPM', 'EffectiveLength', 'NumReads')))
    save_gene_table(counts, MATRICES / 'D0_D7_lengthScaledTPM_counts_all_genes.tsv.gz')
    # 先按计数过滤低表达基因，再用两组独立文库拟合差异模型。
    keep = (counts >= 10).sum(axis=1) >= 4
    design = pd.DataFrame({'condition': np.where(meta.loc[selected, 'day'] == 0, 'D0', 'D7')}, index=selected)
    inference = DefaultInference(n_cpus=2)
    dds = DeseqDataSet(counts=counts.loc[keep].T, metadata=design, design='~condition',
                       refit_cooks=True, inference=inference, quiet=False)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=['condition', 'D7', 'D0'], alpha=0.05,
                      inference=inference, quiet=True)
    stats.summary()
    result = stats.results_df.join(gene_info)
    result['discovery_pass'] = (result.padj < 0.05) & (result.log2FoldChange.abs() >= 1)
    labels = {'baseMean': '平均标准化计数', 'log2FoldChange': 'log2表达倍数变化', 'lfcSE': '效应标准误',
              'stat': 'Wald统计量', 'pvalue': 'P值', 'padj': '校正后P值',
              'symbol': '基因符号', 'biotype': '基因类型', 'discovery_pass': '通过差异筛选'}
    table = result.rename(columns=labels).copy()
    table['通过差异筛选'] = table['通过差异筛选'].map({True: '是', False: '否'})
    save_gene_table(table, TABLES / 'discovery_D7_vs_D0.tsv')
    normalized = pd.DataFrame(dds.layers['normed_counts'].T, index=counts.index[keep], columns=selected)
    save_gene_table(normalized, MATRICES / 'D0_D7_normalized_counts.tsv.gz')
    return result, int(keep.sum())


def main():
    configure_plot_font()
    for directory in (MATRICES, TABLES, FIGURES, REPORTS):
        directory.mkdir(parents=True, exist_ok=True)
    print('读取基因注释与成纤维细胞表达数据。', flush=True)
    genes = annotation()
    meta, matrices, excluded = load_fb()
    print('拟合 D7 对 D0 差异表达模型。', flush=True)
    result, tested = discovery(meta, matrices, genes)
    summary = dict(dataset='GSE252027', cell_type='FB', contrast='D7_vs_D0',
                   discovery_submitted_FB_libraries=len(meta), gene_rows=len(matrices['TPM']),
                   discarded_non_gene_rows=excluded,
                   discovery_D0_n=int((meta.day == 0).sum()), discovery_D7_n=int((meta.day == 7).sum()),
                   discovery_genes_after_count_filter=tested,
                   discovery_finite_padj=int(result.padj.notna().sum()),
                   discovery_pass=int(result.discovery_pass.sum()),
                   annotation='Ensembl GRCm38 release 87',
                   interpretation='GSE252027 成纤维细胞 D7 对 D0 的探索性差异表达结果')
    submitted_library_pca(meta, matrices['TPM'])
    (REPORTS / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(f"{tested:,} 个基因进入检验，{summary['discovery_pass']:,} 个通过差异筛选。", flush=True)

