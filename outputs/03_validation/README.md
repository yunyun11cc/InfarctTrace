# 阶段 03 验证交付

阶段 02 的九个候选均在两个外部矩阵中检出。D7 对各自对照的双侧 Welch 检验按原始 P < 0.05、log₂ 尺度均值差 > 0 判定双队列晋级。 [方法与复现](../../experiments/03_validation/README.md)

| 基因 | GSE111059 均值差 | P | GSE236374 均值差 | P | 晋级 |
| ---     | --- | --- | ---:| ---:| ---|
| C1qtnf3 | +5.86 | 0.00153 | +5.90 | 0.000791 | 是 |
| Cthrc1 | +9.05 | 0.00000589 | +6.54 | 0.0000102 | 是 |
| Postn | +2.61 | 0.00261 | +5.64 | 0.000432 | 是 |
| Wisp1 | +3.29 | 0.00674 | +5.51 | 0.000642 | 是 |
| Cemip | +7.57 | 0.000146 | +6.37 | 0.00138 | 是 |
| Cxadr | +0.14 | 0.259 | −0.74 | 0.0381 | 否 |
| Kif26b | +6.22 | 0.0000469 | +3.52 | 0.00220 | 是 |
| Col24a1 | +1.00 | 0.0495 | +3.36 | 0.0112 | 是 |
| Mrc1 | +2.55 | 0.0848 | +1.39 | 0.0235 | 否 |

双队列原始 P 门槛共晋级 **7 个**。同时要求两队列 BH q < 0.05 时保留 **6 个**，Col24a1 的芯片 q 为 0.0637。芯片对照为未损伤组，全心 RNA-seq 对照为 Sham。

每个外部队列 D7 对对照为 3 对 3。八个基因的表达分数在两个队列观察到 AUROC = 1.00，Cxadr 分别为 0.667 和 0.000。冻结九基因 L2 逻辑回归在各自队列留一验证的 AUROC 与 AUPRC（AP）均为 1.00。完全分离时单基因 Bootstrap 产生的 [1, 1] 区间来自观测分布的边界。晋级七基因的留一结果使用了同批外部数据确定的子集，作为筛后描述。

六个标志物共表达的 54 对比较中，芯片 25 对、全心 RNA-seq 44 对达到 BH q < 0.05。三个 ssGSEA 通路的 27 对比较中，芯片 11 对、全心 RNA-seq 27 对达标。相关系数合并各采样日计算，组别和时间变化可能贡献相关性。

单细胞质控后共保留 **15,186 个细胞**。巨噬细胞标注群中 Mrc1 的检出率为 Sham **43.2%**、D3 **77.9%**、D7 **21.0%**。CD45⁺ 富集样本中观察到一个带 Dcn、Lum、Col1a1、Dpt 和 Pdgfra 表达的 **fibroblast-like** 群。该群检测到 Cthrc1 和 Postn，C1qtnf3 只在 Sham 的少数细胞中检出。这一群的标注限于表达相似性，样本的富集设计和文库重复数限制了对肌成纤维细胞状态及巨噬细胞与成纤维细胞互作的推断。每个时点仅一份公开文库，表中单细胞时间动态使用表达均值和检出率变化。

| 交付项 | 文件 |
|---|---|
| 九基因跨队列汇总、晋级状态和细胞定位 | [交接汇总表](tables/cross_cohort_validation.tsv) |
| 九基因切片与映射 | [GSE111059 切片](tables/GSE111059_frozen_expression.tsv) · [GSE236374 切片](tables/GSE236374_frozen_expression.tsv) · [ID 核对](tables/bulk_gene_mapping.tsv) |
| 外部表达箱线图与全部对比 | [芯片图](figures/GSE111059_expression.png) · [全心 RNA-seq 图](figures/GSE236374_expression.png) · [Welch 与 BH 结果](tables/bulk_expression_tests.tsv) |
| D7 外部表达分布与判别指标 | [九基因箱线图、散点及单基因 AUC](figures/bulk_roc.png) · [单基因和联合模型的 AUROC、区间与 AUPRC](tables/bulk_roc_metrics.tsv) |
| 经典标志物共表达及 ssGSEA | [标志物相关点图](figures/marker_correlations.png) · [通路相关点图](figures/pathway_correlations.png) · [标志物统计](tables/marker_correlations.tsv) · [通路统计](tables/pathway_correlations.tsv) |
| 单细胞定位与时间变化 | [UMAP 特征图](figures/single_cell_umap_features.png) · [气泡图](figures/single_cell_dotplot.png) · [逐文库表达](tables/sc_expression_by_group.tsv) · [对 Sham 的描述性变化](tables/sc_time_deltas.tsv) |
| 单细胞质控与标注 | [文库质控](tables/sc_library_qc.tsv) · [簇注释](tables/sc_cluster_labels.tsv) · [独立标志物表达](tables/sc_cluster_marker_expression.tsv) |

这套结果用于交接表达证据和候选晋级依据。成分与靶点关系需要结合药物成分数据与后续实验判断。
