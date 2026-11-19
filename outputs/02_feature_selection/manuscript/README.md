# 阶段 01_1–02 主图与候选清单

三张图使用 GSE252027 中纳入主分析的 22 份小鼠成纤维细胞文库。Table 1 冻存了下一阶段要验证的九个候选基因及其 Ensembl ID。

| 交付物 | 下载 |
|---|---|
| Figure 1 文库质控、时序变化与候选池 | [PDF](figure_1.pdf) · [SVG](figure_1.svg) · [PNG](figure_1.png) · [TIFF](figure_1.tiff) |
| Figure 2 特征筛选、稳定性与表达轨迹 | [PDF](figure_2.pdf) · [SVG](figure_2.svg) · [PNG](figure_2.png) · [TIFF](figure_2.tiff) |
| Figure 3 通路得分与标志物共表达 | [PDF](figure_3.pdf) · [SVG](figure_3.svg) · [PNG](figure_3.png) · [TIFF](figure_3.tiff) |
| Table 1 冻存的九基因候选清单 | [TSV](table_02_candidate_genes_summary.tsv) |

图宽统一为 183 mm，高度依次为 151、161 和 91 mm。PDF 嵌入 TrueType 字体，SVG 保留可编辑文字，TIFF 为 600 dpi。 [Nature 图件规范](https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/)

## Figure 1 文库质控、时序变化与候选池

**a**，对按原研究记录排除两个 D3 文库后的 22 份文库做 PCA。D0、D1、D3、D7、D14、D28 分别包含 4、4、2、4、4、4 份文库。19,312 个基因的 VST 表达按基因中心化。空心符号对应单份文库，实心符号对应时间点质心，箭头连接相邻时间点质心。PC1 和 PC2 分别解释 31.5% 和 20.8% 的方差。 [样本取舍与处理口径](../../../experiments/01_1_preprocessing_differential/protocol.md) · [PCA 坐标](source_data/figure_1_library_pca.tsv)

**b**，D3 与 D7 分别对 D0 做 Wald 检验，差异门槛为 BH q < 0.05 且 |log₂FC| ≥ 1。六时点 LRT 采用 BH q < 0.05。D3 和 D7 两份差异基因集合的交集含 1,141 个基因，其中 1,139 个同时达到 LRT 门槛，另有 2 个保留在 Wald 交集中。柱高表示互斥交集的基因数。 [交集计数](source_data/figure_1_screen_overlap.tsv) · [逐基因筛选结果](source_data/figure_1_screen_membership.tsv)

**c**，对 1,141 个基因计算每个时间点的文库平均 VST，再对每个基因的六点轨迹做 z 标准化。用欧氏距离和 Ward 层次聚类分为五组，右侧为组内基因的平均 z 轨迹。五组依次为 D3 尖峰、D3 高表达后逐步回落、D7 峰值、D14 峰值及 D0 高表达，基因数依次为 309、205、340、109 和 178。 [逐基因轨迹](source_data/figure_1_candidate_profiles.tsv) · [簇均值](source_data/figure_1_cluster_centroids.tsv)

## Figure 2 特征筛选、稳定性与表达轨迹

**a**，使用 1,000 次分层 Bootstrap 数据重绘。每次有放回抽取 4 份 D7 和 18 份非 D7 文库，比较各方法入选集合和同一方法的全队列入选集合。箱体为四分位间距，中线为中位数，须到 1.5 倍四分位间距内的末端观察值，散点为须外观察值。四种方法的 Jaccard 中位数为 0.600、0.212、0.212 和 0.175。 [逐次 Bootstrap 结果](source_data/figure_2a_bootstrap_jaccard.tsv)

**b**，在固定的 3,907 基因并集候选池及阶段 01_1 VST 上做四折、50 次重复分层交叉验证。每个验证折含一份 D7 文库，曲线表示 200 折 AUROC 的均值。每折只用训练数据做特征排序和标准化，以同一 L2 逻辑回归（C = 1、类别平衡权重）比较四个基选择器和三个聚合策略。在 5–30 的每个特征数上均采用精确 k，排序边界并列时按 Ensembl ID 确定顺序。稀疏选择器得分不足 k 个正值时，剩余名额也依该规则填补。

聚合是对 Ding 等人三种思路的本项目实现。EFS-RA 取四种选择器绝对得分的平均降序名次，并列取平均名次。EFS-SA 取各方法 min–max 归一化后得分的平均值。EFS-IA 按基因同时进入四个排序前缀要求的最小前缀长度排序，完整计入并列块，再用 EFS-RA 打破并列。这三种实现使用原始四个选择器。统一分类器下，Elastic Net 选择器最高为 **0.99075（k = 26）**，EFS-RA 和 EFS-IA 最高均为 **0.98275（k = 25）**，EFS-SA 最高为 **0.97250（k = 23）**。空心星为历史 Elastic Net 筛选及原模型重拟合方案的 Top 20 AUROC **0.95725**。候选池和 VST 来自全队列，折间复用同一批文库。逐折标准差用于描述波动，曲线峰值供队列内比较。 [方法论文](https://doi.org/10.1016/j.patcog.2025.112809) · [均值与样本标准差](source_data/figure_2b_cv_summary.tsv) · [逐折数据](source_data/figure_2b_cv_folds.tsv) · [历史参考值](source_data/figure_2b_legacy_reference.tsv)

**c**，九个待验证基因按全队列 EFS-RA 名次排列。各形状对应 b 中四个基选择器，深色棒棒糖表示四个 Bootstrap 入选率的算术平均。沿用阶段 02 的 Bootstrap 筛选规则，单次抽样保留 Top 20 边界并列并排除零分。`Mrc1` 的平均入选率为 44.15%，在 3,907 个候选中居首，全队列 EFS-RA 名次为第 8。 [九基因入选率](source_data/figure_2c_candidate_frequencies.tsv) · [3,907 基因 EFS-RA 排名](source_data/figure_2c_efs_rank.tsv)

**d**，每个点代表一份文库，线连接各时间点的文库平均 VST。九个小图使用相同纵轴。横轴按采样时间顺序等距排列，阴影标示 D7。 [逐文库 VST](source_data/figure_2d_candidate_vst.tsv)

## Figure 3 通路得分与标志物共表达

**a**，计算九个候选基因的 VST 表达与七条通路的逐文库 ssGSEA NES 之间的 Spearman ρ。使用固定的 KEGG 和 Reactome 小鼠基因集。**b**，计算九基因与 `Nlrp3`、`Casp1`、`Il1b`、`Tgfb1`、`Smad3`、`Col1a1`、`Acta2` 的 VST 表达之间的 Spearman ρ。两图均使用 22 份文库，做 19,999 次双侧文库配对置换，每图分别对 63 个配对做 BH 校正。星号表示 q < 0.05，双星表示 q < 0.01。 [通路相关](source_data/figure_3a_pathway_correlations.tsv) · [逐文库 NES](source_data/figure_3a_pathway_scores.tsv) · [标志物相关](source_data/figure_3b_marker_correlations.tsv) · [逐文库表达](source_data/figure_3_expression.tsv)

通路面板有 35 个配对达到 q < 0.05，标志物面板有 19 个。`Mrc1` 与 ECM 组织和胶原形成分别呈正相关（ρ = 0.636 和 0.637，均为 q = 0.00589）。它与 NLRP3 炎症小体得分的相关系数为 0.470（q = 0.0516）。从评分背景移除九个候选基因后重新计算 NES，九基因的通路分类保持一致。这里合并了六个时间点，共同的时序变化可能贡献相关性。 [通路来源与原分析说明](../../../experiments/02_feature_selection/README.md#核心基因与通路-ssgsea-相关性)

## Table 1 冻存规则

清单沿用 22 文库并集的均值差排序前八基因，并加入四方法平均 Bootstrap 入选率第一的 `Mrc1`。EFS-RA 排名以全池 3,907 个基因为分母。四列入选率采用 0–1 比例。正相关且 BH q < 0.05 的 ssGSEA 配对决定通路分类。炎症和纤维化两组均有符合条件的配对标为“双通路”。按此规则，八个原候选为“双通路”，`Mrc1` 为“纤维化”。外部验证可直接读取 `基因符号` 和 `基因ID` 两列。

```text
SHA-256  7f7254627a79c5db2c664cf073122089d1ee46c9c1b16eec4921b036245fdb5b
```

## 分析与候选清单

本页图件和源数据为归档结果。在仓库根目录运行以下命令，可重算 EFS 比较并使用已保存的源数据生成九基因清单。

```sh
uv run --locked python experiments/02_feature_selection/manuscript_efs.py
uv run --locked python experiments/02_feature_selection/manuscript_table.py
```
