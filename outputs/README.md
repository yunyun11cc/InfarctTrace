# 输出目录

实验结果随代码纳入版本控制。四个阶段与 `experiments/` 下的同名目录对应，每阶段按用途分为四类。

| 子目录 | 格式和用途 |
|---|---|
| `matrices/` | `.tsv.gz`。表达矩阵、计数矩阵及配套基因注释 |
| `tables/` | `.csv`、`.tsv`。样本信息、质控、差异或候选结果。表头使用中文 |
| `figures/` | `.png`。带中文图题、坐标、图例和必要图注的科研图 |
| `reports/` | `.json`。参数、运行摘要和校验记录，保留稳定的字段名以便脚本读取 |

| 阶段 | 状态 | 输出入口 |
|---|---|---|
| `01_preprocessing_differential/` | Completed | [阶段 01 结果与结论](../experiments/01_preprocessing_differential/results.md) |
| `02_feature_selection/` | Completed | [阶段 02 结果](../experiments/02_feature_selection/results.md) |
| `03_validation/` | Pending | [阶段 03 README](../experiments/03_validation/README.md) |
| `04_prioritization/` | Pending | [阶段 04 README](../experiments/04_prioritization/README.md) |

## 阶段 01 图表说明

| 编号与标题 | 文件 | 读表和读图说明 |
|---|---|---|
| 表 1 GSE252027 成纤维细胞提交文库及采样分组 | `tables/discovery_samples.csv` | 每行一个文库，共 24 个，D0 为未损伤对照 |
| 表 2 成纤维细胞文库的表达定量质控统计 | `tables/discovery_qc.csv` | 按文库汇总特征行数、估计片段数和 TPM 分布信息 |
| 图 1 成纤维细胞提交文库的主成分分析 | `figures/01_submitted_library_pca.png` | 24 个提交文库的描述性 PCA，标注待核对的样本 |
| 图 2 主分析文库标准化前后的计数分布 | `figures/02_timecourse_normalization.png` | 对比同一批文库和基因的长度缩放计数与大小因子归一化计数 |
| 图 3 标准化后成纤维细胞文库的主成分分析 | `figures/03_timecourse_pca.png` | 22 个主分析文库的 VST 表达 PCA |
| 图 4 总体时间效应排名前 25 个基因的表达模式 | `figures/04_timecourse_heatmap.png` | 逐行标准化，比较同一基因的各文库表达 |
| 图 5 心梗后各时间点相对未损伤对照的差异表达 | `figures/05_timecourse_volcano.png` | 五个 Wald 对比，虚线标记效应与校正后 P 值门槛 |
| 图 6 心梗后第 7 天的功能与通路富集 | `figures/06_enrichment_D7.png` | GO 生物过程、KEGG 上调过度富集和正负方向 GSEA |

全时段的样本取舍、差异数量、敏感性与富集表如何阅读，见 [阶段 01 结果解读](../experiments/01_preprocessing_differential/results.md)。
