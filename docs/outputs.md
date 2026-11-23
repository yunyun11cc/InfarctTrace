# 实验产物索引

`outputs/` 按 `experiments/` 的阶段编号和名称组织。产物保存在相应工作流的 `matrices/`、`tables/`、`figures/` 和 `reports/` 中。统计方法和结果解读见各阶段文档。

| 阶段 | 输出组织 | 主要产物 |
|---|---|---|
| [01 预处理与差异分析](../experiments/01_preprocessing_differential/docs/README.md) | [`01_preprocessing_differential/`](../outputs/01_preprocessing_differential/) | [VST 矩阵](../outputs/01_preprocessing_differential/matrices/timecourse_vst.tsv.gz)、[Wald 表](../outputs/01_preprocessing_differential/tables/timecourse_wald.tsv.gz)、[D3/D7 交集](../outputs/01_preprocessing_differential/tables/d3_d7_intersection.tsv)、[D3/D7 并集](../outputs/01_preprocessing_differential/tables/d3_d7_union.tsv)、[富集汇总](../outputs/01_preprocessing_differential/tables/enrichment_summary.tsv) |
| [02 特征筛选](../experiments/02_feature_selection/docs/README.md) | [`02_feature_selection/`](../outputs/02_feature_selection/) 按 `intersection_6`、`intersection_22` 和 `union_22` 区分候选池 | [并集排序](../outputs/02_feature_selection/union_22/tables/feature_ranking.tsv)、[逐方法结果](../outputs/02_feature_selection/union_22/methods/)、[重采样](../outputs/02_feature_selection/union_22/resampling/)、[EFS 聚合](../outputs/02_feature_selection/union_22/aggregation/)、[生物学互证](../outputs/02_feature_selection/union_22/bio_validation/) |
| [03 外部队列验证](../experiments/03_external_validation/docs/README.md) | [`03_external_validation/`](../outputs/03_external_validation/) 按清洗、评估、单细胞和汇总分组 | [队列清洗](../outputs/03_external_validation/preprocessing/)、[外部评估](../outputs/03_external_validation/evaluation/)、[单细胞定位](../outputs/03_external_validation/single_cell/)、[跨队列汇总](../outputs/03_external_validation/summary/) |

## 跨阶段输入

- 阶段 02 读取阶段 01 的样本纳入表、VST、Wald 和 LRT。通路分析另读取基因注释与 Ensembl–Entrez 映射。
- 阶段 03 清洗和评估读取阶段 02 并集版的[冻结候选表](../outputs/02_feature_selection/union_22/aggregation/tables/frozen_candidates.tsv)，配合阶段 01 的映射和注释。跨队列汇总结合阶段 01 的 D7 对 D0 结果、阶段 03 的外部评估和单细胞定位。

## 复跑

运行各阶段 `main.py`，通过 `--steps` 复用已有的上游结果。阶段 01 用 `INFARCTTRACE_STAGE1_OUTPUT` 覆盖默认输出目录，阶段 02 和阶段 03 用 `--output-dir`。阶段 02 的 EFS 比较图、候选入选率图与共表达图可用 `--steps analysis-figures` 从已保存的统计表重绘。
