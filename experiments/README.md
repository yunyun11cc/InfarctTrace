# 实验总览

GSE252027 用于候选发现和队列内稳定性分析。GSE111059 和 GSE236374 为阶段 03 提供两套独立的表达证据。运行与数据来源见 [数据说明](../datasets/README.md)。

| 阶段 | 输入与交付 | 状态 |
|---|---|---|
| [01_1 发现队列数据处理与差异分析](01_1_preprocessing_differential/README.md) | 源数据和注释 → 质控、归一化、时间差异和富集结果 | Completed |
| [01_2 外部验证队列清洗](01_2_validation_preprocessing/README.md) | 芯片表达和 RNA-seq 计数 → 样本分组、清洗矩阵、质控与候选覆盖表 | Completed |
| [02 特征筛选](02_feature_selection/README.md) | D3/D7 差异基因交集与并集 → 候选排序、ssGSEA 相关性、200 折验证指标与 1,000 次 Bootstrap 入选频率 | Completed |
| [03 候选验证](03_validation/README.md) | 冻结候选及清洗后的文库级数据 → 队列内敏感性与外部表达方向 | Pending |
| [04 候选优选](04_prioritization/README.md) | 时间效应、稳定性与通路证据 → 排序与候选清单 | Pending |

## 模块结构

每个阶段的 `README.md` 说明目的、输入、运行方式与状态。阶段 01_1 的分析结果见 [阶段 01_1 解读](01_1_preprocessing_differential/results.md)，阶段 02 各版本和重采样分析见 [阶段 02 结果索引](../outputs/02_feature_selection/README.md)。学习材料集中在 [教学笔记](../notebooks/README.md)，计算公式和真实算例见 [阶段 01_1 数学原理](../docs/math/01_1_preprocessing_differential.md)。

## 后续拓展计划

GSE163465 用于后续查看候选在心脏免疫细胞及巨噬细胞亚群中的表达分布，目前为 `Pending`。
