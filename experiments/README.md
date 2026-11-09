# 工程实验

四个阶段均使用 GSE252027 成纤维细胞数据。每个阶段的工程代码和分析说明放在这里，学习材料集中在 [notebooks](../notebooks/README.md)。结果按同名阶段写入 [outputs](../outputs/README.md)。

| 阶段 | 工作与交付 | 当前状态 |
|---|---|---|
| [01 数据处理与差异分析](01_preprocessing_differential/README.md) | 样本质控、计数导入、时间差异及功能富集 | 已完成六组分析、D28 敏感性分析和五个时间对比的富集 |
| [02 特征筛选](02_feature_selection/README.md) | 明确筛选目标，比较候选基因及稳定性 | 待开展 |
| [03 候选验证](03_validation/README.md) | 队列内复现、样本敏感性与时间趋势检查 | 待开展 |
| [04 候选优选](04_prioritization/README.md) | 综合效应、阶段信息、稳定性和通路证据 | 待开展 |

工程代码用中文注释解释样本口径和统计操作。输出表头与图的说明用中文，保留可追溯的基因 ID、GSM 编号和指标缩写。
