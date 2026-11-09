# InfarctTrace

小鼠心肌梗死后成纤维细胞的时序转录组分析。当前四个实验阶段均使用 **GSE252027**，工作流以 Python 运行。

## 快速开始

1. [下载数据](scripts/README.md) — 配置环境，运行下载脚本并校验源文件。
2. [数据说明](datasets/README.md) — 了解数据来源、采样分组、文件格式和分析注释。
3. [实验总览](experiments/README.md) — 查看各阶段的输入、交付物、执行状态和后续计划。
4. [输出索引](outputs/README.md) — 按实验查找矩阵、统计表、图片与报告。

## 实验框架

**数据下载 → 01 数据处理与差异分析 → 02 特征筛选 → 03 队列内验证 → 04 候选优选**

| 阶段 | 模块入口 | 状态 |
|---|---|---|
| 01 | [数据处理与差异分析](experiments/01_preprocessing_differential/README.md) | Completed |
| 02 | [特征筛选](experiments/02_feature_selection/README.md) | Pending |
| 03 | [候选验证](experiments/03_validation/README.md) | Pending |
| 04 | [候选优选](experiments/04_prioritization/README.md) | Pending |

## 参考与学习

- [项目文档](docs/README.md) — 原始实验方案、数据来源和数学原理。
- [教学笔记](notebooks/README.md) — 阅读样本表、表达矩阵和实验图表。
