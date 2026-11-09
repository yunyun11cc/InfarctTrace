# 实验总览

四个阶段均使用 GSE252027 成纤维细胞数据。运行前完成 [数据下载](../scripts/README.md)，并阅读 [数据说明](../datasets/README.md)。当前工程使用 Python，阶段 03 的验证限定在 GSE252027 队列内部。

| 阶段 | 输入与交付 | 状态 |
|---|---|---|
| [01 数据处理与差异分析](01_preprocessing_differential/README.md) | 源数据和注释 → 质控、归一化、时间差异和富集结果 | Completed |
| [02 特征筛选](02_feature_selection/README.md) | 阶段 01 的表达特征与统计证据 → 候选特征与筛选稳定性 | Pending |
| [03 候选验证](03_validation/README.md) | 候选特征与文库级数据 → 队列内稳定性与时间趋势证据 | Pending |
| [04 候选优选](04_prioritization/README.md) | 时间效应、稳定性与通路证据 → 排序与候选清单 | Pending |

## 模块结构

每个阶段的 `README.md` 说明目的、输入、运行方式与状态。完成的阶段从 README 进入 `results.md`，查看结果、解释和结论。待开展的阶段在各自 README 标明 `Pending`。生成文件按同名阶段归入 [输出目录](../outputs/README.md)。学习材料集中在 [教学笔记](../notebooks/README.md)，计算公式和真实算例见 [阶段 01 数学原理](../docs/math/01_preprocessing_differential.md)。

## 后续拓展计划

以下方向状态均为 `Pending`。启用时另行确定分析设计与时间，当前四阶段仍使用 GSE252027。

| 暂缓队列 | 拓展问题 |
|---|---|
| GSE111059 | 在独立的成纤维细胞芯片研究中复核共有时间点的表达方向，届时使用 R 处理芯片数据 |
| GSE236374 | 检查候选信号在 D7 和 D28 的全心组织层面是否同向 |
| GSE163465 | 查看候选在心脏免疫细胞及巨噬细胞亚群中的表达分布 |
