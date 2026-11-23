# 阶段 01 发现队列数据处理与差异分析

**状态：Completed**

本阶段对应原始实验方案的 2.2，使用 GSE252027 成纤维细胞文库，运行流程采用 Python。输入为定量归档、FB 样本表和小鼠基因注释。 [原始方案入口](../../../docs/README.md) · [数据说明](../../../datasets/README.md)

第一次阅读可从 [阶段 01 流程说明](introduction.md) 开始，沿一个真实基因的记录理解输入、质控、差异检验和富集结果。

## 模块与运行顺序

| 步骤 | 模块 | 职责 |
|---|---|---|
| 1 | [导入与质控](../step01_load_and_qc.py) | 导入定量与基因注释，检查矩阵和样本纳入 |
| 2 | [标准化与 VST](../step02_normalization_vst.py) | 长度缩放、模型拟合、归一化矩阵、大小因子和 VST |
| 3 | [差异检验](../step03_differential_expr.py) | 总体 LRT、五个 Wald 对比和 D28r4 排除敏感性分析 |
| 4 | [候选窗口](../step04_window_selection.py) | D3/D7 候选交集、并集及与总体 LRT 的重叠 |
| 5 | [功能富集](../step05_enrichment.py) | 准备固定的注释快照并计算 GO/KEGG ORA 和 KEGG GSEA |
| 6 | [出图](../step06_visualize.py) | 生成质控图、热图、火山图和富集图 |

先按 [下载数据](../../../scripts/README.md) 准备环境和源文件，再从项目根目录运行。

```powershell
uv run --locked python -X utf8 experiments/01_preprocessing_differential/main.py
```

已有上游结果时，可用参数选择步骤。例如只生成窗口表和图时指定 `--steps 4 6`。
步骤 2 和 3 同时运行时，共用一次主模型拟合，并沿用先计算主检验再导出 VST 的顺序。
单独运行任一步时，从步骤 1 的定量矩阵重新拟合所需模型。
需要在其他位置做复算时，可将环境变量 `INFARCTTRACE_STAGE1_OUTPUT` 设为临时输出目录。

## 结果与结论

当前结果在 `outputs/01_preprocessing_differential/` 下按矩阵、表格、图和报告分类。

- [结果与读图](results.md)
- [阶段结论](results.md#6-阶段结论)
- [阶段输出索引](../../../docs/outputs.md)

## 方法参考

[样本和统计方案](protocol.md) · [数学原理与指标计算](../../../docs/math/01_1_preprocessing_differential.md)
