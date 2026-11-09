# 阶段 01 数据处理与差异分析

**状态：Completed**

本阶段对应原始实验方案的 2.2，使用 GSE252027 成纤维细胞文库，运行流程采用 Python。输入为定量归档、FB 样本表和小鼠基因注释。 [原始方案入口](../../docs/README.md) · [数据说明](../../datasets/README.md)

## 模块与运行顺序

| 模块 | 职责 |
|---|---|
| `prepare.py` | 导入定量结果与基因注释，检查矩阵并写出文库纳入表 |
| `differential.py` | 在 22 个文库上完成长度缩放、标准化、总体 LRT 和五个 PyDESeq2 Wald 对比，另做 D28r4 排除敏感性分析 |
| `go_reference.py`、`kegg_reference.py` | 准备小鼠注释和 KEGG 通路快照 |
| `enrichment.py` | GO、KEGG 过度富集检验和 GSEA |
| `visualize.py` | 生成中文质控图、热图、火山图和富集图 |

先按 [下载数据](../../scripts/README.md) 配置环境并恢复源文件，再从项目根目录按顺序执行。已有的注释快照可以复用。

```powershell
uv run --locked python -X utf8 experiments\01_preprocessing_differential\prepare.py
uv run --locked python -X utf8 experiments\01_preprocessing_differential\differential.py
uv run --locked python -X utf8 experiments\01_preprocessing_differential\go_reference.py
uv run --locked python -X utf8 experiments\01_preprocessing_differential\kegg_reference.py
uv run --locked python -X utf8 experiments\01_preprocessing_differential\enrichment.py
uv run --locked python -X utf8 experiments\01_preprocessing_differential\visualize.py
```

## 结果与结论

当前结果在 `outputs/01_preprocessing_differential/` 下按矩阵、表格、图和报告分类。

- [结果与读图](results.md)
- [阶段结论](results.md#6-阶段结论)
- [阶段输出索引](../../outputs/README.md)

## 方法参考

[样本和统计方案](protocol.md) · [数学原理与指标计算](../../docs/math/01_preprocessing_differential.md) · [中文数据阅读笔记](../../notebooks/README.md)
