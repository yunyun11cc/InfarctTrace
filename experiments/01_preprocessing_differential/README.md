# 阶段 01 数据处理与差异分析

本阶段对应 [原始实验方案](../../docs/README.md) 的 2.2，使用 GSE252027 成纤维细胞文库。已完成样本质控、全时段标准化、DESeq2 差异分析和 GO、KEGG、GSEA。主要结果与读图方法见 [结果解读](results.md)，完整口径见 [分析方案](protocol.md)。

## 模块与运行顺序

| 模块 | 职责 |
|---|---|
| `prepare.py` | 导入原始定量与基因注释，检查矩阵并写出文库纳入表 |
| `differential.R` | 在 22 文库上完成长度缩放、标准化、总体 LRT 和五个 Wald 对比，另做 D28r4 排除敏感性分析 |
| `go_reference.py`、`kegg_reference.py` | 准备可追溯的小鼠注释和 KEGG 通路快照 |
| `enrichment.py` | Python 实现 GO、KEGG 过度富集检验和 GSEA |
| `visualize.py` | 输出中文质控图、热图、火山图和富集图 |
| `baseline.py`、`run.py` | 原 D7 对 D0 的 Python 教学基线及薄兼容入口 |

在项目根目录执行。当前工作区的 R 4.5.3 和包位于 `.tools/`。其他机器可用可用的 `Rscript`，配套 R 4.5 与 Bioconductor 3.22 的 `DESeq2`、`tximport`。两份 GO SQLite 数据库由 Python 模块提取，富集阶段用 `uv.lock` 中的 Python 依赖。

```powershell
uv sync --locked
uv run --locked python -X utf8 scripts\download_data.py
uv run --locked python -X utf8 experiments\01_preprocessing_differential\prepare.py
$rscript = '.\.tools\R-4.5.3\bin\Rscript.exe'
& $rscript --vanilla --encoding=UTF-8 experiments\01_preprocessing_differential\differential.R
uv run --locked python -X utf8 experiments\01_preprocessing_differential\go_reference.py
uv run --locked python -X utf8 experiments\01_preprocessing_differential\kegg_reference.py
uv run --locked python -X utf8 experiments\01_preprocessing_differential\enrichment.py
uv run --locked python -X utf8 experiments\01_preprocessing_differential\visualize.py
```

Windows 上如 R 启动提示 `C.UTF-8` 区域设置失败，在当前终端移除 `LANG`、`LC_ALL` 及相应 `LC_*` 覆盖后再运行，让 R 使用系统的 UTF-8 区域设置。运行环境、数据库版本、校验与分析参数分别记录在 `outputs/01_preprocessing_differential/reports/`。

本机 R 4.5.3 在加载部分包后退出时出现 Windows 访问冲突（`0xC0000005`）。两次全时段运行均先完成输出，统计汇总一致，表格维度和结果字段已单独校验。该环境的 R 模块退出状态仍需修复。Python 富集与绘图模块正常结束。

[阶段输出索引](../../outputs/README.md) · [中文数据阅读笔记](../../notebooks/README.md)
