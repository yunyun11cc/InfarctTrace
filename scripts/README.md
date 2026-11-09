# 下载数据

从仓库根目录开始。当前项目使用 GSE252027 的定量归档与小鼠基因注释。项目依赖由 `pyproject.toml` 和 `uv.lock` 管理，Python 版本入口见 `.python-version`。

## 运行

以下是 Windows PowerShell 的运行方式。系统需有 `uv` 与 `curl.exe`。

```powershell
uv sync --locked
uv run --locked python -X utf8 scripts\download_data.py
```

[下载脚本](download_data.py) 读取 [源文件清单](../provenance/download_manifest.csv)，恢复 GSE252027 定量归档、GEO 样本记录、作者代码、可取得的论文 XML 和 GRCm38.87 注释。已经存在的文件按 SHA-256 校验。校验报告写入 `outputs/01_preprocessing_differential/reports/download_verification.json`。源文件大小、恢复范围和具体来源见 [数据下载指南](../docs/guides/data-download.md)。

下载成功后，先看 [数据说明](../datasets/README.md)，再按 [阶段 01 README](../experiments/01_preprocessing_differential/README.md) 进行预处理和差异分析。GO 与 KEGG 注释快照按阶段 01 的模块说明获取。
