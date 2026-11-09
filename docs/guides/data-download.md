# GSE252027 数据下载和恢复

当前分析使用 GSE252027 总 RNA 系列。仓库记录源地址与 SHA-256。大体积定量归档、基因注释和本地论文全文由下载脚本恢复。

```powershell
uv sync --locked
uv run --locked python -X utf8 scripts\download_data.py
```

脚本按 [下载清单](../../provenance/download_manifest.csv) 和 [注释来源记录](../../provenance/annotation_download_2026-11-08.json) 恢复定量归档、family SOFT、作者代码、可获得的论文 XML 与 GRCm38.87 注释。已有文件按 SHA-256 校验，运行记录写到 `outputs/01_preprocessing_differential/reports/download_verification.json`。

每个源文件及其用途见 [GSE252027 数据手册](../../datasets/GSE252027/README.md)。源清单同时保留论文 PDF 等参考来源，实际分析所需文件由脚本筛选。
