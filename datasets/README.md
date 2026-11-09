# 数据说明

当前分析只使用 **GSE252027** 总 RNA 系列，聚焦小鼠心脏成纤维细胞。先运行 [数据下载脚本](../scripts/README.md) 恢复源文件，再用下面的说明核对分析输入。

| 数据 | 用途 |
|---|---|
| [GSE252027 数据手册](GSE252027/README.md) | 采样时间、FB 文库、三个 metadata CSV、定量格式、TPM 与计数、质控状态 |
| [基因注释来源](../provenance/annotation_download_2026-11-08.json) | GRCm38.87 注释的来源、校验值与使用范围 |
| [数据与论文参考](../docs/guides/dataset-reference.md) | GEO 记录、原始论文、作者代码和引用 |

`GSE252027/original/` 保存源归档，`GSE252027/metadata/` 保存样本信息和检查记录，`annotation/` 保存基因注释。预处理生成的表达矩阵和纳入表按阶段写入 [输出目录](../outputs/README.md)。

下一步按 [实验总览](../experiments/README.md) 查看各阶段使用的数据与交付物。
