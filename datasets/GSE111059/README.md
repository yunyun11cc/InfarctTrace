# GSE111059：心梗后心脏成纤维细胞的微阵列表达

更新日期：2026-11-08。当前副本已恢复处理后矩阵、原始 CEL 归档和完整 SOFT，并通过历史 SHA-256 校验。已用于[实验 01](../../experiments/01_fibroblast_screening/tutorial.md)的探索性独立复核；未与其他数据集合并。下文下载内容描述包括历史整理记录，论文 PDF 尚未在本副本恢复。

## 已下载文件

- `original/GSE111059_series_matrix.txt.gz`：GEO 公开处理后表达矩阵，28,846 个探针／转录簇行 × 14 个样本。包含对照探针，不应直接称为 28,846 个不同基因。
- `original/GSE111059_RAW.tar`：16,568,320 字节，内含 14 个原始 CEL.gz 芯片文件；不是测序 reads。本轮保留 TAR，不重复解包占空间。已逐个读取 gzip 至末尾检查完整性。
- `metadata/GSE111059_family.soft.gz`：GEO 完整 SOFT，包含平台注释、系列说明、14 个样本及各自处理后表达。
- `metadata/GSE111059_series.soft.txt`：系列级 GEO 原始文本记录。
- `metadata/GPL23038_annotation.extracted.tsv.gz`：从完整 SOFT 原样提取的平台表（29,129 行），未重新注释或映射基因。
- `metadata/samples.csv`：从 series matrix 的样本字段整理的 14 个样本及分组、取材说明和原始文件地址。
- `metadata/data_inspection.json`、`metadata/raw_archive_contents.json`：矩阵和原始文件初查。
- `download_records.json`：来源、字节数及 SHA-256，路径相对项目根目录。

## 实验设计

小鼠（Mus musculus），TCF21-MerCreMer;R26GFP 谱系示踪，FACS 分选 TCF21 谱系心脏成纤维细胞。平台为 GPL23038，Clariom S Mouse。

| 分组 | GEO 样本编号 | 样本数 |
|---|---|---:|
| 未损伤 | GSM3021170–GSM3021172 | 3 |
| MI 后 3 天 | GSM3021173–GSM3021175 | 3 |
| MI 后 7 天 | GSM3021176–GSM3021178 | 3 |
| MI 后 14 天 | GSM3021179–GSM3021180 | 2 |
| MI 后 28 天 | GSM3021181–GSM3021183 | 3 |

对照取自未损伤左心室，MI 样本说明为梗死区域；不要把未损伤对照自行改为 sham。样本原始名称标注 biological rep；本表忠实记录这些编号，不推断未提供的动物身份。

## 表达值与使用边界

GEO 的处理字段写明使用 Affymetrix Expression Console 生成 `.sst-rma-gene-full.chp`。公开矩阵值范围约 2.256–19.832，属于处理后芯片表达，不能作为 RNA-seq 原始计数送入 DESeq2。正式分析前仍需核对基因映射、控制探针及必要的质量检查。

该数据用于成纤维细胞层面的阶段表达复核，与其他 RNA-seq／单细胞队列按平台分别处理。实验 01 已完成注释映射、D7 对未损伤的 Welch 检验和 BH 校正；正式推进前需用 limma 复核小样本芯片统计。未训练预测模型或确定候选药物靶点。

## 原论文

Fu X, Khalil H, Kanisicak O, et al. **Specialized fibroblast differentiated states underlie scar formation in the infarcted mouse heart.** Journal of Clinical Investigation. 2018;128(5):2127–2143. DOI：[10.1172/JCI98215](https://doi.org/10.1172/JCI98215)；PMID 29664017；PMCID PMC5957472。

[正式论文 PDF](https://pmc-oa-opendata.s3.amazonaws.com/PMC5957472.1/PMC5957472.1.pdf)（17 页）；[补充材料 PDF](https://pmc-oa-opendata.s3.amazonaws.com/PMC5957472.1/jci-128-98215-s001.pdf)（11 页）；XML 全文与书目信息同目录保存。PDF 均来自官方 PMC 开放数据，已核文件头、页数并渲染确认正文首页。

[GEO 原始记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111059) · [PMC 原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC5957472/)
