# GSE252027：主要发现队列

2026-11-08 状态：当前副本已恢复 GEO 官方完整定量包、96 个样本元数据、论文 XML 和作者代码包，源文件通过历史 SHA-256 校验。已完成 FB D0/D7 的探索性差异筛选，见[实验说明](../../experiments/01_fibroblast_screening/tutorial.md)；历史记录中的论文 PDF 尚未在本副本恢复。

## 下载内容

- [GSE252027_RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252027/suppl/GSE252027_RAW.tar)：123,955,200 字节，含 96 个 `.quant.genes.sf.gz` 文件。这是作者提交的 Salmon 表达定量文件，不是 FASTQ。
- [完整样本表](metadata/samples.csv)：CM、EC、FB、HC 四种细胞类型，0/1/3/7/14/28 天，每个组合 4 个测序样本。
- [FB 样本表](metadata/samples_FB.csv)：成纤维细胞 24 个提交样本，是原方案的重点部分。
- [GEO family SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252027/soft/GSE252027_family.soft.gz)、[GEO 系列记录](metadata/GSE252027_series.soft.txt)。
- [检查结果](metadata/data_inspection.json)、[压缩包逐文件清单](metadata/archive_contents.csv)。
- [原作者代码 ZIP](https://zenodo.org/api/records/14929916/files/AMI_revised_submission.zip/content)，来自论文提供的 [Zenodo 14929916](https://zenodo.org/records/14929916)，MD5 与官方记录一致。它不是本项目编写的分析代码。

## 初查与后续需要处理的问题

已完整解压读取 tar 内全部 96 个 gzip 定量文件，GSM 与样本表一一对应；定量数值均有限、非负。每个文件含 `Name`、`Length`、`EffectiveLength`、`TPM`、`NumReads`。`NumReads` 是 Salmon 估计值，可能含小数，不把 TPM 当原始计数。

文件名虽然都写 `quant.genes.sf`，内部特征层级并不统一：CM/EC 各 24 份为 137,734 个转录本；FB 24 份为 51,398 行，其中 51,029 个 ENSMUSG 基因 ID 和 369 个 ENSMUST 转录本 ID；HC 的 24 份中两类格式均存在。因此暂不跨细胞类型拼接，基因／转录本注释需后续核对。

每个测序重复由三只小鼠的细胞池化，不能拆作三个独立测序样本。D0 是未损伤对照，不自动等同 sham。原论文说 total RNA 的 FB D3 剔除两个异常重复，FB 分析量为 22；本次保留全部 24 个源文件。作者代码 `deg_analysis_FB.R` 保留 D3 r1/r4，并另涉及 D28 r4 省略，与正文描述尚需进一步核对；不擅自删除样本或宣布已经完成质控。

## 对应原论文

Maji RK 等，2025，Scientific Data。*A cell type-specific expression atlas of small and total RNA in the heart after myocardial infarction.* DOI：[10.1038/s41597-025-05061-1](https://doi.org/10.1038/s41597-025-05061-1)。

[正式 PDF](https://pmc-oa-opendata.s3.amazonaws.com/PMC12089340.1/PMC12089340.1.pdf)，9 页；[全文 XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12089340/fullTextXML)。该论文覆盖总系列 GSE252028；本项目只下载方案指定的 total RNA 子系列 GSE252027，未另加 small RNA 子系列 GSE252026。

来源：[GEO GSE252027](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252027)。每个源文件的来源和哈希见 [download_records.json](download_records.json)。
