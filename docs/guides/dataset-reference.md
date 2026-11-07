# InfarctTrace：数据与论文参考

更新日期：2026-11-08。

本页由项目原数据准备说明整理而来，集中保留队列、论文和检查细节。当前项目入口见 [README](../../README.md)，当前实验入口见[实验 01](../../experiments/01_fibroblast_screening/README.md)。

**当前阶段：计算机选题预研与首轮基因筛选。** 以导师原方案四个数据集为起点，允许按研究问题补充公开数据；目标是找到有依据、可评价的计算方法问题。首轮方向比较见 [计算机研究选题初筛](../research/topic-screening-2026-11-05.md)，最终计算方法题目尚未确定。

**2026-11-08 已恢复四队列数据并核对历史字节数与 SHA-256。** 首轮比较 GSE252027 FB D7 与 D0，并在 GSE111059 同时间比较中复核：17,622 个基因进入发现检验，1,501 个通过发现门槛，4 个通过预设跨队列门槛（Cthrc1、Ltbp2、Crlf1 上调，Gstm1 下调）。这是探索性表达筛选，尚未训练预测模型，也未证明因果或药物靶点。

开始学习请读 [实验设计与筛选规则](../../experiments/01_fibroblast_screening/protocol.md)、[读懂数据与基因筛选](../../experiments/01_fibroblast_screening/tutorial.md)，再运行 [教学 notebook](../../experiments/01_fibroblast_screening/notebook.ipynb)。下载与分析脚本位于 `scripts/`，结果位于本地 `outputs/01_fibroblast_screening/`。本轮验证记录见 [下载校验](../../provenance/download_verification_2026-11-08.json)和[归档完整性检查](../../provenance/archive_integrity_2026-11-08.json)。

已整理的原方案数据同时包含 bulk RNA、微阵列和单细胞数据，来自以下四个 GEO 队列。新增公开队列需另行核查并记录用途；之前误选的 Kuppe 2022 人类数据不自动纳入本方案。

## 仓库内容与数据获取

仓库同步研究文档、分析代码、教学 notebook、数据集说明、公开样本分组、数据检查结果、论文引用及下载清单。表达矩阵、压缩包、原始芯片文件、下载的论文全文及批量分析结果仅保存在本地，由 `.gitignore` 排除。本轮恢复了前三个 bulk/芯片队列的论文 XML，历史记录中的论文 PDF 尚未在本副本恢复。

**新环境先读 [数据与论文下载指南](data-download.md)。** 其中列出四个队列的 GEO 页面、官方直接下载地址、论文来源与本地保存位置；完整 URL、字节数和校验值在 [download_manifest.csv](../../provenance/download_manifest.csv)。

下文“已下载”和检查结果记录本次本地准备状态，不代表克隆仓库后已有数据。论文表的全文链接指向公开下载源；GSE163465 保留作者预印本与正式出版版本的区别。

## 1 下载状态与各队列用途

| 数据集 | 已下载数据 | 原方案中的用途 | 分组与数量 |
|---|---|---|---|
| [GSE252027](../../datasets/GSE252027/README.md) | 官方 Salmon 定量包，96 个 gzip 文件；完整样本元数据与作者代码 ZIP | 主要发现队列，重点关注 FB 成纤维细胞部分 | CM/EC/FB/HC 四种细胞，各 0/1/3/7/14/28 天 × 4 个提交样本；FB 为 24 个 |
| [GSE111059](../../datasets/GSE111059/README.md) | 28,846 个探针簇 × 14 样本的处理后芯片矩阵；14 个原始 CEL.gz；平台注释 | 独立成纤维细胞队列表达复核 | 未损伤、MI 后 3/7/14/28 天，分别 3/3/3/2/3 个 |
| [GSE236374](../../datasets/GSE236374/README.md) | 24,964 个唯一基因符号 × 9 样本的 FeatureCounts 计数矩阵 | 全心组织层面的补充验证 | sham、MI 后 7/28 天，各 3 个 |
| [GSE163465](../../datasets/GSE163465/README.md) | 7 个文库的 MTX/TSV；每库 31,053 个特征，合计 38,064 个提交条形码；已解包 | 免疫细胞／巨噬细胞表达定位与干预探索 | sham、MI 后 3/7/14 天，以及 3/7/14 天丹参酮 IIA 治疗组 |

已下载的是公开表达数据、样本说明和芯片原文件。GEO 文件名中的 `RAW.tar` 只是归档命名：GSE252027 内是定量结果，GSE163465 内是表达矩阵，不代表已下载 FASTQ 原始测序 reads。本轮未下载测序 reads。

所有队列均为小鼠研究，不能表述为人类临床诊断数据。细胞类型、平台和实验分组不同，不直接拼成一个训练集；后续如研究整合，需要单独定义科学问题、注释映射和评价方式。

## 2 数据集原论文

| 对应数据 | 原始研究论文 | 全文获取版本 |
|---|---|---|
| GSE252027 | Maji 等，2025，Scientific Data：*A cell type-specific expression atlas of small and total RNA in the heart after myocardial infarction.* [DOI](https://doi.org/10.1038/s41597-025-05061-1) | [正式 PDF，9 页](https://pmc-oa-opendata.s3.amazonaws.com/PMC12089340.1/PMC12089340.1.pdf)，另有全文 XML |
| GSE111059 | Fu 等，2018，JCI：*Specialized fibroblast differentiated states underlie scar formation in the infarcted mouse heart.* [DOI](https://doi.org/10.1172/JCI98215) | [正式 PDF，17 页](https://pmc-oa-opendata.s3.amazonaws.com/PMC5957472.1/PMC5957472.1.pdf)，另有补充 PDF、全文 XML |
| GSE236374 | Yu、Song 等，2023，PNAS：*Downregulation of apoptotic repressor AVEN exacerbates cardiac injury after myocardial infarction.* [DOI](https://doi.org/10.1073/pnas.2302482120) | [正式 PDF，11 页](https://pmc-oa-opendata.s3.amazonaws.com/PMC10589712.1/PMC10589712.1.pdf)，另有补充 PDF、全文 XML |
| GSE163465 | Jin 等，2022，Small Methods：*Single-Cell RNA Sequencing Reveals the Temporal Diversity and Dynamics of Cardiac Immunity after Myocardial Infarction.* [正式 DOI](https://doi.org/10.1002/smtd.202100752) | [原作者 2021 年预印本 PDF，47 页](https://www.biorxiv.org/content/10.1101/2021.07.19.452902v1.full.pdf?download=1)，另有原文抽取版和正式书目元数据；**未获取正式出版 PDF** |

GSE252027 的原论文覆盖父系列 GSE252028 的 small RNA 与 total RNA。本项目下载原方案指定的 total RNA 子系列 GSE252027，未另加 GSE252026。GSE163465 的预印本是同一研究的作者版本，不能假定与 2022 正式版逐字一致；[版本及获取说明](../../papers/GSE163465/README.md)已单独记录。

## 3 目录与可直接查看的文件

```text
InfarctTrace/  # 此处仅展示数据与文献部分
├── README.md
├── datasets/
│   ├── GSE252027/
│   ├── GSE111059/
│   ├── GSE236374/
│   └── GSE163465/
├── papers/
│   ├── GSE252027/
│   ├── GSE111059/
│   ├── GSE236374/
│   ├── GSE163465/
│   └── references.bib
├── docs/
└── provenance/
```

每个数据集的 `original/` 保留下载源文件，`metadata/` 保存 GEO 记录、样本表和检查结果，`README.md` 说明数据含义与限制。`original/` 指原样保存，并不统一表示原始 reads。

建议先看：

1. [发现队列 FB 样本表](../../datasets/GSE252027/metadata/samples_FB.csv)，再读该数据集原论文。
2. [GSE111059 样本表](../../datasets/GSE111059/metadata/samples.csv)和[平台注释说明](../../datasets/GSE111059/README.md)。
3. [GSE236374 样本表](../../datasets/GSE236374/metadata/samples.csv)。
4. [GSE163465 样本表](../../datasets/GSE163465/metadata/samples.csv)和[各文库矩阵统计](../../datasets/GSE163465/metadata/matrix_summary.csv)。

## 4 历史初查与处理点

- **GSE252027：特征层级需核对。** 文件名均含 `quant.genes.sf`，但 CM/EC 为转录本 ID；FB 为 51,029 个基因 ID 加 369 个转录本 ID；HC 在不同样本间也存在格式差异。尚未跨细胞类型拼接，也没有把 TPM 当计数使用。
- **GSE252027：重复和质控需保留依据。** 每个文库由三只小鼠细胞池化；原论文在 FB D3 排除两个异常重复，本轮保留全部 24 个 FB 提交文件。作者代码与正文对其他样本的取舍仍需核对，不擅自删样本。
- **GSE111059：是芯片，不是 RNA-seq。** 处理后表达值和原始 CEL 均保留，后续使用适合芯片的预处理与注释方式。
- **GSE236374：是整块心脏组织。** 不能当作纯化成纤维细胞队列；细胞组成变化会影响表达。
- **GSE163465：治疗组单列。** `_T` 对应丹参酮 IIA，不能当作普通 MI 重复，更不能当作芪参还五胶囊处理组。每条件只有一个提交文库，条形码数量不能替代独立动物重复。

上述初查记录来自数据准备阶段。2026-11-08 已开展 D0/D7 的基因级长度缩放、表达过滤、差异分析与独立芯片复核；未进行单细胞注释、跨平台直接拼接或预测模型训练。D3 排除记录争议未影响本轮 D0/D7 主要检验，描述性时间曲线仍保留全部 FB 提交文库。

## 5 来源与校验

[统一下载清单 JSON](../../provenance/download_manifest.json)／[CSV](../../provenance/download_manifest.csv)记录源文件、下载链接、字节数和本地 SHA-256；[SHA256SUMS.txt](../../provenance/SHA256SUMS.txt)可从项目根目录重新校验。解包副本和生成的样本表不冒充另行下载的源数据，逐文件核对结果保存在各数据集 `metadata/data_inspection.json`。

官方来源：[GSE252027](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252027)、[GSE111059](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111059)、[GSE236374](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236374)、[GSE163465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163465)。

[导师方案前期可行性评估](../research/proposal-feasibility-2026-11-02.md)保留为历史论证；当前执行范围以本页为准。
