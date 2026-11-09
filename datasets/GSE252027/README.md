# GSE252027 数据手册

GSE252027 是小鼠心肌梗死后心脏细胞的总 RNA 表达数据。本项目的四个阶段均使用该队列，重点分析成纤维细胞（FB）。目标是识别 D0 至 D28 间表达明显变化的基因，并记录变化时间、方向和统计证据。源文件、样本信息、表达矩阵和差异结果各有不同用途。

## 1 实验设计与分析单位

作者用冠状动脉结扎建立小鼠心梗模型，在指定时间取心脏组织、分离细胞并进行 RNA 测序。D0 为未损伤、未手术对照。总 RNA 子系列为 GSE252027，论文还描述了同一研究的小 RNA 数据。 [GEO 数据页](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252027) · [原始研究](https://doi.org/10.1038/s41597-025-05061-1)

源归档包含 4 类细胞 × 6 个时间点 × 4 个提交文库，共 96 份定量文件。细胞类型为心肌细胞（CM）、内皮细胞（EC）、成纤维细胞（FB）和造血细胞（HC）。项目使用以下 24 个 FB 提交文库。

| 时间点 | 含义 | GSM 编号 | 提交文库数 |
|---|---|---|---:|
| D0 | 未损伤对照 | GSM7992326–GSM7992329 | 4 |
| D1 | 心梗后第 1 天 | GSM7992330–GSM7992333 | 4 |
| D3 | 心梗后第 3 天 | GSM7992334–GSM7992337 | 4 |
| D7 | 心梗后第 7 天 | GSM7992338–GSM7992341 | 4 |
| D14 | 心梗后第 14 天 | GSM7992342–GSM7992345 | 4 |
| D28 | 心梗后第 28 天 | GSM7992346–GSM7992349 | 4 |

28 天是观察窗口，实际有六个不等间隔采样点。不同时间点分别取材，因此按组间比较分析。每个文库通常汇总三只小鼠的材料，统计时一个文库算一次观测。`replicate` 是组内编号。同为重复 1 的 D0 和 D7 文库也按独立组处理。作者文中对重复称呼有不同表述，当前按提交文库记录。

## 2 文件用途

| 相对本目录的文件 | 内容和用途 |
|---|---|
| `original/GSE252027_RAW.tar` | GEO 的 123,955,200 字节归档，包含 96 个 gzip 压缩的定量表。归档内是作者提交的 Salmon 结果 |
| [`metadata/samples.csv`](metadata/samples.csv) | 96 个提交文库的样本说明，用于查看细胞类型和采样分组 |
| [`metadata/samples_FB.csv`](metadata/samples_FB.csv) | 筛出的 24 个 FB 文库。工程代码以 `gsm` 匹配文件，以 `day` 建立分组 |
| [`metadata/archive_contents.csv`](metadata/archive_contents.csv) | 每个归档成员的大小、哈希、行数、ID 层级和数值检查结果 |
| `metadata/GSE252027_series.soft.txt` | GEO 系列级说明，记录整体实验设计与样本清单 |
| `metadata/GSE252027_family.soft.gz` | 压缩的 GEO SOFT 记录，用于回查每个样本的原始说明 |
| `metadata/geo_supplement_filelist.txt` | GEO 发布的归档文件名及大小清单 |
| [`metadata/data_inspection.json`](metadata/data_inspection.json) | 归档结构、样本分组和数值检查的汇总 |
| `metadata/AMI_revised_submission.zip` | 作者发布的分析代码，用于查证样本取舍和原始分析操作 |
| `metadata/zenodo_14929916.json` | 作者代码的 Zenodo 发布信息 |
| [`download_records.json`](download_records.json) | 本队列源文件的地址、大小、哈希和历史记录 |
| `README.md` | 本手册 |

配套注释位于 `datasets/annotation/Mus_musculus.GRCm38.87.gtf.gz`，用于把稳定 Ensembl 基因 ID 映射到基因符号和类型。原始论文的书目信息与可用 XML 位于 `papers/GSE252027/`。当前文件的恢复与校验方式见 [数据下载指南](../../docs/guides/data-download.md)。

### 三个 CSV 里的关键列

- `samples.csv` 和 `samples_FB.csv` 共享同一套列。`gsm` 唯一标识文库，`day` 标识采样组，`replicate` 标识组内重复，`cell_type_abbreviation` 用于选 FB，`supplementary_files` 指向原始定量表。
- `pooled_mice_per_library` 记录池化信息，`paper_qc_status` 提示作者的质控信息。`not_assessed_here` 表示该表未做作者质控判定。D3 的标记表示样本取舍需要结合论文与代码核对。
- `archive_contents.csv` 的 `file` 与 `gsm` 连接归档成员和样本表，`rows` 和两个 Ensembl 计数列说明行的构成，`all_numeric_finite` 和 `all_numeric_nonnegative` 记录数值检查，`sha256` 用于确认文件内容。

三个 CSV 负责说明样本和检查文件。每个基因在每个文库中的数值存放在 `.quant.genes.sf.gz` 中。

## 3 一份定量表的结构

`GSM7992326_ctrl_FB_1.quant.genes.sf.gz` 是 D0 的第 1 个 FB 文库。gzip 解压后是制表符分隔的文本表，文件名中的 `GSM7992326` 在样本表中对应一行。

| 列名 | 含义 | 分析用途 |
|---|---|---|
| `Name` | 定量特征的 Ensembl ID，末尾可能带版本后缀 | 区分基因 `ENSMUSG` 与转录本 `ENSMUST`，建立行索引 |
| `Length` | 文件记录的目标长度 | 查看定量记录。基因汇总后的字段可随样本变化 |
| `EffectiveLength` | 定量模型考虑片段抽样后的有效长度 | 参与生成长度缩放计数 |
| `TPM` | 在每百万归一化丰度尺度上的相对表达量 | 查看文库间分布和时间变化，参与计数尺度转换 |
| `NumReads` | 分配给特征的估计片段数，可有小数 | 提供估计计数信息和文库总量 |

`TPM` 适合先读表达趋势。以 `Cthrc1` 为例，D0 的四个文库约为 0.651、0.574、0.504、0.494 TPM，D7 为 160.870、63.499、54.709、115.817 TPM。四个 D7 文库都较高，D7 组内仍有波动。正式差异结论再查看模型估计的 `log2表达倍数变化` 和 `校正后P值`。 [Salmon 格式说明](https://salmon.readthedocs.io/en/latest/file_formats.html)

每个 FB 文件实际有 51,398 行，包含 51,029 个基因行和 369 个转录本行。当前流程保留 `ENSMUSG` 基因行，去除 ID 版本后缀后核对唯一性，再按基因 ID 对齐所有文库。其他细胞类型文件的行层级有所不同，当前工程以 FB 的实际行类型为准。

## 4 整理后到哪一步

阶段 01 生成了 `FB_TPM.tsv.gz`、`FB_NumReads.tsv.gz`、`FB_EffectiveLength.tsv.gz` 三张 **51,029 个基因 × 24 个文库**的矩阵。它们的基因行与 GSM 列一致，格子分别存相对丰度、估计片段数和有效长度。基因注释文件 `gene_annotation.tsv.gz` 给出基因符号和类型。结果位置见 [输出索引](../../outputs/README.md)。

| 步骤 | 全时段分析进度 |
|---|---|
| 样本和数值质控 | 24 个文库完成文件、ID 和数值检查。主分析排除 D3r2、D3r3，纳入 22 个。D28r4 另做排除敏感性分析 |
| 计数导入和标准化 | 按 22 个主分析文库生成长度缩放计数、归一化计数和 VST 矩阵 |
| 差异分析 | 19,312 个基因进入全时段模型。总体时间效应达门槛 11,489 个，另有五个时间点与 D0 的对比结果 |
| GO、KEGG、GSEA | 已对五个时间对比分别完成 GO、KEGG 过度富集与 KEGG GSEA |

全时段差异模型读取 **长度缩放后的整数计数矩阵和六组样本分组表**。代码使用 TPM、平均有效长度和估计片段总量生成计数，再由模型进行文库大小归一化。`timecourse_length_scaled_counts.tsv.gz` 保存模型输入，`timecourse_normalized_counts.tsv.gz` 保存拟合后的归一化计数。此前 D7 对 D0 的 Python 教学基线仍保留，文件名以 `D0_D7` 开头。TPM 时间图用于阅读相对表达。 [工程分析方案](../../experiments/01_preprocessing_differential/protocol.md) · [全时段结果解读](../../experiments/01_preprocessing_differential/results.md) · [tximport 文档](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)

筛选时先核对样本和基因 ID，再看各组重复、变化幅度、校正后 P 值以及变化发生的时间。完整时间分析应区分早期升高后回落、持续变化和后期变化。阶段 02 的机器学习输入会在明确筛选目标后选择并处理。 [工程阶段总览](../../experiments/README.md)
