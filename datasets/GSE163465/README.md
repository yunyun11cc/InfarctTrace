# GSE163465：心梗后小鼠心脏免疫细胞单细胞转录组

来源：[GEO GSE163465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163465)。获取日期：2026-11-03。

这是导师原方案中用于免疫细胞／巨噬细胞亚群表达定位的数据集；与其他三个系列独立保存，尚未整合或训练模型。

## 数据文件与样本说明

- `original/GSE163465_RAW.tar`：GEO 官方提交的 MTX／TSV 数据包。文件名虽然含 RAW，内容为表达矩阵、特征表和细胞条形码，**不是 FASTQ 原始测序 reads**。
- `original/extracted/`：按 GEO 原始文件名解出的七组 10x 风格三件套。原始 TAR 保留，解压文件不作修改。
- `metadata/GSE163465_family.soft.gz`：包含全部七个 GSM 及实验方案的 GEO family SOFT。
- `metadata/GSE163465_series_matrix.txt.gz`：Series Matrix；本系列的细胞表达数据位于 MTX 文件，Series Matrix 主要提供样本说明，不能把这个小文件当完整表达矩阵。
- `metadata/GSE163465_series.soft.txt`：系列级原始记录。
- `metadata/filelist.txt`：GEO 官方归档及成员文件字节数清单。
- `metadata/samples.csv`／`samples.json`：从 SOFT 提取的分组表，药物处理和疾病／时间分别记录。
- `metadata/data_inspection.json`：归档成员、解压文件、gzip 完整性及矩阵维度核对结果。
- `metadata/matrix_summary.csv`：每个 GSM 的基因数、条形码数和非零元素数。

## 七个文库的分组

| GSM | GEO 标题 | 条件 | 心梗后天数 | 丹参酮 IIA 处理 |
|---|---|---|---|---|
| GSM4985022 | CD45+ Sham | sham | 不适用 | 无 Tan IIA 处理组标记 |
| GSM4985023 | CD45+ Day_3 | MI | 3 | 无 Tan IIA 处理组标记 |
| GSM4985024 | CD45+ Day_7 | MI | 7 | 无 Tan IIA 处理组标记 |
| GSM4985025 | CD45+ Day_14 | MI | 14 | 无 Tan IIA 处理组标记 |
| GSM5513773 | Cd45+ Day_3_T | MI | 3 | 有 |
| GSM5513774 | Cd45+ Day_7_T | MI | 7 | 有 |
| GSM5513775 | Cd45+ Day_14_T | MI | 14 | 有 |

研究对象为小鼠心脏富集的 Cd45+ 免疫细胞，10x Single Cell 3′ v2，mm10 参考基因组。`_T` 为丹参酮 IIA 干预，来自原论文／作者预印本 Fig. 6 及对应 GEO 标题；不能当成额外未处理重复，也不能当成芪参还五胶囊处理组。

每个条件／时间点在 GEO 中只有一个 GSM。作者预印本 Fig. 6 说明同一时间点同时处理五只小鼠，但目前没有核实能拆分到单鼠的文库标识，**不能把五只小鼠直接计为五个独立单细胞测序重复，也不能把单个细胞当独立动物重复**。

这是免疫细胞富集数据，不能当作包含完整成纤维细胞谱系的全心图谱。表达矩阵之外的细胞亚型标签仍需核对作者资料或重新注释；本轮没有产生细胞标签或差异基因结论。

## 原论文

Jin K 等，**Single-Cell RNA Sequencing Reveals the Temporal Diversity and Dynamics of Cardiac Immunity after Myocardial Infarction.** Small Methods，2022，6(3):e2100752。DOI：[10.1002/smtd.202100752](https://doi.org/10.1002/smtd.202100752)。

[论文目录及版本说明](../../papers/GSE163465/README.md)中已保存正式书目元数据、作者预印本 PDF 与原文抽取版；预印本 PDF 不标为正式出版 PDF。

下载来源、字节数和 SHA-256 见本目录的 `download_records.json`。解压文件的来源为该官方 TAR，属于展开副本而非另行下载或分析生成的表达数据。

## 本地文件核验结果

TAR 为 **268,195,840 字节**，与 GEO 官方文件清单一致；21 个成员的文件名和字节数全部吻合。全部 gzip 文件已读至末尾以检查完整性，七份 Matrix Market 头均为 `matrix coordinate integer general`。基因表、条形码表和矩阵维度一致，矩阵数据记录数与头部非零元素数一致。

| GSM | 提交矩阵中的基因数 | 提交的细胞条形码数 |
|---|---:|---:|
| GSM4985022 | 31,053 | 2,962 |
| GSM4985023 | 31,053 | 7,506 |
| GSM4985024 | 31,053 | 5,482 |
| GSM4985025 | 31,053 | 3,942 |
| GSM5513773 | 31,053 | 3,928 |
| GSM5513774 | 31,053 | 6,570 |
| GSM5513775 | 31,053 | 7,674 |
| 合计 | 各矩阵均为 31,053 行 | 38,064 |

这些是本次下载的提交矩阵条形码数，**不是独立小鼠数量，也不保证等于原文每一步质控后的最终细胞数**。本轮未重新过滤、标准化或合并表达矩阵。未逐项验证全部坐标和计数取值；上述检查针对归档完整性、gzip、维度及记录数。
