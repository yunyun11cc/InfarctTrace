# 01_2 外部验证队列清洗

**状态：两套队列已清洗，候选效应检验留到阶段 03。**

| 队列 | 取材与样本 | 清洗后矩阵 |
|---|---|---|
| [GSE111059](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111059) | TCF21 谱系分选心脏成纤维细胞。未损伤、D3、D7、D14、D28 分别有 3、3、3、2、3 份芯片 | [28,846 个转录簇 × 14](../../outputs/01_2_validation_preprocessing/GSE111059/matrices/probe_sst_rma.tsv.gz) · [20,079 个 Entrez 基因 × 14](../../outputs/01_2_validation_preprocessing/GSE111059/matrices/gene_sst_rma.tsv.gz) |
| [GSE236374](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236374) | 小鼠心脏组织 bulk RNA-seq。Sham、MI D7、MI D28 各 3 份文库 | [24,964 个基因符号 × 9 原始计数](../../outputs/01_2_validation_preprocessing/GSE236374/matrices/raw_counts.tsv.gz) · [归一化表达](../../outputs/01_2_validation_preprocessing/GSE236374/matrices/log2_normalized_counts.tsv.gz) |

GSE111059 直接使用作者提交的 SST-RMA 表达值。GPL23038 转录簇按冻结的 Clariom S Mouse 注释映射到小鼠 Entrez ID。保留全量转录簇表用于追溯，基因矩阵使用 20,219 个唯一映射的转录簇，同一基因的多个转录簇逐样本取中位数。另有 31 个对照转录簇、779 个歧义转录簇和 7,817 个未映射转录簇，各自记录在映射表中。 [Bioconductor 注释](https://bioconductor.org/packages/clariomsmousetranscriptcluster.db/) · [逐转录簇映射](../../outputs/01_2_validation_preprocessing/GSE111059/tables/probe_mapping.tsv.gz)

GSE236374 从 GEO 样本标题及处理、时间和组织字段核对分组，以源列名对照 GSM。Sham 单独标识，`Day` 留空。原始计数保留所有基因，包括 5,629 个全零基因。在本队列中用各样本均非零的 14,641 个基因估计 median-of-ratios 大小因子，另存 `log2(normalized_count + 1)` 供质控和画图。后续计数模型读取原始整数计数。 [PyDESeq2 归一化方法](https://pydeseq2.readthedocs.io/en/stable/api/docstrings/pydeseq2.preprocessing.deseq2_norm.html)

两套队列按各自平台保留全部 23 份样本，矩阵列顺序与样本纳入表一致。两张样本质控表汇总芯片信号分布、RNA-seq 文库计数和大小因子。冻结清单的九个基因两套均覆盖，芯片通过 Ensembl → Entrez → 转录簇核对，RNA-seq 按提交的基因符号精确匹配。 [芯片样本表](../../outputs/01_2_validation_preprocessing/GSE111059/tables/sample_inclusion.tsv) · [芯片质控](../../outputs/01_2_validation_preprocessing/GSE111059/tables/sample_qc.tsv) · [RNA-seq 样本表](../../outputs/01_2_validation_preprocessing/GSE236374/tables/sample_inclusion.tsv) · [RNA-seq 质控](../../outputs/01_2_validation_preprocessing/GSE236374/tables/sample_qc.tsv) · [九基因覆盖表](../../outputs/01_2_validation_preprocessing/tables/candidate_coverage.tsv)

从项目根目录运行。首次运行用 `curl` 下载缺失源文件到各队列的 `datasets/<GSE>/original/`，每次均按源清单的 SHA-256 校验。两套输出各自存放。 [代码](prepare.py) · [源清单](sources.json) · [清洗报告](../../outputs/01_2_validation_preprocessing/reports/cleaning_summary.json)

```sh
uv run --locked python experiments/01_2_validation_preprocessing/prepare.py
uv run --locked python -m unittest discover -s experiments/01_2_validation_preprocessing -v
```
