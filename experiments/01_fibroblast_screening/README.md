# Experiment 01: Cross-Cohort Fibroblast Gene Screening

研究问题：心梗后第 7 天相对于未损伤对照，哪些小鼠成纤维细胞基因发生表达变化，并在独立芯片队列获得同向支持？

## 阅读顺序

1. [protocol.md](protocol.md)：先了解样本、比较和筛选门槛。
2. [tutorial.md](tutorial.md)：看懂真实样本表、表达矩阵及候选结果。
3. [notebook.ipynb](notebook.ipynb)：逐格读取数据，理解索引与筛选条件。
4. [run.py](run.py)：完整复现导入、质控、差异分析、外部复核和作图。

## 输入

- `datasets/GSE252027/`：FB D0/D7 各 4 个文库用于主要比较，全部 24 个 FB 文库用于描述性检查。
- `datasets/GSE111059/`：未损伤/D7 各 3 个样本用于独立复核。
- `datasets/annotation/Mus_musculus.GRCm38.87.gtf.gz`：统一 Ensembl 基因、转录本及符号注释。

以上路径均相对于项目根目录。源数据恢复方式见[下载指南](../../docs/guides/data-download.md)。

## 运行

在项目根目录执行：

```powershell
.\.venv\Scripts\python.exe -X utf8 experiments\01_fibroblast_screening\run.py
```

Notebook 可从项目根目录或本实验目录运行，均使用项目 `.venv` 的 Python 内核。它读取已有结果，不重新拟合差异分析模型。

## 输出与首轮结果

输出目录：[outputs/01_fibroblast_screening](../../outputs/01_fibroblast_screening/)。文件只保存在本地，可运行脚本重建。

| 文件 | 内容 |
|---|---|
| `summary.json` | 筛选数量及注释映射统计 |
| `discovery_samples.csv`、`discovery_qc.csv` | FB 样本信息和数值检查 |
| `FB_TPM.tsv.gz`、`FB_NumReads.tsv.gz`、`FB_EffectiveLength.tsv.gz` | 按 GSM 对齐的完整 FB 表达矩阵 |
| `discovery_D7_vs_D0.tsv` | 发现队列差异分析 |
| `validation_D7_vs_D0_Welch.tsv` | 独立芯片探索性复核 |
| `all_evidence.tsv` | 合并证据及每个基因的筛选条件 |
| `cross_cohort_supported_candidates.tsv` | 通过当前全部门槛的候选 |
| `01_*.png`、`02_*.png`、`03_*.png` | 质控、候选时间曲线与 Cthrc1 教学图 |

首轮：17,622 个基因进入发现检验，1,501 个通过发现门槛，4 个通过跨队列门槛。Cthrc1、Ltbp2、Crlf1 上调，Gstm1 下调。候选仅代表探索性表达关联，芯片统计仍需 limma 复核。完整方法与解释边界以本目录的 protocol 和 tutorial 为准。
