# InfarctTrace

**Cross-Cohort Transcriptomics of Post-Infarction Cardiac Repair**

心梗后修复的跨队列转录组研究。以公开小鼠心梗数据为起点，研究表达变化如何随时间和细胞类型变化，以及这些信号能否在独立队列中复现。当前重点是计算机方法选题预研与可重复的分析基线，允许根据研究问题补充公开数据。

## 从哪里开始

- **学习当前实验：** [实验 01 入口](experiments/01_fibroblast_screening/README.md) → [中文教程](experiments/01_fibroblast_screening/tutorial.md) → [可运行 notebook](experiments/01_fibroblast_screening/notebook.ipynb)。
- **理解筛选依据：** [实验设计与筛选规则](experiments/01_fibroblast_screening/protocol.md)。
- **了解研究方向：** [选题初筛](docs/research/topic-screening-2026-11-05.md)。其中的推荐仍是待验证假设，最终方法题目尚未确定。
- **恢复公开数据：** [下载指南](docs/guides/data-download.md)与[数据及论文参考](docs/guides/dataset-reference.md)。

## 项目结构

```text
InfarctTrace/                         # 项目名称；当前工作区目录为 mi-dna-choice
├── README.md                        # 总入口与当前进度
├── datasets/                        # 按 GEO 编号保存源数据、样本表和说明
│   ├── GSE252027/                    # 发现队列：重点使用 FB 成纤维细胞
│   ├── GSE111059/                    # 外部复核：成纤维细胞芯片
│   ├── GSE236374/                    # 补充证据：全心组织 RNA-seq
│   ├── GSE163465/                    # 补充证据：免疫细胞单细胞 RNA-seq
│   └── annotation/                   # 基因及转录本注释源文件
├── experiments/
│   ├── README.md                    # 实验目录与编号约定
│   └── 01_fibroblast_screening/
│       ├── README.md                # 本实验的输入、运行和结果入口
│       ├── protocol.md              # 先于分析确定的设计与筛选规则
│       ├── tutorial.md              # 中文数据阅读与结果讲解
│       ├── notebook.ipynb           # 分步运行的教学笔记
│       └── run.py                   # 完整分析入口
├── scripts/
│   └── download_data.py             # 公共数据下载与完整性校验
├── docs/
│   ├── README.md                    # 文档导航
│   ├── guides/                      # 下载指南和数据参考
│   └── research/                    # 选题与历史可行性论证
├── papers/                          # 论文引用、元数据与本地全文
├── provenance/                      # 下载来源、校验值和操作记录
├── requirements/
│   ├── analysis.txt                 # 直接分析与 notebook 依赖
│   └── analysis-lock.txt            # 当前环境的完整依赖快照
└── outputs/                         # 本地生成结果，不进入 Git
    └── 01_fibroblast_screening/       # 表达矩阵、候选表、QC、图与 summary.json
```

每个实验的设计、代码和学习材料集中存放；输出使用同名目录。源数据按队列保存，公共下载脚本由各实验共用。目前只有一个分析实验，暂不引入额外的包结构或通用训练框架。

## 运行

在项目根目录使用 PowerShell。新环境先准备 Python 3.11 与依赖：

```powershell
uv venv --python 3.11 .venv
uv pip install --python .venv\Scripts\python.exe -r requirements\analysis.txt
```

下载数据并运行当前实验：

```powershell
.\.venv\Scripts\python.exe -X utf8 scripts\download_data.py
.\.venv\Scripts\python.exe -X utf8 experiments\01_fibroblast_screening\run.py
```

已有本地环境可直接运行。下载脚本核对已有文件并保留历史清单；分析脚本重建本实验输出。Notebook 使用项目 `.venv` 中的 Python 内核。

## 当前结果与范围

实验 01 比较 GSE252027 FB 的 D7 与 D0（各 4 个提交文库），再用 GSE111059 的同时间比较（各 3 个样本）复核。首轮结果为：**17,622 个基因进入发现检验 → 1,501 个通过发现门槛 → 4 个通过预设跨队列门槛**。Cthrc1、Ltbp2、Crlf1 上调，Gstm1 下调。

这是探索性表达关联，尚未建立预测模型或证明因果与药物直接作用。芯片复核本轮使用 Welch 检验，后续需以 limma 复核。四个队列的细胞类型、平台和条件不同，按队列分别分析；所有当前队列均为小鼠数据。

## 数据与可追溯性

源表达文件、原始芯片、压缩包、论文全文和批量结果保存在本地，由 `.gitignore` 排除。代码、教学 notebook、样本信息、研究文档和校验记录进入版本管理。

2026-11-08 已校验 14 个分析所需源文件，见[下载校验记录](provenance/download_verification_2026-11-08.json)和[归档完整性检查](provenance/archive_integrity_2026-11-08.json)。历史[下载清单](provenance/download_manifest.csv)还包含未在当前副本恢复的论文文件，因此不能把清单状态直接当作当前文件是否存在的证明。原始数据和历史清单的路径在本次整理中保持一致。
