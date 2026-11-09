# InfarctTrace

**Time-course Transcriptomics of Post-Infarction Cardiac Fibroblasts**

本项目分析小鼠心肌梗死后成纤维细胞的表达变化。自 2026-11-10 起，**数据处理、特征筛选、验证和优选四阶段均只使用 GSE252027**。主要研究对象是成纤维细胞（FB），观察时间覆盖 D0、D1、D3、D7、D14、D28。阶段 03 的验证指 GSE252027 队列内部的验证。

## 从哪里开始

- [GSE252027 数据手册](datasets/GSE252027/README.md) 解释样本、原始文件、三张元数据 CSV、表达量和筛选状态。
- [实验总览](experiments/README.md) 说明四个工程阶段的输入与交付物。
- [教学笔记](notebooks/README.md) 帮助阅读数据、理解筛选基线和解释图表。
- [输出目录索引](outputs/README.md) 说明矩阵、表格、图和报告的位置。
- [数据恢复指南](docs/guides/data-download.md) 记录源文件与校验入口。

## 项目结构

```text
InfarctTrace/
├── datasets/
│   ├── GSE252027/                # 主队列的定量归档、元数据与详细手册
│   └── annotation/              # 小鼠基因注释
├── experiments/
│   ├── 01_preprocessing_differential/  # 数据处理与差异分析
│   ├── 02_feature_selection/          # 特征筛选
│   ├── 03_validation/                 # 主队列内部验证
│   └── 04_prioritization/             # 候选基因优选
├── notebooks/                   # 教学说明和 Jupyter 笔记
├── outputs/                     # 与四个实验阶段同名，分 matrices、tables、figures、reports
├── scripts/                     # 主队列下载和源文件校验
├── docs/                        # 工程文档和数据恢复指南
├── papers/                      # 主队列的论文和引用
├── provenance/                  # 主队列的下载记录和历史校验
├── pyproject.toml               # Python 版本和依赖
└── uv.lock                      # uv 锁定版本
```

## 当前基线与运行方式

主队列已整理 24 个 FB 提交文库、51,029 个基因的 TPM、估计片段数和有效长度矩阵。阶段 01 已完成 22 文库主分析、D28r4 排除敏感性分析和五个时间对比的 GO、KEGG、GSEA。19,312 个基因进入全时段模型，11,489 个达到总体时间效应门槛。阶段 02–04 待开展。 [结果解读](experiments/01_preprocessing_differential/results.md)

阶段 01 的正式分析按 [模块运行说明](experiments/01_preprocessing_differential/README.md) 执行。原 D7 对 D0 的 Python 教学基线仍可通过下面的入口复现。

```powershell
uv sync --locked
uv run --locked python -X utf8 scripts\download_data.py
uv run --locked python -X utf8 experiments\01_preprocessing_differential\run.py
uv run --locked jupyter lab
```

项目使用 `.python-version` 和 `uv.lock`，教学笔记选择项目 `.venv` 的 Python 内核。下载脚本只恢复或校验 GSE252027 及其分析注释。分析脚本生成中文表格和带中文图题、坐标及图注的图片，系统需要有可用的中文字体。源归档和生成结果保存在本地，Git 跟踪源说明、代码、教学资料和输出目录索引。

## 后续计划

当前四阶段均限定 GSE252027。完成主队列分析后，以下三个队列可作为**另行启动的拓展计划**，具体比较和启用时间另行确定。

| 暂缓队列 | 后续可回答的问题 |
|---|---|
| GSE111059 | 在独立成纤维细胞芯片研究的共有时间点复核表达方向 |
| GSE236374 | 检查候选信号在 D7 和 D28 的全心组织层面是否同向 |
| GSE163465 | 查看候选在心脏免疫细胞及巨噬细胞亚群中的表达分布 |
