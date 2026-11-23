# InfarctTrace

InfarctTrace 使用小鼠心梗后的成纤维细胞时序转录组 **GSE252027**，筛选与炎症向纤维化转变相关的候选基因。独立成纤维细胞芯片 **GSE111059** 和心脏组织 RNA-seq **GSE236374** 用于复核表达方向。计算分析交付候选排序和表达证据，靶向干预及纤维化表型验证由医学合作方完成。

## 从数据到结果

1. [下载数据](scripts/README.md) — 配置环境，运行下载脚本并校验源文件。
2. [数据说明](datasets/README.md) — 了解数据来源、采样分组、文件格式和分析注释。
3. [实验总览](experiments/README.md) — 查看各阶段的输入、交付物、执行状态和后续计划。
4. [输出索引](docs/outputs.md) — 按实验查找矩阵、统计表、图片与报告。

## 项目结构

```text
InfarctTrace/
├── scripts/                                # 下载及源文件校验
│   └── download_data.py
├── datasets/                               # 分析输入
│   ├── GSE252027/
│   │   ├── metadata/                       # 样本表、分组与归档检查
│   │   └── original/                       # 下载后在本地生成的定量归档
│   ├── GSE111059/                          # 成纤维细胞芯片验证集
│   ├── GSE236374/                          # 心脏组织 RNA-seq 验证集
│   └── annotation/                         # 本地基因与通路注释
├── experiments/                            # 分阶段组织分析代码和说明
│   ├── 01_preprocessing_differential/      # 质控、标准化、差异、候选窗口与富集
│   │   ├── main.py
│   │   └── docs/
│   ├── 02_feature_selection/               # D3/D7 筛选、重采样和生物学互证
│   │   ├── main.py
│   │   ├── methods/
│   │   ├── bio_validation/
│   │   └── docs/
│   ├── 03_external_validation/             # 外部队列清洗、评估及单细胞补充分析
│   │   ├── main.py
│   │   └── docs/
│   └── 04_prioritization/                  # 待办说明
├── outputs/                                # 各阶段的工程产物
│   ├── 01_preprocessing_differential/      # 矩阵、差异统计、富集与质控图
│   ├── 02_feature_selection/               # 两种交集模式与 22 文库并集模式
│   │   ├── intersection_6/
│   │   ├── intersection_22/
│   │   └── union_22/                       # 筛选、重采样、EFS 聚合和生物学互证
│   └── 03_external_validation/
│       ├── preprocessing/                 # 外部队列清洗
│       ├── evaluation/                    # 表达与 ROC 评估
│       ├── single_cell/                   # 单细胞定位
│       └── summary/                       # 跨队列汇总
├── provenance/                             # 数据来源、下载清单与 SHA-256
├── papers/                                 # 论文元数据和参考文献
├── docs/                                   # 实验方案、数据指南和数学原理
│   ├── guides/
│   ├── math/
│   └── outputs.md                          # 产物索引
├── pyproject.toml                          # Python 版本要求与依赖声明
└── uv.lock                                 # 锁定的 Python 依赖版本
```

上图列出主要目录和入口。数据下载脚本在本地准备源文件，实验模块读取这些输入并写入各阶段配置的输出目录。Git 保存代码、样本说明、下载校验、实验输出和结果解读。

## 实验框架

**01 发现队列预处理与差异分析 → 02 特征筛选 → 03 外部队列验证**。候选优选为待办项。

| 阶段 | 模块入口 | 交付物 | 状态 |
|---|---|---|---|
| 01 | [发现队列数据处理与差异分析](experiments/01_preprocessing_differential/docs/README.md) | 样本质控、归一化矩阵、时序差异和通路富集 | Completed |
| 02 | [特征筛选](experiments/02_feature_selection/docs/README.md) | D3/D7 交集和并集排序、ssGSEA 相关性、交叉验证与 Bootstrap 稳定性 | Completed |
| 03 | [外部队列验证](experiments/03_external_validation/docs/README.md) | 两套外部队列的清洗、ROC、通路和共表达；可选单细胞定位 | Completed |
| 04 | [候选优选](experiments/04_prioritization/README.md) | 候选排序、入选依据和证据摘要 | Pending |

## 阶段 02 结果

22 文库并集版使用 D3 vs D0 和 D7 vs D0 的 **3,907 个差异基因**，比较 D7 的 4 份文库与其余五个时点的 18 份文库。交叉验证采用 4 折 × 50 次，每折在训练文库中筛选后重新拟合。Elastic Net 的平均 AUROC 为 **0.957 ± 0.108**。候选池与 VST 来自全队列，这组数值描述固定输入下的队列内判别表现。 [阶段 02 输出索引](experiments/02_feature_selection/docs/results.md) · [逐折性能汇总](outputs/02_feature_selection/union_22/resampling/tables/cv_summary.tsv)

| 分析 | 输出 |
|---|---|
| 并集特征筛选 | [3,907 个基因的排序](outputs/02_feature_selection/union_22/tables/feature_ranking.tsv) · [八次文库扰动明细](outputs/02_feature_selection/union_22/tables/perturbation_selections.tsv) |
| ssGSEA | [前八基因与七条通路的相关性热图](outputs/02_feature_selection/union_22/bio_validation/figures/ssgsea_top_candidates.png) · [相关系数、P 值和校正值](outputs/02_feature_selection/union_22/bio_validation/tables/ssgsea_correlations.tsv) |
| 200 折交叉验证 | [每折性能](outputs/02_feature_selection/union_22/resampling/tables/cv_fold_metrics.tsv) · [各方法的基因入选频率](outputs/02_feature_selection/union_22/resampling/tables/cv_selection_frequency.tsv) |
| 1,000 次分层 Bootstrap | [Jaccard 箱线图和棒棒糖图](outputs/02_feature_selection/union_22/resampling/figures/03_selection_stability.png) · [逐次 Jaccard](outputs/02_feature_selection/union_22/resampling/tables/bootstrap_jaccard.tsv) · [各方法及综合入选频率](outputs/02_feature_selection/union_22/resampling/tables/bootstrap_selection_probability.tsv) |

交集版使用两个时点都达标的 **1,141 个基因**。六文库版比较 D7 与 D3，22 文库版比较 D7 与其余时点。 [22 文库交集排序](outputs/02_feature_selection/intersection_22/tables/feature_ranking.tsv) · [六文库图表与结果解读](experiments/02_feature_selection/docs/results.md)

### 复现 22 文库并集版

在项目根目录运行以下命令。第一步读取已有的阶段 01 输出。若已有并集排序表，可直接运行所需的下游脚本。 [阶段 02 运行说明](experiments/02_feature_selection/docs/README.md)

```sh
uv run --locked python experiments/02_feature_selection/main.py --steps selection resampling ssgsea figures
```

## 参考与学习

- [项目文档](docs/README.md) — 原始实验方案、数据来源和数学原理。
- [阶段 01 导读](experiments/01_preprocessing_differential/docs/introduction.md) — 阅读样本表、表达矩阵和实验图表。
