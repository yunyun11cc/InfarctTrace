# InfarctTrace

InfarctTrace 使用小鼠心梗后的成纤维细胞时序转录组 **GSE252027**，筛选与炎症向纤维化转变相关的候选基因。计算分析交付候选排序和表达证据，靶向干预及纤维化表型验证由医学合作方完成。

## 从数据到结果

1. [下载数据](scripts/README.md) — 配置环境，运行下载脚本并校验源文件。
2. [数据说明](datasets/README.md) — 了解数据来源、采样分组、文件格式和分析注释。
3. [实验总览](experiments/README.md) — 查看各阶段的输入、交付物、执行状态和后续计划。
4. [输出索引](outputs/README.md) — 按实验查找矩阵、统计表、图片与报告。

## 项目结构

```text
InfarctTrace/
├── scripts/                                # 下载及源文件校验
│   └── download_data.py
├── datasets/                               # 分析输入
│   ├── GSE252027/
│   │   ├── metadata/                       # 样本表、分组与归档检查
│   │   └── original/                       # 下载后在本地生成的定量归档
│   └── annotation/                         # 本地基因与通路注释
├── experiments/                            # 分阶段组织分析代码和说明
│   ├── 01_preprocessing_differential/       # 质控、归一化、差异和富集
│   │   ├── prepare.py
│   │   ├── differential.py
│   │   ├── enrichment.py
│   │   ├── visualize.py
│   │   ├── protocol.md                     # 样本纳入与统计口径
│   │   └── results.md                      # 结果、读图与阶段结论
│   ├── 02_feature_selection/               # D3/D7 候选排序、稳定性分析与绘图
│   ├── 03_validation/
│   └── 04_prioritization/
├── outputs/                                # 与实验阶段同名的结果目录
│   ├── README.md                           # 图表和文件索引
│   ├── 01_preprocessing_differential/      # 表达矩阵、差异统计表、图片与报告
│   └── 02_feature_selection/               # 阶段 02 全部结果
│       ├── intersection_6/                # 原始六文库交集筛选
│       ├── intersection_22/               # 全部主分析文库的交集筛选
│       └── union_22/                       # 并集筛选与 ssGSEA 相关性
├── provenance/                             # 数据来源、下载清单与 SHA-256
├── papers/                                 # 论文元数据和参考文献
├── docs/                                   # 实验方案、数据指南和数学原理
│   ├── guides/
│   └── math/
├── notebooks/                              # 教学笔记与数据阅读示例
├── pyproject.toml                          # Python 版本要求与依赖声明
└── uv.lock                                 # 锁定的 Python 依赖版本
```

上图列出主要目录和入口。数据下载脚本在本地准备源文件，实验模块读取这些输入并生成同名阶段的输出。Git 保存代码、样本说明、下载校验、实验输出和结果解读。

## 实验框架

**数据下载 → 01 数据处理与差异分析 → 02 特征筛选 → 03 队列内验证 → 04 候选优选**

| 阶段 | 模块入口 | 交付物 | 状态 |
|---|---|---|---|
| 01 | [数据处理与差异分析](experiments/01_preprocessing_differential/README.md) | 样本质控、归一化矩阵、时序差异和通路富集 | Completed |
| 02 | [特征筛选](experiments/02_feature_selection/README.md) | D3/D7 差异基因交集、候选排序与文库扰动稳定性 | Completed |
| 03 | [候选验证](experiments/03_validation/README.md) | GSE252027 队列内的稳定性与时间趋势证据 | Pending |
| 04 | [候选优选](experiments/04_prioritization/README.md) | 候选排序、入选依据和证据摘要 | Pending |

## 阶段 02 结果

三版结果集中在同一目录。22 文库并集版另附前八基因与七条炎症及纤维化通路的 ssGSEA 相关性。 [阶段 02 输出索引](outputs/02_feature_selection/README.md) · [ssGSEA 热图](outputs/02_feature_selection/union_22/figures/04_ssgsea_correlation.png)

22 文库并集版还提供 4 折 × 50 次分层交叉验证，以及 1,000 次分层 Bootstrap 的筛选稳定性。 [交叉验证与 Bootstrap 结果](outputs/02_feature_selection/README.md#22-文库并集的重复交叉验证与-bootstrap)

22 文库版比较 D7 的 4 份文库与其他五个时点的 18 份文库。以 D3 vs D0 和 D7 vs D0 的差异基因生成候选池，**交集含 1,141 个基因，并集含 3,907 个基因**。分别有 15 个和 14 个基因在全部 8 次扰动中进入 VST 均值差绝对值前 20。原始六文库版比较 D7 与 D3，有 8 个基因在全部扰动中进入前 20。 [版本说明与输出](experiments/02_feature_selection/README.md) · [22 文库交集排序](outputs/02_feature_selection/intersection_22/tables/feature_ranking.tsv) · [22 文库并集排序](outputs/02_feature_selection/union_22/tables/feature_ranking.tsv) · [六文库结果解读](experiments/02_feature_selection/results.md)

- [差异基因交集](outputs/02_feature_selection/intersection_6/figures/01_candidate_intersection.png)
- [前八基因的逐文库表达](outputs/02_feature_selection/intersection_6/figures/02_candidate_expression.png)
- [四种方法的筛选稳定性](outputs/02_feature_selection/intersection_6/figures/03_selection_stability.png)
- [完整基因排序](outputs/02_feature_selection/intersection_6/tables/feature_ranking.tsv) · [八次扰动明细](outputs/02_feature_selection/intersection_6/tables/perturbation_selections.tsv) · [运行摘要](outputs/02_feature_selection/intersection_6/reports/selection_summary.json)

上方图表为原始六文库版。已有阶段 01 输出时，可在项目根目录重跑 22 文库交集版。两个命令都加 `--union` 可复跑并集版。绘图需要中文字体。

```sh
uv run --locked python experiments/02_feature_selection/select_features.py
uv run --locked python experiments/02_feature_selection/visualize.py
```

## 参考与学习

- [项目文档](docs/README.md) — 原始实验方案、数据来源和数学原理。
- [教学笔记](notebooks/README.md) — 阅读样本表、表达矩阵和实验图表。
