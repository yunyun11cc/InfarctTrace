# 阶段 02 特征筛选

**状态：Completed（固定候选池的探索性筛选）**

取阶段 01 中 D3 对 D0 与 D7 对 D0 都通过差异筛选的基因，以 Ensembl ID 求交集。项目差异门槛为 `padj < 0.05` 且 `|log2FC| ≥ 1`。交集内保留同向和反向基因。筛选目标是比较 D3 和 D7 的表达模式，观察基因排序对文库取舍的敏感性。输入是阶段 01 的 Wald 表、主分析文库纳入表和 VST 矩阵。 [阶段 01 入口](../01_preprocessing_differential/README.md)

## 运行

在项目根目录使用项目 Python 环境。

```sh
uv run --locked python experiments/02_feature_selection/select_features.py
uv run --locked python experiments/02_feature_selection/visualize.py
uv run --locked python -m unittest discover -s experiments/02_feature_selection -v
```

默认将文件写入 `outputs/02_feature_selection/`，也可用 `--output-dir` 指定目录。绘图需要中文字体，例如 Noto Sans CJK SC 或 macOS 的 PingFang SC。

## 筛选与稳定性

1. 在两份 D3 和四份 D7 文库上计算 **D7 减 D3 的 VST 均值差**。均值差为正表示 D7 表达较高，绝对值用于排序。
2. 在同一基因池内拟合三种探索性模型。LASSO 和 Elastic Net 均使用带类别平衡权重的逻辑回归，每次拟合独立标准化特征，取系数绝对值作为分数。固定 `C=1`，`l1_ratio` 分别为 `1` 和 `0.5`。随机森林使用 512 棵树及 Gini 不纯度下降的重要性。随机种子固定为 `252027`。
3. 遍历 8 种各移除一份 D3 和 D7 文库的组合，每次用剩下的一份 D3 和三份 D7 重新计算各方法分数。每次按分数绝对值选前 20，保留边界并列，排除零分数。入选频率为入选次数除以 8。方向一致率记录扰动后的均值差与全样本均值差同向的比例。主表按均值差入选频率降序，再按全样本均值差绝对值降序，以基因 ID 打破并列。

候选池与 VST 均来自阶段 01 的全队列。八种扰动组合重用同一批文库，因此频率只描述**固定候选池内的队列内稳定性**。表中的 VST 均值差是描述性效应量。 [特征筛选与评估偏差](https://pmc.ncbi.nlm.nih.gov/articles/PMC124442/)

## 输出

| 文件 | 内容 |
|---|---|
| [完整基因排序](../../outputs/02_feature_selection/tables/feature_ranking.tsv) | 1,141 个基因的排名、两个相对 D0 的 Wald 结果、D7−D3 的 VST 均值差、方向一致率和四种方法的分数与频率 |
| [八次扰动明细](../../outputs/02_feature_selection/tables/perturbation_selections.tsv) | 每个文库扰动组合下各方法入选的基因及分数 |
| [运行摘要](../../outputs/02_feature_selection/reports/selection_summary.json) | 实际文库编号、各模型参数、每次入选数量以及与全样本名单的 Jaccard 指数 |
| [差异基因交集](../../outputs/02_feature_selection/figures/01_candidate_intersection.png) | 显示交集规模和两个相对 D0 的表达变化方向 |
| [前八基因逐文库表达](../../outputs/02_feature_selection/figures/02_candidate_expression.png) | 显示 D3、D7 的文库表达值与组均值 |
| [筛选稳定性](../../outputs/02_feature_selection/figures/03_selection_stability.png) | 比较四种方法在八次扰动中的名单重合度和实际入选数 |

[阶段 02 结果解读](results.md)
