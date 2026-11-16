# 阶段 02 特征筛选

**状态：Completed（固定候选池的探索性筛选）**

用阶段 01 的 D3 对 D0 与 D7 对 D0 的差异筛选结果，以 Ensembl ID 建立候选池。项目差异门槛为 `padj < 0.05` 且 `|log2FC| ≥ 1`。交集要求两次对比均达标，共 1,141 个基因。并集收纳至少一次对比达标的基因，共 3,907 个。三版均读取阶段 01 的 Wald 表、主分析文库纳入表和 VST 矩阵。 [阶段 01 入口](../01_preprocessing_differential/README.md)

## 运行

在项目根目录使用项目 Python 环境。

```sh
uv run --locked python experiments/02_feature_selection/select_features.py
uv run --locked python experiments/02_feature_selection/visualize.py
uv run --locked python -m unittest discover -s experiments/02_feature_selection -v
```

默认使用 **22 × 1,141** 的交集输入，并写入 `outputs/02_feature_selection/intersection_22/`。用 `--union` 运行 **22 × 3,907** 的并集版，结果写入 `outputs/02_feature_selection/union_22/`。

```sh
uv run --locked python experiments/02_feature_selection/select_features.py --union
uv run --locked python experiments/02_feature_selection/visualize.py --union
```

原始六文库结果保存在 `outputs/02_feature_selection/intersection_6/`。可用以下命令将六文库版复跑到另一个目录。

```sh
uv run --locked python experiments/02_feature_selection/select_features.py --d3-d7-only --output-dir outputs/02_feature_selection/intersection_6_rerun
uv run --locked python experiments/02_feature_selection/visualize.py --output-dir outputs/02_feature_selection/intersection_6_rerun
```

绘图需要中文字体，例如 Noto Sans CJK SC 或 macOS 的 PingFang SC。

## 筛选与稳定性

1. `fit_rankers()` 以 D7 为正类，以其余输入文库为参照。22 文库版比较 **D7 的 4 份 vs 非 D7 的 18 份**，六文库版比较 **D7 的 4 份 vs D3 的 2 份**。组内按文库求 VST 均值，绝对差用于排序。
2. 在同一基因池内拟合三种探索性模型。LASSO 和 Elastic Net 均使用带类别平衡权重的逻辑回归，每次拟合独立标准化特征，取系数绝对值作为分数。固定 `C=1`，`l1_ratio` 分别为 `1` 和 `0.5`。随机森林使用 512 棵树及 Gini 不纯度下降的重要性。随机种子固定为 `252027`。
3. 遍历 8 种各移除一份 D3 和 D7 文库的组合。22 文库版每次在剩余 20 份文库上拟合，六文库版每次使用剩余 4 份文库。每次按分数绝对值选前 20，保留边界并列，排除零分数。入选频率为入选次数除以 8。方向一致率记录扰动后的均值差与全样本均值差同向的比例。主表按均值差入选频率降序，再按全样本均值差绝对值降序，以基因 ID 打破并列。

候选池与 VST 均来自阶段 01 的全队列。八种扰动组合重用同一批文库，因此频率只描述**固定候选池内的队列内稳定性**。表中的 VST 均值差是描述性效应量。 [特征筛选与评估偏差](https://pmc.ncbi.nlm.nih.gov/articles/PMC124442/)

## 22 文库并集输出

并集包含仅 D3 达标的 2,244 个、两个时点均达标的 1,141 个，以及仅 D7 达标的 522 个基因。表中 `D3通过差异筛选` 和 `D7通过差异筛选` 记录各自的门槛判定。例如 `仅D3上调` 表示 D3 通过门槛且方向为上调。两个时点的 Wald 效应估计均可在表中查看。

按均值差排序，14 个基因在全部 8 次扰动中进入前 20，均属于原交集。稳定入选的前四是 `C1qtnf3`、`Cthrc1`、`Kif26b` 和 `Col24a1`。

- [完整基因排序](../../outputs/02_feature_selection/union_22/tables/feature_ranking.tsv) · [八次扰动明细](../../outputs/02_feature_selection/union_22/tables/perturbation_selections.tsv) · [运行摘要](../../outputs/02_feature_selection/union_22/reports/selection_summary.json)
- [差异基因并集](../../outputs/02_feature_selection/union_22/figures/01_candidate_union.png) · [六时点逐文库表达](../../outputs/02_feature_selection/union_22/figures/02_candidate_expression.png) · [筛选稳定性](../../outputs/02_feature_selection/union_22/figures/03_selection_stability.png)

## 核心基因与通路 ssGSEA 相关性

从 `union_22` 主排序取前八基因。对阶段 01 的 **19,312 个 VST 基因 × 22 份文库**逐文库打分，使用 GSEApy `ssgsea` 的 NES。基因集为 KEGG 小鼠 NF-κB、TNF、细胞因子–受体相互作用和 TGF-β，以及 Reactome 小鼠 NLRP3 炎症小体、细胞外基质组织和胶原形成。使用阶段 01 的一对一 Entrez→Ensembl 映射。每条通路匹配到 9–213 个基因，完整成员和下载来源随代码保存。 [GSEApy](https://doi.org/10.1093/bioinformatics/btac757) · [KEGG 小鼠通路](https://www.kegg.jp/kegg-bin/show_organism?org=mmu) · [Reactome 下载与映射](https://reactome.org/download-data/)

每格为基因表达与通路得分在 22 份文库上的 Spearman ρ。星号用 19,999 次双侧文库配对置换的原始 P 值，`*` 为 P < 0.05，`**` 为 P < 0.01。对全部 56 对检验计算 BH 校正值。`Rho_without_core` 列来自移除全部八个核心基因后重新打分的敏感性分析。候选筛选和相关分析使用同一批文库，P 和 q 用于探索性解读。这里汇总了六个时间点，相关性可能反映共同的时序变化。 [Spearman 置换检验说明](https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.stats.spearmanr.html)

```sh
uv run --locked python experiments/02_feature_selection/pathway_activity.py
```

[热图、打分矩阵及统计表](../../outputs/02_feature_selection/README.md#22-文库并集的-ssgsea)

## 重复分层交叉验证与 Bootstrap

22 文库并集版使用固定的 3,907 个基因，D7 为正类。采用 4 折分层交叉验证，重复 50 次。每折验证集有 1 个 D7 文库，训练集有 3 个 D7 和 13–14 个非 D7 文库。每折只用训练文库拟合标准化和模型，从训练模型取前 20 个特征，再用这些特征重新拟合后为验证文库打分。LASSO 和 Elastic Net 按绝对系数排序，随机森林按特征重要性排序。保留边界并列并排除零分数。汇报 200 折 AUROC、AUPRC（Average Precision）和准确率的均值及样本标准差。CV 入选频率为该基因在相应方法的 200 个训练折中入选的比例。

Bootstrap 对 D7 和非 D7 分层，每次分别有放回抽 4 和 18 份文库，重复 1,000 次。每次重新运行 VST 均值差、LASSO、Elastic Net 和随机森林。按每种方法统计入选频率，综合选中概率取四种方法频率的算术平均。箱线图展示每次抽样的名单与对应全样本名单的 Jaccard，相邻的棒棒糖图展示综合平均概率最高的 20 个基因。

候选池和阶段 01 的 VST 都由全队列生成。交叉验证表反映固定候选池和输入矩阵下的队列内表现。不同折重复使用同一批文库，因此折间标准差用于描述波动。 [预筛选对交叉验证的影响](https://pmc.ncbi.nlm.nih.gov/articles/PMC124442/)

```sh
uv run --locked python -u experiments/02_feature_selection/resampling.py
uv run --locked python experiments/02_feature_selection/visualize_resampling.py
```

[性能表、概率表和稳定性图](../../outputs/02_feature_selection/README.md#22-文库并集的重复交叉验证与-bootstrap)

## 22 文库交集输出

均值差排序中有 15 个基因在每次扰动中都进入前 20。排名前四为 `C1qtnf3`、`Cthrc1`、`Kif26b` 和 `Col24a1`。

- [完整基因排序](../../outputs/02_feature_selection/intersection_22/tables/feature_ranking.tsv) · [八次扰动明细](../../outputs/02_feature_selection/intersection_22/tables/perturbation_selections.tsv) · [运行摘要](../../outputs/02_feature_selection/intersection_22/reports/selection_summary.json)
- [差异基因交集](../../outputs/02_feature_selection/intersection_22/figures/01_candidate_intersection.png) · [六时点逐文库表达](../../outputs/02_feature_selection/intersection_22/figures/02_candidate_expression.png) · [筛选稳定性](../../outputs/02_feature_selection/intersection_22/figures/03_selection_stability.png)

## 六文库输出

| 文件 | 内容 |
|---|---|
| [完整基因排序](../../outputs/02_feature_selection/intersection_6/tables/feature_ranking.tsv) | 1,141 个基因的排名、两个相对 D0 的 Wald 结果、D7−D3 的 VST 均值差、方向一致率和四种方法的分数与频率 |
| [八次扰动明细](../../outputs/02_feature_selection/intersection_6/tables/perturbation_selections.tsv) | 每个文库扰动组合下各方法入选的基因及分数 |
| [运行摘要](../../outputs/02_feature_selection/intersection_6/reports/selection_summary.json) | 实际文库编号、各模型参数、每次入选数量以及与全样本名单的 Jaccard 指数 |
| [差异基因交集](../../outputs/02_feature_selection/intersection_6/figures/01_candidate_intersection.png) | 显示交集规模和两个相对 D0 的表达变化方向 |
| [前八基因逐文库表达](../../outputs/02_feature_selection/intersection_6/figures/02_candidate_expression.png) | 显示 D3、D7 的文库表达值与组均值 |
| [筛选稳定性](../../outputs/02_feature_selection/intersection_6/figures/03_selection_stability.png) | 比较四种方法在八次扰动中的名单重合度和实际入选数 |

[阶段 02 结果解读](results.md)
