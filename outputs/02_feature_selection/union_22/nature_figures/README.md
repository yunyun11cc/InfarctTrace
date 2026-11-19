# 22 文库并集的四张单图

四张图均由 Python/Matplotlib 直接读取已保存的结果。每张只有一个数据面板，统一宽 183 mm，输出 600 dpi PNG。

| 图 | 展示内容 | 源数据 |
|---|---|---|
| [01 差异效应散点图](01_candidate_union.png) | 每个点代表一个基因。横纵轴分别为 D3、D7 对 D0 的 Wald log₂FC。图例显示仅 D3 达标的 2,244 个、仅 D7 达标的 522 个，以及两天均达标的 1,141 个基因的调控方向。达标定义为 BH q < 0.05 且 &#124;log₂FC&#124; ≥ 1。 | [3,907 个基因的筛选标记和效应](../tables/feature_ranking.tsv) |
| [02 逐文库表达矩阵](02_candidate_expression.png) | 八个基因共用一套 VST 色标。六个时点内各列为单份文库，带边框的末列为该时点均值。右侧数字给出 D7 与非 D7 的 VST 均值差和八次扰动中的入选次数。 | [阶段 01_1 VST](../../../01_1_preprocessing_differential/matrices/timecourse_vst.tsv.gz) · [文库分组](../reports/selection_summary.json) · [基因排名和入选次数](../tables/feature_ranking.tsv) |
| [03 八次扰动的稳定性](03_selection_stability.png) | 横轴为相对全队列名单的 Jaccard，纵轴为实际入选基因数。每次分别移除一份 D3 和一份 D7 文库。圆点大小和倍数标记表示坐标相同的扰动次数。LASSO 和 Elastic Net 的点在 20 两侧略微错开，实际均为 20。并列基因保留。 | [四种方法的逐次 Jaccard 和名单大小](../reports/selection_summary.json) · [原 Bootstrap 稳定性图](../resampling/figures/03_selection_stability.png) |
| [04 通路相关性点阵](04_ssgsea_correlation.png) | 八个核心基因 × 七条通路。点面积随 &#124;ρ&#124; 增大，颜色区分 Spearman ρ 的正负。星号对应 19,999 次双侧置换的原始 P 值（* P < 0.05，** P < 0.01）。分析合并了全部 22 份文库，相关性可能反映共同的时序变化。 | [相关系数、P 值和 BH q 值](../tables/ssgsea_correlations.tsv) · [逐文库 ssGSEA NES](../tables/ssgsea_scores.tsv) |

在仓库根目录运行：

```sh
uv run --locked python experiments/02_feature_selection/nature_union_figures.py
```

[原始四张图](../figures) · [绘图脚本](../../../../experiments/02_feature_selection/nature_union_figures.py)
