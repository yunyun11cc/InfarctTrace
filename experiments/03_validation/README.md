# 阶段 03 冻结候选的外部验证

在两套独立的心梗数据中复核表达与判别表现，并在单细胞数据中描述细胞来源。候选使用阶段 02 冻结的九个 Ensembl ID 与基因符号。 [冻结清单](../../outputs/02_feature_selection/manuscript/table_02_candidate_genes_summary.tsv) · [外部队列清洗](../01_2_validation_preprocessing/README.md)

| 数据 | 主要比较 | 输入与取材 |
|---|---|---|
| GSE111059 | D7 对未损伤组，各 3 份 | 作者已标准化的 SST-RMA 芯片信号，样本为 TCF21 谱系分选的心脏成纤维细胞 |
| GSE236374 | D7 对 Sham，各 3 份 | FeatureCounts 原始计数换算为 log2(CPM+1)，样本为心脏组织 |
| GSE163465 | Sham、D3、D7，每组一份公开文库 | 心脏 CD45⁺ 免疫细胞富集。三份未给药 10x 矩阵 |

[GSE111059 原研究](https://www.jci.org/articles/view/98215) · [GSE236374 原研究](https://doi.org/10.1073/pnas.2302482120) · [GSE163465 原研究](https://doi.org/10.1002/smtd.202100752)

## 方法与判定

1. **差异表达**。D7 是两外部队列共同采样日，也是阶段 02 的目标时点。晋级要求两队列的双侧 Welch 原始 P < 0.05 且 log₂ 尺度均值差为正。芯片均值差来自已归一化的 log₂ 信号，RNA-seq 来自 log₂(CPM+1)，单位差异需在解读时保留。每个对比对九基因另做 BH 校正。其他时点与合并 MI 相对本队列对照的比较作为补充，单独列在结果表。
2. **判别表现**。单基因分数直接用外部表达值，方向预先固定为表达越高越偏向 MI。95% 区间用 10,000 次分层百分位 Bootstrap。每组仅 3 份，完全分离的样本会产生 [1, 1] 的退化区间，区间无法反映新样本的不确定性。联合模型用冻结九基因，L2 逻辑回归 C=1，标准化与拟合均在各留一训练折内完成。另列晋级七基因子集的留一结果，因该子集由同一批外部数据判定，只作筛后描述。发现集 AUROC 是同队列 D7 对 D0 描述量。
3. **共表达及通路得分**。两队列各自合并其全部采样日。报告九基因与六个标志物的 Spearman ρ、双侧 P 及每个队列的 BH q。使用阶段 02 固定的 Reactome 小鼠 NLRP3 inflammasome、Extracellular matrix organization 和 Collagen formation 成员，在外部全基因矩阵计算 ssGSEA NES，再与九基因求相关。相关性可能部分来自共同的时间或组别变化，通路成员重叠也可能影响得分相关。
4. **单细胞定位**。每细胞保留 300–6,500 个检出基因且线粒体 UMI 比例 < 20%。归一到每万 UMI 后 log1p；按 GSM 选取高变基因并排除九候选和线粒体基因，用 30 个 PC、15 个近邻、Leiden 0.65 和 UMAP。宽类群按独立标志物审查。报告每个公开文库内各类群的检出率和含零平均表达。每时点一份公开文库，跨时点 P 值和 AUROC 留空。CD45⁺ 富集的取样范围限制了成纤维细胞定位及细胞间互作推断。 [单细胞生物学重复与推断](https://doi.org/10.1038/s41467-021-25960-2)

## 运行

先按阶段 01_2 准备外部矩阵。阶段 01_1 的小鼠基因注释缓存位于 \`.tools/org.Mm.eg.sqlite\`，GSE163465 使用 \`datasets/GSE163465/original/GSE163465_RAW.tar\`。从项目根目录执行。

\`\`\`sh
uv run --locked python experiments/03_validation/bulk.py
uv run --locked --with 'scanpy[leiden]==1.11.5' python experiments/03_validation/single_cell.py
uv run --locked python experiments/03_validation/summarize.py
\`\`\`

[结果和交接文件](../../outputs/03_validation/README.md)
