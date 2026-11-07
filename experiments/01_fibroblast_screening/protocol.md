# 实验 01：心梗后成纤维细胞的表达变化与跨队列复核

制定日期：2026-11-08。以下规则在运行差异分析、查看候选名单之前确定。本轮是探索性计算基线，不作为已经成立的方法创新。

## 问题与主要比较

主要问题：心梗后第 7 天，相对于未损伤对照，哪些小鼠成纤维细胞基因出现足够大且可信的表达变化，并能在另一项研究的相同时间比较中复现？

- 发现队列：GSE252027 的 FB D7 与 D0，各 4 个提交文库。每文库来源于三只小鼠池化，统计单位为提交文库；原文对 technical/biological replicate 的用词不一致，独立性仍是证据限制。
- 复核队列：GSE111059 的 D7 与 uninjured，各 3 个生物学重复标签样本。平台为芯片，按芯片表达数据单独分析。
- 选择 D7 是因为两队列有相同时间点，并避开尚有剔除争议的发现集 D3；这不是疾病严重度、药物疗效或未来结局预测。
- 所有 24 个 FB 文库用于数据结构、表达分布和 PCA 的描述性检查。D3 的两个异常重复及作者代码额外省略 D28 r4 的问题单列，不据本次候选结果决定删样本。

## 预处理与统计规则

1. 原始归档只读，下载核对历史字节数和 SHA-256。基因 ID 与样本 ID 均保留。
2. 发现集仅保留 ENSMUSG 基因行。混入的 ENSMUST 转录本行单独记录，不与基因行叠加，避免重复计数。将原论文指定的 Ensembl GRCm38.87 用于稳定 ID、转录本 ID 与符号注释；本轮不重新定量。
3. 导入 Salmon 的 TPM、EffectiveLength、NumReads。按 tximport 官方 `lengthScaledTPM` 算法，将 TPM 乘各基因跨所分析样本的平均有效长度，再按原估计计数的列总量缩放，四舍五入后用于 PyDESeq2。该步骤不是把原始 TPM 直接当计数。均值长度只在 D0/D7 的 8 个文库中计算；该路线仍需日后用 R tximport/DESeq2 复核实现一致性。
4. 在发现集 8 个文库中，保留至少 4 个文库具有转换计数 ≥10 的基因。按 `~ condition` 拟合负二项模型，进行 D7 对 D0 的 Wald 检验和 BH 多重检验校正。不因同为 replicate 1 就假定跨时间配对。
5. 发现候选门槛：`padj < 0.05` 且 `abs(log2FoldChange) >= 1`。这是点估计效应量门槛，不是已证明真实效应超过两倍的阈值检验；LFC 暂不收缩。
6. 芯片使用公开处理后的 SST-RMA 表达值；不重复 log2，不送入计数模型。用 GPL23038 注释中的 ENSMUST，经 Ensembl GRCm38.87 映射到 ENSMUSG；多基因歧义探针簇剔除。多个唯一映射探针簇按基因取每样本中位数。映射丢失和歧义均计数。
7. 芯片使用 D7 与未损伤的双侧 Welch t 检验，在所有成功映射且具有有限检验值的基因上 BH 校正。它是小样本探索性复核，未采用 limma 的经验贝叶斯方差收缩，不将其等同于正式 limma 结果。
8. 跨队列支持门槛：通过发现门槛、芯片 `padj < 0.05`、芯片 `abs(mean_D7 - mean_D0) >= 0.5`，且两队列效应方向一致。芯片差值是其 log2 表达尺度上的效应。
9. 候选按两队列较大的校正 p 值由小到大展示；这个排序值不是联合 p 值或新的 FDR。报告上调与下调，允许没有通过复核的候选。

## 交付和解释边界

交付样本表、原始行示例、基因计数/TPM矩阵、排除与映射记录、PCA、主要差异结果、外部复核结果、筛选漏斗和候选表达曲线。每个候选保留完整证据列，而不是只给一个基因名字。

“跨队列支持”仅表示在两个相关实验中的同向表达变化；不能据此称为促进纤维化的因果基因、药物直接靶点或人类临床诊断标志物。两队列均用于本轮探索，未来若做预测模型，需要另外保留测试队列，特征筛选必须限定训练数据。

依据：[tximport 官方导入说明](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)、[长度缩放官方实现](https://github.com/thelovelab/tximport/blob/master/R/helper.R)、[PyDESeq2 官方示例](https://pydeseq2.readthedocs.io/en/stable/auto_examples/plot_minimal_pydeseq2_pipeline.html)、[发现集原论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC12089340/)、[复核集原论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC5957472/)。
