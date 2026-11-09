# GSE252027 预处理与差异分析方案

## 现有探索性基线

- 主要对比为 D7 对未损伤 D0，每组 4 个提交文库。统计单位为文库，每个文库通常汇总三只小鼠。组内重复编号用于区分文库。
- 24 个 FB 文库生成描述性表达矩阵、样本质控表和 PCA。D3 和 D28 的样本纳入在全时段分析前另行核查。
- 定量表保留 ENSMUSG 基因行，按去掉版本后缀的 Ensembl ID 对齐。GRCm38.87 提供基因符号和类型。
- D0/D7 的 TPM 乘以这些文库中该基因的平均有效长度，再按各文库估计片段总量缩放和取整。以至少 4 个文库计数大于或等于 10 为低表达过滤条件。
- 使用 PyDESeq2 的 `~condition` 模型和 Wald 检验，模型估计文库大小因子。候选门槛为校正后 P 值小于 0.05 且 log2 表达倍数变化绝对值大于或等于 1。

## 全时段主分析

1. 主分析纳入 D0、D1、D3、D7、D14、D28 分别 4、4、2、4、4、4 个文库。论文记录 FB D3 排除两个重复，公开代码保留 D3r1、D3r4，故排除 D3r2、D3r3。主分析保留 D28r4。由于公开脚本的 FB 片段只列 D28r1–r3，再单独用 21 个文库复算总体时间效应与 D28 对 D0，报告基因集重合程度。
2. 使用与 `lengthScaledTPM` 一致的长度缩放公式，按纳入文库重建整数尺度计数。保留至少一个时间组内有两个文库计数达到 10 的基因，PyDESeq2 估计文库大小因子和离散度。VST 用于 PCA 和热图。
3. 时间点作为分类因子。总体 LRT 固定 PyDESeq2 在全模型中估计的离散度和大小因子，比较六个组均值与一个共同均值的负二项最大似然。使用 5 个自由度，排除 Cook 异常后对有效 P 值进行 BH 校正。PyDESeq2 的 Wald 检验用于 D1、D3、D7、D14、D28 与 D0 的五个比较。每个比较保留 PyDESeq2 的独立过滤，并使用 `|log2FC| ≥ 1` 和 `padj < 0.05` 的门槛。

## 功能与通路富集

- 注释采用 Bioconductor 3.22 的 `org.Mm.eg.db` 和 `GO.db` 3.22.0，以及固定到本次快照的 KEGG 小鼠 `mmu` 基因集。对本次检验基因做稳定 Ensembl ID 与 Entrez ID 的一对一映射，排除两侧有歧义的编号。
- GO 对 BP、CC、MF 分开分析。GO 与 KEGG 的过度富集检验均用超几何分布，按各对比的上调、下调分别处理。背景是该对比有有限校正 P 值且能对应到该注释集的基因。背景内集合大小范围为 10–500，按对比、方向和注释集分别做 BH 校正。
- GSEA 用所有有有限 Wald 统计量和原始 P 值的映射基因，以带符号的 Wald 统计量排序。采用 GSEApy 的 multilevel 方法，KEGG 集合大小为 15–500，权重为 1，各对比使用固定种子 252028–252032。名义 P 值估计下限为 1e−10，每个对比内进行 BH 校正。并列分数按 Entrez ID 排序，运行时个别对比出现不超过 0.05% 的并列。

[主队列数据手册](../../datasets/GSE252027/README.md) · [DESeq2 的 LRT 方法说明](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) · [PyDESeq2 文档](https://pydeseq2.readthedocs.io/en/stable/api/docstrings/pydeseq2.dds.DeseqDataSet.html) · [GSEApy 文档](https://gseapy.readthedocs.io/en/latest/run.html)
