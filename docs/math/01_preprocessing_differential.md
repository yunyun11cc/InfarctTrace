# 实验 01 数学原理与指标计算

这份说明对应 InfarctTrace 阶段 01 的 **GSE252027 成纤维细胞全时段 Python 分析**。用同一条数据流理解各个数字：

**定量结果 → 长度缩放计数 → 文库大小归一化 → 总体时间检验与逐时间点比较 → GO、KEGG、GSEA。**

主分析为 D0、D1、D3、D7、D14、D28 六组，共 22 个文库。D28r4 排除敏感性分析用 21 个文库。早期 D7 对 D0 教学基线单独使用 8 个文库。以下实际数值以全时段 Python 输出为准。 [实验方案](../../experiments/01_preprocessing_differential/protocol.md) · [阶段 01 结果解读](../../experiments/01_preprocessing_differential/results.md)

## 1. TPM、计数与归一化

用 $i$ 表示基因，用 $j$ 表示文库。$R_{ij}$ 是定量文件的 `NumReads`，$L_{ij}$ 是 `EffectiveLength`，$T_{ij}$ 是 `TPM`。一列对应一个测序文库。

### 1.1 Salmon 提供的 TPM

用简化的定量关系理解 TPM。先按有效长度校正估计片段数，再除以该样本中所有定量目标的丰度之和。

$$
T_{ij}=10^6\frac{R_{ij}/L_{ij}}{\sum_h R_{hj}/L_{hj}}
$$

较长的目标通常有更多机会被测到，所以先除以有效长度。再乘一百万后，结果表示在该文库的**相对丰度**。例如 10 TPM 表示每百万份长度校正丰度中占约 10 份。分母对应定量时的目标集合。项目读取提交的 TPM，筛出基因行以后不重新计算这一列。`NumReads` 则保留 Salmon 的片段分配估计，可能出现小数。 [Salmon 定量格式](https://salmon.readthedocs.io/en/latest/file_formats.html)

### 1.2 将 TPM 转回计数尺度

阶段 01 用与 `lengthScaledTPM` 相同的公式。设 $S$ 为**本次纳入**的文库集合，$G$ 为保留的基因行集合。

$$
\overline L_i=\frac{1}{|S|}\sum_{j\in S}L_{ij},\qquad
U_{ij}=T_{ij}\overline L_i
$$

$$
C_{ij}=\operatorname{round}\left(
U_{ij}\cdot\frac{\sum_{h\in G}R_{hj}}{\sum_{h\in G}U_{hj}}
\right).
$$

第一步给同一个基因使用跨文库的平均有效长度。第二步用每列原始估计片段总量缩放。取整前该列总量与保留基因行的估计片段总量一致。这里的 $C_{ij}$ 是**长度缩放后的整数计数**，供负二项模型使用。22 文库主分析和 21 文库敏感性分析分别计算 $\overline L_i$。

以 Cthrc1（`ENSMUSG00000054196`）在 `GSM7992326` 为例。该文库的 TPM 为 $0.651286$，原估计片段数为 $28$。22 文库中这个基因的平均有效长度约为 $472.3783$，保留基因行的本列估计片段总量约为 $25{,}863{,}604.909$，$\sum_h U_{hj}$ 约为 $312{,}645{,}303.596$。

$$
0.651286\times472.3783\times
\frac{25{,}863{,}604.909}{312{,}645{,}303.596}
\approx25.4506\ \longrightarrow\ C_{ij}=25.
$$

[项目转换代码](../../experiments/01_preprocessing_differential/prepare.py#L95) · [tximport 方法](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)

### 1.3 低表达过滤和文库大小因子

全时段分析的过滤规则是：存在至少一个时间组，该基因在组内至少两个文库中 $C_{ij}\ge 10$。例如只在 D3 有足够表达的基因，也有机会通过。这个过滤把 51,029 个基因缩小到 19,312 个模型输入基因。

PyDESeq2 用中位数比值法估计文库大小因子 $s_j$。在用于估计的各文库均为正计数的基因中，先对每个基因求参考几何均值 $g_i$，再为每个文库求对数比值的中位数。

$$
g_i=\exp\!\left(\frac{1}{|S|}\sum_{j\in S}\ln C_{ij}\right),\quad
s_j=\exp\!\left[\operatorname{median}_i(\ln C_{ij}-\ln g_i)\right],\quad
C_{ij}^{\mathrm{norm}}=C_{ij}/s_j.
$$

实际拟合时，$s_j$ 作为模型里的已知偏移量使用。归一化矩阵用于查看表达大小等。结果表中的 `平均标准化计数`（`baseMean`）是这个基因在纳入文库中的 $C_{ij}^{\mathrm{norm}}$ 均值，和某一组的均值是两个不同的量。

示例文库 `GSM7992326` 的大小因子约为 $0.99936774$，Cthrc1 的标准化计数约为 $25/0.99936774=25.01582$。Cthrc1 在 22 个文库的 `baseMean` 约为 $2386.59$。 [项目模型入口](../../experiments/01_preprocessing_differential/differential.py#L22) · [PyDESeq2 源码](https://github.com/scverse/PyDESeq2/blob/v0.5.4/pydeseq2/preprocessing.py)

## 2. 负二项模型和效应量

对每个基因，长度缩放后的计数 $C_{ij}$ 服从负二项模型。$q_{i,d(j)}$ 是文库 $j$ 对应时间组的表达参数，$\alpha_i$ 是基因的离散度。

$$
C_{ij}\sim\mathrm{NB}(\mu_{ij},\alpha_i),\quad
\mu_{ij}=s_jq_{i,d(j)},\quad
\operatorname{Var}(C_{ij})=\mu_{ij}+\alpha_i\mu_{ij}^{2}.
$$

离散度描述同一时间组内的额外波动。PyDESeq2 先估计各基因的离散度，再拟合离散度与平均表达的趋势，并借助该趋势稳定逐基因估计。模型在对数链接上写成 $\ln\mu_{ij}=\ln s_j+x_j^\mathsf{T}\beta_i$。$x_j$ 记录该文库所属时间组，D0 为对比基准。 [DESeq2 负二项模型](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

### 2.1 log2FoldChange

以 D7 对 D0 为例，$c^\mathsf{T}\hat\beta_i$ 取出两组在自然对数尺度的差，除以 $\ln 2$ 就得到结果表里的效应量。

$$
\widehat{\log_2FC}_i
=\frac{c^\mathsf{T}\hat\beta_i}{\ln 2}
=\log_2(\hat q_{i,\mathrm{D7}}/\hat q_{i,\mathrm{D0}}).
$$

$+1$ 对应模型估计为 2 倍，$-1$ 对应 0.5 倍。Cthrc1 在 D7 对 D0 的 `log2表达倍数变化` 为 $7.8289$，对应模型估计约 $2^{7.8289}=227.37$ 倍。这是模型效应的点估计。若直接用两组平均 TPM 相除，得到的是另一种描述量。

### 2.2 标准误、Wald 统计量和原始 P 值

`效应标准误` 衡量 $\widehat{\log_2FC}$ 的估计不确定性。PyDESeq2 使用模型的加权信息矩阵。$X$ 是设计矩阵，$W$ 的对角元素是 $\mu_{ij}/(1+\alpha_i\mu_{ij})$，$M=X^\mathsf{T}WX$，$H=(M+\Lambda)^{-1}$，$\Lambda$ 是实现中的微小岭稳定项。

$$
SE_{\log_2FC}
=\frac{\sqrt{(Hc)^\mathsf{T}M(Hc)}}{\ln 2},\qquad
Z=\frac{\widehat{\log_2FC}}{SE_{\log_2FC}}.
$$

逐时点比较的零假设为**该时间点和 D0 的模型表达相同**，即 $\log_2FC=0$。在零假设和模型近似下，$Z$ 用标准正态分布计算双侧尾概率。

$$
p_{\mathrm{Wald}}=\Pr_{H_0}(|Z_{\mathrm{null}}|\ge |Z_{\mathrm{obs}}|)
=2[1-\Phi(|Z_{\mathrm{obs}}|)].
$$

因此 P 值的概率对象是：**假设零效应成立，得到当前这么极端或更极端的统计量的概率**。同样的变化幅度，如果文库间波动大、标准误大，$|Z|$ 就会变小。Cthrc1 的 D7 结果为 $7.8289/0.304924\approx25.67494$，原始 P 值约为 $2.23\times10^{-145}$。实际计算直接用正态分布的尾概率函数，避免极小值在 `1-CDF` 运算中因精度丢失变成零。 [PyDESeq2 Wald 实现](https://github.com/scverse/PyDESeq2/blob/v0.5.4/pydeseq2/utils.py)

## 3. 总体时间检验 LRT

对同一基因，比较两个嵌套模型：

- **简化模型**使用一个共同表达均值，六组共用。
- **完整模型**允许六个采样组各有自己的表达均值。

两者沿用完整模型估计的离散度与文库大小因子。对一组文库设 $\mu_j=s_j q$，固定 $\alpha$ 时，对数似然中随 $q$ 变化的部分为：

$$
\ell(q)=\sum_j\left[
C_j\ln\mu_j-(C_j+1/\alpha)\ln(1+\alpha\mu_j)
\right].
$$

与 $q$ 无关的常数在两模型相减时抵消。代码用二分法求一组文库的最大似然，找使 $\sum_j (C_j-\mu_j)/(1+\alpha\mu_j)=0$ 的组均值。全零组按零均值的极限处理。完整模型的最大对数似然是六组最大值之和。

$$
D=2(\ell_{\mathrm{full}}-\ell_{\mathrm{reduced}}),\qquad
p_{\mathrm{LRT}}=\Pr(\chi^2_{5}\ge D).
$$

自由度为 $6-1=5$。LRT 的零假设是**六组表达均值相等**。结果表的 `LRT统计量` 表示增加时间组参数带来的拟合改善。若总体检验显著，再读 Wald 表定位哪一天变化及其方向。

Cthrc1 的总体 LRT 统计量为 $657.9594$，原始 P 值约为 $6.03\times10^{-140}$，校正后约为 $5.77\times10^{-136}$。D7 的 Wald P 值和总体 P 值回答不同问题。 [项目 LRT 实现](../../experiments/01_preprocessing_differential/differential.py#L37) · [时间因子与 LRT](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#likelihood-ratio-test)

## 4. Cook 距离、校正后 P 值与筛选门槛

### 4.1 Cook 距离

Cook 距离衡量某个文库对某个基因的模型拟合影响。PyDESeq2 的计算可概括为：

$$
D_{ij}^{\mathrm{Cook}}=
\frac{(C_{ij}-\hat\mu_{ij})^2}{p(\hat\mu_{ij}+\tilde\alpha_i\hat\mu_{ij}^{2})}
\cdot\frac{h_{ij}}{(1-h_{ij})^{2}}.
$$

这里 $p$ 是模型参数数，$h_{ij}$ 是拟合杠杆值，$\tilde\alpha_i$ 是专供 Cook 计算的稳健离散度。大的残差叠加大的杠杆会推高距离。包用 $F_{p,n-p}$ 的 99% 分位作门槛，寻找最大 Cook 时考虑至少有三个重复的组，并进一步检查计数。被标记的基因 P 值在对应检验中留空。本项目关闭 Cook 异常计数的自动替换。文库纳入表中的 D3 取舍属于样本级决定。 [PyDESeq2 Cook 实现](https://github.com/scverse/PyDESeq2/blob/v0.5.4/pydeseq2/dds.py)

### 4.2 BH 校正

一次检验成千上万个基因时，要在一个检验家族内处理多重比较。把可校正的 $m$ 个 P 值从小到大排序为 $p_{(1)},\dots,p_{(m)}$，第 $i$ 个的 BH 校正结果为：

$$
q_{(i)}=\min\left(1,\min_{k\ge i}\frac{m}{k}p_{(k)}\right).
$$

例如 3 个 P 值 $0.010,0.040,0.045$，先乘 $m/k$ 得 $0.030,0.060,0.045$，再从右向左取最小值，得到 $0.030,0.045,0.045$。最后映射回原基因。FDR 指重复执行筛选时，发现名单中假发现比例的期望值。`padj` 是针对多重检验计算的校正值。

阶段 01 使用以下校正口径：

| 输出 | 哪些 P 值一起校正 | 如何判定 |
|---|---|---|
| 总体 LRT | 全部通过低表达过滤、且有有限 P 值的基因 | `校正后P值 < 0.05` |
| 各时间点 Wald | 每个比较单独校正，保留 PyDESeq2 的独立过滤 | `校正后P值 < 0.05` 且 `\|log2FC\| ≥ 1` |
| GO、KEGG 过度富集 | 每个比较、方向、注释集各自的可检验条目 | `校正后P值 < 0.05` |
| KEGG GSEA | 每个比较内有有效名义 P 值的通路 | `校正后P值 < 0.05` |

Wald 的独立过滤使用 `baseMean`。PyDESeq2 在 50 个候选分位阈值上分别计算 BH 结果，观察在给定 $\alpha=0.05$ 下的拒绝数，再按平滑趋势选择阈值。落在阈值以下的基因可以有原始 P 值，而校正列留空。Cthrc1 在 D7 对 D0 的 `padj` 约为 $4.27\times10^{-141}$，和 $7.8289$ 的效应量一起满足两个预设门槛。这个门槛组合是零效应检验加效应量点估计过滤。 [SciPy 的 BH 接口](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html) · [PyDESeq2 独立过滤](https://github.com/scverse/PyDESeq2/blob/v0.5.4/pydeseq2/ds.py)

### 4.3 D28 敏感性与 Jaccard 指数

设 $A$ 是 22 文库主分析的入选基因，$B$ 是排除 D28r4 后的入选基因。

$$
J(A,B)=\frac{|A\cap B|}{|A\cup B|}.
$$

D28 对 D0 的两份名单分别有 866 和 609 个基因，共同入选 524 个，所以 $J=524/(866+609-524)=0.551$。总体时间效应的 $J$ 约为 0.869。这个数衡量名单重合程度。要看某个候选是否受 D28r4 影响，还需同时查看两次拟合的效应量与校正后 P 值。 [敏感性结果说明](../../experiments/01_preprocessing_differential/results.md)

## 5. 标准化图、PCA、热图和火山图

- **图 1 提交文库 PCA**先算 $\log_2(TPM+1)$，按基因跨文库方差取前 1,000 个，使用全部 24 个提交文库。**图 3 主分析 PCA**使用 22 个文库的 VST，按方差取前 500 个基因。
- **图 2 箱线图**把同一组基因的长度缩放计数、归一化计数各变为 $\log_2(x+1)$，每箱代表一份文库。箱体覆盖第 25 至 75 百分位，横线是中位数。两图比较的是同一批文库。
- **图 4 热图**按总体 LRT 校正 P 值选前 25 个显著基因，对每个基因在 22 个文库中的 VST 取行内标准分数 $z_{ij}=(v_{ij}-\bar v_i)/s_i$，其中 $s_i$ 为样本标准差。颜色是同一行内的相对高低，色标在 −2.5 和 2.5 处截断显示。
- **图 5 火山图**横坐标是 $\log_2FC$，纵坐标是 $-\log_{10}(padj)$。横虚线为 $-\log_{10}(0.05)\approx1.301$，竖虚线为 $-1$ 和 $+1$。例如 $padj=0.001$ 对应纵坐标 3。

PCA 把每个文库当成一个多维点，所选基因是坐标。对各基因在文库间中心化得到 $X_c$，对 $X_c^\mathsf{T}X_c/(n-1)$ 做特征分解。文库沿特征向量的投影 $X_cv_k$ 是 PC 得分，轴括号里的百分比是 $\lambda_k/\sum_j\lambda_j$。样本点距离用于描述所选表达空间中的接近程度。

VST 是方差稳定变换，用来减弱表达越高、波动越大的现象。本轮采用 PyDESeq2 参数化趋势 $\alpha(x)=a_0+a_1/x$，对归一化计数 $x$ 使用：

$$
v(x)=\log_2\!\left(
\frac{1+a_1+2a_0x+2\sqrt{a_0x(1+a_1+a_0x)}}{4a_0}
\right).
$$

例如 Cthrc1 在 `GSM7992326` 的归一化计数为 25.01582，当前 VST 值约为 6.53641。VST 用于图 3、图 4。 [项目绘图实现](../../experiments/01_preprocessing_differential/visualize.py) · [PyDESeq2 VST 实现](https://github.com/scverse/PyDESeq2/blob/v0.5.4/pydeseq2/dds.py)

## 6. GO 与 KEGG 过度富集

这两个表使用同一个**超几何检验**。不同的是基因集来源。GO 以生物过程（BP）、细胞组分（CC）、分子功能（MF）组织基因，KEGG 以通路组织基因。先做 Ensembl 到 Entrez 的双向一对一映射。

在某个比较、方向、注释集合下定义：

| 记号 | 对应结果列 | 含义 |
|---|---|---|
| $N$ | `背景基因数` | 该比较有有效校正 P 值并能在对应注释库找到的基因数 |
| $K$ | `背景条目大小` | 背景中属于当前条目的基因数 |
| $n$ | `候选数` | 背景中通过差异门槛且属于当前方向的基因数 |
| $k$ | `命中数` | 候选中也属于当前条目的基因数 |

零假设下，从 $N$ 个背景基因中等概率、不放回地选 $n$ 个。条目命中数 $X$ 服从超几何分布。P 值取观察到 $k$ 个或更多命中的右尾。

$$
E[X]=n\frac{K}{N},\qquad
\text{富集倍数}=\frac{k/n}{K/N}=\frac{k}{E[X]},\qquad
p=\Pr(X\ge k)=\sum_{x=k}^{\min(K,n)}
\frac{\binom Kx\binom{N-K}{n-x}}{\binom Nn}.
$$

代码用 `hypergeom.sf(k - 1, N, K, n)`。仅检验背景内大小在 10–500 的条目，零命中条目仍计入相应的 BH 家族。富集倍数衡量命中比例相对背景的提升，统计证据要再看校正 P 值。

| D7 上调例子 | $N$ | $K$ | $n$ | $k$ | 期望命中 | 富集倍数 | 原始 P 值 | 校正后 P 值 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| GO 细胞外基质组织 `GO:0030198` | 14,884 | 284 | 1,169 | 54 | 22.31 | 2.42 | $8.38\times10^{-10}$ | $1.34\times10^{-7}$ |
| KEGG 细胞外基质与受体相互作用 `mmu04512` | 6,881 | 82 | 600 | 26 | 7.15 | 3.64 | $2.57\times10^{-9}$ | $8.96\times10^{-7}$ |

同一个时间比较在 GO 和 KEGG 中有各自的背景与候选数，因为两个数据库的覆盖范围不同。 [项目 ORA 实现](../../experiments/01_preprocessing_differential/enrichment.py#L55) · [SciPy 超几何分布](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.hypergeom.html)

## 7. GSEA 的 ES、NES、P 值与核心基因

GSEA 使用**该对比所有具备有效 Wald 统计量和原始 P 值的一对一映射基因**。先按带符号的 Wald 统计量从大到小排序。较大的正值在上调端，较大的负值在下调端。这个列表比已经通过差异门槛的名单更广。

设排序中有 $N$ 个基因，通路集合为 $S$，与排序重合的基因数为 $m$，第 $j$ 名的 Wald 分数为 $r_j$。本次权重指数为 1。沿列表从头到尾走，命中通路基因时按 $|r_j|$ 增加，其他位置等量扣减。走到第 $i$ 个位置的累计分数为：

$$
R(i)=
\sum_{\substack{j\le i\\g_j\in S}}\frac{|r_j|}{\sum_{g_h\in S}|r_h|}
-\sum_{\substack{j\le i\\g_j\notin S}}\frac{1}{N-m}.
$$

取累计曲线离 0 最远的**带符号**峰值为 ES。正 ES 指通路成员倾向排序前端，负 ES 指倾向后端。正向的核心基因是通路成员中位于正峰之前的那部分，负向对应负峰之后的成员。

零分布来自随机基因集的富集分数。NES 将实际 ES 除以对应方向零分布 ES 的平均绝对大小，使不同集合的分数更可比。名义 P 值衡量同方向下随机基因集得到至少同样极端 ES 的尾概率。本轮使用 GSEApy multilevel 估计该概率，集合重合大小限制为 15–500，估计下限为 $10^{-10}$。代码再按时间比较对名义 P 值作 BH 校正。

例如 D7 对 D0 的 `mmu04512` 共有 **16,368 个基因参与排序**，其中 **82 个属于该通路**。`ES` 为 $0.58339$，`NES` 为 $1.97520$，说明它偏向上调端。该例名义 P 值约为 $9.01\times10^{-7}$，校正后约为 $1.68\times10^{-5}$。NES 衡量整条基因集在排序中的偏向程度。 [项目 GSEA 实现](../../experiments/01_preprocessing_differential/enrichment.py#L80) · [GSEApy prerank 文档](https://gseapy.readthedocs.io/en/latest/run.html#gseapy.prerank)

## 8. 指标与结果文件索引

以下文件位于 `outputs/01_preprocessing_differential/`。

| 文件 | 主要用来回答什么 |
|---|---|
| `tables/sample_inclusion.csv` | 哪些文库进主分析和敏感性分析 |
| `tables/discovery_qc.csv` | 每个文库有多少特征、估计片段总量、TPM 总和和 TPM 达到 1 的基因数 |
| `matrices/FB_TPM.tsv.gz` | 某基因在某文库的相对丰度 |
| `matrices/timecourse_length_scaled_counts.tsv.gz` | 负二项模型使用的整数尺度计数 |
| `tables/timecourse_size_factors.tsv`、`matrices/timecourse_normalized_counts.tsv.gz` | 文库大小因子和归一化计数 |
| `matrices/timecourse_vst.tsv.gz` | 主分析 PCA 和热图使用的变换值 |
| `tables/timecourse_lrt.tsv.gz` | 全时间窗口内是否存在表达差异 |
| `tables/timecourse_wald.tsv.gz` | 哪一天变化、方向如何、效应和统计证据有多大 |
| `tables/timecourse_sensitivity_summary.tsv` | 更换 D28r4 纳入规则后，名单有多少重合 |
| `tables/enrichment_go_ora.tsv.gz`、`tables/enrichment_kegg_ora.tsv.gz` | 已筛出的差异基因中，哪些功能和通路命中比例偏高 |
| `tables/enrichment_kegg_gsea.tsv.gz` | 完整排序中，通路成员偏向哪一端，哪些成员推动富集 |

[输出图表索引](../../outputs/README.md)
