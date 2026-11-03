# GSE236374：全心组织补充验证队列

状态：已下载官方计数矩阵、9 个样本的元数据、正式原论文和补充 PDF；未作标准化、差异分析或模型训练。

## 下载内容与检查

- [FeatureCounts 计数矩阵](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE236nnn/GSE236374/suppl/GSE236374_MI_FeatureCounts_matrix.txt.gz)：379,384 字节，24,964 个唯一基因符号 × 9 个样本。
- [样本表](metadata/samples.csv)：sham、MI 后 7 天、28 天各 3 个样本，记录 GSM 与矩阵列名对应关系。
- [GEO family SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE236nnn/GSE236374/soft/GSE236374_family.soft.gz)、[系列记录](metadata/GSE236374_series.soft.txt)。
- [矩阵检查结果](metadata/data_inspection.json)：全部数值有限、非负、为整数，9 列与元数据匹配，无重复基因符号。

来源为小鼠整块心脏组织 bulk RNA-seq，永久 LAD 结扎模型；不能当成纯化成纤维细胞队列。表达差异可能同时受细胞比例改变影响。

## 对应原论文

Yu P、Song S 等，2023，PNAS。*Downregulation of apoptotic repressor AVEN exacerbates cardiac injury after myocardial infarction.* DOI：[10.1073/pnas.2302482120](https://doi.org/10.1073/pnas.2302482120)。

- [正式 PDF](https://pmc-oa-opendata.s3.amazonaws.com/PMC10589712.1/PMC10589712.1.pdf)，11 页。
- [全文 XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10589712/fullTextXML)。
- [原论文补充材料 PDF](https://pmc-oa-opendata.s3.amazonaws.com/PMC10589712.1/pnas.2302482120.sapp.pdf)。

来源：[GEO GSE236374](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236374)。每个源文件的来源和哈希见 [download_records.json](download_records.json)。未下载原始测序 reads。
