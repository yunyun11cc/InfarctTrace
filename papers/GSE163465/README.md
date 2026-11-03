# GSE163465 原论文与作者预印本

核对及获取日期：2026-11-03。

## 正式发表论文

Jin K, Gao S, Yang P, Guo R, Li D, Zhang Y, Lu X, Fan G, Fan X. **Single-Cell RNA Sequencing Reveals the Temporal Diversity and Dynamics of Cardiac Immunity after Myocardial Infarction.** Small Methods. 2022;6(3):e2100752. DOI：[10.1002/smtd.202100752](https://doi.org/10.1002/smtd.202100752)。PMID：[35023642](https://pubmed.ncbi.nlm.nih.gov/35023642/)。

GEO 的 `Series_pubmed_id` 为 35023642，与 Europe PMC 和 Crossref 元数据一致。`paper_metadata_europepmc.json`、`paper_metadata_crossref.json` 保存了核对来源。Europe PMC 标注正式版需订阅、未给出 PMCID；OpenAlex 给出的开放版本指向作者预印本。

**正式版 PDF 尚未下载成功。** Wiley 的 PDF、pdfdirect、epdf 请求均返回 HTTP 403；通过 Exa 获取正式全文页超时。没有把错误页面、摘要或生成摘要保存为正式原论文 PDF。

## 已获取的可读原文

[作者预印本 PDF](https://www.biorxiv.org/content/10.1101/2021.07.19.452902v1.full.pdf?download=1)：47 页，1,154,846 字节，真实 PDF 1.4；已核对标题、作者、首页 DOI 与 2021-07-19 版本日期，并渲染预览首页。这是作者预印本，不能标为正式出版 PDF。

[作者预印本官方全文](https://www.biorxiv.org/content/10.1101/2021.07.19.452902v1.full)可在线阅读；原文抽取版仅在本地保留，不随仓库上传。

- 预印本 DOI：[10.1101/2021.07.19.452902](https://doi.org/10.1101/2021.07.19.452902)，版本 v1，2021-07-19。
- `preprint_metadata_biorxiv.json` 中 `published` 字段明确对应正式 DOI `10.1002/smtd.202100752`，作者和题目一致。
- 这是使用 Exa `web_fetch_exa` 从官方页面提取的原文 Markdown，保留正文各节、图注和参考文献，**不是我们撰写的摘要，也不是正式排版 PDF**。
- 网页抽取未保证完整保留图片、数学公式和补充材料；阅读图像与公式应回到原文页面或正式版。
- bioRxiv 官方 PDF 初次请求返回 HTTP 403/429；稍后通过带 `?download=1` 的同一官方 PDF 地址成功获取，下载来源记录在 `download_records.json`。

下载处理后表达矩阵、GEO 样本说明及每组数据量见 [数据目录](../../datasets/GSE163465/README.md)。
