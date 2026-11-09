# 历史 R 对照。当前 GSE252027 流程使用 differential.py。
# 此脚本会写入同名正式输出，重现旧版时应在独立工作副本执行。
source('experiments/01_preprocessing_differential/common.R', encoding = 'UTF-8')
suppressPackageStartupMessages(library(DESeq2))

admitted <- read.csv(out_dir('tables', 'sample_inclusion.csv'), check.names = FALSE, fileEncoding = 'UTF-8-BOM')
quant <- lapply(c('TPM', 'EffectiveLength', 'NumReads'), function(field) {
  as.matrix(read_tsv(out_dir('matrices', paste0('FB_', field, '.tsv.gz')), row.names = 1))
})
stopifnot(all(vapply(quant, function(mat) identical(dimnames(mat), dimnames(quant[[1]])), logical(1))),
          identical(colnames(quant[[1]]), admitted[['GSM编号']]))
annotation <- read_tsv(out_dir('matrices', 'gene_annotation.tsv.gz'), row.names = 1)
times <- paste0('D', c(0, 1, 3, 7, 14, 28))

fit_time_model <- function(keep_sample) {
  subset <- admitted[keep_sample, ]
  ids <- subset[['GSM编号']]
  coldata <- data.frame(day = factor(paste0('D', subset[['时间点（天）']]), levels = times), row.names = ids)
  stopifnot(all(table(coldata$day) >= 2))
  # 按纳入文库计算平均有效长度。这里导入的是长度缩放后、文库总量尺度的计数。
  counts <- round(tximport::makeCountsFromAbundance(
    quant[[3]][, ids], quant[[1]][, ids], quant[[2]][, ids], 'lengthScaledTPM'))
  # 任一时间组内至少两个文库达到 10，保留只在一个时间段表达的基因。
  expressed <- Reduce(`|`, lapply(times, function(day) {
    rowSums(counts[, coldata$day == day, drop = FALSE] >= 10) >= 2
  }))
  dds <- DESeqDataSetFromMatrix(countData = counts[expressed, , drop = FALSE], colData = coldata, design = ~ day)
  DESeq(dds, test = 'LRT', reduced = ~ 1, minReplicatesForReplace = Inf, quiet = TRUE)
}

attach_gene <- function(r) {
  ids <- rownames(r)
  data.frame('基因ID' = ids, '基因符号' = annotation[ids, '基因符号'], check.names = FALSE)
}

hits <- function(r) is.finite(r$padj) & r$padj < 0.05 & is.finite(r$log2FoldChange) & abs(r$log2FoldChange) >= 1
matrix_table <- function(x) data.frame('基因ID' = rownames(x), x, check.names = FALSE)

cat('拟合 22 个文库的全时段主模型。\n')
main <- fit_time_model(admitted[['主分析纳入']] == '是')
lrt <- results(main, alpha = 0.05)
lrt_pass <- is.finite(lrt$padj) & lrt$padj < 0.05
# LRT 同时检验五个时间系数，因此总体结果只报统计量与 P 值。
overall <- cbind(attach_gene(lrt), data.frame('平均标准化计数' = lrt$baseMean, 'LRT统计量' = lrt$stat,
                 'P值' = lrt$pvalue, '校正后P值' = lrt$padj, '总体时间效应显著' = ifelse(lrt_pass, '是', '否'), check.names = FALSE))
write_tsv(overall[order(overall[['校正后P值']]), ], out_dir('tables', 'timecourse_lrt.tsv.gz'))

# 在同一批计数、大小因子与离散度上拟合 Wald 系数，得到有方向的对比。
wald <- nbinomWaldTest(main, betaPrior = FALSE, quiet = TRUE)
comparisons <- lapply(times[-1], function(day) results(wald, contrast = c('day', day, 'D0'), alpha = 0.05))
names(comparisons) <- times[-1]
contrast_rows <- lapply(times[-1], function(day) {
  r <- comparisons[[day]]
  cbind(attach_gene(r), data.frame('比较' = paste0(day, '_vs_D0'), '平均标准化计数' = r$baseMean,
        'log2表达倍数变化' = r$log2FoldChange, '效应标准误' = r$lfcSE, 'Wald统计量' = r$stat,
        'P值' = r$pvalue, '校正后P值' = r$padj, '通过差异筛选' = ifelse(hits(r), '是', '否'), check.names = FALSE))
})
write_tsv(do.call(rbind, contrast_rows), out_dir('tables', 'timecourse_wald.tsv.gz'))
summary <- do.call(rbind, lapply(times[-1], function(day) {
  r <- comparisons[[day]]
  data.frame('比较' = paste0(day, '_vs_D0'), '文库数' = sum(main$day == day),
             '有限校正P值基因数' = sum(is.finite(r$padj)),
             '上调基因数' = sum(hits(r) & r$log2FoldChange > 0),
             '下调基因数' = sum(hits(r) & r$log2FoldChange < 0), check.names = FALSE)
}))
write_tsv(summary, out_dir('tables', 'timecourse_contrast_summary.tsv'))
write_tsv(matrix_table(counts(main)), out_dir('matrices', 'timecourse_length_scaled_counts.tsv.gz'))
write_tsv(matrix_table(counts(main, normalized = TRUE)), out_dir('matrices', 'timecourse_normalized_counts.tsv.gz'))
vst <- varianceStabilizingTransformation(main, blind = FALSE)
write_tsv(matrix_table(assay(vst)), out_dir('matrices', 'timecourse_vst.tsv.gz'))
write_tsv(data.frame('GSM编号' = colnames(main), '时间组' = main$day, '文库大小因子' = sizeFactors(main),
                     check.names = FALSE), out_dir('tables', 'timecourse_size_factors.tsv'))

# D28r4 在作者脚本与正文之间有取舍差异，用排除后的结果量化敏感性。
cat('拟合排除 D28r4 的 21 文库敏感性模型。\n')
sensitivity <- fit_time_model(admitted[['敏感性分析纳入']] == '是')
sens_lrt <- results(sensitivity, alpha = 0.05)
sens_wald <- nbinomWaldTest(sensitivity, betaPrior = FALSE, quiet = TRUE)
sens_d28 <- results(sens_wald, contrast = c('day', 'D28', 'D0'), alpha = 0.05)
primary_d28 <- comparisons[['D28']]
shared <- intersect(rownames(primary_d28), rownames(sens_d28))
names(lrt_pass) <- rownames(lrt)
sens_lrt_pass <- is.finite(sens_lrt$padj) & sens_lrt$padj < 0.05
names(sens_lrt_pass) <- rownames(sens_lrt)
write_tsv(data.frame('基因ID' = shared,
                     '主分析总体校正后P值' = lrt[shared, 'padj'],
                     '敏感性总体校正后P值' = sens_lrt[shared, 'padj'],
                     '主分析D28效应' = primary_d28[shared, 'log2FoldChange'],
                     '敏感性D28效应' = sens_d28[shared, 'log2FoldChange'],
                     '主分析D28校正后P值' = primary_d28[shared, 'padj'],
                     '敏感性D28校正后P值' = sens_d28[shared, 'padj'], check.names = FALSE),
          out_dir('tables', 'timecourse_d28_sensitivity.tsv.gz'))
jaccard <- function(a, b) {
  combined <- union(a, b)
  if (length(combined)) length(intersect(a, b)) / length(combined) else NA_real_
}
main_lrt_ids <- names(lrt_pass)[lrt_pass]
sens_lrt_ids <- names(sens_lrt_pass)[sens_lrt_pass]
main_d28_ids <- rownames(primary_d28)[hits(primary_d28)]
sens_d28_ids <- rownames(sens_d28)[hits(sens_d28)]
write_tsv(data.frame('检验' = c('总体时间效应', 'D28_vs_D0'),
                     '主分析入选基因数' = c(length(main_lrt_ids), length(main_d28_ids)),
                     '敏感性分析入选基因数' = c(length(sens_lrt_ids), length(sens_d28_ids)),
                     '共同入选基因数' = c(length(intersect(main_lrt_ids, sens_lrt_ids)), length(intersect(main_d28_ids, sens_d28_ids))),
                     'Jaccard指数' = c(jaccard(main_lrt_ids, sens_lrt_ids), jaccard(main_d28_ids, sens_d28_ids)),
                     check.names = FALSE), out_dir('tables', 'timecourse_sensitivity_summary.tsv'))
writeLines(capture.output(sessionInfo()), out_dir('reports', 'differential_r_session.txt'))
cat('全时段进入模型的基因数', nrow(main), '总体时间效应显著数', sum(lrt_pass), '\n')
print(summary, row.names = FALSE)
