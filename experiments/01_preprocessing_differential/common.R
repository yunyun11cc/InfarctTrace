# 阶段 01 的 R 模块从项目根目录运行，共用输出位置与 UTF-8 表格格式。
stage_out <- file.path('outputs', '01_preprocessing_differential')
out_dir <- function(kind, name) file.path(stage_out, kind, name)
for (kind in c('matrices', 'tables', 'reports')) dir.create(file.path(stage_out, kind), recursive = TRUE, showWarnings = FALSE)
local_lib <- file.path('.tools', 'R-library')
if (dir.exists(local_lib)) .libPaths(c(normalizePath(local_lib), .libPaths()))

read_tsv <- function(path, ...) {
  read.delim(path, check.names = FALSE, fileEncoding = 'UTF-8-BOM', ...)
}

write_tsv <- function(x, path) {
  con <- if (endsWith(path, '.gz')) gzfile(path, 'wt', encoding = 'UTF-8') else file(path, 'wt', encoding = 'UTF-8')
  on.exit(close(con))
  write.table(x, con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
}
