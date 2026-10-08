zscore <- function(x) {
  s <- apply(x, 1, stats::sd); s[!is.finite(s) | s == 0] <- 1
  sweep(sweep(x, 1, rowMeans(x), '-'), 1, s, '/')
}
validate_rlog <- function(x, max_cells) {
  if (ncol(x) > max_cells) stop('rlog exceeds --rlog-max-cells; DESeq2 rlog is expensive. Use arcsinh or explicitly raise the limit.')
  if (any(x != floor(x)) || any(x > .Machine$integer.max)) stop('rlog requires integer counts; mean intensities cannot be silently rounded. Use integer summed counts or arcsinh.')
  if (ncol(x) < 3 || any(colSums(x) == 0)) stop('rlog requires at least 3 cells and no all-zero cells')
}
transform_counts <- function(x, method, cfg) {
  base <- sub('_zscore$', '', method)
  if (base == 'arcsinh') y <- asinh(x/cfg$cofactor) else {
    validate_rlog(x, cfg$rlog_max_cells)
    dds <- DESeq2::DESeqDataSetFromMatrix(round(x), data.frame(row.names=colnames(x)), ~1)
    dds <- DESeq2::estimateSizeFactors(dds, type='poscounts')
    y <- SummarizedExperiment::assay(DESeq2::rlog(dds, blind=TRUE, fitType='mean'))
  }
  if (grepl('_zscore$', method)) y <- zscore(y)
  if (any(!is.finite(y))) stop('Transformation produced non-finite values')
  y
}
