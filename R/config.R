default_config <- function() list(input = NULL, samples = NULL, output = NULL,
  intensities = 'intensities', regionprops = 'regionprops', images = 'img', masks = 'masks',
  transforms = 'arcsinh,arcsinh_zscore', cofactor = 1, min_area = 5,
  pixel_size = 1, seed = 220225, plot_cells = 5000, heatmap_cells = 500,
  embedding_cells = 0, snr_cells = 30000, neighbors = 15, perplexity = 30,
  rlog_max_cells = 1000, low_snr = 2, high_snr = 10, min_signal = 2,
  exclude = '', image_channels = '', skip_images = FALSE, skip_embeddings = FALSE)

help_text <- function() {
  cat('CyTOF QC\nUsage: Rscript bin/cytof-qc.R --input DIR --samples mapping.csv --output DIR [options]\n\n')
  cat('Mapping requires roi_id,sample_name; optional columns are retained.\n')
  for (n in names(default_config())) cat(sprintf('  --%-20s %s\n', gsub('_', '-', n), default_config()[[n]] %||% '(required)'))
  cat('\nBoolean switches take no value. Transform choices: arcsinh,arcsinh_zscore,rlog,rlog_zscore.\nSee README.md for input layout, rlog restrictions and output schemas.\n')
}
`%||%` <- function(x, y) if (is.null(x)) y else x
split_list <- function(x) { if (!nzchar(x)) character() else trimws(strsplit(x, ',', fixed=TRUE)[[1]]) }
parse_args <- function(args) {
  cfg <- default_config(); i <- 1L
  while (i <= length(args)) {
    key <- gsub('-', '_', sub('^--', '', args[i]))
    if (!startsWith(args[i], '--') || !key %in% names(cfg)) stop('Unknown option: ', args[i])
    if (is.logical(cfg[[key]])) cfg[[key]] <- TRUE else {
      i <- i + 1L
      if (i > length(args) || startsWith(args[i], '--')) stop('Missing value for ', key)
      cfg[[key]] <- if (is.numeric(cfg[[key]])) suppressWarnings(as.numeric(args[i])) else args[i]
    }
    i <- i + 1L
  }
  for (key in c('input','samples','output')) if (is.null(cfg[[key]])) stop('Required: --', key)
  nums <- names(cfg)[vapply(cfg, is.numeric, logical(1))]
  for (key in nums) if (!is.finite(cfg[[key]]) || cfg[[key]] < 0) stop('Invalid nonnegative number: ', key)
  for (key in c('cofactor','pixel_size','perplexity')) if (cfg[[key]] <= 0) stop(key, ' must be positive')
  for (key in c('seed','plot_cells','heatmap_cells','embedding_cells','snr_cells','neighbors','rlog_max_cells'))
    if (cfg[[key]] != floor(cfg[[key]]) || (!key %in% c('seed','embedding_cells') && cfg[[key]] < 2)) stop('Invalid integer: ', key)
  cfg$transforms <- unique(split_list(cfg$transforms))
  if (!length(cfg$transforms) || any(!cfg$transforms %in% c('arcsinh','arcsinh_zscore','rlog','rlog_zscore'))) stop('Invalid transforms')
  cfg
}
check_dependencies <- function(cfg) {
  packages <- c('mclust','dittoSeq','SingleCellExperiment','S4Vectors','SummarizedExperiment',
    'ggplot2','patchwork','ggrastr','viridis','RColorBrewer')
  if (!nzchar(Sys.which('pdfunite'))) stop('Missing pdfunite: install Poppler utilities (poppler-utils on Ubuntu).')
  if (!cfg$skip_embeddings) packages <- c(packages, 'uwot', 'Rtsne', 'scater')
  if (!cfg$skip_images) packages <- c(packages, 'EBImage', 'tiff', 'cytomapper')
  if (any(grepl('rlog', cfg$transforms))) packages <- c(packages, 'DESeq2', 'SummarizedExperiment')
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing)) stop('Missing packages: ', paste(missing, collapse=', '), '. Run Rscript scripts/install.R')
}
