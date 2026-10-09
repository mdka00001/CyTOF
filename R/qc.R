snr_partition <- function(x, positive) {
  if (!any(positive) || all(positive)) return(c(signal=NA, noise=NA, snr=NA))
  signal <- mean(x[positive]); noise <- mean(x[!positive])
  c(signal=signal, noise=noise, snr=if (noise == 0) Inf else signal/noise)
}
cell_snr <- function(x, cfg) {
  idx <- sort(sample.int(ncol(x), min(ncol(x), cfg$snr_cells)))
  do.call(rbind, lapply(seq_len(nrow(x)), function(i) {
    v <- x[i,idx]; status <- 'ok'; result <- c(signal=NA, noise=NA, snr=NA)
    if (length(unique(v)) < 3) status <- 'insufficient_variation' else {
      fit <- tryCatch(mclust::Mclust(asinh(v/cfg$cofactor), G=2, verbose=FALSE), error=function(e) NULL)
      if (is.null(fit) || length(unique(fit$classification)) != 2) status <- 'mixture_failed' else {
        means <- tapply(v, fit$classification, mean)
        result <- snr_partition(v, fit$classification == as.integer(names(which.max(means))))
      }
    }
    data.frame(marker=rownames(x)[i], t(result), status=status, cells_used=length(idx))
  }))
}
marker_review <- function(snr, cfg) {
  snr$review <- vapply(seq_len(nrow(snr)), function(i) {
    r <- snr[i,]
    if (is.na(r$snr)) return('Review: SNR unavailable or constant marker')
    if (!is.finite(r$snr)) return('Review: zero background; infinite SNR is not proof of specificity')
    if (r$signal <= cfg$min_signal) return('Potential exclusion: low positive signal; inspect images')
    if (r$snr < cfg$low_snr) return('Potential exclusion: low signal/background separation')
    if (r$snr >= cfg$high_snr) return('High SNR: potentially informative; confirm localization')
    'Intermediate SNR: assess biology and staining'
  }, character(1))
  snr
}
roi_metrics <- function(data, cfg) {
  do.call(rbind, lapply(seq_len(nrow(data$images)), function(i) {
    roi <- data$images$roi_id[i]; d <- data$cells[data$cells$roi_id == roi,,drop=FALSE]
    pixels <- data$images$width_px[i]*data$images$height_px[i]
    kept <- d$area >= cfg$min_area
    data.frame(roi_id=roi, cells_before=nrow(d), cells_after=sum(kept),
      coverage=sum(d$area)/pixels, cells_per_mm2=sum(kept)/(pixels*cfg$pixel_size^2/1e6))
  }))
}
compute_embeddings <- function(x, use, cfg, idx) {
  selected <- use & apply(x, 1, stats::sd) > 0
  if (sum(selected) < 2 || length(idx) < 5) stop('Embeddings require at least 2 variable selected markers and 5 cells; use --skip-embeddings')
  object <- SingleCellExperiment::SingleCellExperiment(assays=list(exprs=x[,idx,drop=FALSE]))
  n <- length(idx)
  perplexity <- min(cfg$perplexity, (n-2)/3)
  set.seed(cfg$seed)
  object <- scater::runUMAP(object, subset_row=selected, exprs_values='exprs',
    n_neighbors=min(cfg$neighbors,n-1), n_threads=1, n_sgd_threads=1)
  set.seed(cfg$seed)
  object <- scater::runTSNE(object, subset_row=selected, exprs_values='exprs',
    perplexity=perplexity, num_threads=1)
  list(UMAP=SingleCellExperiment::reducedDim(object,'UMAP'),
       TSNE=SingleCellExperiment::reducedDim(object,'TSNE'), cell_id=colnames(x)[idx],
       markers=rownames(x)[selected], perplexity=perplexity, neighbors=min(cfg$neighbors,n-1),
       engine='scater::runUMAP / scater::runTSNE')
}
