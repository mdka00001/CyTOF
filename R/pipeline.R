run_qc <- function(cfg) {
  check_dependencies(cfg)
  if (dir.exists(cfg$output) && length(list.files(cfg$output,all.files=TRUE,no..=TRUE))) stop('Output directory must be empty; choose a new directory')
  data <- read_input(cfg)
  keep <- data$cells$area >= cfg$min_area
  if (sum(keep)<2) stop('Fewer than two cells remain after area filtering')
  x <- data$counts[,keep,drop=FALSE]; cells <- data$cells[keep,,drop=FALSE]
  if (any(grepl('rlog',cfg$transforms))) validate_rlog(x,cfg$rlog_max_cells)
  for (d in c('', 'plots','tables','matrices','objects','coordinates')) dir.create(file.path(cfg$output,d),recursive=TRUE,showWarnings=FALSE)
  out <- normalizePath(cfg$output)
  writeLines('running',file.path(out,'STATUS'))
  success <- FALSE
  on.exit(if (!success) writeLines('failed: outputs may be incomplete; see CLI error',file.path(out,'STATUS')),add=TRUE)
  set.seed(cfg$seed)
  saveRDS(cfg,file.path(out,'objects','config.rds'))
  capture.output(dput(cfg),file=file.path(out,'config.txt'))
  report <- new_report(out)
  on.exit(if (report$pdf %in% grDevices::dev.list()) grDevices::dev.off(report$pdf),add=TRUE)
  report_text(report,'Run overview',c(sprintf('%d markers; %d input cells; %d retained cells; %d ROIs.',nrow(x),ncol(data$counts),ncol(x),nrow(data$images)),
    paste('Transformations:',paste(cfg$transforms,collapse=', ')),
    sprintf('Cells with area < %g pixels excluded. Pixel size: %g micrometers. Seed: %d.',cfg$min_area,cfg$pixel_size,cfg$seed),
    'SNR flags are review suggestions, not automatic marker exclusions. High SNR can be inflated by near-zero background. Biological absence and technical failure cannot be distinguished by SNR alone.',
    'Workflow reference: https://bodenmillergroup.github.io/IMCDataAnalysis/image-and-cell-level-quality-control.html'))
  data$cells$kept <- keep
  write_csv(data$cells,file.path(out,'tables','cells_all.csv'))
  write_csv(data$panel,file.path(out,'tables','panel.csv'))
  write_matrix(data$counts,file.path(out,'matrices','counts_all.csv.gz'))
  write_matrix(x,file.path(out,'matrices','counts_retained.csv.gz'))
  metrics <- roi_metrics(data,cfg); write_csv(metrics,file.path(out,'tables','roi_qc.csv'))
  for (metric in c('coverage','cells_per_mm2','cells_before','cells_after')) {
    report_plot(report,paste('ROI',metric),switch(metric,
      coverage='Sum of segmented cell areas divided by image pixel area, before filtering. Low coverage may indicate sparse tissue or segmentation failure; values above one indicate inconsistent measurements.',
      cells_per_mm2='Retained cell count divided by ROI area in square millimeters using the configured pixel size. Density also depends on tissue and cell composition.',
      'Cell counts per ROI before or after the area filter; ROIs with no retained cells remain listed.'),function() {
        graphics::barplot(metrics[[metric]],names.arg=metrics$roi_id,las=2,cex.names=.5,ylab=metric,main=paste('ROI',metric))
      })
  }
  report_plot(report,'Cell area distributions','Pre-filter cell areas by ROI. Differences may reflect segmentation bias or biological cell-size differences. The red line is the configured minimum area.',function() {
    graphics::boxplot(area~roi_id,data=data$cells,las=2,cex.axis=.5,outline=FALSE,ylab='Area (pixels)',main='Cell area by ROI')
    graphics::abline(h=cfg$min_area,col='red',lty=2)
  })
  snr <- marker_review(cell_snr(data$counts,cfg),cfg)
  write_csv(snr,file.path(out,'tables','cell_snr_marker_review.csv'))
  report_plot(report,'Cell-level SNR','A two-component Gaussian mixture is fitted to arcsinh intensities in a seeded cell subset before filtering. Signal/background means and their ratio use raw intensities. Constant markers and failed mixtures are flagged. Infinite or undefined values are omitted from the scatter and retained in the table.',function() plot_snr(snr,'Cell-level SNR'))
  for (start in seq.int(1,nrow(snr),by=12)) {
    d <- snr[start:min(start+11,nrow(snr)),]
    report_text(report,'Marker review',sprintf('%s: SNR=%s; signal=%s. %s',d$marker,format(d$snr,digits=3),format(d$signal,digits=3),d$review))
  }
  if (!cfg$skip_images) {
    pixel <- image_qc(data,cfg,report); write_csv(pixel,file.path(out,'tables','pixel_snr.csv'))
    for (filtered in c(FALSE,TRUE)) {
      p <- if (filtered) pixel[!is.na(pixel$signal) & pixel$signal>cfg$min_signal,,drop=FALSE] else pixel
      avg <- do.call(rbind,lapply(split(p,p$marker),function(d) data.frame(marker=d$marker[1],signal=mean(d$signal,na.rm=TRUE),snr=mean(d$snr,na.rm=TRUE))))
      if (is.null(avg)) avg <- data.frame(marker=character(),signal=numeric(),snr=numeric())
      write_csv(avg,file.path(out,'tables',paste0('pixel_snr_summary_',filtered,'.csv')))
      report_plot(report,paste('Pixel-level SNR; signal filter',filtered),paste('Otsu threshold separates signal and background per image/channel. Arithmetic ROI means are plotted on log2 axes. Signal filtering removes ROI-marker pairs with positive signal <=',cfg$min_signal,'. This helps expose inflated SNR from very weak staining. Undefined/infinite points remain in CSV tables.'),function() plot_snr(avg,'Pixel-level SNR'))
    }
  } else report_text(report,'Image QC skipped','Image/mask previews and pixel-level SNR were explicitly disabled.')
  idx <- sort(sample.int(ncol(x),min(ncol(x),cfg$embedding_cells)))
  hidx <- sort(sample.int(ncol(x),min(ncol(x),cfg$heatmap_cells)))
  pidx <- sort(sample.int(ncol(x),min(ncol(x),cfg$plot_cells)))
  saveRDS(list(counts=x,cells=cells,panel=data$panel,selected_markers=data$panel$name[data$use],roi_qc=metrics,cell_snr=snr),file.path(out,'objects','qc_data.rds'))
  for (method in cfg$transforms) {
    message('Transforming: ',method)
    y <- transform_counts(x,method,cfg)
    write_matrix(y,file.path(out,'matrices',paste0(method,'.csv.gz')))
    saveRDS(y,file.path(out,'objects',paste0(method,'.rds')))
    report_plot(report,paste(method,'cell heatmap'),sprintf('Seeded subset of %d retained cells; selected markers only. Rows and columns are clustered within this bounded subset. Values are not additionally scaled. Z-scores describe relative expression across retained cells, not absolute signal.',length(hidx)),function() plot_heatmap(y[data$use,hidx,drop=FALSE],paste(method,'cell heatmap')))
    groups <- split(seq_len(ncol(y)),cells$roi_id)
    means <- vapply(groups,function(j) rowMeans(y[,j,drop=FALSE]),numeric(nrow(y)))
    dimnames(means) <- list(rownames(y),names(groups))
    write_csv(data.frame(marker=rownames(means),means,check.names=FALSE),file.path(out,'tables',paste0(method,'_roi_means.csv')))
    report_plot(report,paste(method,'ROI means'),'Mean of transformed retained-cell intensities per ROI (not transform of raw means). Differences may arise from biology, tissue composition or staining. Empty post-filter ROIs are absent here but present in the ROI QC table.',function() plot_heatmap(means[data$use,,drop=FALSE],paste(method,'ROI means'),TRUE))
    for (m in which(data$use)) {
      report_plot(report,paste(method,'distribution:',rownames(y)[m]),'Expression distributions by mapped sample in a bounded, seeded subset of retained cells. Shifts can indicate staining differences or different biological cell populations. Boxplots omit outlier points for readability.',function() {
        graphics::boxplot(split(y[m,pidx],cells$sample_name[pidx]),las=2,outline=FALSE,ylab=method,main=rownames(y)[m])
      })
    }
    if (!cfg$skip_embeddings) {
      emb <- compute_embeddings(y,data$use,cfg,idx)
      saveRDS(emb,file.path(out,'objects',paste0(method,'_embeddings.rds')))
      eidx <- seq_len(min(length(idx),cfg$plot_cells))
      for (kind in c('UMAP','TSNE')) {
        coords <- emb[[kind]]; colnames(coords) <- paste0(kind,1:2)
        write_csv(data.frame(cell_id=emb$cell_id,roi_id=cells$roi_id[idx],sample_name=cells$sample_name[idx],coords),file.path(out,'coordinates',paste0(method,'_',kind,'.csv')))
        for (field in intersect(c('sample_name','roi_id','patient_id','indication','batch'),names(cells))) {
          report_plot(report,paste(method,kind,field),sprintf('Embedding of %d retained cells using selected, variable markers and up to 30 PCs. Plot shows %d cells. Neighbor count=%d; t-SNE perplexity=%g. Separation may reflect biology or batch effects; distances and clusters are exploratory. Unsampled cells have no inferred coordinates.',length(idx),length(eidx),emb$neighbors,emb$perplexity),function() plot_embedding(coords[eidx,,drop=FALSE],cells[[field]][idx[eidx]],paste(method,kind,field)))
        }
        for (m in seq_len(nrow(y))) report_plot(report,paste(method,kind,rownames(y)[m]),'The same embedding is colored by transformed marker expression, including markers excluded from embedding construction. The legend shows minimum and maximum plotted values.',function() plot_embedding(coords[eidx,,drop=FALSE],y[m,idx[eidx]],paste(method,kind,rownames(y)[m]),TRUE))
      }
    }
    rm(y); invisible(gc(FALSE))
  }
  if (cfg$skip_embeddings) report_text(report,'Embeddings skipped','UMAP and t-SNE were explicitly disabled.')
  report_text(report,'Reproducibility and interpretation',c('See config.txt and sessionInfo.txt for parameters and package versions. Matrix CSV files have one cell per row and marker columns; cell_id joins all outputs. qc_data.rds contains retained raw data and QC metadata; transformations and embeddings are stored separately to limit memory duplication.',
    'Arcsinh uses the configured cofactor. Marker-wise z-scoring uses all retained cells and maps constant markers to zero. DESeq2 rlog, when requested, uses integer counts, poscounts size factors, a mean dispersion fit and blind=TRUE; its RNA-seq model is not validated as an IMC-specific normalization. No log1p approximation or intensity rounding is used.',
    'HTML references the plots directory; distribute both together. PDF contains the same explanations and figures. No automatic marker exclusion or batch correction is performed.'))
  finish_report(report)
  capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
  files <- list.files(out,recursive=TRUE,full.names=TRUE)
  files <- files[basename(files)!='STATUS']
  write_csv(data.frame(file=substring(files,nchar(out)+2),bytes=file.info(files)$size,md5=unname(tools::md5sum(files))),file.path(out,'manifest.csv'))
  writeLines('complete',file.path(out,'STATUS')); success <- TRUE
  message('Complete: ',out)
  invisible(out)
}
