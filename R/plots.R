# The plotting stack and layouts follow 1.qc_script_1.R / 1.qc_script_2.R.
plot_metadata <- function(cells) {
  cells$sample_id <- cells$roi_id
  if (!'patient_id' %in% names(cells)) cells$patient_id <- cells$sample_name
  if (!'ROI' %in% names(cells)) {
    ids <- unique(cells$roi_id); labels <- sub('.*[_.-]', '', ids)
    if (anyDuplicated(labels)) labels <- ids
    cells$ROI <- labels[match(cells$roi_id, ids)]
  }
  for (field in intersect(c('patient_id','ROI','sample_id','sample_name','indication','batch'),names(cells)))
    cells[[field]] <- factor(cells[[field]], levels=unique(as.character(cells[[field]])))
  cells
}
plot_colors <- function(cells) {
  fields <- intersect(c('patient_id','ROI','sample_id','sample_name','indication','batch'),names(cells))
  setNames(lapply(fields, function(field) {
    lev <- levels(factor(cells[[field]]))
    scheme <- if (field=='ROI') 'BrBG' else if (field=='indication') 'Set2' else 'Set1'
    base <- RColorBrewer::brewer.pal(if (scheme=='BrBG') 11 else if (scheme=='Set2') 8 else 9, scheme)
    colors <- if (length(lev)<=length(base)) base[seq_along(lev)] else grDevices::colorRampPalette(base)(length(lev))
    setNames(colors,lev)
  }), fields)
}
expression_object <- function(x, cells) {
  stopifnot(identical(colnames(x),as.character(cells$cell_id)))
  rownames(cells) <- colnames(x)
  SingleCellExperiment::SingleCellExperiment(assays=list(exprs=x),colData=S4Vectors::DataFrame(cells))
}
script_heatmap <- function(object, markers, title, colors, aggregate=FALSE) {
  if (!length(markers)) { graphics::plot.new(); graphics::title('No selected markers'); return(invisible(NULL)) }
  fields <- if (aggregate) intersect(c('indication','patient_id','ROI'),names(SummarizedExperiment::colData(object))) else 'patient_id'
  args <- list(object=object, genes=markers, assay='exprs', scale='none',
    cluster_cols=ncol(object)>1, cluster_rows=length(markers)>1,
    annot.by=fields, annotation_colors=colors[fields],
    heatmap.colors=if (aggregate) viridis::viridis(100) else grDevices::colorRampPalette(c('blue','white','red'))(100),
    show_colnames=aggregate, fontsize=10, fontsize_row=10, fontsize_col=7,
    main=title, silent=TRUE)
  if (!aggregate) {
    args$breaks <- seq(-4,4,length.out=101)
    args$legend_breaks <- c(-4,-2,0,2,4)
  }
  heat <- do.call(dittoSeq::dittoHeatmap,args)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(width=.9,height=.9))
  grid::grid.draw(heat$gtable)
  grid::popViewport()
  invisible(heat)
}
marker_distributions <- function(object, markers, colors, type) {
  plots <- dittoSeq::multi_dittoPlot(object, vars=markers, group.by='patient_id',
    plots=type, assay='exprs', color.panel=unname(colors$patient_id),
    list.out=TRUE, legend.show=FALSE)
  plots <- lapply(plots, function(p) p + ggplot2::theme_minimal(base_size=11) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1)))
  patchwork::wrap_plots(plots, ncol=min(4L,length(markers)))
}
embedding_metadata_grid <- function(object, colors, fields) {
  plots <- unlist(lapply(fields, function(field) lapply(c('UMAP','TSNE'), function(kind) {
    ggrastr::rasterise(dittoSeq::dittoDimPlot(object,var=field,reduction.use=kind,size=.2), layers='Point',dpi=150) +
      ggplot2::scale_color_manual(values=colors[[field]],drop=FALSE) +
      ggplot2::ggtitle(paste(field,kind)) + ggplot2::theme_minimal(base_size=11)
  })),recursive=FALSE)
  patchwork::wrap_plots(plots,ncol=2,guides='keep')
}
embedding_marker_grid <- function(object, kind) {
  plots <- lapply(rownames(object),function(marker) {
    ggrastr::rasterise(dittoSeq::dittoDimPlot(object,var=marker,reduction.use=kind,
      assay='exprs',size=.2),layers='Point',dpi=150) +
      ggplot2::scale_color_viridis_c(name=marker) + ggplot2::ggtitle(marker) +
      ggplot2::theme_minimal(base_size=11) +
      ggplot2::theme(plot.title=ggplot2::element_text(size=12,face='bold',hjust=.5),
        plot.margin=ggplot2::margin(14,18,14,18,unit='pt'),panel.grid=ggplot2::element_blank(),
        legend.position='right',legend.title=ggplot2::element_text(size=9),
        legend.text=ggplot2::element_text(size=8),legend.key.height=grid::unit(1.2,'cm'))
  })
  patchwork::wrap_plots(plots,ncol=min(4L,length(plots)),guides='keep')
}
report_expression <- function(report, y, cells, use, method, cfg, hidx, embeddings=NULL) {
  cells <- plot_metadata(cells); colors <- plot_colors(cells)
  object <- expression_object(y,cells); markers <- rownames(y)[use]
  report_plot(report,paste(method,'cell heatmap'),sprintf('dittoHeatmap with patient annotations, blue-white-red colors, and the -4 to 4 display range used in the reference scripts. Values outside that range saturate. %d seeded cells; selected markers only; scale=none.',length(hidx)),
    function() script_heatmap(object[,hidx],markers,paste(method,'cell heatmap'),colors),width=12,height=10)
  groups <- split(seq_len(ncol(y)),cells$roi_id)
  means <- vapply(groups,function(j) rowMeans(y[,j,drop=FALSE]),numeric(nrow(y)))
  dimnames(means) <- list(rownames(y),names(groups))
  write_csv(data.frame(marker=rownames(means),means,check.names=FALSE),file.path(report$out,'tables',paste0(method,'_roi_means.csv')))
  meta <- cells[vapply(groups,function(j) j[1],integer(1)),,drop=FALSE]; meta$cell_id <- colnames(means)
  aggregated <- expression_object(means,meta)
  report_plot(report,paste(method,'ROI means'),'dittoHeatmap of mean transformed expression by ROI, with viridis intensity legend and patient/ROI annotations, as in the reference scripts. Biological composition and staining can both produce differences.',
    function() script_heatmap(aggregated,markers,paste(method,'ROI means'),colors,TRUE),width=16,height=12)
  if (length(markers)) for (type in c('boxplot','ridgeplot')) {
    p <- marker_distributions(object,markers,colors,type)
    report_plot(report,paste(method,'all-marker',type),'All selected markers grouped by patient_id (sample_name is the fallback if patient_id is absent). Uses dittoSeq multi_dittoPlot. All retained cells are included; each panel keeps its own expression axis. Differences can reflect biology or staining.',
      function() print(p),width=5*min(4,length(markers)),height=4*ceiling(length(markers)/4))
  }
  if (!is.null(embeddings)) {
    idx <- match(embeddings$cell_id,cells$cell_id)
    if (anyNA(idx)) stop('Embedding cell IDs do not match retained cells')
    embedded <- object[,idx]
    for (kind in c('UMAP','TSNE')) SingleCellExperiment::reducedDim(embedded,kind) <- embeddings[[kind]]
    # Plot ALL embedded cells: taking the first N of sorted IDs hid later patients.
    fields <- intersect(c('patient_id','ROI','indication','batch'),names(cells))
    p <- embedding_metadata_grid(embedded,colors,fields)
    report_plot(report,paste(method,'UMAP and t-SNE by metadata'),sprintf('dittoDimPlot with paired UMAP/t-SNE panels, as in the reference scripts. All %d embedded cells are drawn; no prefix truncation. Separation can reflect biology or technical effects.',ncol(embedded)),function() print(p),width=18,height=6*length(fields))
    write_csv(as.data.frame(table(patient_id=cells$patient_id[idx])),file.path(report$out,'tables',paste0(method,'_embedding_patient_counts.csv')))
    for (kind in c('UMAP','TSNE')) {
      coords <- embeddings[[kind]]; colnames(coords) <- paste0(kind,1:2)
      write_csv(data.frame(cell_id=embeddings$cell_id,roi_id=cells$roi_id[idx],sample_name=cells$sample_name[idx],patient_id=cells$patient_id[idx],coords),file.path(report$out,'coordinates',paste0(method,'_',kind,'.csv')))
      p <- embedding_marker_grid(embedded,kind)
      report_plot(report,paste(method,kind,'all markers'),'All markers on one page using dittoDimPlot and patchwork, with rasterized points and separate continuous viridis colorbars. All embedded cells are plotted, including markers not used to construct the embedding.',
        function() print(p),width=5*min(4,nrow(y)),height=4*ceiling(nrow(y)/4))
    }
  }
}
