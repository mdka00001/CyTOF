fcs_help <- function() {
  cat(paste0(
    'Usage: cytof-qc fcs --qc-dir DIR --output DIR [options]\n\n',
    'Build FCS files from a completed cytof-qc run.\n',
    '  --base arcsinh              Normalized assay used as the source\n',
    '  --scales minmax,zscore      FCS expression scales to export\n',
    '  --area-min 5                Inclusive cell-area lower bound\n',
    '  --area-max 250              Inclusive cell-area upper bound\n',
    '  --split-by time_group       Also write one FCS per group (if present)\n',
    '  --help                      Show this help\n\n',
    'Every scale also produces one combined FCS. Event metadata is in matching CSV files.\n'
  ))
}
parse_fcs_args <- function(args) {
  cfg <- list(qc_dir=NULL,output=NULL,base='arcsinh',scales='minmax,zscore',
              area_min=5,area_max=250,split_by='time_group')
  i <- 1L
  while (i<=length(args)) {
    key <- gsub('-','_',sub('^--','',args[i]))
    if (key=='help') { fcs_help(); quit(status=0) }
    if (!key %in% names(cfg)) stop('Unknown fcs option: ',args[i])
    i <- i+1L
    if (i>length(args) || startsWith(args[i],'--')) stop('Missing value for --',gsub('_','-',key))
    cfg[[key]] <- if (key %in% c('area_min','area_max')) suppressWarnings(as.numeric(args[i])) else args[i]
    i <- i+1L
  }
  if (is.null(cfg$qc_dir) || is.null(cfg$output)) stop('Required: cytof-qc fcs --qc-dir DIR --output DIR')
  if (!is.finite(cfg$area_min) || !is.finite(cfg$area_max) || cfg$area_min<0 || cfg$area_max<cfg$area_min) stop('Area bounds must satisfy 0 <= --area-min <= --area-max')
  cfg$scales <- unique(split_list(cfg$scales))
  if (!length(cfg$scales) || any(!cfg$scales %in% c('minmax','zscore'))) stop('--scales accepts minmax,zscore')
  if (!cfg$base %in% c('arcsinh','rlog')) stop('--base accepts arcsinh or rlog')
  cfg
}
fcs_channel_names <- function(x) {
  names <- make.names(x,unique=TRUE,allow_=TRUE)
  names <- substr(names,1,64)
  if (anyDuplicated(toupper(names))) stop('FCS parameter names collide after normalization')
  names
}
minmax_markers <- function(x) {
  lo <- apply(x,1,min); span <- apply(x,1,function(v) max(v)-min(v))
  span[!is.finite(span) | span==0] <- 1
  sweep(sweep(x,1,lo,'-'),1,span,'/')
}
fcs_event_data <- function(expr, cells) {
  # FCS stores numeric channels. Preserve original categorical labels in the
  # companion event CSV; add numeric event channels for each metadata field.
  channels <- list()
  for (field in names(cells)) {
    if (field %in% c('cell_id','kept')) next
    value <- cells[[field]]
    if (is.numeric(value) || is.integer(value) || is.logical(value)) {
      number <- as.numeric(value)
      if (any(!is.finite(number))) stop('FCS metadata contains missing/non-finite values in ',field)
    } else {
      labels <- as.character(value)
      if (anyNA(labels)) stop('FCS metadata contains missing values in ',field)
      number <- match(labels,unique(labels))
    }
    channels[[field]] <- number
  }
  channel_matrix <- do.call(cbind,channels)
  if (is.null(dim(channel_matrix))) channel_matrix <- matrix(channel_matrix,ncol=1)
  colnames(channel_matrix) <- names(channels)
  values <- cbind(t(expr),channel_matrix)
  if (any(!is.finite(values))) stop('FCS event values must all be finite')
  storage.mode(values) <- 'double'
  values
}
write_fcs_file <- function(values, cells, path) {
  original_names <- colnames(values)
  channel_names <- fcs_channel_names(original_names)
  colnames(values) <- channel_names
  min_values <- apply(values,2,min)
  max_values <- apply(values,2,max)
  parameters <- data.frame(name=channel_names,desc=original_names,
                           range=as.integer(pmax(1,ceiling(max_values-min_values+1))),
                           minRange=min_values,maxRange=max_values,
                           row.names=channel_names,check.names=FALSE)
  frame <- flowCore::flowFrame(exprs=values,
    parameters=Biobase::AnnotatedDataFrame(parameters))
  flowCore::write.FCS(frame,filename=path,what='numeric')
  sidecar <- data.frame(event_index=seq_len(nrow(values)),
    cell_id=cells$cell_id,
    cells[,setdiff(names(cells),c('cell_id','kept')),drop=FALSE],
    check.names=FALSE)
  write_csv(sidecar,sub('\\.fcs$','_events.csv',path,ignore.case=TRUE))
}
safe_group_name <- function(x) {
  out <- gsub('[^A-Za-z0-9._-]+','_',as.character(x))
  out[!nzchar(out)] <- 'empty'
  out
}
run_fcs <- function(cfg) {
  if (!requireNamespace('flowCore',quietly=TRUE)) stop('Missing package flowCore. Install it with Rscript scripts/install.R.')
  input <- normalizePath(cfg$qc_dir,mustWork=TRUE)
  if (!file.exists(file.path(input,'STATUS')) || readLines(file.path(input,'STATUS'),n=1L,warn=FALSE)!='complete') stop('--qc-dir must be a completed cytof-qc output directory')
  object_path <- file.path(input,'objects','qc_data.rds')
  if (!file.exists(object_path)) stop('Missing objects/qc_data.rds in --qc-dir')
  data <- readRDS(object_path)
  cells <- as.data.frame(data$cells,stringsAsFactors=FALSE)
  if (exists('plot_metadata',mode='function')) cells <- plot_metadata(cells)
  if (!all(c('cell_id','area') %in% names(cells))) stop('QC object must contain cell_id and area')
  if (!is.numeric(cells$area) || any(!is.finite(cells$area))) stop('Cell area must be finite and numeric')
  keep <- cells$area>=cfg$area_min & cells$area<=cfg$area_max
  if (!any(keep)) stop('No cells remain in the requested area range')
  cells <- cells[keep,,drop=FALSE]
  if (!cfg$split_by %in% c('','none') && !cfg$split_by %in% names(cells)) {
    warning('Split column ',cfg$split_by,' is absent; writing combined FCS files only.')
    cfg$split_by <- ''
  }
  if (dir.exists(cfg$output) && length(list.files(cfg$output,all.files=TRUE,no..=TRUE))) stop('FCS output directory must be empty; choose a new directory')
  dir.create(cfg$output,recursive=TRUE,showWarnings=FALSE)
  out <- normalizePath(cfg$output)
  base_path <- file.path(input,'objects',paste0(cfg$base,'.rds'))
  if (!file.exists(base_path)) stop('Missing normalized matrix: objects/',basename(base_path))
  base <- readRDS(base_path)
  if (!identical(colnames(base),as.character(data$cells$cell_id))) stop('Normalized assay cell IDs do not match qc_data.rds')
  base <- base[,keep,drop=FALSE]
  for (scale in cfg$scales) {
    expr <- if (scale=='minmax') minmax_markers(base) else zscore(base)
    values <- fcs_event_data(expr,cells)
    destinations <- list(all=seq_len(nrow(cells)))
    if (nzchar(cfg$split_by)) {
      group <- as.character(cells[[cfg$split_by]])
      groups <- split(seq_along(group),factor(group,levels=unique(group)))
      safe <- safe_group_name(names(groups))
      if (anyDuplicated(safe)) stop('Group labels collide after converting to safe FCS filenames')
      for (j in seq_along(groups)) destinations[[paste0(cfg$split_by,'_',safe[j])]] <- groups[[j]]
    }
    for (name in names(destinations)) {
      idx <- destinations[[name]]
      stem <- paste0(name,'_',scale)
      write_fcs_file(values[idx,,drop=FALSE],cells[idx,,drop=FALSE],file.path(out,paste0(stem,'.fcs')))
    }
    message('Wrote ',length(destinations),' ',scale,' FCS file(s) with ',nrow(base),' markers and ',nrow(cells),' filtered cells.')
  }
  write_csv(data.frame(cell_id=cells$cell_id,area=cells$area,kept=TRUE),file.path(out,'included_cells.csv'))
  writeLines(c('FCS export from cytof-qc',paste('QC run:',input),paste('Source assay:',cfg$base),
    paste('Cell area range:',cfg$area_min,'to',cfg$area_max),
    paste('Scales:',paste(cfg$scales,collapse=', ')),
    paste('Group split:',if(nzchar(cfg$split_by)) cfg$split_by else 'none'),
    'Each .fcs has a matching _events.csv file. Categorical metadata are encoded as numeric FCS channels; original text values are preserved in the event CSV.'),
    file.path(out,'README.txt'))
  message('FCS export complete: ',out)
  invisible(out)
}
