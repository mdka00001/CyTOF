read_table <- function(path) {
  if (!file.exists(path)) stop('Missing file: ', path)
  first <- readLines(path, n=1, warn=FALSE)
  sep <- if (grepl('\t', first, fixed=TRUE)) '\t' else ','
  utils::read.table(path, header=TRUE, sep=sep, quote='"', comment.char='',
                    check.names=FALSE, stringsAsFactors=FALSE)
}
assert_unique <- function(x, label) {
  if (anyNA(x) || any(!nzchar(as.character(x))) || anyDuplicated(x)) stop(label, ' must be nonempty and unique')
}
read_input <- function(cfg) {
  panel <- read_table(file.path(cfg$input, 'panel.csv'))
  if (!all(c('channel','name') %in% names(panel))) stop('panel.csv requires channel,name')
  if ('keep' %in% names(panel)) {
    k <- tolower(as.character(panel$keep))
    if (anyNA(k) || any(!k %in% c('0','1','true','false'))) stop('panel keep must be 0/1 or TRUE/FALSE')
    panel <- panel[k %in% c('1','true'),,drop=FALSE]
  }
  if (!nrow(panel)) stop('No retained panel channels')
  markers <- as.character(panel$name); assert_unique(markers, 'Panel marker names')
  assert_unique(panel$channel, 'Panel channels')
  use <- rep(TRUE, nrow(panel))
  if ('use_channel' %in% names(panel)) {
    u <- tolower(as.character(panel$use_channel))
    if (anyNA(u) || any(!u %in% c('true','false','1','0'))) stop('use_channel must be TRUE/FALSE or 1/0')
    use <- u %in% c('true','1')
  }
  excluded <- split_list(cfg$exclude)
  if (any(!excluded %in% markers)) stop('Unknown excluded marker')
  use <- use & !markers %in% excluded
  images <- read_table(file.path(cfg$input, 'images.csv'))
  if (!all(c('image','width_px','height_px') %in% names(images))) stop('images.csv requires image,width_px,height_px')
  if (!nrow(images)) stop('No images listed')
  images$roi_id <- tools::file_path_sans_ext(basename(images$image))
  assert_unique(images$roi_id, 'Image ROI IDs')
  for (k in c('width_px','height_px')) if (!is.numeric(images[[k]]) || any(!is.finite(images[[k]]) | images[[k]] <= 0)) stop('Invalid image dimensions')
  mapping <- read_table(cfg$samples)
  if (!all(c('roi_id','sample_name') %in% names(mapping))) stop('Mapping requires roi_id,sample_name')
  assert_unique(mapping$roi_id, 'Mapping ROI IDs')
  if (anyNA(mapping$sample_name) || any(!nzchar(mapping$sample_name))) stop('Empty sample name')
  reserved <- c('cell_id','Object','area','width_px','height_px','image','kept')
  if (any(names(mapping) %in% reserved)) stop('Mapping uses a reserved metadata column')
  mi <- match(images$roi_id, mapping$roi_id)
  if (anyNA(mi)) stop('Unmapped ROI IDs: ', paste(images$roi_id[is.na(mi)], collapse=', '))
  if (any(!mapping$roi_id %in% images$roi_id)) warning('Mapping contains unused ROI IDs')
  mats <- metas <- vector('list', nrow(images))
  for (i in seq_len(nrow(images))) {
    roi <- images$roi_id[i]
    x <- read_table(file.path(cfg$input, cfg$intensities, paste0(roi, '.csv')))
    props <- read_table(file.path(cfg$input, cfg$regionprops, paste0(roi, '.csv')))
    if (!'Object' %in% names(x) || !all(c('Object','area') %in% names(props))) stop('Tables require Object, and regionprops requires area: ', roi)
    assert_unique(x$Object, paste(roi, 'intensity Object')); assert_unique(props$Object, paste(roi, 'regionprops Object'))
    idx <- match(x$Object, props$Object)
    if (anyNA(idx) || nrow(props) != nrow(x)) stop('Intensity/regionprops objects do not match: ', roi)
    cols <- if (all(markers %in% names(x))) markers else as.character(panel$channel)
    if (!all(cols %in% names(x))) stop('Intensity columns must match panel names or channels: ', roi)
    if (!all(vapply(x[cols], is.numeric, logical(1)))) stop('Non-numeric intensities: ', roi)
    mat <- t(as.matrix(x[cols])); rownames(mat) <- markers
    if (any(!is.finite(mat)) || any(mat < 0)) stop('Intensities must be finite and nonnegative')
    if (any(setdiff(names(mapping), 'roi_id') %in% names(props))) stop('Mapping columns collide with region properties')
    meta <- props[idx,,drop=FALSE]
    if (!is.numeric(meta$area) || any(!is.finite(meta$area) | meta$area < 0)) stop('Invalid cell areas')
    meta$cell_id <- paste(roi, x$Object, sep='::')
    for (k in names(mapping)) meta[[k]] <- rep(mapping[[k]][mi[i]], nrow(meta))
    meta$width_px <- rep(images$width_px[i], nrow(meta)); meta$height_px <- rep(images$height_px[i], nrow(meta))
    colnames(mat) <- meta$cell_id
    mats[[i]] <- mat; metas[[i]] <- meta
  }
  counts <- do.call(cbind, mats); meta <- do.call(rbind, metas)
  if (!ncol(counts)) stop('No cells found')
  assert_unique(meta$cell_id, 'Cell IDs')
  rownames(meta) <- meta$cell_id
  list(counts=counts, cells=meta, panel=panel, use=use, images=images)
}
write_csv <- function(x, path) utils::write.csv(x, path, row.names=FALSE, na='NA')
write_matrix <- function(x, path, chunk=1000L) {
  con <- gzfile(path, 'wt'); on.exit(close(con))
  for (start in seq.int(1L, ncol(x), by=chunk)) {
    idx <- start:min(ncol(x), start+chunk-1L)
    d <- data.frame(cell_id=colnames(x)[idx], t(x[,idx,drop=FALSE]), check.names=FALSE)
    utils::write.table(d, con, sep=',', row.names=FALSE, col.names=start==1L, append=start!=1L)
  }
}
