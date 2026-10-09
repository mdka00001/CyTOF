# Preserve raw TIFF values, especially 16-bit segmentation labels.
read_raw_tiff <- function(path) {
  info <- tiff::readTIFF(path, all=TRUE, payload=FALSE)
  if (any(!info$bits.per.sample %in% c(8,16,32))) stop('Expected 8/16-bit integer or 32-bit float TIFF: ',path)
  pages <- lapply(seq_len(nrow(info)), function(i)
    tiff::readTIFF(path, all=i, as.is=info$bits.per.sample[i] <= 16)[[1]])
  if (!length(pages) || any(vapply(pages,function(p) length(dim(p))!=2,logical(1)))) stop('Expected grayscale TIFF pages: ',path)
  shape <- dim(pages[[1]])
  if (any(vapply(pages,function(p) !identical(dim(p),shape),logical(1)))) stop('Inconsistent TIFF page sizes: ',path)
  if (length(pages)==1) return(t(pages[[1]]))
  a <- array(0,c(rev(shape),length(pages)))
  for (i in seq_along(pages)) a[,,i] <- t(pages[[i]])
  a
}

image_qc <- function(data, cfg, report) {
  rows <- list(); k <- 0L
  channels <- split_list(cfg$image_channels)
  if (!length(channels)) channels <- head(as.character(data$panel$name)[data$use],3)
  if (!length(channels) || length(channels)>3 || any(!channels %in% data$panel$name)) stop('--image-channels requires 1 to 3 panel marker names')
  for (i in seq_len(nrow(data$images))) {
    roi <- data$images$roi_id[i]
    path <- file.path(cfg$input,cfg$images,basename(data$images$image[i]))
    maskpath <- file.path(cfg$input,cfg$masks,basename(data$images$image[i]))
    if (!file.exists(path) || !file.exists(maskpath)) stop('Missing image/mask for ',roi,'; configure --images/--masks or --skip-images')
    img <- read_raw_tiff(path); mask <- read_raw_tiff(maskpath)
    if (length(dim(img))==2) dim(img) <- c(dim(img),1)
    if (length(dim(img))!=3 || dim(img)[3]!=nrow(data$panel)) stop('Image channels differ from panel order: ',roi)
    if (length(dim(mask))!=2 || any(dim(mask)!=dim(img)[1:2])) stop('Image/mask dimensions differ: ',roi,'; masks are never silently resized')
    if (any(dim(img)[1:2] != c(data$images$width_px[i],data$images$height_px[i]))) stop('Image dimensions differ from images.csv: ',roi)
    if (any(!is.finite(img)) || any(img<0) || any(!is.finite(mask)) || any(mask<0) || any(mask!=floor(mask))) stop('Invalid image intensities or mask labels: ',roi)
    ids <- data$cells$Object[data$cells$roi_id==roi]
    if (!setequal(unique(as.vector(mask[mask>0])),ids)) stop('Mask labels differ from cell Object IDs: ',roi)
    for (j in seq_len(dim(img)[3])) {
      v <- as.vector(img[,,j]); result <- c(signal=NA,noise=NA,snr=NA); threshold <- NA_real_
      if (diff(range(v))>0) {
        threshold <- EBImage::otsu(EBImage::Image(img[,,j]),range=range(v),levels=65536)
        result <- snr_partition(v,v>threshold)
      }
      k <- k+1L; rows[[k]] <- data.frame(roi_id=roi,marker=data$panel$name[j],t(result),threshold=threshold)
    }
    nr <- nrow(mask); nc <- ncol(mask); edge <- matrix(FALSE,nr,nc)
    if (nr>1) { d <- mask[-1,,drop=FALSE]!=mask[-nr,,drop=FALSE]; edge[-1,] <- edge[-1,]|d; edge[-nr,] <- edge[-nr,]|d }
    if (nc>1) { d <- mask[,-1,drop=FALSE]!=mask[,-nc,drop=FALSE]; edge[,-1] <- edge[,-1]|d; edge[,-nc] <- edge[,-nc]|d }
    edge <- edge & mask>0
    # Original script's plotPixels approach, including its white outline channel.
    display <- array(0,c(nr,nc,length(channels)+1L))
    for (j in seq_along(channels)) display[,,j] <- img[,,match(channels[j],data$panel$name)]
    display[,,length(channels)+1L] <- edge*.03
    preview <- cytomapper::CytoImageList(setNames(list(EBImage::Image(display)),roi))
    cytomapper::channelNames(preview) <- c(channels,'CellOutline')
    channel_colors <- setNames(lapply(c('red','blue','green')[seq_along(channels)],function(color) c('black',color)),channels)
    channel_colors$CellOutline <- c('black','white')
    report_plot(report,paste('Segmentation:',roi),paste('cytomapper plotPixels, following the reference script. Channels:',paste(channels,collapse=', '),
      '. The displayed color key identifies every channel and the white cell outlines. Image channels are normalized for display by plotPixels; raw values used for SNR are unchanged.'),function() {
        cytomapper::plotPixels(preview,colour_by=c(channels,'CellOutline'),colour=channel_colors,
          missing_colour='white',bcg=list(CellOutline=c(0,.3,1)),
          image_title=list(text=paste('ROI:',roi),cex=.8),
          legend=list(colour_by.title.cex=.8,colour_by.labels.cex=.8))
      },width=10,height=10)

  }
  do.call(rbind,rows)
}
