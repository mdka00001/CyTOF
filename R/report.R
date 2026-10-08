escape_html <- function(x) {
  x <- gsub('&', '&amp;', x, fixed=TRUE); x <- gsub('<','&lt;',x,fixed=TRUE)
  gsub('>', '&gt;', x, fixed=TRUE)
}
new_report <- function(out) {
  r <- new.env(parent=emptyenv()); r$out <- out; r$sections <- list()
  grDevices::pdf(file.path(out,'report.pdf'), width=11, height=8.5, onefile=TRUE)
  r$pdf <- grDevices::dev.cur(); r
}
report_text <- function(r, title, text) {
  grDevices::dev.set(r$pdf); graphics::par(mar=c(1,1,3,1)); graphics::plot.new()
  graphics::title(main=title)
  lines <- unlist(lapply(text, function(t) c(strwrap(t, width=110), '')))
  graphics::text(0, 1, paste(lines,collapse='\n'), adj=c(0,1), cex=.85)
  r$sections[[length(r$sections)+1L]] <- list(title=title,text=paste(text,collapse='\n'),file=NULL)
}
report_plot <- function(r, title, explanation, draw) {
  name <- sprintf('plot_%04d.png', length(r$sections)+1L)
  render <- function() {
    graphics::par(mar=c(7,5,4,2)); draw()
  }
  grDevices::dev.set(r$pdf); render()
  # Explanation on a separate PDF page avoids clipping long text or plot labels.
  grDevices::png(file.path(r$out,'plots',name), width=1400,height=950,res=140)
  tryCatch(render(), finally=grDevices::dev.off())
  report_text(r, title, explanation)
  r$sections[[length(r$sections)]]$file <- paste0('plots/',name)
}
finish_report <- function(r) {
  grDevices::dev.off(r$pdf)
  parts <- c('<!doctype html><html lang="en"><meta charset="utf-8"><title>CyTOF QC</title>',
    '<style>body{font:16px system-ui;max-width:1100px;margin:40px auto;padding:20px}img{max-width:100%}pre{white-space:pre-wrap}section{border-top:1px solid #ddd;margin-top:30px}</style><body><h1>CyTOF QC report</h1>')
  for (s in r$sections) parts <- c(parts, paste0('<section><h2>',escape_html(s$title),'</h2><pre>',escape_html(s$text),'</pre>'),
    if (!is.null(s$file)) paste0('<img alt="',escape_html(s$title),'" src="',s$file,'">'), '</section>')
  writeLines(c(parts,'</body></html>'), file.path(r$out,'report.html'))
}
plot_snr <- function(d, title) {
  ok <- is.finite(d$snr) & d$snr > 0 & is.finite(d$signal) & d$signal > 0
  if (!any(ok)) { graphics::plot.new(); graphics::title(paste(title,'(no finite values)')); return(invisible(NULL)) }
  graphics::plot(log2(d$signal[ok]),log2(d$snr[ok]),pch=19,xlab='Positive signal (log2)',ylab='SNR (log2)',main=title)
  graphics::text(log2(d$signal[ok]),log2(d$snr[ok]),labels=d$marker[ok],pos=3,cex=.65)
}
plot_heatmap <- function(x, title, labels=FALSE) {
  if (!length(x)) { graphics::plot.new(); graphics::title('No selected markers'); return() }
  if (nrow(x)>1) x <- x[stats::hclust(stats::dist(x))$order,,drop=FALSE]
  if (ncol(x)>1) x <- x[,stats::hclust(stats::dist(t(x)))$order,drop=FALSE]
  graphics::image(seq_len(ncol(x)),seq_len(nrow(x)),t(x),col=grDevices::hcl.colors(64,'Viridis'),axes=FALSE,xlab='Cells / ROIs',ylab='',main=title)
  graphics::axis(2,at=seq_len(nrow(x)),labels=rownames(x),las=2,cex.axis=.6)
  if (labels) graphics::axis(1,at=seq_len(ncol(x)),labels=colnames(x),las=2,cex.axis=.5)
}
plot_embedding <- function(coords, values, title, continuous=FALSE) {
  if (continuous) {
    breaks <- range(values); bins <- if (diff(breaks)==0) rep(32,length(values)) else 1+floor(63*(values-breaks[1])/diff(breaks))
    colors <- grDevices::hcl.colors(64,'Viridis')[bins]
  } else { f <- factor(values); pal <- grDevices::hcl.colors(nlevels(f),'Dark 3'); colors <- pal[f] }
  graphics::plot(coords,pch=16,cex=.3,col=colors,xlab='Dimension 1',ylab='Dimension 2',main=title)
  if (!continuous) graphics::legend('topright',legend=levels(f),col=pal,pch=16,cex=.5)
  else graphics::legend('topright',legend=format(range(values),digits=3),col=grDevices::hcl.colors(64,'Viridis')[c(1,64)],pch=16,cex=.7)
}
