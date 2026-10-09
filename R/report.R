escape_html <- function(x) {
  x <- gsub('&', '&amp;', x, fixed=TRUE); x <- gsub('<','&lt;',x,fixed=TRUE)
  gsub('>', '&gt;', x, fixed=TRUE)
}
new_report <- function(out) {
  if (!nzchar(Sys.which('pdfunite'))) stop('PDF assembly requires pdfunite (Poppler utilities).')
  r <- new.env(parent=emptyenv()); r$out <- out; r$sections <- list(); r$pages <- character()
  r$page_dir <- file.path(out, '.report-pages'); dir.create(r$page_dir, showWarnings=FALSE)
  r
}
report_page <- function(r, draw, width=11, height=8.5) {
  path <- file.path(r$page_dir, sprintf('page_%04d.pdf', length(r$pages)+1L))
  grDevices::pdf(path, width=width, height=height, onefile=TRUE, useDingbats=FALSE)
  tryCatch(draw(), finally=grDevices::dev.off())
  r$pages <- c(r$pages, path)
}
report_text <- function(r, title, text) {
  report_page(r, function() {
    graphics::par(mar=c(1,1,3,1)); graphics::plot.new(); graphics::title(main=title)
    lines <- unlist(lapply(text, function(t) c(strwrap(t, width=110), '')))
    graphics::text(0, 1, paste(lines,collapse='\n'), adj=c(0,1), cex=.85)
  })
  r$sections[[length(r$sections)+1L]] <- list(title=title,text=paste(text,collapse='\n'),file=NULL)
}
report_plot <- function(r, title, explanation, draw, width=11, height=8.5) {
  name <- sprintf('plot_%04d.png', length(r$sections)+1L)
  render <- function() { graphics::par(mar=c(7,5,4,2)); draw() }
  report_page(r, render, width, height)
  # Match PDF page proportions so a multi-panel page also stays together in HTML.
  grDevices::png(file.path(r$out,'plots',name), width=width, height=height, units='in', res=100)
  tryCatch(render(), finally=grDevices::dev.off())
  report_text(r, title, explanation)
  r$sections[[length(r$sections)]]$file <- paste0('plots/',name)
}
finish_report <- function(r) {
  status <- system2(Sys.which('pdfunite'), shQuote(c(r$pages, file.path(r$out,'report.pdf'))))
  if (status != 0L) stop('Unable to assemble report.pdf')
  unlink(r$page_dir, recursive=TRUE)
  parts <- c('<!doctype html><html lang="en"><meta charset="utf-8"><title>CyTOF QC</title>',
    '<style>body{font:16px system-ui;max-width:1100px;margin:40px auto;padding:20px}img{max-width:100%}pre{white-space:pre-wrap}section{border-top:1px solid #ddd;margin-top:30px}</style><body><h1>CyTOF QC report</h1>')
  for (s in r$sections) parts <- c(parts, paste0('<section><h2>',escape_html(s$title),'</h2><pre>',escape_html(s$text),'</pre>'),
    if (!is.null(s$file)) paste0('<a href="',s$file,'"><img alt="',escape_html(s$title),'" src="',s$file,'"></a>'), '</section>')
  writeLines(c(parts,'</body></html>'), file.path(r$out,'report.html'))
}
plot_snr <- function(d, title) {
  ok <- is.finite(d$snr) & d$snr > 0 & is.finite(d$signal) & d$signal > 0
  if (!any(ok)) { graphics::plot.new(); graphics::title(paste(title,'(no finite values)')); return(invisible(NULL)) }
  graphics::plot(log2(d$signal[ok]),log2(d$snr[ok]),pch=19,xlab='Positive signal (log2)',ylab='SNR (log2)',main=title)
  graphics::text(log2(d$signal[ok]),log2(d$snr[ok]),labels=d$marker[ok],pos=3,cex=.65)
}
