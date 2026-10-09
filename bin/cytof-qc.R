#!/usr/bin/env Rscript
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
# R encodes spaces in --file arguments as ~+~. Decode before resolving modules.
script <- gsub('~+~', ' ', script, fixed=TRUE)
root <- dirname(dirname(normalizePath(script, mustWork=TRUE)))
for (f in c('config','io','transforms','qc','report','plots','images','pipeline','fcs')) source(file.path(root,'R',paste0(f,'.R')))
args <- commandArgs(TRUE)
if (length(args) && args[1]=='fcs') {
  fcs_args <- args[-1]
  if (!length(fcs_args) || '--help' %in% fcs_args || '-h' %in% fcs_args) { fcs_help(); quit(status=0) }
  tryCatch(run_fcs(parse_fcs_args(fcs_args)),error=function(e) { message('ERROR: ',conditionMessage(e)); quit(status=1) })
} else {
  if (!length(args) || '--help' %in% args || '-h' %in% args) { help_text(); quit(status=0) }
  tryCatch(run_qc(parse_args(args)),error=function(e) { message('ERROR: ',conditionMessage(e)); quit(status=1) })
}
