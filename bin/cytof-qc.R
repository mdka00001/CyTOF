#!/usr/bin/env Rscript
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
# R encodes spaces in --file arguments as ~+~. Decode before resolving modules.
script <- gsub('~+~', ' ', script, fixed=TRUE)
root <- dirname(dirname(normalizePath(script, mustWork=TRUE)))
for (f in c('config','io','transforms','qc','report','images','pipeline')) source(file.path(root,'R',paste0(f,'.R')))
args <- commandArgs(TRUE)
if (!length(args) || '--help' %in% args || '-h' %in% args) { help_text(); quit(status=0) }
tryCatch(run_qc(parse_args(args)),error=function(e) { message('ERROR: ',conditionMessage(e)); quit(status=1) })
