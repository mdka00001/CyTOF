#!/usr/bin/env Rscript
options(repos=c(CRAN='https://cloud.r-project.org'))
cran <- c('mclust','uwot','Rtsne','tiff','BiocManager','ggplot2','patchwork','ggrastr','viridis','RColorBrewer')
missing <- cran[!vapply(cran,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) install.packages(missing)
bioc <- c('EBImage','DESeq2','dittoSeq','scater','cytomapper','flowCore')
missing <- bioc[!vapply(bioc,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) BiocManager::install(missing,ask=FALSE,update=FALSE)
