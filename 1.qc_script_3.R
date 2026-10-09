###Save fcs file

library(Spectre)

Spectre::package.check(type = 'spatial')

Spectre::package.load(type = 'spatial')

 

#Extact single cell data and save as .fcs using Spectre

 

#First check that the cell names of the assay and metadata are exactly the same

identical(rownames(t(assay(spe, "rescaled"))), rownames(colData(spe)))

 

#Extract data, localization  and metadata, and merge. Finally, save as .fcs

              data <- t(assay(spe, "rescaled"))

             data <- data.frame(data)

 

              metadata <- colData(spe)

              metadata <- data.frame(metadata)

 

              loc <- spatialCoords (spe)

              loc <- data.frame(loc)

 

 

              all <- cbind(data,metadata,loc)

 

              setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

              dir.create("data/fcs files 23.01.2025")

              write.files(all, file.prefix = 'data/fcs files 23.01.2025/_all_rescaled', divide.by = "time_group",write.csv = FALSE, write.fcs = TRUE)

              write.files(all, file.prefix = 'data/fcs files 23.01.2025/all_rescaled', write.csv = FALSE, write.fcs = TRUE)

 

 

 

#I use this for the z-score:

# Arcsinh transformation

assay(spe, "asinh") <- asinh(counts(spe)/1)

 

# Min-max rescaling

assay(spe, "rescaled") <- t(apply(assay(spe, "asinh"), 1, function(x){

  (x - min(x))/(max(x)-min(x))}))

 

# Z-score rescaling

assay(spe, "z-score") <- t(apply(assay(spe, "asinh"), 1, function(x){

  (x - mean(x))/(sd(x))}))

 

 

#And I also remove the cells that are too small or too big (here treshhold 5 – 250 µm, but check in your data):

colData(spe) %>%

    as.data.frame() %>%

    group_by(ROI) %>%

    ggplot() +

        geom_boxplot(aes(ROI, area)) +

        theme_minimal(base_size = 15) +

        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +

        ylab("Cell area") + xlab("")+

  ggtitle("before size cutoff")

 

# # size cut-off by percentile

# summary(spe$area)

# size.cutoff.lower <- quantile(spe$area, 0.01)

# cat(paste0("The 1st percentile of all object sizes is ", size.cutoff.lower, " µm^2. \n The objects smaller than this will be excluded.\n"))

# size.cutoff.upper <- quantile(spe$area, 0.99)

# cat(paste0("The 99th percentile of all object sizes is ", size.cutoff.upper, " µm^2. \n The objects larger than this will be excluded.\n"))

 

# manual size cut-off

size.cutoff.lower <- 5

size.cutoff.upper <- 250

cat(paste0("The lower and upper cut-off of object size are set to ",size.cutoff.lower, " and ",

           size.cutoff.upper, " µm^2 \nThe object not within this range will be excluded!"))

n.exclude <- sum(sum(spe$area < size.cutoff.lower), sum(spe$area > size.cutoff.upper))

 

spe$size.filter[between(spe$area, size.cutoff.lower, size.cutoff.upper)] <- "include"

spe$size.filter[!between(spe$area,size.cutoff.lower, size.cutoff.upper)] <- "exclude"

spe.no.filter <- spe

spe <- spe[,spe$size.filter == "include"]

 