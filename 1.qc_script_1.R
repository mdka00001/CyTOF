library(imcRtools)
library(cytomapper)
library(dplyr)
library(dittoSeq)
library(RColorBrewer)
library(stringr)
############################## Read steinbock output #################################
spe <- read_steinbock("/media/usb/Nagham_cytof_2nd_3rd/20082026_synovium")
spe

#unique(colnames(spe))
### check the counts matrix ###
counts(spe)[1:5,1:5]


### check the metadata ###
head(colData(spe))


### check spatial coords ###
head(spatialCoords(spe))


## cell cell interaction neighbours ###
colPair(spe, "neighborhood")


### channel info ###
head(rowData(spe))



############################## Fix Metadata ######################
 ###### 1. add sample metadata ########
 ###### steinbock output #########
mapping_samples <- c("Synovium_B1223485_B1235864_B1236739_B1238297_20082026_001" = "B12_35864",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_002" = "B12_35864",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_003" = "B12_35864",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_004" = "B12_36739",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_005" = "B12_36739",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_006" = "B12_36739",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_007" = "B12_36739",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_008" = "B12_36739",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_009" = "B12_23485",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_010" = "B12_23485",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_011" = "B12_23485",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_012" = "B12_23485",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_013" = "B12_23485",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_014" = "B12_23485",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_015" = "B12_23485",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_016" = "B12_23485",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_017" = "B12_23485",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_018" = "B12_23485",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_019" = "B12_23485",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_020" = "B12_23485",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_021" = "B12_23485",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_022" = "B12_38297",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_023" = "B12_38297",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_024" = "B12_38297",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_025" = "B12_38297",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_026" = "B12_38297",
                     "Synovium_B1223485_B1235864_B1236739_B1238297_20082026_027" = "B12_38297"

)

spe$patient_id <- unname(mapping_samples[as.character(spe$sample_id)])

###################################################################
###################################################################


### prepare metadata ###

extract_roi <- function(sample_ids) {
  # Split each sample_id by standard delimiters (_, -, or .)
  tokens_list <- strsplit(sample_ids, "[_\\.\\-]")
  
  # Determine token lengths and align matrix representation
  max_len <- max(sapply(tokens_list, length))
  padded_tokens <- t(sapply(tokens_list, function(x) {
    c(x, rep(NA, max_len - length(x)))
  }))
  
  # Identify fields based on variance across samples
  is_numeric <- apply(padded_tokens, 2, function(col) {
    all(grepl("^\\d+$", na.omit(col)))
  })
  
  unique_counts <- apply(padded_tokens, 2, function(col) length(unique(na.omit(col))))
  
  # Strategy: Find the RIGHTMOST numeric field with HIGH variance
  # This captures the finest-resolution ROI identifier
  numeric_high_var <- which(is_numeric & unique_counts > 1)
  
  if (length(numeric_high_var) > 0) {
    # Use the rightmost high-variance numeric field
    roi_idx <- max(numeric_high_var)
  } else {
    # No high-variance numeric field found
    # Try to find any numeric field; if none, use the last field
    numeric_fields <- which(is_numeric)
    if (length(numeric_fields) > 0) {
      roi_idx <- max(numeric_fields)
    } else {
      # No numeric fields at all - use the last field
      roi_idx <- max_len
    }
  }
  
  # Extract and return only ROI
  roi_vector <- padded_tokens[, roi_idx]
  return(roi_vector)
}

spe$ROI <- extract_roi(as.character(spe$sample_id))
unique(spe$ROI)

unique(spe$sample_id)
head(colData(spe))

###################################################################
###################################################################

############################ arcsinh transformation ###############

assay(spe, "exprs") <- asinh(counts(spe)/1)


###################################################################
###################################################################


################# define color scheme #################

### spe ###
color_vectors <- list()

ROI <- setNames(colorRampPalette(brewer.pal(11, "BrBG"))(length(unique(spe$ROI))), 
                unique(spe$ROI))
patient_id <- setNames(brewer.pal(length(unique(spe$patient_id)), name = "Set1"), 
                unique(spe$patient_id))
sample_id <- setNames(colorRampPalette(c(brewer.pal(6, "YlOrRd")[3:5],
                                           brewer.pal(6, "PuBu")[3:6],
                                           brewer.pal(6, "YlGn")[3:5],
                                           brewer.pal(6, "BuPu")[3:6]))(length(unique(spe$sample_id))),
                unique(spe$sample_id))
indication <- setNames(brewer.pal(length(unique(spe$indication)), name = "Set2"), 
                unique(spe$indication))

color_vectors$ROI <- ROI
color_vectors$patient_id <- patient_id
color_vectors$sample_id <- sample_id
color_vectors$indication <- indication

metadata(spe)$color_vectors <- color_vectors

#####################################################################
#####################################################################

######################### Read Images ###############################


images <- loadImages("/media/usb/Nagham_cytof_2nd_3rd/20082026_synovium/img")

masks <- loadImages("/media/usb/Nagham_cytof_2nd_3rd/20082026_synovium/masks_deepcell/", as.is = TRUE)

channelNames(images) <- rownames(spe)
images
masks

all.equal(names(images), names(masks))


# Extract patient id from image name using existing mapping
patient_id <- unname(mapping_samples[names(images)])

patient_id

# Store patient and image level information in elementMetadata
mcols(images) <- mcols(masks) <- DataFrame(sample_id = names(images),
                                           patient_id = patient_id)

### spe ###

##########################################################################

#### save files ###
saveRDS(spe, "/media/usb/Nagham_cytof_2nd_3rd/20082026_synovium/steinbock_r_objects/spe.rds")
saveRDS(images, "/media/usb/Nagham_cytof_2nd_3rd/20082026_synovium/steinbock_r_objects/images.rds")
saveRDS(masks, "/media/usb/Nagham_cytof_2nd_3rd/20082026_synovium/steinbock_r_objects/masks.rds")


#####################################################################################################




