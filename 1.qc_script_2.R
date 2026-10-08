library(imcRtools)
library(cytomapper)
library(EBImage)
library(abind)
library(grDevices)
library(dittoSeq)
library(viridis)
library(pheatmap)
library(cytoviewer)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(scuttle)
library(mclust)
library(scater)
library(patchwork)
library(ggrastr)
#################### spe ###############################################
#################### Image and cell evel quality QC ####################

images <- readRDS("/media/usb/Nagham_cytof_2nd_3rd/20082026_synovium/steinbock_r_objects/images.rds")
masks <- readRDS("/media/usb/Nagham_cytof_2nd_3rd/20082026_synovium/steinbock_r_objects/masks.rds")
spe <- readRDS("/media/usb/Nagham_cytof_2nd_3rd/20082026_synovium/steinbock_r_objects/spe.rds")

#########################################################################


##########################
### zscore on top of arcsinh norm ######
###########################


# Start from arcsinh-transformed expression values
exprs <- as.matrix(assay(spe, "exprs"))

# Z-score each marker across cells
marker_mean <- rowMeans(exprs)
marker_sd <- apply(exprs, 1, sd)

# Avoid division by zero for constant markers
marker_sd[marker_sd == 0] <- 1

assay(spe, "zscore") <- sweep(
  sweep(exprs, 1, marker_mean, "-"),
  1, marker_sd, "/"
)

###########verify#######################
z <- as.matrix(assay(spe, "zscore"))

range(z, na.rm = TRUE)
quantile(z, c(0, 0.01, 0.5, 0.99, 1), na.rm = TRUE)

# Across ALL cells used to calculate Z-scores:
# means should be approximately 0
range(rowMeans(z))

# SDs should be approximately 1 (0 for constant markers)
range(apply(z, 1, sd))

# Identify markers with the largest positive Z-scores
head(sort(apply(z, 1, max), decreasing = TRUE), 10)

####################################################


######################### antibody specificity ####################################
#assay(spe, "exprs") <- asinh(counts(spe)/1)
cur_cells <- sample(seq_len(ncol(spe)), 30000)

# 1. Open 300 DPI PNG graphics device
png(
  filename = "ditto_heatmap_highres_zscore_deepcell.png",
  width = 10,        # Width in inches
  height = 8,        # Height in inches
  units = "in",
  res = 300          # 300 DPI resolution
)

# 2. Render heatmap with global font size set to 14
dittoHeatmap(
  spe[, cur_cells], 
  genes = rownames(spe)[rowData(spe)$use_channel],
  assay = "zscore", 
  cluster_cols = TRUE, 
  scale = "none",
  heatmap.colors = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-4, 4, length.out = 101),
  legend_breaks = c(-4, -2, 0, 2, 4), 
  annot.by = "patient_id",
  annotation_colors = list(patient_id = metadata(spe)$color_vectors$patient_id),
  fontsize = 14,             # Sets overall font size to 14
  fontsize_row = 14,         # Ensures gene row labels are font 14
  fontsize_col = 14          # Ensures cell column labels are font 14
)

# 3. Close device to write file
dev.off()

######################### antibody specificity asinh ####################################
#assay(spe, "exprs") <- asinh(counts(spe)/1)
cur_cells <- sample(seq_len(ncol(spe)), 30000)

# 1. Open 300 DPI PNG graphics device
png(
  filename = "ditto_heatmap_highres_asinh_deepcell.png",
  width = 10,        # Width in inches
  height = 8,        # Height in inches
  units = "in",
  res = 300          # 300 DPI resolution
)

# 2. Render heatmap with global font size set to 14
dittoHeatmap(
  spe[, cur_cells], 
  genes = rownames(spe)[rowData(spe)$use_channel],
  assay = "exprs", 
  cluster_cols = TRUE, 
  scale = "none",
  heatmap.colors = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-4, 4, length.out = 101),
  legend_breaks = c(-4, -2, 0, 2, 4), 
  annot.by = "patient_id",
  annotation_colors = list(patient_id = metadata(spe)$color_vectors$patient_id),
  fontsize = 14,             # Sets overall font size to 14
  fontsize_row = 14,         # Ensures gene row labels are font 14
  fontsize_col = 14          # Ensures cell column labels are font 14
)

# 3. Close device to write file
dev.off()

###########################################################
###################################


# ------------------------------------------------------------------------------
# 1. Configuration & Parameters
# ------------------------------------------------------------------------------
# Output file path
pdf_output_path <- "multi_channel_roi_plots_deepcell.pdf"

# Target channels to visualize (Select 2 to 3 channels present in channelNames)
target_channels <- c("ICSK1", "DNA1", "ICSK2")

# Channel colors mapping
channel_colors <- list(
  ICSK1       = c("black", "red"),
  DNA1        = c("black", "blue"),
  ICSK2       = c("black", "green"),
  CellOutline = c("black", "white")
)

# ------------------------------------------------------------------------------
# 2. Rescale ALL Masks to Match Image Dimensions
# ------------------------------------------------------------------------------
masks_spe2_scaled <- CytoImageList(lapply(names(images), function(id) {
  img <- images[[id]]
  msk <- masks[[id]]
  
  EBImage::resize(
    msk, 
    w = dim(img)[1], 
    h = dim(img)[2], 
    filter = "none" # Nearest-neighbor interpolation preserves discrete cell IDs
  )
}))


# ------------------------------------------------------------------------------
# 4. Extract 1-Pixel Perimeters & Stack Outline Channel
# ------------------------------------------------------------------------------
old_channels <- channelNames(images)

images_spe2_with_outline <- CytoImageList(lapply(names(images), function(id) {
  img <- images[[id]]
  msk <- masks[[id]]
  
  m_mat <- as.matrix(msk)
  nr <- nrow(m_mat)
  nc <- ncol(m_mat)
  
  # 4-neighbor difference on integer cell IDs
  diff_r <- m_mat[-1, ] != m_mat[-nr, ]
  diff_c <- m_mat[, -1] != m_mat[, -nc]
  
  perim_mat <- matrix(FALSE, nrow = nr, ncol = nc)
  perim_mat[-nr, ] <- perim_mat[-nr, ] | diff_r
  perim_mat[-1, ]  <- perim_mat[-1, ]  | diff_r
  perim_mat[, -nc] <- perim_mat[, -nc] | diff_c
  perim_mat[, -1]  <- perim_mat[, -1]  | diff_c
  
  perim <- perim_mat & (m_mat > 0)
  
  # Stack low-intensity perimeter channel (0.03 for a delicate outline)
  img_mat <- as.array(img)
  outline_mat <- array(as.numeric(perim) * 0.03, dim = c(nr, nc, 1))
  
  combined_mat <- abind::abind(img_mat, outline_mat, along = 3)
  return(Image(combined_mat))
}))

names(images_spe2_with_outline) <- names(images)
mcols(images_spe2_with_outline) <- mcols(images)
channelNames(images_spe2_with_outline) <- c(old_channels, "CellOutline")

# ------------------------------------------------------------------------------
# 5. Export All ROIs to a Single PDF File
# ------------------------------------------------------------------------------
roi_names <- names(images_spe2_with_outline)

# Open PDF device (A4 paper size layout)
pdf(file = pdf_output_path, width = 8, height = 10)

for (roi in roi_names) {
  # Subset single ROI
  cur_img <- images_spe2_with_outline[roi]
  
  # Render plot page
  plotPixels(
    image = cur_img,
    missing_colour = "white",
    colour_by = c(target_channels, "CellOutline"),
    colour = channel_colors,
    bcg = list(
      CellOutline = c(0, 0.3, 1) # Soft outline gain control
    ),
    image_title = list(
      text = paste("ROI:", roi),
      cex = 1.0
    ),
    legend = list(
      colour_by.title.cex = 0.8,
      colour_by.labels.cex = 0.8
    )
  )
}

# Close graphic device
dev.off()

cat("Successfully exported", length(roi_names), "ROIs to", pdf_output_path, "\n")

################################################################################
#################### Image and cell evel quality QC ####################

cur_snr <- lapply(names(images), function(x){
    img <- images[[x]]
    mat <- apply(img, 3, function(ch){
        # Otsu threshold
        thres <- otsu(ch, range = c(min(ch), max(ch)), levels = 65536)
        # Signal-to-noise ratio
        snr <- mean(ch[ch > thres]) / mean(ch[ch <= thres])
        # Signal intensity
        ps <- mean(ch[ch > thres])
        
        return(c(snr = snr, ps = ps))
    })
    t(mat) %>% as.data.frame() %>% 
        mutate(image = x,
               marker = colnames(mat)) %>% 
        pivot_longer(cols = c(snr, ps))
})

cur_snr <- do.call(rbind, cur_snr)


# 1. Open 300 DPI PNG graphics device
png(
  filename = "snr_img_level_deepcell.png",
  width = 10,        # Width in inches
  height = 8,        # Height in inches
  units = "in",
  res = 300          # 300 DPI resolution
)


cur_snr %>% 
    group_by(marker, name) %>%
    summarize(log_mean = log2(mean(value))) %>%
    pivot_wider(names_from = name, values_from = log_mean) %>%
    ggplot() +
    geom_point(aes(ps, snr)) +
    geom_label_repel(aes(ps, snr, label = marker)) +
    theme_minimal(base_size = 15) + ylab("Signal-to-noise ratio [log2]") +
    xlab("Signal intensity [log2]")
dev.off()


cur_snr <- cur_snr %>% 
    pivot_wider(names_from = name, values_from = value) %>%
    filter(ps > 2) %>%
    pivot_longer(cols = c(snr, ps))


# 1. Open 300 DPI PNG graphics device
png(
  filename = "snr_img_level_deepcell_filtered.png",
  width = 10,        # Width in inches
  height = 8,        # Height in inches
  units = "in",
  res = 300          # 300 DPI resolution
)


cur_snr %>% 
    group_by(marker, name) %>%
    summarize(log_mean = log2(mean(value))) %>%
    pivot_wider(names_from = name, values_from = log_mean) %>%
    ggplot() +
    geom_point(aes(ps, snr)) +
    geom_label_repel(aes(ps, snr, label = marker)) +
    theme_minimal(base_size = 15) + ylab("Signal-to-noise ratio [log2]") +
    xlab("Signal intensity [log2]")

dev.off()


########## qc by image area ################

cell_density <- colData(spe) %>%
    as.data.frame() %>%
    group_by(sample_id) %>%
    # Compute the number of pixels covered by cells and 
    # the total number of pixels
    summarize(cell_area = sum(area),
              no_pixels = mean(width_px) * mean(height_px)) %>%
    # Divide the total number of pixels 
    # by the number of pixels covered by cells
    mutate(covered_area = cell_area / no_pixels)


png(
  filename = "cell_coverage_per_sample.png",
  width = 10,        # Width in inches
  height = 8,        # Height in inches
  units = "in",
  res = 300          # 300 DPI resolution
)

# Visualize the image area covered by cells per image
ggplot(cell_density) +
        geom_point(aes(reorder(sample_id,covered_area), covered_area)) + 
        theme_minimal(base_size = 15) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15)) +
        ylim(c(0, 1)) +
        ylab("% covered area") + xlab("")

dev.off()


############ heatmap cells aggregated asinh ###############################

image_mean <- aggregateAcrossCells(spe, 
                                   ids = spe$sample_id, 
                                   statistics="mean",
                                   use.assay.type = "exprs")
#assay(image_mean, "exprs") <- asinh(counts(image_mean))

png(
  filename = "heatmap_cells_aggregated.png",
  width = 10,        # Width in inches
  height = 12,        # Height in inches
  units = "in",
  res = 300          # 300 DPI resolution
)

dittoHeatmap(image_mean, genes = rownames(spe)[rowData(spe)$use_channel],
             assay = "exprs", cluster_cols = TRUE, scale = "none",
             heatmap.colors = viridis(100), 
             annot.by = c( "patient_id", "ROI"),
             annotation_colors = list(indication = metadata(spe)$color_vectors$indication,
                                      patient_id = metadata(spe)$color_vectors$patient_id,
                                      ROI = metadata(spe)$color_vectors$ROI),
             show_colnames = TRUE)

dev.off()


############ heatmap cells aggregated zscore ###############################

image_mean <- aggregateAcrossCells(spe, 
                                   ids = spe$sample_id, 
                                   statistics="mean",
                                   use.assay.type = "zscore")
#assay(image_mean, "exprs") <- asinh(counts(image_mean))

png(
  filename = "heatmap_cells_aggregated_zscore.png",
  width = 10,        # Width in inches
  height = 12,        # Height in inches
  units = "in",
  res = 300          # 300 DPI resolution
)

dittoHeatmap(image_mean, genes = rownames(spe)[rowData(spe)$use_channel],
             assay = "zscore", cluster_cols = TRUE, scale = "none",
             heatmap.colors = viridis(100), 
             annot.by = c( "patient_id", "ROI"),
             annotation_colors = list(indication = metadata(spe)$color_vectors$indication,
                                      patient_id = metadata(spe)$color_vectors$patient_id,
                                      ROI = metadata(spe)$color_vectors$ROI),
             show_colnames = TRUE)

dev.off()



##################################################################
################ cell level QC ########################################
###################################################################


set.seed(220224)
mat <- sapply(seq_len(nrow(spe)), function(x){
    cur_exprs <- assay(spe, "exprs")[x,]
    cur_counts <- assay(spe, "counts")[x,]
    
    cur_model <- Mclust(cur_exprs, G = 2)
    mean1 <- mean(cur_counts[cur_model$classification == 1])
    mean2 <- mean(cur_counts[cur_model$classification == 2])
    
    signal <- ifelse(mean1 > mean2, mean1, mean2)
    noise <- ifelse(mean1 > mean2, mean2, mean1)
    
    return(c(snr = signal/noise, ps = signal))
})
    
cur_snr <- t(mat) %>% as.data.frame() %>% 
        mutate(marker = rownames(spe))


png(
  filename = "snr_cell_deepcell.png",
  width = 10,        # Width in inches
  height = 8,        # Height in inches
  units = "in",
  res = 300          # 300 DPI resolution
)
cur_snr %>% ggplot() +
    geom_point(aes(log2(ps), log2(snr))) +
    geom_label_repel(aes(log2(ps), log2(snr), label = marker)) +
    theme_minimal(base_size = 15) + ylab("Signal-to-noise ratio [log2]") +
    xlab("Signal intensity [log2]")

dev.off()



### Next, we observe the distributions of cell size across the individual images.

png(
  filename = "cell_size_dist_deepcell.png",
  width = 20,        # Width in inches
  height = 8,        # Height in inches
  units = "in",
  res = 300          # 300 DPI resolution
)

dittoPlot(spe, var = "area", 
          group.by = "sample_id", 
          plots = "boxplot") +
        ylab("Cell area") + xlab("")
dev.off()


summary(spe$area)

sum(spe$area < 5)

spe <- spe[,spe$area >= 5]


######### cell density per mm2

cell_density <- colData(spe) %>%
    as.data.frame() %>%
    group_by(sample_id) %>%
    summarize(cell_count = n(),
           no_pixels = mean(width_px) * mean(height_px)) %>%
    mutate(cells_per_mm2 = cell_count/(no_pixels/1000000))


png(
  filename = "cell_density_per_mm2_deepcell.png",
  width = 12,        # Width in inches
  height = 8,        # Height in inches
  units = "in",
  res = 300          # 300 DPI resolution
)

ggplot(cell_density) +
    geom_point(aes(reorder(sample_id,cells_per_mm2), cells_per_mm2)) + 
    theme_minimal(base_size = 15) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
    ylab("Cells per mm2") + xlab("")

dev.off()


######### technical/batch effects QC per channel

png(
  filename = "batch_signal_channel_deepcell.png",
  width = 18,        # Width in inches
  height = 18,        # Height in inches
  units = "in",
  res = 300          # 300 DPI resolution
)

multi_dittoPlot(spe, vars = rownames(spe)[rowData(spe)$use_channel],
               group.by = "patient_id", plots = "ridgeplot", 
               assay = "exprs", 
               color.panel = metadata(spe)$color_vectors$patient_id)

dev.off()



############################################################################
######################## Dim. Red. TSNE UMAP (arcsh) ###############################
############################################################################

set.seed(220225)
spe <- runUMAP(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs") 
spe <- runTSNE(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs") 

png(
  filename = "tsne_umap_roi_patient_arsh.png",
  width = 18,        # Width in inches
  height = 18,        # Height in inches
  units = "in",
  res = 300          # 300 DPI resolution
)


# visualize patient id 
p1 <- dittoDimPlot(spe, var = "patient_id", reduction.use = "UMAP", size = 0.2) + 
    scale_color_manual(values = metadata(spe)$color_vectors$patient_id) +
    ggtitle("Patient ID on UMAP")
p2 <- dittoDimPlot(spe, var = "patient_id", reduction.use = "TSNE", size = 0.2) + 
    scale_color_manual(values = metadata(spe)$color_vectors$patient_id) +
    ggtitle("Patient ID on TSNE")

# visualize region of interest id
p3 <- dittoDimPlot(spe, var = "ROI", reduction.use = "UMAP", size = 0.2) + 
    scale_color_manual(values = metadata(spe)$color_vectors$ROI) +
    ggtitle("ROI ID on UMAP")
p4 <- dittoDimPlot(spe, var = "ROI", reduction.use = "TSNE", size = 0.2) + 
    scale_color_manual(values = metadata(spe)$color_vectors$ROI) +
    ggtitle("ROI ID on TSNE")



(p1 + p2) / (p3 + p4) 

dev.off()


   ########## plot marker expression ###################
######################## UMAP expression for all markers ########################
# Include every marker in spe, regardless of use_channel; reuse the existing UMAP.
umap_markers <- rownames(spe)
stopifnot(length(umap_markers) > 0L)
umap_ncol <- min(4L, length(umap_markers))
umap_nrow <- ceiling(length(umap_markers) / umap_ncol)

umap_marker_plots <- lapply(umap_markers, function(marker) {
    ggrastr::rasterise(dittoDimPlot(
        spe, var = marker, reduction.use = "UMAP", assay = "exprs",
        size = 0.2, raster = FALSE
    ), layers = "Point", dpi = 150) +
        ggtitle(marker) +
        theme_minimal(base_size = 11) +
        theme(
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
            plot.margin = margin(14, 18, 14, 18, unit = "pt"),
            panel.grid = element_blank(),
            legend.position = "right",
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 8),
            legend.key.height = grid::unit(1.2, "cm")
        )
})

umap_all_markers <- patchwork::wrap_plots(
    umap_marker_plots, ncol = umap_ncol, guides = "keep"
)

# A single page with 5 x 4 inch panels and separate expression legends.
# Vector points and text remain sharp at any zoom; PDF does not need a raster DPI.
ggsave(
    filename = "umap_all_markers_deepcell.pdf",
    plot = umap_all_markers,
    device = grDevices::pdf,
    width = 5 * umap_ncol, height = 4 * umap_nrow, units = "in",
    limitsize = FALSE, bg = "white", useDingbats = FALSE
)




############################################################################
######################## Dim. Red. TSNE UMAP (zscore) ###############################
############################################################################

set.seed(220225)
spe <- runUMAP(spe, subset_row = rowData(spe)$use_channel, exprs_values = "zscore") 
spe <- runTSNE(spe, subset_row = rowData(spe)$use_channel, exprs_values = "zscore") 

png(
  filename = "tsne_umap_roi_patient_zscore.png",
  width = 18,        # Width in inches
  height = 18,        # Height in inches
  units = "in",
  res = 300          # 300 DPI resolution
)


# visualize patient id 
p1 <- dittoDimPlot(spe, var = "patient_id", reduction.use = "UMAP", size = 0.2) + 
    scale_color_manual(values = metadata(spe)$color_vectors$patient_id) +
    ggtitle("Patient ID on UMAP")
p2 <- dittoDimPlot(spe, var = "patient_id", reduction.use = "TSNE", size = 0.2) + 
    scale_color_manual(values = metadata(spe)$color_vectors$patient_id) +
    ggtitle("Patient ID on TSNE")

# visualize region of interest id
p3 <- dittoDimPlot(spe, var = "ROI", reduction.use = "UMAP", size = 0.2) + 
    scale_color_manual(values = metadata(spe)$color_vectors$ROI) +
    ggtitle("ROI ID on UMAP")
p4 <- dittoDimPlot(spe, var = "ROI", reduction.use = "TSNE", size = 0.2) + 
    scale_color_manual(values = metadata(spe)$color_vectors$ROI) +
    ggtitle("ROI ID on TSNE")



(p1 + p2) / (p3 + p4) 

dev.off()


   ########## plot marker expression ###################
######################## UMAP expression for all markers ########################
# Include every marker in spe, regardless of use_channel; reuse the existing UMAP.
umap_markers <- rownames(spe)
stopifnot(length(umap_markers) > 0L)
umap_ncol <- min(4L, length(umap_markers))
umap_nrow <- ceiling(length(umap_markers) / umap_ncol)

umap_marker_plots <- lapply(umap_markers, function(marker) {
    ggrastr::rasterise(dittoDimPlot(
        spe, var = marker, reduction.use = "UMAP", assay = "zscore",
        size = 0.2, raster = FALSE
    ), layers = "Point", dpi = 150) +
        ggtitle(marker) +
        theme_minimal(base_size = 11) +
        theme(
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
            plot.margin = margin(14, 18, 14, 18, unit = "pt"),
            panel.grid = element_blank(),
            legend.position = "right",
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 8),
            legend.key.height = grid::unit(1.2, "cm")
        )
})

umap_all_markers <- patchwork::wrap_plots(
    umap_marker_plots, ncol = umap_ncol, guides = "keep"
)

# A single page with 5 x 4 inch panels and separate expression legends.
# Vector points and text remain sharp at any zoom; PDF does not need a raster DPI.
ggsave(
    filename = "umap_all_markers_deepcell_zscore.pdf",
    plot = umap_all_markers,
    device = grDevices::pdf,
    width = 5 * umap_ncol, height = 4 * umap_nrow, units = "in",
    limitsize = FALSE, bg = "white", useDingbats = FALSE
)
