## Load the regions bed files and perform the enrichment analysis vs remap2022

library(dplyr)
library(ggplot2)
library(regioneReloaded)
# library(BSgenome.Hsapiens.UCSC.hg38)
output_directory = "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/"
input_rds_path = "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/sampling_1000_ngn2_with_background_full_remap_TF_chip_increasing_decreasing_CrosswiseMatrix.rds"
undiffWTC11_input_rds_path = "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/sampling_1000_undiffWTC11_with_background_full_remap_TF_chip_increasing_decreasing_CrosswiseMatrix.rds"

ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix <- readRDS(input_rds_path)
undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix <- readRDS(undiffWTC11_input_rds_path)


# plot a subset of the data in the matrix form (working)

# Get the matrix of Z-scores (or transformed values)
matrix_data <- regioneReloaded::getMatrix(ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix)
matrix_data <- regioneReloaded::getMatrix(undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix)

dim(matrix_data)

colnames(matrix_data) # These should be your TF names
rownames(matrix_data) # These should be your Alist names

# Get the enrichment values for sequence set 1
enrichment_set1 <- matrix_data["Neuron_upregulating_dCREs", ]
head(enrichment_set1)
# Get the enrichment values for sequence set 2
enrichment_set2 <- matrix_data["Neuron_downregulating_dCREs", ]
head(enrichment_set2)

# Order TFs by enrichment for set 1 (descending)
sorted_tfs_set1 <- names(sort(enrichment_set1, decreasing = TRUE))

# Order TFs by enrichment for set 2 (descending)
sorted_tfs_set2 <- names(sort(enrichment_set2, decreasing = TRUE))

# Select the top 15 enriched TFs for each set
top15_set1_tfs <- head(sorted_tfs_set1, 15)
top15_set2_tfs <- head(sorted_tfs_set2, 15)


# Combine and unique the TF names to avoid duplicates if a TF is in both top 15 lists
selected_tfs <- unique(c(top15_set1_tfs, top15_set2_tfs))

# You might want to ensure that these selected_tfs are indeed in your matrix_data colnames
selected_tfs <- selected_tfs[selected_tfs %in% colnames(matrix_data)]

# Subset the matrix data by selected TFs (columns)
subset_matrix_data <- matrix_data[, selected_tfs]

# Check dimensions of the new matrix
dim(subset_matrix_data)

# Make a copy of the original genoMatriXeR object
ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix_subset <- ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix
undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix_subset <- undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix

# Replace the 'matrix' slot with the subsetted matrix data
# The @matrix slot itself is a list, and the actual matrix is usually under 'GMat' or 'LZM'
# Let's inspect the structure of the @matrix slot first
# str(ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix@matrix)

# Assuming the matrix is stored in the 'GMat' element of the @matrix slot
ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix_subset@matrix$GMat <- subset_matrix_data
undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix_subset@matrix$GMat <- subset_matrix_data

# If you also had clustering information stored in @matrix (e.g., FitCol for column clustering),
# you'll need to either remove it (to force plotCrosswiseMatrix to re-cluster if you set clusterize_cols=TRUE)
# or update it if you manually re-cluster. For simplicity, let's remove it and let plotCrosswiseMatrix handle it.

# Remove existing column clustering information if it exists, to allow new clustering on subset
# This depends on how the clustering information is stored. It's often in FitCol/FitRow.
if (!is.null(ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix_subset@matrix$FitCol)) {
    ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix_subset@matrix$FitCol <- NULL
}
if (!is.null(undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix_subset@matrix$FitCol)) {
    undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix_subset@matrix$FitCol <- NULL
}

subset_plot <- plotCrosswiseMatrix(ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix_subset)
ggsave("/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/ngn2_subset_TFBS_crosswise_matrix_plot.pdf", plot = subset_plot, width = 10, height = 8)

# subset_plot <- plotCrosswiseMatrix(undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix_subset)
# ggsave("/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/undiffWTC11_subset_TFBS_crosswise_matrix_plot.pdf", plot = subset_plot, width = 10, height = 8)


# Extract the two datasets
up_data <- ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix@multiOverlaps$Neuron_upregulating_dCREs %>%
    mutate(regulation = "increasing transcription") # Add a column to label the dataset
down_data <- ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix@multiOverlaps$Neuron_downregulating_dCREs %>%
    mutate(regulation = "decreasing transcription") # Add a column to label the dataset
# up_data <- undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix@multiOverlaps$Neuron_upregulating_dCREs %>%
#     mutate(regulation = "increasing transcription") # Add a column to label the dataset
# down_data <- undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix@multiOverlaps$Neuron_downregulating_dCREs %>%
#     mutate(regulation = "decreasing transcription") # Add a column to label the dataset

head(up_data)

# Select top x other enriched experiments (not in specific_experiments) for each dataset
n_top = 15
sig_level = 0.1
top5_up <- up_data %>%
    filter(adj.p_value < sig_level & is.finite(norm_zscore)) %>%
    arrange(desc(norm_zscore)) %>%
    slice_head(n = n_top) %>%
    pull(name)

top5_down <- down_data %>%
    filter(adj.p_value < sig_level & is.finite(norm_zscore)) %>%
    arrange(desc(norm_zscore)) %>%
    slice_head(n = n_top) %>%
    pull(name)

head(top5_up)
head(top5_down)

# Filter for specific experiments
filtered_up <- up_data %>%
    filter(name %in% c(top5_up))
filtered_down <- down_data %>%
    filter(name %in% c(top5_down))

# Combine all data
combined_data <- bind_rows(filtered_up, filtered_down)
head(combined_data)
x_limits <- range(c(filtered_up$norm_zscore, filtered_down$norm_zscore), na.rm = TRUE)


# Plot the data
pdf_name_enhancer <- "ngn2_dotplot_enhancer_remap2022_TF_enrichment_sampling.pdf"
pdf_name_enhancer <- "undiffWTC11_dotplot_enhancer_remap2022_TF_enrichment_sampling.pdf"
pdf(file.path(output_directory, pdf_name_enhancer), width = 8, height = 6) # Adjust width and height
ggplot(filtered_up, aes(x = norm_zscore, y = reorder(name, norm_zscore), fill = regulation)) +
    geom_point(aes(size = -log10(adj.p_value)), shape = 21, alpha = 0.6, stroke = NA) + # Ensure fill is mapped and color is NA
    labs(
        x = "Enrichment (normalized z-score)",
        y = "Experiment",
        title = "Enrichment of CREs"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        legend.position = "top",
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 10, margin = margin(r = 5)), # Reduce y-axis text size for density
        axis.text.x = element_text(size = 12) # Adjust x-axis text size
    ) +
    scale_fill_manual(
        values = c("Upregulating" = "#33a02c"), # , "Downregulating" = "#e31a1c"),  # Green for Enhancer, Red for Silencing
        labels = c("Upregulating" = "Enhancer"), # , "Downregulating" = "Silencer"),  # Rename labels
        name = "Regulation"
    ) +
    scale_size_continuous(
        range = c(2, 5), # Adjust dot size range for better visibility
        name = "Significance (-log10 p-value)"
    ) +
    guides(
        fill = guide_legend(override.aes = list(size = 6, shape = 21)) # Increase legend dot size
    ) +
    scale_x_continuous(limits = x_limits) +
    coord_cartesian(clip = "off") # Ensure no points are cut off
dev.off()

# Plot the data
pdf_name_silencer <- "ngn2_dotplot_silencer_remap2022_TF_enrichment_sampling.pdf"
pdf_name_silencer <- "undiffWTC11_dotplot_silencer_remap2022_TF_enrichment_sampling.pdf"
pdf(file.path(output_directory, pdf_name_silencer), width = 8, height = 6) # Adjust width and height
ggplot(filtered_down, aes(x = norm_zscore, y = reorder(name, norm_zscore), fill = regulation)) +
    geom_point(aes(size = -log10(adj.p_value)), shape = 21, alpha = 0.6, stroke = NA) + # Ensure fill is mapped and color is NA
    labs(
        x = "Enrichment (normalized z-score)",
        y = "Experiment",
        title = "Enrichment of CREs"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        legend.position = "top",
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 10, margin = margin(r = 5)), # Reduce y-axis text size for density
        axis.text.x = element_text(size = 12) # Adjust x-axis text size
    ) +
    scale_fill_manual(
        values = c("Downregulating" = "#e31a1c"), # Green for Enhancer, Red for Silencing
        labels = c("Downregulating" = "Silencer"), # Rename labels
        name = "Regulation"
    ) +
    scale_size_continuous(
        range = c(2, 5), # Adjust dot size range for better visibility
        name = "Significance (-log10 p-value)"
    ) +
    guides(
        fill = guide_legend(override.aes = list(size = 6, shape = 21)) # Increase legend dot size
    ) +
    scale_x_continuous(limits = x_limits) +
    coord_cartesian(clip = "off") # Ensure no points are cut off
dev.off()
