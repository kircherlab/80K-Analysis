## Load the regions bed files and perform the enrichment analysis vs remap2022
library(dplyr)
library(ggplot2)
library(regioneReloaded)
library(BSgenome.Hsapiens.UCSC.hg38)

### NGN2 focus: undiffWTC11 dCREs are not filtered out + very small subset as background set


input_path <- "/home/kisa/coding/80K_MPRA/remap_chipseq_pia" # no tailing slash allowed
input_path <- "/sc-projects/sc-proj-bih-reg-seqs/users/kisa11/projects/80k_region_enrichment_regioneReloaded" # no tailing slash allowed
output_directory <- "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results" # no trailing slash allowed
output_directory <- "/sc-projects/sc-proj-bih-reg-seqs/users/kisa11/projects/80k_region_enrichment_regioneReloaded/80k_results" # no trailing slash allowed
n_threads <- 64
# remap2022:
# zcat remap2022_nr_macs2_hg38_v1_0.bed.gz | head -n 100000 | gzip >remap2022_nr_macs2_hg38_v1_0_subset_100k.bed.gz
# remap2022 <- file.path(input_path, "remap2022_nr_macs2_hg38_v1_0_subset_100k.bed.gz")
# remap2022 <- file.path(input_path, "remap2022_nr_macs2_hg38_v1_0_subset_10k.bed.gz")
# remap2022 <- file.path(input_path, "remap2022_nr_macs2_hg38_v1_0.bed.gz")
# 68.655.741 regions
remap2022 <- file.path(input_path, "remap2022_nr_macs2_hg38_v1_0_target_column.bed.gz")

# background sequences: (why 0.01? Stringent definition of around zero activity; Why 0.04 a bigger set of sequences without statistical significant activity different to scrambled sequences. looked for a threshold which uses the most regions but still has no overlap with the considered active regions)
# around_zero_activity_both_cell_types_001 <- file.path(input_path, "around_zero_activity_both_cell_types_2_duplicated_regions_not_stand_aware.bed")
around_zero_activity_both_cell_types_001 <- file.path(input_path, "around_zero_activity_both_cell_types_no_duplicated_regions_001_1000_subsampled.bed")

# NGN2 focus
upbed <- file.path(input_path, "ngn2_increasing_regions_duplicates_without_considering_strand.bed")
downbed <- file.path(input_path, "ngn2_decreasing_regions_duplicates_without_considering_strand.bed")

data_up <- read.table(upbed, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
data_down <- read.table(downbed, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))

head(data_up)

ngn2_increasing_activity_regions <- GenomicRanges::makeGRangesFromDataFrame(data_up, seqnames.field = "chr", start.field = "start", end.field = "end")
ngn2_decreasing_activity_regions <- GenomicRanges::makeGRangesFromDataFrame(data_down, seqnames.field = "chr", start.field = "start", end.field = "end")
# ngn2_dCREs <- list(Neuron_upregulating_dCREs = ngn2_increasing_activity_regions, Neuron_downregulating_dCREs = ngn2_decreasing_activity_regions)

# # background sequences
background_conservative <- read.table(around_zero_activity_both_cell_types_001, header = F, sep = "\t", col.names = c("chr", "start", "end", "name", "strand", "score"))
both_cell_types_around_zero_activity_001 <- GenomicRanges::makeGRangesFromDataFrame(background_conservative, seqnames.field = "chr", start.field = "start", end.field = "end")
ngn2_dCREs_and_background <- list(Neuron_upregulating_dCREs = ngn2_increasing_activity_regions, Neuron_downregulating_dCREs = ngn2_decreasing_activity_regions, background_conservative = both_cell_types_around_zero_activity_001)

# max: use library(data.table) and fread() for faster reading of large files
remap_df <- read.table(remap2022)
TFs <- GenomicRanges::makeGRangesFromDataFrame(remap_df, seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = T)

split_granges <- split(TFs, GenomicRanges::mcols(TFs)$V10)
split_granges <- lapply(split_granges, function(gr) {
    GenomicRanges::mcols(gr)$V4 <- NULL
    return(gr)
})
tflist <- c(split_granges)

# wrtie tflist

# now we performe the crosswise permutation test, using 5000 sampling of the data and 100 permutations.
# The function that we use for randomization if "randomizeRegions" (NOTE: It can take a while, but you
# # can do it in parallel using the argument mc.cores) and the evalutation function is "numOveralps".

# print to the user: running crosswisePermTest
cat("Running crosswisePermTest with sampling...\n")

set.seed(42)
# WITH background sequences
ngn2_remap_chip_increasing_decreasing <- crosswisePermTest(
    Alist = ngn2_dCREs_and_background,
    Blist = tflist,
    sampling = TRUE,
    fraction = 0.15, # fraction of the regions to sample
    min_sampling = 500, # if the fraction of the regions is smaller than the minimum sampling
    ranFUN = "resampleRegions", # "randomizeRegions", # we might want to use resampleRegions() as well and compare
    evFUN = "numOverlaps",
    ntimes = 1000, # might want to increase to 5000
    genome = "hg38",
    mc.cores = n_threads
)
cat("Finished crosswisePermTest with sampling.\n")

# # WITHOUT background sequences
# ngn2_remap_chip_increasing_decreasing <- crosswisePermTest(
#     Alist = ngn2_dCREs,
#     Blist = tflist,
#     sampling = FALSE,
#     ranFUN = "randomizeRegions",
#     evFUN = "numOverlaps",
#     ntimes = 1000,
#     genome = "hg38",
#     mc.cores = n_threads,
# )

# # WITH background sequences
ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix <- makeCrosswiseMatrix(ngn2_remap_chip_increasing_decreasing, clusterize = TRUE, transform =TRUE) # clusterize = F, transform = T) # , transform=T) #, clusterize=T)

# # WITHOUT background sequences
# # set clusterize to FALSE if only 2 sequence sets => here more then 2 sequence sets
# ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix <- makeCrosswiseMatrix(ngn2_remap_chip_increasing_decreasing, clusterize = TRUE, transform = TRUE) # clusterize = F, transform = T) # , transform=T) #, clusterize=T)

output_name <- "sampling_1000_ngn2_with_background_full_remap_TF_chip_increasing_decreasing_CrosswiseMatrix.rds"

# Use file.path() to combine directory and file name
output_path <- file.path(output_directory, output_name)
saveRDS(ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix, file = output_path)

# Extract the two datasets
up_data <- ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix@multiOverlaps$Neuron_upregulating_dCREs %>%
    mutate(regulation = "increasing transcription") # Add a column to label the dataset
down_data <- ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix@multiOverlaps$Neuron_downregulating_dCREs %>%
    mutate(regulation = "decreasing transcription") # Add a column to label the dataset

# Select top x other enriched experiments (not in specific_experiments) for each dataset
n_top = 15
sig_level = 0.05
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

# Filter for specific experiments
filtered_up <- up_data %>%
    filter(name %in% c(top5_up))
filtered_down <- down_data %>%
    filter(name %in% c(top5_down))

# Combine all data
combined_data <- bind_rows(filtered_up, filtered_down)

x_limits <- range(c(filtered_up$norm_zscore, filtered_down$norm_zscore), na.rm = TRUE)


# Plot the data
pdf_name_enhancer <- "sampling_1000_dotplot_enhancer_with_background_full_remap2022_TF_enrichment.pdf"
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
pdf_name_silencer <- "sampling_1000_dotplot_silencer_with_background_full_remap2022_TF_enrichment.pdf"
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