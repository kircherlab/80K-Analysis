# install.packages("dplyr")
# install.packages("devtools")
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager", repos = 'https://cran.uni-muenster.de/')
# BiocManager::install("mpra")
# install.packages("dplyr")
# install.packages("devtools")
# install.packages("tidyr")
# # install ggplot
# install.packages("ggplot2")
# library(mpra)
library(ggplot2)
library(dplyr)
library(devtools)
library(tidyr)
library(arrow)

load_all("/home/kisa/coding/80K_MPRA/bc_mpralm")

# set output
output_dir <- "/home/kisa/coding/80K_MPRA/bc_MPRAlm_results/results/"
output_name <- "80k_cell_type_differences_BCalm_2025_05_normalized_counts.tsv"

# WTC11: Normalized counts using mpralib
input <- "/home/kisa/coding/80K_MPRA/element_analysis_output/202502_element_analysis/80k_combined_celltype_element_data_2025_05.tsv.gz"
input <- "/home/kisa/coding/80K_MPRA/element_analysis_output/202502_element_analysis/80k_combined_celltype_element_data_2025_05.tsv.gz"

# start preprocessing
df <- read.table(file = gzfile(input), sep = "\t", header = TRUE, comment.char = "")
colnames(df) <- c(
    "barcode", "name", "dna_count_1", "rna_count_1", "dna_count_2",
    "rna_count_2", "dna_count_3", "rna_count_3", "cell_type",
    "name_without_cell_type"
)

# filtering normalized counts:
df <- df %>% filter(if_all(matches("dna", ignore.case = TRUE), ~ . > 0))
df <- df %>% filter(if_all(matches("rna", ignore.case = TRUE), ~ . > 0))

print("Number of unique names loaded after dna filtering: ")
print(length(unique(df$name)))


# variant map
element_map = "/home/kisa/coding/80K_MPRA/element_analysis_output/202502_element_analysis/80k_combined_cell_types_matching_names.tsv.gz"
map_df <- read.table(file = gzfile(element_map), sep = "\t", header = TRUE, comment.char = "")
map_df <- map_df %>% select(ID, REF, ALT)
head(map_df)
head(df)
print("Create var df")

# var_df <- create_var_df(df, map_df)

# Merge on REF
df_ref <- merge(df, map_df, by.x = "name", by.y = "REF", all.x = FALSE)
head(df_ref)
df_ref$cell_type <- "NGN2_neurons"
df_ref$ALT <- NULL
print("Number of different oligos after left join with variant map for NGN2 Neurons")
print(length(unique(df_ref$name)))
# Merge on ALT
df_alt <- merge(df, map_df, by.x = "name", by.y = "ALT", all.x = FALSE)
df_alt$cell_type <- "undiff_WTC11"
df_alt$REF <- NULL
print("Number of different oligos after left join with variant map for undiff WTC11")
print(length(unique(df_alt$name)))
# Combine the results
df_combined <- rbind(df_ref, df_alt)

# Select and rename columns as necessary
var_df <- df_combined %>% select(variant_id = ID, cell_type, barcode, matches("count"))

# Note: ref is NGN2_neuron and alt is undiff_WTC11
print("Number of unique names loaded after create_var_df: ")
print(length(unique(var_df$variant_id)))

# Filter for oligos with at least min_bc barcodes in both alleles
min_bc <- 10
df_filt <- var_df %>%
    group_by(variant_id, cell_type) %>%
    filter(n() >= min_bc) %>%
    ungroup()
df_filt <- df_filt %>%
    group_by(variant_id) %>%
    filter(n_distinct(cell_type) == 2) %>%
    ungroup()

print("Number of unique variants loaded after filtering for barcodes and allele existence: ")
print(length(unique(df_filt$variant_id)))
head(df_filt)
# write output
# output_name_filtered_element_between_cell_types = "80k_NGN2_undiff_WTC11_tested_within_both_celltypes.tsv"
# write.table(df_filt, file = paste(output_dir, output_name_filtered_element_between_cell_types, sep = ""), sep = "\t", row.names = FALSE)

# 53764 elements within both cell-types tested
# downsample barcode df
df_smaller <- downsample_barcodes(df_filt, id_column_name = 'variant_id')
dna <- create_dna_df(df_smaller, id_column_name = "variant_id", allele_column_name="cell_type")
rna <- create_rna_df(df_smaller, id_column_name = "variant_id", allele_column_name = "cell_type")

n_reps <- 3

n_sample <- n_reps

# number of n_samples is number of columns, divided by (nr of samples * number of alleles)
bcs <- ncol(dna) / (n_sample * 2)

print("Call MPRASet ...")
mpra <- MPRASet(DNA = dna, RNA = rna, eid = row.names(dna), barcode = NULL)
design <- data.frame(intcpt = 1, alt = grepl("undiff_WTC11", colnames(mpra)))
print("Finished with MPRASet")

# the replicate where each barcode belongs to is a blocking factor, indicated by the block_vector.
block_vector <- rep(1:n_sample, each = bcs * 2)

print("Start with BCalm...")
start <- Sys.time()

# if you normalized before: # NOTE because of different cell-types => independent groups
mpralm_fit <- mpralm(object = mpra, design = design, aggregate = "none", normalize = FALSE, model_type = "indep_groups", plot = FALSE, block = block_vector)
cat("Finished ... running time: ", Sys.time() - start, "\n")


# TODO: Write output

# write output:
print("Start writing output...")

# Finding significant differences
toptab_allele_bc <- topTable(mpralm_fit, coef = 2, number = Inf)
toptab_allele_bc$variant_id <- row.names(toptab_allele_bc)
write.table(toptab_allele_bc, file = paste(output_dir, output_name, sep = ""), sep = "\t", row.names = FALSE)

print("done")
