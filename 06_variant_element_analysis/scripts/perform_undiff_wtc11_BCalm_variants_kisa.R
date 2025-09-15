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


# define filtering and other variables
bc_threshold <- 10
n_reps <- 3

# input path:
# NGN2: normalized counts using mpralib
input <- "/home/kisa/coding/80K_MPRA/element_analysis_output/202502_element_analysis/processing_mpralib/80k_NGN2_normalized_counts.tsv"
input <- "/home/kisa/coding/80K_MPRA/element_analysis_output/202502_element_analysis/processing_mpralib/80k_NGN2_normalized_counts_renamed_ID2barcode.tsv.gz"
# undifferentiated WTC11: normalized counts using mpralib
input <- "/home/kisa/coding/80K_MPRA/element_analysis_output/202502_element_analysis/processing_mpralib/80k_WTC11_resequencing_normalized_counts.tsv.gz"

# sanity checked variant map for the normalized counts from mpralib (ICCB notebook)
variant_map <- "/home/kisa/coding/80K_MPRA/server_results/80k_counts_after_metadatafile/unique_id_renamed_variant_region_map_mprasnakeflow_input.tsv.gz"
variant_map <- "/home/kisa/coding/80K_MPRA/server_results/80k_counts_after_metadatafile/unique_id_renamed_variant_region_map_mprasnakeflow_input_tilde2comma.tsv.gz"
# does not contain ","

# output paths:
# tested sequences vs negative control
output_dir <- "/home/kisa/coding/80K_MPRA/bc_MPRAlm_results/results/"
output_name <- "NGN2_variant_bcalm_80k_bbmap_no_outlier_removal_normalized_counts_removed_additional_normalization_2025_05_mpralib_design_variant_map_unique_id_all_rep_DNA_RNA.tsv"
output_name <- "testing_undiffWTC11_variant_bcalm_80k_bbmap_no_outlier_removal_normalized_counts_2025_05_mpralib_design_variant_map_unique_id_all_rep_DNA_RNA.tsv"
output_name <- "undiffWTC11_variant_bcalm_80k_bbmap_no_outlier_removal_normalized_counts_2025_05_mpralib_design_variant_map_unique_id_all_rep_DNA_RNA.tsv"


# start preprocessing
df <- read.table(file = gzfile(input), sep = "\t", header = TRUE, comment.char = "")
colnames(df) <- c(
    "barcode", "name", "dna_count_1", "rna_count_1", "dna_count_2",
    "rna_count_2", "dna_count_3", "rna_count_3"
)
# # filtering: raw counts:
# df <- df %>% filter(if_all(matches("DNA", ignore.case = TRUE), ~ . >= 1))
# df <- df %>% filter(if_all(matches("RNA", ignore.case = TRUE), ~ . >= 1))
# filtering normalized counts:
df <- df %>% filter(if_all(matches("DNA", ignore.case = TRUE), ~ . > 0))
df <- df %>% filter(if_all(matches("RNA", ignore.case = TRUE), ~ . > 0))

print("Number of unique names loaded after dna and rna filtering: ")
print(length(unique(df$name))) # undiffWTC11 76213
colnames(df)
# variant map
map_df <- read.table(file = gzfile(variant_map), sep = "\t", header = TRUE, comment.char = "")
map_df <- map_df %>% select(ID, REF, ALT)

print("Create var df")
# var_df <- create_var_df(df, map_df)

# Merge on REF
df_ref <- merge(df, map_df, by.x = "name", by.y = "REF", all.x = FALSE)
head(df_ref)
df_ref$allele <- "ref"
df_ref$ALT <- NULL
print("Number of different oligos after left join with variant map for REF")
print(length(unique(df_ref$name))) # 18047
# Merge on ALT
df_alt <- merge(df, map_df, by.x = "name", by.y = "ALT", all.x = FALSE)
df_alt$allele <- "alt"
df_alt$REF <- NULL
print("Number of different oligos after left join with variant map for ALT")
print(length(unique(df_alt$name))) # 44637
# Combine the results
df_combined <- rbind(df_ref, df_alt)

# Select and rename columns as necessary
var_df <- df_combined %>% select(variant_id = ID, allele, barcode, matches("count"))


print("Number of unique names loaded after create_var_df: ")
print(length(unique(var_df$variant_id))) # 45429

# Filter for oligos with at least min_bc barcodes in both alleles
min_bc <- 10
df_filt <- var_df %>%
    group_by(variant_id, allele) %>%
    filter(n() >= min_bc) %>%
    ungroup()

print("Number of unique variants loaded after filtering for barcodes: ")
print(length(unique(df_filt$variant_id))) # 35230 (lost 10k variants)

df_filt <- df_filt %>%
    group_by(variant_id) %>%
    filter(n_distinct(allele) == 2) %>%
    ungroup()

print("Number of unique variants loaded after filtering for barcodes and allele existence: ")
print(length(unique(df_filt$variant_id))) # wrong names: 26721; correct names: 27384 (lost another 7k variants because of two alleles)

# downsample barcode df
df_smaller <- downsample_barcodes(df_filt, id_column_name = "variant_id")
dna <- create_dna_df(df_smaller, id_column_name = "variant_id")
rna <- create_rna_df(df_smaller, id_column_name = "variant_id")

n_sample <- n_reps

# number of n_samples is number of columns, divided by (nr of samples * number of alleles)
bcs <- ncol(dna) / (n_sample * 2)

print("Call MPRASet ...")
mpra <- MPRASet(DNA = dna, RNA = rna, eid = row.names(dna), barcode = NULL)
design <- data.frame(intcpt = 1, alt = grepl("alt", colnames(mpra)))
print("Finished with MPRASet")

# the replicate where each barcode belongs to is a blocking factor, indicated by the block_vector.
block_vector <- rep(1:n_sample, each = bcs * 2)

print("Start with BCalm...")
start <- Sys.time()
# if you normalized before:
mpralm_fit <- mpralm(object = mpra, design = design, aggregate = "none", normalize = FALSE, model_type = "corr_groups", plot = FALSE, block = block_vector)
cat("Finished ... running time: ", Sys.time() - start, "\n")


# write output:
print("Start writing output...")
# Finding significant variants
toptab_allele_bc <- topTable(mpralm_fit, coef = 2, number = Inf)
toptab_allele_bc$variant_id <- row.names(toptab_allele_bc)
write.table(toptab_allele_bc, file = paste(output_dir, output_name, sep = ""), sep = "\t", row.names = FALSE)
print("done")
