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

load_all('/home/kisa/coding/80K_MPRA/bc_mpralm')

# read input
# input = 'data/ProxProm/raw_counts/HepG2_allreps_merged_barcode_assigned_counts.tsv.gz'
# input = '/data/cephfs-1/home/users/kisa11_c/unmirrored/projects/MPRA/IGVF_Y1_design/experiment/final_resequencing/results/experiments/standard_bwa/assigned_counts/assignmentFixDuplicates/standardConfig/NGN2_allreps_merged_barcode_assigned_counts.tsv.gz' # snakemake
# input = '/data/cephfs-1/scratch/groups/kircher/MPRA/IGVF_Y1_design/experiment/09102024_80K_development_MPRAsnakeflow/results/experiments/bbmapStandardMapq35/assigned_counts/assignmentFixDuplicates/standardConfig/NGN2_allreps_merged_barcode_assigned_counts.tsv.gz' # snakemake
input = "/data/cephfs-1/scratch/groups/kircher/MPRA/IGVF_Y1_design/experiment/09102024_80K_development_MPRAsnakeflow/unique_variant_current_development_bbmap_std_mapq35_80K_mpra_final_resequencing/experiments/bbmapStandardMapq35/assigned_counts/assignmentFixDuplicates/standardConfig/NGN2_allreps_merged_barcode_assigned_counts.tsv.gz" # snakemake # unique_variant_bbmap 35: 70500 => testable: 38362
input = "/home/kisa/coding/80K_MPRA/server_results/MPRAsnakeflow/bbmap_standard_mapq10NoLength_unique_variant_id/NGN2_allreps_merged_barcode_assigned_counts.tsv.gz" # unique_variant_bbmap 10: 70798 => testable: 38589
input = "/home/kisa/coding/80K_MPRA/server_results/MPRAsnakeflow/bbmap_standard_mapq10NoLength_unique_variant_id/testing_ngn2_MRPAsnakeflow_counts_302067.tsv.gz" # unique_variant_bbmap 10:  => testable: 323
input = "/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/experiment/20241113_80K_MPRAsnakeflow/results/experiments/mpra80KNeuronbbmapmapq30BC10DNA1RNA1/assigned_counts/assignmentFixDuplicates/default/NGN2_allreps_merged_barcode_assigned_counts.tsv.gz"

# removed hashtag rows (comment char problems and not important for the variant results)
input <- "/home/kisa/coding/80K_MPRA/server_results/80k_counts_after_metadatafile/NGN2_allreps_merged_barcode_assigned_counts.tsv.gz"
input = "/home/kisa/coding/80K_MPRA/server_results/80k_counts_after_metadatafile/NGN2_allreps_merged_barcode_assigned_counts_hashtag_expanded.tsv.gz"
input <- "/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/experiment/20241113_80K_MPRAsnakeflow/results/experiments/mpra80KNeuronbbmapmapq30BC10DNA1RNA1/assigned_counts/assignmentFixDuplicates/default/NGN2_allreps_merged_barcode_assigned_counts.tsv.gz"
input <- "/home/kisa/coding/80K_MPRA/server_results/80k_counts_after_metadatafile/NGN2_allreps_merged_barcode_assigned_counts_hashtag_expanded.tsv.gz"
input <- "/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/experiment/80K_collision_skip_with_adapter/results/experiments/standard_bwa/assigned_counts/assignmentFixDuplicates/standardConfig/NGN2_allreps_merged_barcode_assigned_counts.tsv.gz"
input <- "/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/experiment/strand_sensitive_80K_NGN2/results/experiments/mpra80KNGN2bbmapmapq30BC10DNA1RNA1/assigned_counts/assignmentFixDuplicates/default/NGN2_allreps_merged_barcode_assigned_counts.tsv.gz"
# Finished ... running time:  38.09017
input <- "/home/kisa/coding/80K_MPRA/server_results/80k_counts_after_metadatafile/bwa_finest_NGN2_allreps_merged_barcode_assigned_counts_outlier_removed.tsv.gz"
input <- "/home/kisa/coding/80K_MPRA/server_results/80k_counts_after_metadatafile/bwa_finest_NGN2_allreps_merged_barcode_assigned_counts_hashtag_expanded.tsv.gz"
input <- "/home/kisa/coding/80K_MPRA/server_results/80k_counts_after_metadatafile/NGN2_allreps_merged_barcode_assigned_counts_hashtag_expanded.tsv.gz"

# wtc11 standard design strand sensitive + outlier detection
input <- "/home/kisa/coding/80K_MPRA/server_results/80k_counts_after_metadatafile/strand_sensitive_WTC11_allreps_merged_barcode_assigned_counts_outlier_removed.tsv.gz"


# NGN2: normalized counts using mpralib
input <- "/home/kisa/coding/80K_MPRA/element_analysis_output/202502_element_analysis/processing_mpralib/80k_NGN2_normalized_counts.tsv"
input <- "/home/kisa/coding/80K_MPRA/element_analysis_output/202502_element_analysis/processing_mpralib/80k_NGN2_normalized_counts_renamed_ID2barcode.tsv.gz"
input <- "/home/kisa/coding/80K_MPRA/element_analysis_output/202502_element_analysis/processing_mpralib/80k_NGN2_resequencing_SelfmadeStrandSensitive_normalized_counts.tsv.gz"


variant_map = "/home/kisa/coding/80K_MPRA/server_results/80k_counts_after_metadatafile/renamed_variant_region_map_mprasnakeflow_input.tsv.gz"

# sanity checked variant map for the normalized counts from mpralib (ICCB notebook)
variant_map = "/home/kisa/coding/80K_MPRA/server_results/80k_counts_after_metadatafile/unique_id_renamed_variant_region_map_mprasnakeflow_input.tsv.gz"
# tilde2comma variant map (needed for undiffWTC11)
variant_map <- "/home/kisa/coding/80K_MPRA/server_results/80k_counts_after_metadatafile/unique_id_renamed_variant_region_map_mprasnakeflow_input_tilde2comma.tsv.gz"
n_reps <- 3 # as.integer(args[5]) #

# set output
output_dir = '/home/kisa/coding/80K_MPRA/bc_MPRAlm_results/results/'
output_name = "bcalm_test_toptable_bc_mpralm_80K_bbmap_std_mapq35.tsv"
output_name = "bcalm_test_toptable_bc_mpralm_80K_bbmap_std_mapq10.tsv"
output_name = "testing_documentation_tmp_variant"
output_name = "variant_bcalm_80K_bbmap_std_mapq30_no_hashtag.tsv"
output_name <- "variant_bcalm_80K_bbmap_std_mapq30_no_hashtag_variant_map_unique_variant_id.tsv"
output_name <- "variant_bcalm_80K_bbmap_std_mapq30_no_hashtag_no_outlier_removal_variant_map_unique_variant_id.tsv"
output_name <- "variant_bcalm_80K_bwa_no_outlier_removal_no_hashtag_exists_variant_map_unique_variant_id.tsv"
output_name <- "variant_bcalm_80K_bbmap_std_mapq30_no_outlier_removal_no_hashtag_exists_variant_map_unique_variant_id.tsv"
output_name <- "variant_bcalm_80k_bwa_finest_outlier_removal_no_hashtag_exists_variant_map_unique_id.tsv"
output_name <- "variant_bcalm_80k_bwa_finest_withouth_outlier_removal_no_hashtag_exists_variant_map_unique_id.tsv"
output_name <- "variant_bcalm_80k_bbmap_outlier_removal_no_hashtag_exists_variant_map_unique_id.tsv"
output_name <- "WTC11_variant_bcalm_80k_bbmap_outlier_removal_default_design_variant_map_unique_id.tsv"
output_name <- "NGN2_variant_bcalm_80k_bbmap_no_outlier_removal_normalized_counts_mpralib_design_variant_map_unique_id_all_rep_DNA_RNA.tsv"
output_name <- "NGN2_variant_bcalm_80k_bbmap_no_outlier_removal_normalized_counts_removed_additional_normalization_2025_05_mpralib_design_variant_map_unique_id_all_rep_DNA_RNA.tsv"
output_name <- "NGN2_variant_bcalm_80k_bbmap_selfmadeStrandSensitive_no_outlier_removal_normalized_counts_removed_additional_normalization_2025_05_mpralib_design_variant_map_unique_id_all_rep_DNA_RNA.tsv"


# start preprocessing
df <- read.table(file=gzfile(input), sep='\t', header=TRUE, comment.char = "")
colnames(df)
colnames(df) <- c(
    "barcode", "name", "dna_count_1", "rna_count_1", "dna_count_2",
    "rna_count_2", "dna_count_3", "rna_count_3"
)
# # filtering:
# df <- df %>% filter(if_all(matches("DNA", ignore.case = TRUE), ~ . >= 1))
# df <- df %>% filter(if_all(matches("RNA", ignore.case = TRUE), ~ . >= 1))
# filtering normalized counts:
df <- df %>% filter(if_all(matches("DNA", ignore.case = TRUE), ~ . > 0))
df <- df %>% filter(if_all(matches("RNA", ignore.case = TRUE), ~ . > 0))

print("Number of unique names loaded after dna filtering: ")
print(length(unique(df$name)))

# variant map
map_df <- read.table(file = gzfile(variant_map), sep = "\t", header = TRUE, comment.char = "")
map_df <- map_df %>% select(ID, REF, ALT)
# map_df[1,3]
colnames(df)
print("Create var df")
# var_df <- create_var_df(df, map_df)

# Merge on REF
df_ref <- merge(df, map_df, by.x = "name", by.y = "REF", all.x = FALSE)
df_ref$allele <- "ref"
df_ref$ALT <- NULL
print("Number of different oligos after left join with variant map for REF")
print(length(unique(df_ref$name))) # NGN2: 18365; undiffWTC11: 18047

# Merge on ALT
df_alt <- merge(df, map_df, by.x = "name", by.y = "ALT", all.x = FALSE)
df_alt$allele <- "alt"
df_alt$REF <- NULL
print("Number of different oligos after left join with variant map for ALT")
print(length(unique(df_alt$name))) # NGN2: 45456; undiffWTC11: 44637
# Combine the results
df_combined <- rbind(df_ref, df_alt)

# Select and rename columns as necessary
var_df <- df_combined %>% select(variant_id = ID, allele, barcode, matches("count"))

print("Number of unique names loaded after create_var_df: ")
print(length(unique(var_df$variant_id))) # NGN2: 45847; undiffWTC11: 45429

colnames(var_df)
# Filter for oligos with at least min_bc barcodes in both alleles
min_bc <- 10
df_filt <- var_df %>%
    group_by(variant_id, allele) %>%
    filter(n() >= min_bc) %>%
    ungroup()
df_filt <- df_filt %>%
    group_by(variant_id) %>%
    filter(n_distinct(allele) == 2) %>%
    ungroup()

print("Number of unique variants loaded after filtering for barcodes and allele existence: ")
print(length(unique(df_filt$variant_id))) # NGN2: 39498; undiffWTC11: correct names: 27384

# downsample barcode df
df_smaller <- downsample_barcodes(df_filt, id_column_name = 'variant_id')
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
# mpralm_fit <- mpralm(object = mpra, design = design, aggregate = "none", normalize = TRUE, model_type = "corr_groups", plot = FALSE, block = block_vector)
# if you normalized before:
mpralm_fit <- mpralm(object = mpra, design = design, aggregate = "none", normalize = FALSE, model_type = "corr_groups", plot = FALSE, block = block_vector)
cat("Finished ... running time: ", Sys.time() - start, "\n") # NGN2 35.38816

# write output:
print("Start writing output...")
# Finding significant variants
toptab_allele_bc <- topTable(mpralm_fit, coef = 2, number = Inf)
toptab_allele_bc$variant_id <- row.names(toptab_allele_bc)
write.table(toptab_allele_bc, file = paste(output_dir, output_name, sep = ""), sep = "\t", row.names = FALSE)
# write_feather(toptab_allele_bc, paste(output_dir, output_name_feather, sep=""))

print("done")