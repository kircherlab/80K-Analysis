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

variant_map = "/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/resources/80k_Ngn2/renamed_variant_region_map_mprasnakeflow_input.tsv.gz"
n_reps <- 3 # as.integer(args[5]) #

# set output
output_dir = '/home/kisa/coding/80K_MPRA/bc_MPRAlm_results/results/'
output_name = "bcalm_test_toptable_bc_mpralm_80K_bbmap_std_mapq35.tsv"
output_name = "bcalm_test_toptable_bc_mpralm_80K_bbmap_std_mapq10.tsv"
output_name = "testing_documentation_tmp_variant"

# start preprocessing
df <- read.table(file=gzfile(input), sep='\t', header=TRUE)

# filtering:
df <- df %>% filter(if_all(matches("DNA", ignore.case = TRUE), ~ . >= 1))
df <- df %>% filter(if_all(matches("RNA", ignore.case = TRUE), ~ . >= 1))

print("Number of unique names loaded after dna filtering: ")
print(length(unique(df$name)))

# variant map
map_df <- read.table(file = gzfile(variant_map), sep = "\t", header = TRUE)
map_df <- map_df %>% select(ID, REF, ALT)
# map_df[1,3]
var_df <- create_var_df(df, map_df)
colnames(df)
colnames(map_df)
print("Number of unique names loaded after create_var_df: ")
print(length(unique(var_df$variant_id)))
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
print(length(unique(df_filt$variant_id)))

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
mpralm_fit <- mpralm(object = mpra, design = design, aggregate = "none", normalize = TRUE, model_type = "corr_groups", plot = FALSE, block = block_vector)
# weights <- mpralm(object = mpra, design = design, aggregate = "none", normalize = TRUE, model_type = "corr_groups", plot = TRUE, block = block_vector, return_weights=TRUE)
cat("Finished ... running time: ", Sys.time() - start, "\n")

# write output:
print("Start writing output...")
# Finding significant variants
toptab_allele_bc <- topTable(mpralm_fit, coef = 2, number = Inf)
toptab_allele_bc$variant_id <- row.names(toptab_allele_bc)
write.table(toptab_allele_bc, file = paste(output_dir, output_name, sep = ""), sep = "\t", row.names = FALSE)
# write_feather(toptab_allele_bc, paste(output_dir, output_name_feather, sep=""))

print("done")