# install.packages("dplyr")
# install.packages("devtools")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = 'https://cran.uni-muenster.de/')
BiocManager::install("mpra")
install.packages("dplyr")
install.packages("devtools")
install.packages("tidyr")
# install ggplot
install.packages("ggplot2")
library(ggplot2)
library(mpra)
library(dplyr)
library(devtools)
library(tidyr)

load_all('/home/kisa/coding/80K_MPRA/bc_mpralm')

# read input
input = 'data/ProxProm/raw_counts/HepG2_allreps_merged_barcode_assigned_counts.tsv.gz'
input = '/data/cephfs-1/home/users/kisa11_c/unmirrored/projects/MPRA/IGVF_Y1_design/experiment/final_resequencing/results/experiments/standard_bwa/assigned_counts/assignmentFixDuplicates/standardConfig/NGN2_allreps_merged_barcode_assigned_counts.tsv.gz' # snakemake
input = '/data/cephfs-1/scratch/groups/kircher/MPRA/IGVF_Y1_design/experiment/09102024_80K_development_MPRAsnakeflow/results/experiments/bbmapStandardMapq35/assigned_counts/assignmentFixDuplicates/standardConfig/NGN2_allreps_merged_barcode_assigned_counts.tsv.gz' # snakemake
variant_map = "/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/resources/80k_Ngn2/renamed_variant_region_map_mprasnakeflow_input.tsv.gz"

# set output
output_dir = '/home/kisa/coding/80K_MPRA/bc_MPRAlm_results/results/'
output_name = 'toptable_bc_mpralm_80K_bbmap_std_mapq35.tsv'

# start preprocessing
df <- read.table(file=gzfile(input), sep='\t', header=TRUE)

# downsample barcode df
df_smaller <- downsample_barcodes(df)
dna <- create_dna_df(df_smaller, id_column_name="name")
rna <- create_rna_df(df_smaller, id_column_name="name")

map_df <- read.table(file=gzfile(variant_map), sep='\t', header=TRUE)
map_df <- map_df %>% select(ID, REF, ALT)
var_df <- create_var_df(df_smaller, map_df)

# TODO: add optional filtering
# ...



# Create MPRASet
# number of barcodes is number of columns, divided by (nr of samples * number of alleles)
n_sample = 3
bcs <- ncol(dna) / (n_sample*2)
print('Call MPRASet ...')
mpra_bc_set <- MPRASet(DNA = dna, RNA = rna, eid = row.names(dna), barcode = NULL)
design <- data.frame(intcpt = 1, alt = grepl("alt", colnames(mpra_bc_set)))
print('Finished with MPRASet')

# the replicate where each barcode belongs to is a blocking factor, indicated by the block_vector.
block_vector <- rep(1:n_sample, each=bcs*2)

# Fit mpralm (using barcodes)
print('Start with mpralm...')
mpralm_fit_bc <- mpralm(object = mpra_bc_set, design = design, aggregate = "none", normalize = TRUE, block = block_vector, model_type = "corr_groups", plot = TRUE)
#weights <- mpralm(object = mpra_bc_set, design = design, aggregate = "none", normalize = TRUE, model_type = "corr_groups", plot = TRUE, block = block_vector, return_weights=TRUE)
print('Finished')
print('Start writing output...')
# Finding significant variants
toptab_allele_bc <- topTable(mpralm_fit_bc, coef = 2, number = Inf)
toptab_allele_bc$variant_id <- row.names(toptab_allele_bc)
write.table(toptab_allele_bc, file=paste(output_dir, output_name, sep=""), sep="\t", row.names = FALSE)
dev.off()

toptab_allele_bc_not_na = toptab_allele_bc[!is.na(toptab_allele_bc$logFC),]
toptab_allele_bc_not_na
print('done')

#### starting of pias code:
library(devtools)
suppressPackageStartupMessages(load_all("BCalm"))
# suppressPackageStartupMessages(library(mpralm_bc))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
name <- "all" # args[1] #
output <- "/data/humangen_kircherlab/MPRA/tt_mpra_analysis/results/HepG2_outputs/mpralm_bc/" # args[4] #
nr_reps <- 3 # as.integer(args[5]) #


# input = "/data/humangen_kircherlab/MPRA/tt_mpra_analysis/data/ProxProm/mpralm/mpralm_no_outliers_input_HepG2.tsv.gz"
# df <- read.table(file=gzfile(input), sep='\t', header=TRUE)


input <- "data/ProxProm/raw_counts/HepG2_allreps_merged_barcode_assigned_counts.tsv.gz"
df <- read.table(file = gzfile(input), sep = "\t", header = TRUE)
df <- df %>% filter(if_all(matches("DNA", ignore.case = TRUE), ~ . >= 1))
df <- df %>% filter(if_all(matches("RNA", ignore.case = TRUE), ~ . >= 1))

variant_map <- "data/ProxProm/HepG2_variantTable.tsv.gz"
map_df <- read.table(file = gzfile(variant_map), sep = "\t", header = TRUE)
map_df <- map_df %>% select(ID, REF, ALT)
var_df <- create_var_df(df, map_df)

# Function to detect outliers using the z-score method
is_outlier <- function(x) {
    z_scores <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    abs(z_scores) > 3
}
# Filter out rows where any RNA count is an outlier
# var_df <- var_df %>% filter(!if_any(all_of(grep("rna", colnames(df), value = TRUE)), is_outlier))

# Filter for oligos with at least min_bc barcodes in both alleles
df_filt <- var_df %>%
    group_by(variant_id, allele) %>%
    filter(n() >= 10) %>%
    ungroup()
df_filt <- df_filt %>%
    group_by(variant_id) %>%
    filter(n_distinct(allele) == 2) %>%
    ungroup()


df_smaller <- downsample_barcodes(df_filt, "variant_id")
dna <- create_dna_df(df_smaller)
rna <- create_rna_df(df_smaller)

s <- nr_reps

# number of barcodes is number of columns, divided by (nr of samples * number of alleles)
bcs <- ncol(dna) / (s * 2)

mpra <- MPRASet(DNA = dna, RNA = rna, eid = row.names(dna), barcode = NULL)
design <- data.frame(intcpt = 1, alt = grepl("alt", colnames(mpra)))
# the replicate where each barcode belongs to is a blocking factor, indicated by the block_vector.
block_vector <- rep(1:s, each = bcs * 2)

start <- Sys.time()

mpralm_fit <- mpralm(object = mpra, design = design, aggregate = "none", normalize = TRUE, model_type = "corr_groups", plot = FALSE, block = block_vector)
# weights <- mpralm(object = mpra, design = design, aggregate = "none", normalize = TRUE, model_type = "corr_groups", plot = TRUE, block = block_vector, return_weights=TRUE)
cat("running time: ", Sys.time() - start, "\n")

# Finding significant variants
toptab_allele <- topTable(mpralm_fit, coef = 2, number = Inf)
toptab_allele$variant_id <- row.names(toptab_allele)
# dev.off()
# write_feather(toptab_allele, paste(output,"toptable_",name,"_HepG2.feather", sep=""))