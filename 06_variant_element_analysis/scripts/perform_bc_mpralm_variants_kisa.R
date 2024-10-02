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
