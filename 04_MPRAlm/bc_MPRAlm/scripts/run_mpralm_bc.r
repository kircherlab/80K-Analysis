#!/usr/bin/env Rscript

library(mpra)
library(dplyr)
library(arrow)
# install.packages("arrow")


# old input
# input_rna = "/data/humangen_kircherlab/MPRA/tt_mpra_analysis/data/ProxProm/mpralm/mpralm_bc_rna_input_HepG2.tsv" #args[1]
# input_dna = "/data/humangen_kircherlab/MPRA/tt_mpra_analysis/data/ProxProm/mpralm/mpralm_bc_dna_input_HepG2.tsv"
# weights_file = "results/HepG2_outputs/mpralm/weights_all_HepG2.feather"
# png(file = "results/HepG2_outputs/mpralm_bc/voom_all.png")
# # args = commandArgs(trailingOnly=TRUE)

# # general 80K input
# input_dir = "/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/"
# output_dir = "/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/bc_MPRAlm/"
# dir.create(file.path(output_dir))

# # input_rna = paste0(input_dir, "preprocess/mpralm_bc_rna_input_NGN2.tsv") #args[1]
# # input_dna = paste0(input_dir, "preprocess/mpralm_bc_dna_input_NGN2.tsv") # args[2]
# name = "lowConfig" # args[4]
# name = "standard" # args[4]
# input_rna = paste0(input_dir, "preprocess/", name,  "_mpralm_bc_rna_input_NGN2.tsv") #args[1]
# input_dna = paste0(input_dir, "preprocess/", name, "_mpralm_bc_dna_input_NGN2.tsv") # args[2]
# weights_file = paste0(input_dir, "MPRAlm/weights_", name, "_NGN2.feather") # args[3]
# png(file = file.path(output_dir, paste0("voom_", name, ".png")))
# output_name = paste0("toptable_", name, "NGN2.feather")

## 80K with controls test (only GC_Selvarajan):
# name = "GC_Selvarajan_standard"

## 80K with some variant controls and barcode filter of 10
name = "standard_with_all_controls"
name = "standard_with_all_variant_controls"
name = "standard_with_all_variant_controls_mendelian_no_downsampling"
name = "different_filtering_variant_controls_mendelian"

input_dir = "/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/"
output_dir = "/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/bc_MPRAlm/no_downsampling/"
input_dir = "/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/preprocess/resequencing/"
output_dir = "/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/bc_MPRAlm/resequencing/"

png(file = file.path(output_dir, paste0("voom_", name, ".png")))

input_dna = paste0(input_dir, name, "_bc_10_mpralm_bc_dna_input_NGN2.tsv")
input_rna = paste0(input_dir, name, "_bc_10_mpralm_bc_rna_input_NGN2.tsv")
# weights_file = paste0(input_dir, "MPRAlm/test_controls_concat/weights_", name, "_NGN2.feather") # args[3]
# output_name = paste0("toptable_", name, "NGN2.feather")
# # local: 
# name = "standard"
# input_dir = "/home/kisa/coding/80K_MPRA/80K-Analysis/04_MPRAlm/results/preprocess/"
# output_dir = "/home/kisa/coding/80K_MPRA/80K-Analysis/04_MPRAlm/results/bc_MPRAlm/locally_computed/"

# # local windows: 
# name = "standard"
# input_dir = "C:\\Users\\user\\OneDrive - Charité - Universitätsmedizin Berlin\\Documents\\coding\\80K_analysis\\bc_mpralm\\results\\preprocess\\"

# output_dir = "C:\\Users\\user\\OneDrive - Charité - Universitätsmedizin Berlin\\Documents\\coding\\80K_analysis\\bc_mpralm\\results\\bc_MPRAlm\\locally_computed\\"


# input_dna = paste0(input_dir, name, "_mpralm_bc_dna_input_NGN2.tsv")
# input_rna = paste0(input_dir, name, "_mpralm_bc_rna_input_NGN2.tsv")
# png(file = file.path(output_dir, paste0("voom_", name, ".png")))
output_name_feather = paste0("toptable_", name, "_NGN2.feather")
output_name = paste0("toptable_", name, "_NGN2.tsv")

# Reading the barcode level input data
rna <- read.table(file=input_rna, sep='\t', header=TRUE)
dna <- read.table(file=input_dna, sep='\t', header=TRUE)
s <- 3

row.names(rna) <- rna$variant_id
rna[,1] <- NULL

row.names(dna) <- dna$variant_id
dna[,1] <- NULL

# number of barcodes is number of columns, divided by (nr of samples * number of alleles)
bcs <- ncol(dna) / (s*2)

mpra_bc_set <- MPRASet(DNA = dna, RNA = rna, eid = row.names(dna), barcode = NULL)
design <- data.frame(intcpt = 1, alt = grepl("alt", colnames(mpra_bc_set)))

# to do element level, use only intercept: design <- data.frame(intcpt = 1)

# the replicate where each barcode belongs to is a blocking factor, indicated by the block_vector. 
block_vector <- rep(1:s, each=bcs*2)

compute_logratio2 <- function(object, aggregate = c("mean", "sum", "none")) {
	dna <- getDNA(object, aggregate = FALSE)
	rna <- getRNA(object, aggregate = FALSE)
	logr <- log2(rna + 1) - log2(dna + 1)
}

# instead of calculating the precision weights, let the mpra package return our calculated weights (optional, I find that it actually makes the model slightly worse)
# assignInNamespace("get_precision_weights", get_precision_weights2, ns="mpra")

# logratio calculation in MPRAlm package turns NAs into 0s, exactly what we don't want.
assignInNamespace("compute_logratio", compute_logratio2, ns="mpra")

mpralm_fit_bc <- mpralm(object = mpra_bc_set, design = design, aggregate = "none", normalize = TRUE, model_type = "corr_groups", plot = TRUE, block = block_vector)
#weights <- mpralm(object = mpra_bc_set, design = design, aggregate = "none", normalize = TRUE, model_type = "corr_groups", plot = TRUE, block = block_vector, return_weights=TRUE)

# Finding significant variants
toptab_allele_bc <- topTable(mpralm_fit_bc, coef = 2, number = Inf)
toptab_allele_bc$variant_id <- row.names(toptab_allele_bc)
write.table(toptab_allele_bc, file=paste(output_dir, output_name, sep=""), sep="\t", row.names = FALSE)
#write_feather(toptab_allele_bc, paste(output_dir, output_name_feather, sep=""))
dev.off()

# For element levels, use treat instead of toptable:
# define user-defined-threshold by plotting the logFC of the negative controls. Find the logFC of the 95th percentile of the negative controls, and set this as the user defined thrshold.
# With the following code, all elements above this threshold will be seen as active.

# tr <- treat(fit, lfc=<user-defined-threshold>)
# mpra_result <- topTreat(tr, coef = 1, number = Inf)

