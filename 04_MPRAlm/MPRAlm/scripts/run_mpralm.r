#!/usr/bin/env Rscript

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("mpra")
library(mpra)
library(dplyr)
library(arrow)

args = commandArgs(trailingOnly=TRUE)
# # input = "/data/humangen_kircherlab/MPRA/tt_mpra_analysis/data/ProxProm/mpralm/mpralm_input_HepG2.tsv" #args[1]
# # old assignment
# # input = "/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/preprocess/mpralm_input_NGN2.tsv" #args[1]
# # output_dir = "/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/MPRAlm/"
# # name = "all" #args[2]
# # new assignment
# name = "lowConfig" #args[2]
# name = "standard" #args[2]
# input = paste0("/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/preprocess/", name, "_mpralm_input_NGN2.tsv") #args[1]
# output_dir = "/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/MPRAlm/"



## 80K with some variant controls and barcode filter of 10
name = "standard_with_variant_controls_no_downsampling"

# newest version: matching problem of mendelian solved
name = "standard_with_all_variant_controls_mendelian"

# newest version: matching problem of mendelian solved + no downsampling
name = "standard_with_all_variant_controls_mendelian_no_downsampling"

input_dir = "/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/preprocess/no_downsampling/"
output_dir = "/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/MPRAlm/no_downsampling/"

input = paste0(input_dir, name, "_bc_10_mpralm_input_NGN2.tsv") #args[1]

# # control example + cardiac_neuro_cava_random
# name = "GC_Selvarajan_standard"
# input_dir = "/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/preprocess/test_controls_concat/"
# output_dir = "/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/MPRAlm/test_controls_concat/"
# input = paste0(input_dir, name, "_mpralm_input_NGN2.tsv") #args[1]


png(file=file.path(output_dir, paste0("voom_mpralm_",name,".png")))
output_name_feather = paste0("toptable_", name, "_mpralm_NGN2.feather")
output_name = paste0("toptable_", name, "_mpralm_NGN2.tsv")

# Reading the barcode level input data
raw_data <- read.table(file=input, sep='\t', header=TRUE)
s <- 3

# Summing the barcode level data per oligo to create mpralm input
rna <- raw_data[, grep("variant_id|allele|RNA", colnames(raw_data))] %>% group_by(variant_id, allele) %>% summarise(across(everything(), sum),.groups = 'drop')
rna <- cbind(filter(rna,allele=="ref")[,-2], filter(rna,allele=="alt")[,3:(s+2)])
names(rna)[-1] <- paste0("sample", c(1:s), "_", rep(c("ref", "alt"), each = s))
row.names(rna) <- rna$variant_id
rna[,1] <- NULL

dna <- raw_data[, grep("variant_id|allele|DNA", colnames(raw_data))] %>% group_by(variant_id, allele) %>% summarise(across(everything(), sum),.groups = 'drop')
dna <- cbind(filter(dna,allele=="ref")[,-2], filter(dna,allele=="alt")[,3:(s+2)])
names(dna)[-1] <- paste0("sample", c(1:s), "_", rep(c("ref", "alt"), each = s))
row.names(dna) <- dna[,1]
dna[,1] <- NULL

# Running MPRAlm
mpralm_set <- MPRASet(DNA = dna, RNA = rna, eid = row.names(dna), barcode = NULL)
design <- data.frame(intcpt = 1, alt = grepl("alt", colnames(mpralm_set)))
block_vector <- rep(1:s, 2)
# for 80K: 21102 rows 
weights_mpralm <- mpralm(object = mpralm_set, design = design, aggregate = "none", normalize = TRUE, block = block_vector, model_type = "corr_groups", return_weights=TRUE)

# writing cannot open local file (no such file) if relative paths
write_feather(data.frame(weights_mpralm), file.path(output_dir, paste("weights_", name, "_mpralm_NGN2.feather", sep="")))

mpralm_fit_mpralm <- mpralm(object = mpralm_set, design = design, aggregate = "none", normalize = TRUE, block = block_vector, model_type = "corr_groups", plot = TRUE)

# Finding significant variants
toptab_allele_mpralm <- topTable(mpralm_fit_mpralm, coef = 2, number = Inf)
toptab_allele_mpralm$variant_id <- row.names(toptab_allele_mpralm)
write.table(toptab_allele_mpralm, file=paste(output_dir, output_name, sep=""), sep="\t", row.names = FALSE)
#write_feather(toptab_allele_mpralm, file.path(output_dir, output_name))

dev.off()





