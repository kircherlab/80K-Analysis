#!/usr/bin/env Rscript

library(mpra)
library(dplyr)
library(arrow)

args = commandArgs(trailingOnly=TRUE)
input = "/data/humangen_kircherlab/MPRA/tt_mpra_analysis/data/ProxProm/mpralm/mpralm_input_HepG2.tsv" #args[1]
png(file = "results/HepG2_outputs/mpralm/voom_all.png")
name = "all" #args[2]

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
mpra <- MPRASet(DNA = dna, RNA = rna, eid = row.names(dna), barcode = NULL)
design <- data.frame(intcpt = 1, alt = grepl("alt", colnames(mpra)))
block_vector <- rep(1:s, 2)
weights <- mpralm(object = mpra, design = design, aggregate = "none", normalize = TRUE, block = block_vector, model_type = "corr_groups", return_weights=TRUE)
write_feather(data.frame(weights), paste("./results/HepG2_outputs/mpralm/weights_",name,"_HepG2.feather", sep=""))

mpralm_fit <- mpralm(object = mpra, design = design, aggregate = "none", normalize = TRUE, block = block_vector, model_type = "corr_groups", plot = TRUE)

# Finding significant variants
toptab_allele <- topTable(mpralm_fit, coef = 2, number = Inf)
toptab_allele$variant_id <- row.names(toptab_allele)
write_feather(toptab_allele, paste("./results/HepG2_outputs/mpralm/toptable_",name,"_HepG2.feather", sep=""))
dev.off()






