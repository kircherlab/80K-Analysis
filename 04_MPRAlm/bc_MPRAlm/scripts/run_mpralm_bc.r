#!/usr/bin/env Rscript

library(mpra)
library(dplyr)
library(arrow)


args = commandArgs(trailingOnly=TRUE)
input_rna = "/data/humangen_kircherlab/MPRA/tt_mpra_analysis/data/ProxProm/mpralm/mpralm_bc_rna_input_HepG2.tsv" #args[1]
input_dna = "/data/humangen_kircherlab/MPRA/tt_mpra_analysis/data/ProxProm/mpralm/mpralm_bc_dna_input_HepG2.tsv"
weights_file = "results/HepG2_outputs/mpralm/weights_all_HepG2.feather"
png(file = "results/HepG2_outputs/mpralm_bc/voom_all.png")
name = args[2]

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

mpra <- MPRASet(DNA = dna, RNA = rna, eid = row.names(dna), barcode = NULL)
design <- data.frame(intcpt = 1, alt = grepl("alt", colnames(mpra)))
# the replicate where each barcode belongs to is a blocking factor, indicated by the block_vector. 
block_vector <- rep(1:s, each=bcs*2)

# Get the weights as calculated from the aggregated data, and manipulate to fit to the structure of the barcode data. 
# This is a bit hard coded using the column names
agg_weights <- read_feather(weights_file, as_data_frame = TRUE, mmap = TRUE)
col_order <- gsub("_bc[0-9]*","",colnames(mpra))
weights <- agg_weights[col_order]

get_precision_weights2 <- function(...) {return(weights)}

compute_logratio2 <- function(object, aggregate = c("mean", "sum", "none")) {
	dna <- getDNA(object, aggregate = FALSE)
	rna <- getRNA(object, aggregate = FALSE)
	logr <- log2(rna + 1) - log2(dna + 1)
}

# instead of calculating the precision weights, let the mpra package return our calculated weights (optional, I find that it actually makes the model slightly worse)
# assignInNamespace("get_precision_weights", get_precision_weights2, ns="mpra")

# logratio calculation in MPRAlm package turns NAs into 0s, exactly what we don't want.
assignInNamespace("compute_logratio", compute_logratio2, ns="mpra")

mpralm_fit <- mpralm(object = mpra, design = design, aggregate = "none", normalize = TRUE, model_type = "corr_groups", plot = TRUE, block = block_vector)
#weights <- mpralm(object = mpra, design = design, aggregate = "none", normalize = TRUE, model_type = "corr_groups", plot = TRUE, block = block_vector, return_weights=TRUE)

# Finding significant variants
toptab_allele <- topTable(mpralm_fit, coef = 2, number = Inf)
toptab_allele$variant_id <- row.names(toptab_allele)
write_feather(toptab_allele, paste("./results/HepG2_outputs/mpralm_bc/toptable_all_HepG2.feather", sep=""))
dev.off()

