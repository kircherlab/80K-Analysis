# install.packages("dplyr")
# install.packages("devtools")
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager", repos = 'https://cran.uni-muenster.de/')
# BiocManager::install("mpra")
# install.packages("dplyr")
# install.packages("devtools")
# install.packages("tidyr")
# install ggplot
# library(mpra)
# install.packages("ggplot2")
library(ggplot2)
library(dplyr)
library(devtools)
library(tidyr)

# bcalm environment: conda-forge::r-dplyr conda-forge::r-devtools conda-forge::r-tidyr conda-forge::r-ggplot2
# mamba create -n BCalm conda-forge::r-dplyr conda-forge::r-devtools conda-forge::r-tidyr conda-forge::r-ggplot2 -y
# BiocGenerics, limma, SummarizedExperiment
# bioconda::bioconductor-biocgenerics bioconda::bioconductor-limma bioconda::bioconductor-summarizedexperiment
load_all("/home/kisa/coding/80K_MPRA/bc_mpralm")


# define filtering and other variables
bc_threshold <- 10
nr_reps <- 3
# compare this dataset only against the scrambled control
# group names for the group comparison
test_label = "cardiac_neuro_cava_random"
# negative_control = "negative_neuron_ctrl"
scramble_control = "scramble_ctrl"

# input path:
# undifferentiated WTC11: normalized counts using mpralib
input <- "/home/kisa/coding/80K_MPRA/element_analysis_output/202502_element_analysis/processing_mpralib/80k_WTC11_resequencing_normalized_counts.tsv.gz"
# label file for group comparisons
labelfile <- "/home/kisa/coding/80K_MPRA/element_analysis_output/202502_element_analysis/processing_mpralib/80k_metadata_scramble_control_label.tsv.gz" # created from: /home/kisa/coding/80K_MPRA/80K-Analysis/06_variant_element_analysis/notebooks/qc_of_80k_mpra.ipynb

# output paths:
# # tested sequences vs negative control
# topTreat_result_diff_negative_control <- "/home/kisa/coding/80K_MPRA/element_analysis_output/undiffWTC11_202505_80k_element_bbmap_no_outlier_removal_bc10_toptreat_results_negative_ctrl_result_095_both_sied.tsv"
# # tested sequences vs scrambled controls
# topTreat_result_diff_scrambled <- "/home/kisa/coding/80K_MPRA/element_analysis_output/undiffWTC11_202505_80k_element_bbmap_no_outlier_removal_bc10_toptreat_results_scrambled_NP_MK_result_095_both_sided.tsv"
topTreat_result_diff_scrambled <- "/home/kisa/coding/80K_MPRA/element_analysis_output/undiffWTC11_202505_80k_element_bbmap_no_outlier_removal_bc10_toptreat_results_scrambled_NP_MK_result_075_both_sided.tsv"


# read count input
df <- read.table(file = gzfile(input), sep = "\t", header = TRUE, comment.char = "")
print("Number of unique names loaded: ") # undiffWTC11: 77876
print(length(unique(df$name)))

# filter counts for minimum number of barcodes
df_filt <- df %>%
    group_by(name) %>%
    filter(n() >= bc_threshold) %>%
    ungroup()
nrow(df_filt)
# undiffWTC11: barcodes: 5508573

# number of oligos in file remaining:
print("Number of unique names after barcode filtering: ")
print(length(unique(df_filt$name)))
# undiffWTC11: oligos: 72897

# subsample the count data for easier computations:
df_smaller <- downsample_barcodes(df_filt)
dna <- create_dna_df(df_smaller, id_column_name = "name")
rna <- create_rna_df(df_smaller, id_column_name = "name")


# read label file for the different groups and filter them
# NOTE I assume no header given:
labels <- read.table(gzfile(labelfile), header = FALSE, sep = "\t", col.names = c("name", "label"), comment.char = "")

labels_vec <- as.vector(labels$label)
names(labels_vec) <- labels$name
# Use only these labels of the sequences that remained after filtering
labels_vec <- labels_vec[rownames(dna)]

# Start processing the count data:
# number of barcodes is number of columns, divided by (nr of samples * number of alleles)
bcs <- ncol(dna) / nr_reps
print("Start with MPRASet...")
mpra <- MPRASet(DNA = dna, RNA = rna, eid = row.names(dna), barcode = NULL, label = labels_vec)

# the replicate where each barcode belongs to is a blocking factor, indicated by the block_vector.
block_vector <- rep(1:nr_reps, each = bcs)

print("Start with BCalm...")
start <- Sys.time()
# NOTE: no normalization
mpralm_fit <- fit_elements(object = mpra, normalize = FALSE, block = block_vector, plot = FALSE)
print("Finnished BCalm ...") # 22 minuten
cat("running time: ", Sys.time() - start, "\n")

# testing the sequences for group differences compared to a negative control and a scrambled control

# # VS: negative controls
# testing_element_negative_control <- mpra_treat(mpralm_fit, 0.95, neg_label = negative_control)
# # View(testing_element_scrambled)

# element_toptreat <- topTreat(testing_element_negative_control, coef = 1, number = Inf)

# element_toptreat$name <- row.names(element_toptreat)
# write.table(element_toptreat, file = topTreat_result_diff_negative_control, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# test_element_plot <- plot_groups(mpralm_fit, 0.95, neg_label = negative_control, test_label = test_label)
# print(test_element_plot)

# VS: scrambled controls
testing_element_scrambled <- mpra_treat(mpralm_fit, 0.75, neg_label = scramble_control)
element_toptreat <- topTreat(testing_element_scrambled, coef = 1, number = Inf)
element_toptreat$name <- row.names(element_toptreat)
write.table(element_toptreat, file = topTreat_result_diff_scrambled, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

test_element_plot <- plot_groups(mpralm_fit, 0.75, neg_label = scramble_control, test_label = test_label)
print(test_element_plot)
