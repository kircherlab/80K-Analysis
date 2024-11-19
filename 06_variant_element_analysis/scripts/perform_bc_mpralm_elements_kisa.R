# input = 'data/ProxProm/raw_counts/HepG2_allreps_merged_barcode_assigned_counts.tsv.gz' # snakemake
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
load_all('/home/kisa/coding/80K_MPRA/bc_mpralm')

input = '/data/cephfs-1/home/users/kisa11_c/unmirrored/projects/MPRA/IGVF_Y1_design/experiment/final_resequencing/results/experiments/standard_bwa/assigned_counts/assignmentFixDuplicates/standardConfig/NGN2_allreps_merged_barcode_assigned_counts.tsv.gz' # snakemake
input = '/home/kisa/coding/80K_MPRA/element_analysis_output/element_assigned_barcodes.tsv.gz' # modified snakemake output (see: /home/kisa/coding/80K_MPRA/80K-Analysis/07_quality_control/notebooks/state_of_data.ipynb)
input = '/home/kisa/coding/80K_MPRA/element_analysis_output/scrambled_vs_tested_element_assigned_barcodes.tsv.gz' # modified snakemake output (see: /home/kisa/coding/80K_MPRA/80K-Analysis/07_quality_control/notebooks/state_of_data.ipynb)
input = '/home/kisa/coding/80K_MPRA/element_analysis_output/negative_neuron_NP_vs_tested_element_assigned_barcodes_no_alt.tsv.gz' # modified snakemake output (see: /home/kisa/coding/80K_MPRA/80K-Analysis/07_quality_control/notebooks/state_of_data.ipynb)
input = '/home/kisa/coding/80K_MPRA/element_analysis_output/negative_neuron_CTRLs_vs_tested_element_assigned_barcodes_no_alt.tsv.gz' # modified snakemake output (see: /home/kisa/coding/80K_MPRA/80K-Analysis/07_quality_control/notebooks/state_of_data.ipynb)
input = "/home/kisa/coding/80K_MPRA/server_results/MPRAsnakeflow/bbmap_standard_mapq10NoLength_unique_variant_id/NGN2_allreps_merged_barcode_assigned_counts.tsv.gz"
# input = "/home/kisa/coding/80K_MPRA/server_results/MPRAsnakeflow/bbmap_standard_mapq10NoLength_unique_variant_id/testing_ngn2_MRPAsnakeflow_counts_302067.tsv.gz"

bc_threshold <- 50
bc_threshold <- 10
dna_count_threshold <- 1
rna_count_threshold <- 1
nr_reps <- 3 # if element for variant see below
df <- read.table(file=gzfile(input), sep='\t', header=TRUE)
print("Number of unique names loaded: ")
print(length(unique(df$name)))
# Filter all rows where DNA and RNA are both >= 1 in all replicates
df <- df %>% filter(if_all(matches("dna", ignore.case=TRUE), ~ . >= dna_count_threshold))
print("Number of unique names loaded after dna filtering: ")
print(length(unique(df$name)))
df <- df %>% filter(if_all(matches("rna", ignore.case=TRUE), ~ . >= rna_count_threshold))
print("Number of unique names loaded after rna filtering: ")
print(length(unique(df$name)))
# Filter for oligos with at least min_bc barcodes per oligo
df_filt <- df %>% group_by(name) %>% filter(n() >= bc_threshold) %>% ungroup()
nrow(df_filt)

# number of oligos in file remaining:
print("Number of unique names after barcode filtering: ")
print(length(unique(df_filt$name)))


# element assigned: 27295
# Filter: removing "ALT_" rows (not working)
#df_filt_alt <- df_filt %>%
#  filter(grepl(":ATL", name))
#typeof(df_filt)

# write df_filt to file
# write.table(df_filt, file = "/home/kisa/coding/80K_MPRA/element_analysis_output/negative_neuron_NP_vs_tested_element_assigned_barcodes_no_alt_filtered.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

df_smaller <- downsample_barcodes(df_filt)
dna <- create_dna_df(df_smaller, id_column_name="name")
rna <- create_rna_df(df_smaller, id_column_name="name")

# labelfile = 'data/ProxProm/labels.tsv'
labelfile = '/data/cephfs-1/home/users/kisa11_c/unmirrored/projects/MPRA/IGVF_Y1_design/resources/association_data/design_no_duplicates_sequence_and_header_label.tsv'
labelfile = '/home/kisa/coding/80K_MPRA/design_data/removed_brackets_design_no_duplicates_sequence_and_header_label.fa'
labelfile = '/home/kisa/coding/80K_MPRA/element_analysis_output/element_assigned_barcodes_label.tsv'
labelfile = '/home/kisa/coding/80K_MPRA/element_analysis_output/negative_neuron_CTRLs_vs_tested_element_assigned_barcodes_no_alt_label_file.tsv'
labelfile = "/home/kisa/coding/80K_MPRA/design_data/design_info/renamed_design_no_duplicates_sequence_and_header_with_adapter_no_brackets_no_collisions_label.tsv"

labels <- read.table(labelfile, header=TRUE, sep='\t', col.names=c('name', 'label'))
labels_vec <- as.vector(labels$label)
names(labels_vec) <- labels$name
# Use only these labels of the sequences that remained after filtering
labels_vec <- labels_vec[rownames(dna)]

length(unique(labels$label))
unique(labels$label)

# For now I just remove all alt alleles
# TODO: remove ALT
# number of barcodes is number of columns, divided by (nr of samples * number of alleles)
bcs <- ncol(dna) / nr_reps

mpra <- MPRASet(DNA = dna, RNA = rna, eid = row.names(dna), barcode = NULL, label=labels_vec)

# the replicate where each barcode belongs to is a blocking factor, indicated by the block_vector.
block_vector <- rep(1:nr_reps, each=bcs)
start <- Sys.time()
mpralm_fit <- fit_elements(object = mpra, normalize=TRUE, block = block_vector)
cat("running time: ", Sys.time() - start, "\n")

# filter mpralm_fit by only having tested sequences and C_negative_neuron_NP
# Right? (or are we interested in the performance of the controls?)

# how many rows with logratio below -0.6
# filtered_df <- mpralm_fit[mpralm_fit$logratio >= -0.6, ]
# to_test_test <- mpralm_fit[mpralm_fit$label == "cardiac_neuro_cava_random", ]
# to_negative_test <- mpralm_fit[mpralm_fit$label == "C_negative_neuron_NP", ]
# higher than -0.6 ( which is the threshol 0.95 percentile)
# sum(to_test_test$logratio >= -0.6, na.rm = TRUE)

# sum(to_test_test$logratio <= -0.6, na.rm = TRUE)

# length(to_negative_test$p.value)

# test_filtered_df <- to_test_test[to_test_test$logratio >= -0.6, ]
# nrow(filtered_df)
# # double dipping test: I filter for sequences which have high values: logratio > -0.6 and use treat again
# filtered_df <- mpralm_fit[mpralm_fit$logratio > -0.6, ]

# neu_negative_np_mean <- mean(mpralm_fit[mpralm_fit$label == "C_negative_neuron_NP", ]$logratio)
# neu_negative_np_mean
# mpralm_fit[mpralm_fit$label == "C_negative_neuron_NP", ]$label
# # mpralm_fit$neu_negative_np_normalized_logratio <- mpralm_fit$

# test_element <- mpralm_fit
# test_element$neu_negative_np_normalized_logratio <- TRUE

testing_element <- mpra_treat(mpralm_fit, 0.95, neg_label = "C_negative_neuron_NP", test_label = NULL, side = "right") # 68297 expected 69257 - 92
View(testing_element)
testing_element_both_sided <- mpra_treat(mpralm_fit, 0.95, neg_label = "C_negative_neuron_NP", test_label = NULL, side = "both") # 69165 expected 69257

write.table(testing_element, file = "/home/kisa/coding/80K_MPRA/element_analysis_output/testing_element_treat_result.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
toptab_allele_bc <- topTable(mpralm_fit, coef = 2, number = Inf)
toptab_allele_bc$name <- row.names(toptab_allele_bc)
write.table(toptab_allele_bc, file = "/home/kisa/coding/80K_MPRA/element_analysis_output/testing_element_toptable_results.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
test_element_plot <- plot_groups(mpralm_fit, percentile = 0.95, neg_label = "C_negative_neuron_NP", test_label = "cardiac_neuro_cava_random")
print(test_element_plot)

# # define negative label: C_negative_neuron_NP
# test_mpra_treat_negative_np <- mpra_treat(mpralm_fit, 0.95, "C_negative_neuron_NP", "cardiac_neuro_cava_random", side="right")
# test_mpra_treat_positive_NP <- mpra_treat(mpralm_fit, 0.95, "C_positive_neuron_NP", "cardiac_neuro_cava_random", side="left")
# test_mpra_treat_mk_scrambled_only_two_groups <- mpra_treat(mpralm_fit, 0.95, "MK_scrambled", "cardiac_neuro_cava_random", side="both")
# test_mpra_treat_negative_neuron_NP_only_two_groups <- mpra_treat(mpralm_fit, 0.95, "C_negative_neuron_NP", "cardiac_neuro_cava_random", side="right")
# test_mpra_treat_negative_neuron_combined_controls_only_two_groups <- mpra_treat(mpralm_fit, 0.95, "C_combined_negative_neuron_MK_NP", "cardiac_neuro_cava_random", side="right")

# View(test_mpra_treat)
# View(test_mpra_treat_mk_scrambled_only_two_groups)
# View(test_mpra_treat_negative_neuron_NP_only_two_groups)
# View(test_mpra_treat_negative_neuron_combined_controls_only_two_groups)
# # with our without row names
# # most important: log; pvalue...
# # write.table(test_mpra_treat, file = "/home/kisa/coding/80K_MPRA/element_analysis_output/bc_mpralm_elements_with_alt_all_controls.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(test_mpra_treat_negative_np, file = "/home/kisa/coding/80K_MPRA/element_analysis_output/bc_mpralm_elements_with_subset_controls_NP_neuro_negative_right_sided_tested.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(test_mpra_treat_positive_NP, file = "/home/kisa/coding/80K_MPRA/element_analysis_output/bc_mpralm_elements_with_subset_controls_NP_neuro_positive_left_sided_tested.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(test_mpra_treat_mk_scrambled, file = "/home/kisa/coding/80K_MPRA/element_analysis_output/bc_mpralm_elements_with_subset_controls_MK_scrambled_both_sided_tested.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(test_mpra_treat_mk_scrambled_only_two_groups, file = "/home/kisa/coding/80K_MPRA/element_analysis_output/bc_mpralm_MK_scrambled_vs_tested_elements_both_sided_no_alt.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(test_mpra_treat_negative_neuron_NP_only_two_groups, file = "/home/kisa/coding/80K_MPRA/element_analysis_output/bc_mpralm_negative_neuron_NP_vs_tested_elements_right_side_no_alt_unique.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(test_mpra_treat_negative_neuron_NP_only_two_groups, file = "/home/kisa/coding/80K_MPRA/element_analysis_output/bc_mpralm_negative_neuron_NP_vs_all_tested_elements_right_side_unique.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(test_mpra_treat_negative_neuron_combined_controls_only_two_groups, file = "/home/kisa/coding/80K_MPRA/element_analysis_output/bc_mpralm_negative_neuron_combined_vs_tested_elements_right_side_no_alt.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# # plot groups
# test_plotting_groups_negative_neuron_np <- plot_groups(mpralm_fit, percentile = 0.95, neg_label= "C_negative_neuron_NP", test_label = "cardiac_neuro_cava_random")
# test_plotting_groups_positive_neuron_np <- plot_groups(mpralm_fit, percentile = 0.95, neg_label= "C_positive_neuron_NP", test_label = "cardiac_neuro_cava_random")
# test_plotting_groups_mk_scrambled <- plot_groups(mpralm_fit, percentile = 0.95, neg_label= "MK_scrambled", test_label = "cardiac_neuro_cava_random")
# test_plotting_groups_mk_scrambled_two_groups <- plot_groups(mpralm_fit, percentile = 0.95, neg_label= "MK_scrambled", test_label = "cardiac_neuro_cava_random")
# test_mpra_treat_negative_neuron_NP_only_two_groups_two_groups <- plot_groups(mpralm_fit, percentile = 0.95, neg_label= "C_negative_neuron_NP", test_label = "cardiac_neuro_cava_random")
# test_plotting_negative_neuron_NP_MK_combined <- plot_groups(mpralm_fit, percentile = 0.95, neg_label= "C_combined_negative_neuron_MK_NP", test_label = "cardiac_neuro_cava_random")

# print(test_plotting_groups_negative_neuron_np)
# print(test_plotting_groups_positive_neuron_np)
# print(test_plotting_groups_mk_scrambled)
# print(test_plotting_groups_mk_scrambled_two_groups)
# print(test_plotting_groups_mk_scrambled_two_groups)
# print(test_mpra_treat_negative_neuron_NP_only_two_groups_two_groups)
# print(test_plotting_negative_neuron_NP_MK_combined)


# # plot all labels: leave the neg_label
# test_plotting_all_groups <- plot_groups(mpralm_fit, percentile = 0.95, test_label = "cardiac_neuro_cava_random")
# print(test_plotting_all_groups)

# # plot all important labels:
# mpralm_fit_important_groups <- mpralm_fit %>%
#   filter(!grepl('GC', label))

# test_plotting_important_groups <- plot_groups(mpralm_fit, percentile = 0.95, test_label = "cardiac_neuro_cava_random")
