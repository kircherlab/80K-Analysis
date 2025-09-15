# server results: /data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/experiment/20241113_80K_MPRAsnakeflow/results/experiments/mpra80KNeuronbbmapmapq30BC10DNA1RNA1SelfmadeStrandsensitivity/assigned_counts/assignmentFixDuplicates/default/NGN2_allreps_merged_barcode_assigned_counts.tsv.gz
# generate BCalm inputs elements:
# mpralib sequence-design get-counts --input  /home/kisa/coding/80K_MPRA/ismb/NGN2_neurons/reporter_experiment.barcode.NGN2.defaultAssignment.default.all.tsv.gz  --sequence-design /home/kisa/coding/80K_MPRA/80K-Analysis/07_quality_control/notebooks/control_metadata/MPRA_80215_April_server.tsv.gz --barcodes --all-oligos --normalized-counts --bc-threshold 1 --output  /home/kisa/coding/80K_MPRA/80K-Analysis/06_variant_element_analysis/scripts/test_mpralib_output/80k_NGN2_normalized_counts.tsv
# NGN2
mpralib sequence-design get-counts --input  /home/kisa/coding/80K_MPRA/ismb/NGN2_neurons/renamed_reporter_experiment.barcode.NGN2.defaultAssignment.default.all.tsv.gz  --sequence-design /home/kisa/coding/80K_MPRA/80K-Analysis/07_quality_control/notebooks/control_metadata/MPRA_80215_April_server.tsv.gz --barcodes --all-oligos --normalized-counts --bc-threshold 1 --output  /home/kisa/coding/80K_MPRA/80K-Analysis/06_variant_element_analysis/scripts/test_mpralib_output/80k_NGN2_normalized_counts.tsv

# undiff WTC11
mpralib sequence-design get-counts --input  /home/kisa/coding/80K_MPRA/ismb/undiff_WTC11/renamed_reporter_experiment.barcode.WTC11.defaultAssignment.default.all.tsv.gz  --sequence-design /home/kisa/coding/80K_MPRA/80K-Analysis/07_quality_control/notebooks/control_metadata/MPRA_80215_April_server.tsv.gz --barcodes --all-oligos --normalized-counts --bc-threshold 1 --output  /home/kisa/coding/80K_MPRA/80K-Analysis/06_variant_element_analysis/scripts/test_mpralib_output/80k_WTC11_normalized_counts.tsv


# running bcalm: (this script does not use BCalm (use /home/kisa/coding/80K_MPRA/80K-Analysis/06_variant_element_analysis/scripts/perform_bc_mpralm_elements_kisa.R))
Rscript /home/kisa/coding/80K_MPRA/80K-Analysis/06_variant_element_analysis/scripts/call_bcalm_elements.R \
--counts /home/kisa/coding/80K_MPRA/80K-Analysis/06_variant_element_analysis/scripts/test_mpralib_output/80k_NGN2_normalized_counts.tsv \
--labels /home/kisa/coding/80K_MPRA/80K-Analysis/06_variant_element_analysis/scripts/test_mpralib_output/80k_metadata_scramble_control_label.tsv \
--test-label cardiac_neuro_cava_random \
--control-label scrambled_control \
--percentile 0.975 \
--output /home/kisa/coding/80K_MPRA/80K-Analysis/06_variant_element_analysis/scripts/test_mpralib_output/NGN2_bcalm_tested_elements_vs_scrambled.tsv \
--output-vulcano-plot /home/kisa/coding/80K_MPRA/80K-Analysis/06_variant_element_analysis/scripts/test_mpralib_output/NGN2_bcalm_tested_elements_vs_scrambled_volcano.png \
--output-density-plot /home/kisa/coding/80K_MPRA/80K-Analysis/06_variant_element_analysis/scripts/test_mpralib_output/NGN2_bcalm_tested_elements_vs_scrambled_density.png \
--normalize FALSE


parser <- ArgumentParser(description = "Process BCALM element data")
parser$add_argument("--counts", type = "character", required = TRUE, help = "Path to the counts file")
parser$add_argument("--labels", type = "character", required = TRUE, help = "Path to the labels file")
parser$add_argument("--test-label", type = "character", required = TRUE, help = "Name of the test group")
parser$add_argument("--control-label", type = "character", required = TRUE, help = "Name of the control group")
parser$add_argument("--percentile",
    type = "double", default = 0.975,
    help = "Percentile of control to test on. Default is 0.975"
)
parser$add_argument("--output", type = "character", required = TRUE, help = "Path to the output file")
parser$add_argument("--output-vulcano-plot", type = "character", required = FALSE, help = "Path to store the vulcano plot")
parser$add_argument("--output-density-plot", type = "character", required = FALSE, help = "Path to store the density plot")
parser$add_argument("--normalize", type = "logical", default = TRUE, help = "Whether to normalize the data (TRUE or FALSE)")



# generate BCalm input variants
mpralib sequence-design get-variant-counts --input  <bc_file>--sequence-design <metadata file> \
--barcodes --normalized-counts --bc-threshold <mabe just 1?> --output  <output-file>

# optional: compute variant map
