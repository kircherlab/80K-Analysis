import pandas as pd
import argparse
import os

# read the tsv file from command line
parser = argparse.ArgumentParser(description='Process AlphaGenome TSV file.')
parser.add_argument(
    'alphaGenome_tsv_path',
    type=str,
    help='Path to the input AlphaGenome TSV file (e.g., ~/results/my_scores.tsv.gz)'
)
parser.add_argument(
    '-o', '--output_directory',
    type=str,
    default='.', # Default to current directory if not specified
    help='Directory to save output files. Creates if it does not exist.'
)

# Parse the arguments from the command line
args = parser.parse_args()

# Access the path provided by the user
alphaGenome_tsv_path = args.alphaGenome_tsv_path
output_directory = args.output_directory

# Path Expansion and Directory Creation
alphaGenome_tsv_path = os.path.expanduser(alphaGenome_tsv_path)
output_directory = os.path.expanduser(output_directory) # Also expand output directory path

if not os.path.exists(output_directory):
    os.makedirs(output_directory)
    print(f"Created output directory: {output_directory}")

alphaGenome_tsv = pd.read_csv(alphaGenome_tsv_path, sep="\t", low_memory=False)

# Ensure required columns exist before filtering
required_filter_columns = ['output_type', 'variant_scorer', 'variant_id', 'raw_score']
for col in required_filter_columns:
    if col not in alphaGenome_tsv.columns:
        raise ValueError(f"Input TSV is missing required column for filtering: '{col}'. "
                         "Please check your AlphaGenome output structure.")


atac_diff_log2_sum_df_chunk = alphaGenome_tsv[
            (alphaGenome_tsv['output_type'] == 'ATAC') &
            (alphaGenome_tsv["variant_scorer"] == "CenterMaskScorer(requested_output=ATAC, width=501, aggregation_type=DIFF_LOG2_SUM)")
        ]

expected_variants_in_chunk = 100
if atac_diff_log2_sum_df_chunk['variant_id'].nunique() != expected_variants_in_chunk:
    print("The number of unique variant IDs in the AlphaGenome TSV does not match the expected count of 100.")

# compute the mean absolute difference over the cell-types
aggregated_atac_diff_log2_alt_ref = atac_diff_log2_sum_df_chunk.groupby('variant_id')[
    'raw_score'
].apply(lambda x: x.abs().mean()).reset_index()
aggregated_atac_diff_log2_alt_ref.rename(
    columns={'raw_score': 'mean_abs_atac_diff_log2_alt_ref_alphaGenome'}, inplace=True)

# compute the max absolute difference over the cell-types
max_abs_atac_diff_log2_alt_ref = atac_diff_log2_sum_df_chunk.groupby('variant_id')[
    'raw_score'
].apply(lambda x: x.abs().max()).reset_index()
max_abs_atac_diff_log2_alt_ref.rename(
    columns={'raw_score': 'max_abs_atac_diff_log2_alt_ref_alphaGenome'}, inplace=True
)

# combine both
alphaGenome_chunk_abs_max_mean_atac_diff_pred = aggregated_atac_diff_log2_alt_ref.merge(
    max_abs_atac_diff_log2_alt_ref, on='variant_id', how='inner'
)

# Save the processed scores to a tsv file
output_name = alphaGenome_tsv_path.split("/")[-1].split(".")[0] + "_processed_scores.tsv.gz"
output_file = os.path.join(output_directory, output_name)
alphaGenome_chunk_abs_max_mean_atac_diff_pred.to_csv(output_file, sep="\t", index=False)
