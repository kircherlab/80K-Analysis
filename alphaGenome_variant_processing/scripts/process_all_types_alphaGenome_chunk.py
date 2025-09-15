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

# prediction_types = [
#     'ATAC', 'CAGE', 'CHIP_HISTONE', 'CHIP_TF', 'DNASE', 'PROCAP',
#     'RNA_SEQ', 'SPLICE_JUNCTIONS', 'SPLICE_SITES', 'SPLICE_SITE_USAGE'
# ]

prediction_types = list(alphaGenome_tsv['output_type'].unique())

# Final merged result
all_scores_df = None

for prediction_type in prediction_types:
    pred_df = alphaGenome_tsv[alphaGenome_tsv['output_type'] == prediction_type].copy()

    # Compute mean, median, max of abs(raw_score)
    grouped = pred_df.groupby('variant_id')['raw_score']

    stats_df = pd.DataFrame({
        'variant_id': grouped.apply(lambda x: x.name),
        f'{prediction_type}_mean_abs': grouped.apply(lambda x: x.abs().mean()),
        f'{prediction_type}_median_abs': grouped.apply(lambda x: x.abs().median()),
        f'{prediction_type}_max_abs': grouped.apply(lambda x: x.abs().max()),
        f'{prediction_type}_max_diff': grouped.apply(lambda x: x.max()),
    }).reset_index(drop=True)

    # Merge all prediction types on variant_id
    if all_scores_df is None:
        all_scores_df = stats_df
    else:
        all_scores_df = all_scores_df.merge(stats_df, on='variant_id', how='outer')


# Save the processed scores to a tsv file
output_name = alphaGenome_tsv_path.split("/")[-1].split(".")[0] + "_all_types_processed_scores.tsv.gz"
output_file = os.path.join(output_directory, output_name)
# Write to file
all_scores_df.to_csv(output_file, sep="\t", index=False)


