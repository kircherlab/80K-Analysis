# process_gwas.py

import pandas as pd
import os
import sys
import gzip

def process_single_gwas_file(gwas_filepath, variant_info_filepath, output_dir):
    """
    Processes a single gzipped GWAS result file, merges it with variant information,
    filters for overlapping variants, adds a source file column, and saves the result.

    Args:
        gwas_filepath (str): Path to the gzipped GWAS results file.
        variant_info_filepath (str): Path to the variant information TSV file.
        output_dir (str): Directory to save the processed output files.
    """
    try:
        # Extract the base name of the GWAS file without extension for output naming
        gwas_filename = os.path.basename(gwas_filepath)
        # Remove all extensions (.tsv.gz) to get the clean base name
        base_filename_no_ext = gwas_filename.replace('.tsv.gz', '')

        print(f"Processing file: {gwas_filepath}")

        # Load the GWAS results file
        # Using 'engine="c"' for faster reading of TSV files
        gwas_results_df = pd.read_csv(gwas_filepath, sep='\t', compression='gzip', engine='c')
        print(f"Loaded GWAS file: {gwas_filepath} with {len(gwas_results_df)} rows.")

        # Load the variant information file
        variant_info_df = pd.read_csv(variant_info_filepath, sep='\t', engine='c')
        print(f"Loaded variant info file: {variant_info_filepath} with {len(variant_info_df)} rows.")

        # Perform the merge operation as specified by the user
        # Assumption: 'SPDI_elem_chrom' in variant_info_df and 'variant' in coffee_consumed_df
        # contain comparable identifiers for merging. If 'SPDI_elem_chrom' is just a chromosome,
        # this merge will likely not yield specific variant overlaps.
        variant_info_filtered_spdi_tested_coffee_consumed = variant_info_df.merge(
            gwas_results_df[['variant', 'minor_AF', 'tstat', 'pval']],
            left_on='SPDI_elem_chrom',  # Key from variant_info_df
            right_on='variant',         # Key from coffee_consumed_df
            how='left'
        ).copy()
        print(f"Merged data. Resulting shape: {variant_info_filtered_spdi_tested_coffee_consumed.shape}")

        # Filter out rows where 'pval' is NaN (i.e., no overlap found)
        variant_info_filtered_spdi_tested_coffee_consumed_annotated = \
            variant_info_filtered_spdi_tested_coffee_consumed.loc[
                ~variant_info_filtered_spdi_tested_coffee_consumed['pval'].isna()
            ].copy()
        print(f"Filtered for overlapping variants. Remaining rows: {len(variant_info_filtered_spdi_tested_coffee_consumed_annotated)}")

        # Check if the resulting table is not empty
        if not variant_info_filtered_spdi_tested_coffee_consumed_annotated.empty:
            # Add the 'source_file' column
            variant_info_filtered_spdi_tested_coffee_consumed_annotated['source_file'] = gwas_filename

            # Define the output filename based on the input file and number of rows
            num_rows = len(variant_info_filtered_spdi_tested_coffee_consumed_annotated)
            output_filename = f"{base_filename_no_ext}_overlap_{num_rows}.tsv.gz"
            output_filepath = os.path.join(output_dir, output_filename)

            # Ensure the output directory exists
            os.makedirs(output_dir, exist_ok=True)

            # Write the annotated DataFrame to a gzipped TSV file
            variant_info_filtered_spdi_tested_coffee_consumed_annotated.to_csv(
                output_filepath,
                sep='\t',
                compression='gzip',
                index=False
            )
            print(f"Successfully saved {num_rows} overlapping variants to: {output_filepath}")
        else:
            print(f"No overlapping variants found for {gwas_filepath}. No output file generated.")

    except FileNotFoundError:
        print(f"Error: One of the input files not found. GWAS: {gwas_filepath}, Variant Info: {variant_info_filepath}")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: GWAS file {gwas_filepath} is empty or malformed.")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while processing {gwas_filepath}: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python process_gwas.py <gwas_filepath> <variant_info_filepath> <output_directory>")
        sys.exit(1)

    gwas_filepath = sys.argv[1]
    variant_info_filepath = sys.argv[2]
    output_dir = sys.argv[3]

    process_single_gwas_file(gwas_filepath, variant_info_filepath, output_dir)
