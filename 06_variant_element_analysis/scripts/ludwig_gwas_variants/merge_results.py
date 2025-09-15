# merge_results.py

import pandas as pd
import os
import sys
import glob
import gzip

def merge_all_results(input_dir, output_merged_filepath):
    """
    Finds all processed GWAS overlap files in a directory and merges them
    into a single gzipped TSV file. It also checks for duplicate variant
    entries in the final merged dataset using 'SPDI_elem_chrom'.

    Args:
        input_dir (str): Directory containing the individual processed
                         _overlap_*.tsv.gz files.
        output_merged_filepath (str): Path to save the final merged gzipped TSV file.
    """
    print(f"Starting merge process from directory: {input_dir}")

    # Use glob to find all files matching the pattern
    # The pattern matches files like '22601_11232861__significant_gwas_overlap_123.tsv.gz'
    file_pattern = os.path.join(input_dir, "*_overlap_*.tsv.gz")
    processed_files = glob.glob(file_pattern)

    if not processed_files:
        print(f"No processed files found matching pattern '{file_pattern}'. Nothing to merge.")
        return

    print(f"Found {len(processed_files)} files to merge.")

    all_dataframes = []
    for i, filepath in enumerate(processed_files):
        try:
            # Read each gzipped TSV file
            df = pd.read_csv(filepath, sep='\t', compression='gzip', engine='c')
            all_dataframes.append(df)
            if (i + 1) % 100 == 0: # Print progress every 100 files
                print(f"Read {i+1}/{len(processed_files)} files...")
        except pd.errors.EmptyDataError:
            print(f"Warning: File {filepath} is empty or malformed and will be skipped.")
            continue
        except Exception as e:
            print(f"Error reading file {filepath}: {e}. Skipping this file.")
            continue

    if not all_dataframes:
        print("No valid dataframes to concatenate. Exiting.")
        return

    # Concatenate all dataframes into a single one
    merged_df = pd.concat(all_dataframes, ignore_index=True)
    print(f"Successfully concatenated all data. Total rows: {len(merged_df)}")

    # Check for duplicate variant entries after concatenation using 'SPDI_elem_chrom'
    if 'SPDI_elem_chrom' in merged_df.columns:
        num_unique_variants = merged_df['SPDI_elem_chrom'].nunique()
        if len(merged_df) != num_unique_variants:
            num_duplicate_entries = len(merged_df) - num_unique_variants
            print(f"Detected {num_duplicate_entries} Variants which occur more than once (variants identified by 'SPDI_elem_chrom') in the merged dataset. All instances are kept.")
            print("Number of unique variants with a GWAS associated hit:", num_unique_variants)

        else:
            print("No duplicate variant entries found in the merged dataset.")
    else:
        print("Warning: 'SPDI_elem_chrom' column not found in the merged dataframe, cannot check for duplicate variants.")

    # Ensure the output directory for the merged file exists, handling cases where only a filename is given
    output_dir = os.path.dirname(output_merged_filepath)
    if output_dir: # Only create directory if output_merged_filepath contains a directory component
        os.makedirs(output_dir, exist_ok=True)

    # Save the final merged DataFrame to a gzipped TSV file
    merged_df.to_csv(output_merged_filepath, sep='\t', compression='gzip', index=False)
    print(f"Merged results saved to: {output_merged_filepath}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python merge_results.py <input_directory> <output_merged_filepath>")
        sys.exit(1)

    input_directory = sys.argv[1]
    output_file = sys.argv[2]

    merge_all_results(input_directory, output_file)
