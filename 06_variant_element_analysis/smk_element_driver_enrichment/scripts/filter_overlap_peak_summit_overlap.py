import pandas as pd
import argparse

# example usage: python filter_overlap_peak_summit_overlap.py --input_path all_SCREE_CREs_remap2022_overlap.bed.gz --output_path all_SCREE_CREs_remap2022_filtered_summit_overlap_no_border_overlap_min_10bp.bed.gz --writing
def filter_remap2022_CRE_overlap(remap2022_CRE_overlap_path, output_path="all_SCREE_CREs_remap2022_filtered_summit_overlap_no_border_overlap_min_10bp.bed.gz"):
    remap2022_CRE_overlap = pd.read_csv(remap2022_CRE_overlap_path, sep="\t", header=None, low_memory=False)
    remap2022_CRE_overlap.columns = ["chrom", "start", "end", "name", "score", "strand", "chrom_B", "start_B", "end_B", "name_remap", "score_B", "strand_B", "peak_start", "peak_end", "spectrum", "TF_name", "overlap_number"]


    remap2022_CRE_overlap["overlap_number"].hist()

    remap2022_CRE_overlap_filtered_df = remap2022_CRE_overlap[
        (remap2022_CRE_overlap['peak_start'] >= remap2022_CRE_overlap['start'] + 5) &
        (remap2022_CRE_overlap['peak_end'] <= remap2022_CRE_overlap['end'] - 5)  &
        (remap2022_CRE_overlap['overlap_number'] >= 10)
    ]

    print(f"Number of interactions after summit filtering: {len(remap2022_CRE_overlap_filtered_df)}")
    print(f"Proportion of interactions after summit filtering: {len(remap2022_CRE_overlap_filtered_df)/len(remap2022_CRE_overlap)}")
    remap2022_CRE_overlap_filtered_df[["chrom", "start", "end", "name", "score", "strand", "name_remap", "TF_name", "overlap_number"]].to_csv(output_path, sep="\t", index=False)
    print(f"Filtered file saved to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter remap2022 CRE overlap file based on summit overlap and minimum overlap number.")
    parser.add_argument("--input_path", type=str, help="Path to the remap2022 CRE overlap file (BED format).")
    parser.add_argument("--output_path", type=str, default="all_SCREE_CREs_remap2022_filtered_summit_overlap_no_border_overlap_min_10bp.bed.gz", help="Output path for the filtered BED file.")
    args = parser.parse_args()

    filter_remap2022_CRE_overlap(args.input_path, args.output_path)
