import pandas as pd
import matplotlib.pyplot as plt

remap2022_CRE_overlap_path = "./tested_CREs_remap2022.bed.gz"
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

# write the filtered remap overlap:
output_path = "tested_80k_MPRA_CREs_remap2022_filtered_summit_overlap_no_border_overlap_min_10bp.bed.gz"
writing = True
if writing:
    remap2022_CRE_overlap_filtered_df[["chrom", "start", "end", "name", "score", "strand", "name_remap", "TF_name", "overlap_number"]].to_csv(output_path, sep="\t", index=False)

