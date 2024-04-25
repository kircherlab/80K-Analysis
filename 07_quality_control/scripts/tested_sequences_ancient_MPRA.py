# This script shows the way to the filtered tested sequences for the ancient MPRA design of Thorben in 2024
"""

Filter oligo counts for all replicates of MPRAsnakeflow in oder to get x, y and z oligos of positive, negative and random sampled mean log2 rations 

:Author: Kilian Salomon
:Contact: kilian.salomon@bih-charite.de
:Date: *25.04.2024
:Type: tool
:Input: BAM
:Output: BAM

"""

# qulity measures
min_barcodes = 50
min_dna_count = 10
min_rna_count = 10

top_x = 100
bottom_y = 100
random_z = 200
seed = 123

# input
input_path = "/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/experiments/standard_bwa/assigned_counts/assignmentFixDuplicates/standardConfig/NGN2_allreps_minThreshold_merged.tsv.gz"

# output
output_path = "/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/07_quality_control/results/ancient_mpra/filtered_tested_sequences.tsv"


def combine_replicates(df_allreps, total_dna_counts, total_rna_counts):
    df_allreps = df_allreps.groupby(by=["condition", "name"]).aggregate(
        {
            "replicate": "count",
            "dna_counts": ["sum", "mean"],
            "rna_counts": ["sum", "mean"],
            "dna_normalized": "mean",
            "rna_normalized": "mean",
            "ratio": "mean",
            "log2": "mean",
            "n_obs_bc": ["sum", "mean"],
        }
    )
    
    df_allreps = df_allreps.reset_index()
    df_out = df_allreps.iloc[:, 0:2]
    df_out.columns = ["condition", "name"]

    df_out["replicates"] = df_allreps.replicate["count"]

    scaling = 10**6

    df_out["dna_counts"] = df_allreps.dna_counts["sum"]
    df_out["rna_counts"] = df_allreps.rna_counts["sum"]


    df_out["dna_normalized"] = df_out["dna_counts"] / total_dna_counts * scaling
    df_out["rna_normalized"] = df_out["rna_counts"] / total_rna_counts * scaling

    df_out["ratio"] = df_out["rna_normalized"] / df_out["dna_normalized"]
    df_out["log2"] = np.log2(df_out.ratio)

    df_out["mean_dna_counts"] = df_allreps.dna_counts["mean"]
    df_out["mean_rna_counts"] = df_allreps.rna_counts["mean"]
    df_out["mean_dna_normalized"] = df_allreps.dna_normalized["mean"]
    df_out["mean_rna_normalized"] = df_allreps.rna_normalized["mean"]
    df_out["mean_ratio"] = df_allreps.ratio["mean"]
    df_out["mean_log2"] = df_allreps.log2["mean"]
    
    df_out["mean_n_obs_bc"] = df_allreps.n_obs_bc["mean"].apply(int)
    return df_out


def add_log2_name(row):
    return row['name'] + "__mean_log2_ratio" + str(round(row['mean_log2'], 4))


merged_tsv_path = input_path
merged_tsv = pd.read_csv(merged_tsv_path, sep='\t')

# get total read counts for normalization
total_dna_counts = merged_tsv[["dna_counts"]].sum().iloc[0]
total_rna_counts = merged_tsv[["rna_counts"]].sum().iloc[0]

# only tested sequences
tested_merged_tsv = merged_tsv.loc[merged_tsv['name'].str.startswith('cardiac_neuro_cava_random')]
print('Number of sequences form the group: cardiac_neuro_cava_random: ', tested_merged_tsv.['name'].nunique())

# quality filter 


# barcode threshold
filtered_tested_merged_tsv = tested_merged_tsv.loc[tested_merged_tsv['n_obs_bc'] > min_barcodes]
print('barcodes: ', filtered_tested_merged_tsv['name'].nunique())

# min dna and rna count
filtered_tested_merged_tsv = filtered_tested_merged_tsv[filtered_tested_merged_tsv['dna_counts'] > min_dna_count]
filtered_tested_merged_tsv = filtered_tested_merged_tsv[filtered_tested_merged_tsv['rna_counts'] > min_rna_count]
print('dna rna counts: ', filtered_tested_merged_tsv['name'].nunique())

# occuring in all 3 replicates: count each name and keep only those with count 3
filtered_tested_merged_tsv['count'] = filtered_tested_merged_tsv.groupby('name')['name'].transform('count')
filtered_tested_merged_tsv = filtered_tested_merged_tsv.loc[filtered_tested_merged_tsv['count'] == 3]
print('replicates: ', filtered_tested_merged_tsv['name'].nunique())

# print shape 
print('After filtering: ', filtered_tested_merged_tsv['name'].nunique())



# merge the sequences over the replicates
combined_filtered_df = combine_replicates(filtered_tested_merged_tsv, total_dna_counts, total_rna_counts)

# add the mean_log2 to the header
combined_filtered_df['name_mean_log2'] = combined_filtered_df.apply(add_log2_name, axis = 1)

# drop unused columns: ('dna_counts', 'rna_counts', 'dna_normalized', 'rna_normalized', 'ratio') 
combined_filtered_final_df = combined_filtered_df.drop(columns=['dna_counts', 'rna_counts', 'dna_normalized', 'rna_normalized', 'ratio'])

# sort by log2
sorted_combined_filtered_df = combined_filtered_final_df.sort_values('mean_log2', ascending=False)

# make 3 sets of variants, then merge

top = sorted_combined_filtered_df.head(top_x)
random_inbetween = sorted_combined_filtered_df.head(-bottom_y).tail(-top_x).sample(n=random_z, random_state=123)
bottom = sorted_combined_filtered_df.tail(bottom_y)

tested_sequences_filtered_df = pd.concat([top, random_inbetween, bottom])

# write to output file
tested_sequences_filtered_df.to_csv(output_path, sep="\t", index=False)