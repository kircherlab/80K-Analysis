import pandas as pd
import polars as pl
import numpy as np
import polars.selectors as cs

# read a polars dataframe
barcode_counts_path = "/data/cephfs-2/unmirrored/groups/ag-kircher/MPRA/IGVF_Y1_design/experiment/final_resequencing/results/experiments/standard_bwa/assigned_counts/assignmentFixDuplicates/standardConfig/NGN2_allreps_merged_barcode_assigned_counts.tsv.gz"
barcode_counts_path = "/data/cephfs-2/unmirrored/groups/ag-kircher/MPRA/IGVF_Y1_design/experiment/resequencing_results/new_bam_filtering/results/experiments/standard_bwa/assigned_counts/assignmentFixDuplicates/standardConfig/NGN2_allreps_merged_barcode_assigned_counts.tsv.gz"
df = pl.read_csv(barcode_counts_path, has_header=True, separator="\t")
# Simply summing all the barcodes per oligo:
df_summed = df.group_by("name").agg(cs.contains("count").sum())
print(df_summed.head())
# Helper functions:

def calculate_logrs(df):
	df_norm = df.with_columns(cs.contains("count") / cs.contains("count").quantile(0.75)  * 1e6)
	df_logrs = df_norm.with_columns(logr1 = np.log2(pl.col("rna_count_replicate_1")/pl.col("dna_count_replicate_1")),
									logr2 = np.log2(pl.col("rna_count_replicate_2")/pl.col("dna_count_replicate_2")),
									logr3 = np.log2(pl.col("rna_count_replicate_3")/pl.col("dna_count_replicate_3"))
									)
	return df_logrs

def print_logr_correlation(df_logrs):
	print("corr btwn rep 1 and 3: "+str(df_logrs.fill_nan(None).drop_nulls().select(pl.corr("logr1", "logr3")).item()))
	print("corr btwn rep 1 and 2: "+str(df_logrs.fill_nan(None).drop_nulls().select(pl.corr("logr1", "logr2")).item()))
	print("corr btwn rep 2 and 3: "+str(df_logrs.fill_nan(None).drop_nulls().select(pl.corr("logr2", "logr3")).item()))

# print before filtering
print_logr_correlation(calculate_logrs(df_summed))


# Filter barcodes with different effect directions:
same_sign = df.filter((np.sign(pl.col("rna_count_replicate_1") - pl.col("dna_count_replicate_1")) ==
							  np.sign(pl.col("rna_count_replicate_2") - pl.col("dna_count_replicate_2"))) & (
								   np.sign(pl.col("rna_count_replicate_1") - pl.col("dna_count_replicate_1")) ==
							  np.sign(pl.col("rna_count_replicate_3") - pl.col("dna_count_replicate_3"))))
same_sign_summed = same_sign.group_by("name").agg(cs.contains("count").sum(), bcs=pl.len())
df_logrs = calculate_logrs(same_sign_summed.filter(pl.col("bcs") >= 5))
print("Percentage of oligos after removing barcodes with different effect directions: " + str(df_logrs.height/df_summed.height))
print_logr_correlation(df_logrs)

