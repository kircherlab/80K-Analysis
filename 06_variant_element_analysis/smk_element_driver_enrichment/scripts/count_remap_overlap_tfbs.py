
import pandas as pd

def count_tfbs_overlaps(filtered_overlap_path):
	"""
	Reads a filtered overlap BED file and returns a DataFrame:
	- Rows: sequences (chrom, start, end, name)
	- Columns: TF_name
	- Values: number of overlaps of each TF with each sequence
	"""
	df = pd.read_csv(filtered_overlap_path, sep="\t", header=None, low_memory=False)
	# Assumes columns: chrom, start, end, name, ... TF_name, overlap_number
	# Find TF_name and overlap_number columns
	colnames = ["chrom", "start", "end", "name", "score", "strand", "chrom_B", "start_B", "end_B", "name_remap", "score_B", "strand_B", "peak_start", "peak_end", "spectrum", "TF_name", "overlap_number"]
	# Find TF_name column index (last or second last)
	tf_col = None
	for i in range(len(df.columns)):
		if df.iloc[:,i].astype(str).str.contains("TF", case=False).any():
			tf_col = i
	if tf_col is None:
		tf_col = -2  # fallback
	overlap_col = -1
	# Add column names
	df.columns = colnames + [f"col_{i}" for i in range(len(df.columns)-4)]
	df = df.rename(columns={df.columns[tf_col]: "TF_name", df.columns[overlap_col]: "overlap_number"})
	# Group by sequence and TF_name, count overlaps
	overlap_table = df.groupby(["chrom", "start", "end", "name", "TF_name"]).size().reset_index(name="overlap_count")
	# Pivot to wide format: sequences as rows, TFs as columns
	big_table = overlap_table.pivot_table(index=["chrom", "start", "end", "name"], columns="TF_name", values="overlap_count", fill_value=0)
	return big_table

def tfbs_contingency_table(filtered_overlap_path, output_path):
	"""
	For each TF, count:
	- n_overlap: number of sequences with at least one overlap
	- n_nonoverlap: number of sequences with zero overlap
	Save as DataFrame to output_path.
	"""
	big_table = count_tfbs_overlaps(filtered_overlap_path)
	result = []
	n_sequences = big_table.shape[0]
	for tf in big_table.columns:
		n_overlap = (big_table[tf] > 0).sum()
		n_nonoverlap = n_sequences - n_overlap
		result.append({"TFBS_name": tf, "n_overlap": n_overlap, "n_nonoverlap": n_nonoverlap})
	df_out = pd.DataFrame(result)
	df_out.to_csv(output_path, sep="\t", index=False)
	return df_out

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Count TFBS overlaps and generate contingency table.")
	parser.add_argument("filtered_overlap_path", type=str, help="Path to filtered overlap BED file.")
	parser.add_argument("output_path", type=str, help="Path to save contingency table TSV.")
	args = parser.parse_args()
	tfbs_contingency_table(args.filtered_overlap_path, args.output_path)

