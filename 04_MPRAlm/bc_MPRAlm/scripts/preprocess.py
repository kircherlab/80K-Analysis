## File to put data into correct format for a barcode level MPRAlm, also creates the input for the aggregated (normal) MPRAlm
## to make sure the same sequences are used in both. 

import polars as pl
import polars.selectors as cs

input = "/data/humangen_kircherlab/MPRA/tt_mpra_analysis/data/ProxProm/malacoda/malacoda_input_HepG2.tsv.gz"
output_rna = "/data/humangen_kircherlab/MPRA/tt_mpra_analysis/data/ProxProm/mpralm/mpralm_bc_rna_input_HepG2.tsv"
output_dna = "/data/humangen_kircherlab/MPRA/tt_mpra_analysis/data/ProxProm/mpralm/mpralm_bc_dna_input_HepG2.tsv"
output_both = "/data/humangen_kircherlab/MPRA/tt_mpra_analysis/data/ProxProm/mpralm/mpralm_input_HepG2.tsv"

df = pl.read_csv(input, separator="\t")

# filter for barcodes with a DNA count of at least 2
df = df.filter(pl.all_horizontal(pl.col("^DNA.*$") > 1))
df = df.filter(pl.all_horizontal(pl.col("^RNA.*$") > 1))
# filter for oligos with at least 10 barcodes in both alleles
df = df.filter(pl.count().over(["variant_id", "allele"]) >= 10)
df = df.filter(pl.n_unique("allele").over("variant_id") == 2)

# the max number of barcodes used per oligo is the number of barcodes of the 0.95th quantile of the sequences.
max_bc = df.group_by(["variant_id", "allele"]).count().select("count").quantile(0.95).item()

# downsample barcodes
df = df.filter(pl.int_range(0, pl.count()).shuffle().over(["variant_id", "allele"]) < max_bc)

# input for normal mpralm
df.write_csv(output_both, separator="\t")

# pivot dataframes so that each BC is a new column
df = df.with_columns(pl.int_range(1, pl.col("barcode").count()+1).over(["variant_id", "allele"]).alias("bc"))
df = df.with_columns(new_idx = pl.col("bc").cast(pl.Utf8) + "_" + pl.col("allele"))
df_rna = df.pivot(values=cs.contains("RNA"), index="variant_id", columns="new_idx", sort_columns=True)
df_dna = df.pivot(values=cs.contains("DNA"), index="variant_id", columns="new_idx", sort_columns=True)
df_rna = df_rna.select(
    pl.all().name.map(lambda col_name: col_name.replace('RNA', 'sample').replace('new_idx_', 'bc'))
).sort("variant_id")
df_dna = df_dna.select(
    pl.all().name.map(lambda col_name: col_name.replace('DNA', 'sample').replace('new_idx_', 'bc'))
).sort("variant_id")

df_rna.write_csv(output_rna, separator="\t")
df_dna.write_csv(output_dna, separator="\t")