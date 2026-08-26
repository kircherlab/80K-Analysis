#!/usr/bin/env python3
import gzip
import pandas as pd
import pyranges as pr

# =============================
# INPUT FILES
# =============================
gtex_file = "data/gtex_prefiltered_by_mpra_1_based.tsv.gz"
sig_variants_file = "data/sig_variants.bed"

output_file = "data/GTEx_prefiltered_sig_variants_only.tsv.gz"

# =============================
# SETTINGS
# =============================
CHUNKSIZE = 1_000_000


# Define brain tissues: simple + transparent rule.
# You can expand later (e.g. add Spinal_cord, Pituitary, etc.)
def is_brain_tissue(series: pd.Series) -> pd.Series:
    # GTEx brain tissues typically start with "Brain_"
    return series.astype(str).str.startswith("Brain_")


# =============================
# LOAD SIGNIFICANT VARIANTS BED
# =============================
print("Loading significant MPRA variants...")

sig_vars = pd.read_csv(sig_variants_file, sep="\t", header=None, usecols=[0, 1, 2])
sig_vars.columns = ["Chromosome", "Start", "End"]

# Ensure types
sig_vars["Chromosome"] = sig_vars["Chromosome"].astype(str)
sig_vars["Start"] = sig_vars["Start"].astype(int)
sig_vars["End"] = sig_vars["End"].astype(int)

sig_pr = pr.PyRanges(sig_vars)

print(f"Significant variants loaded: {len(sig_vars):,}")

# =============================
# STREAM GTEx + FILTER
# =============================
print("Streaming GTEx file and filtering to significant variants...")

first_chunk = True
kept_rows_total = 0
seen_rows_total = 0

# write gz output incrementally
out_handle = gzip.open(output_file, "wt")

for chunk in pd.read_csv(gtex_file, sep="\t", chunksize=CHUNKSIZE):
    seen_rows_total += len(chunk)

    # Required columns for overlap: chrom + pos (hg38)
    # Your file has both "chromosome/start/end" (often legacy coords) and "chrom/pos" (hg38).
    # We'll use hg38: chunk["chrom"] and chunk["pos"].
    if "chrom" not in chunk.columns or "pos" not in chunk.columns:
        raise ValueError(
            "GTEx file must contain 'chrom' and 'pos' columns (hg38 position)."
        )

    # Clean coords
    chunk["chrom"] = chunk["chrom"].astype(str)
    # pos may be float if read strangely; coerce safely
    chunk["pos"] = pd.to_numeric(chunk["pos"], errors="coerce").astype("Int64")

    chunk = chunk.dropna(subset=["chrom", "pos"])
    chunk["pos"] = chunk["pos"].astype(int)

    # Create 1bp interval for overlap (BED-like: [pos-1, pos])
    var_df = pd.DataFrame(
        {
            "Chromosome": chunk["chrom"].values,
            "Start": (chunk["pos"].values - 1).astype(int),
            "End": chunk["pos"].values.astype(int),
        }
    )

    var_pr = pr.PyRanges(var_df)

    # Find overlaps to significant variants
    ov = var_pr.join(sig_pr)
    if len(ov) == 0:
        continue

    # Unique overlapping variant coordinates
    ov_coords = ov.df[["Chromosome", "Start", "End"]].drop_duplicates()

    # Merge back to chunk using the constructed coordinates
    chunk["_ov_chr"] = var_df["Chromosome"].values
    chunk["_ov_start"] = var_df["Start"].values
    chunk["_ov_end"] = var_df["End"].values

    filtered = chunk.merge(
        ov_coords,
        left_on=["_ov_chr", "_ov_start", "_ov_end"],
        right_on=["Chromosome", "Start", "End"],
        how="inner",
    )

    # Drop helper columns
    filtered = filtered.drop(
        columns=["_ov_chr", "_ov_start", "_ov_end", "Chromosome", "Start", "End"]
    )

    if len(filtered) == 0:
        continue

    # Add brain_flag based on tissue
    if "tissue" not in filtered.columns:
        raise ValueError("GTEx file must contain a 'tissue' column.")
    filtered["brain_flag"] = is_brain_tissue(filtered["tissue"])

    # Strip Ensembl version for easier joins later (ENSG... .15 -> ENSG...)
    if "gene" in filtered.columns:
        filtered["gene_id"] = (
            filtered["gene"].astype(str).str.replace(r"\..*$", "", regex=True)
        )
    else:
        filtered["gene_id"] = pd.NA

    kept_rows_total += len(filtered)

    # Write
    if first_chunk:
        filtered.to_csv(out_handle, sep="\t", index=False, header=True)
        first_chunk = False
    else:
        filtered.to_csv(out_handle, sep="\t", index=False, header=False)

    print(f"Processed {seen_rows_total:,} rows; kept {kept_rows_total:,} rows")

out_handle.close()
print("Done.")
print(f"Output written to: {output_file}")