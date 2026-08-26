import pandas as pd
import numpy as np
import re


# NOTE: the fine-mapped GTEx data for 49 tissues was downloaded from 25.05.2025 https://www.dropbox.com/scl/fo/bjp6o8hgixt5occ6ggq2o/ADTTMnyH-4rzDFS6g7oIcEE?rlkey=7558jp42yvmyjhgmcilbhlnu7&e=1&dl=0
eqtl_path = "data/GTEx_49tissues_release1.tsv.bgz"

def parse_variant_hg38(v):
    """
    Returns (chrom, pos1, ref, alt) from strings like:
        chr1_13550_G_A_b38
        chr1:13550:G:A
        chr1_13550_G_A
    Safe against whitespace and trailing characters.
    """
    if pd.isna(v):
        return None

    v = str(v).strip().replace(":", "_")
    v = re.sub(r"_b38$", "", v)

    parts = v.split("_")
    if len(parts) < 4:
        return None

    chrom, pos_str, ref, alt = parts[:4]

    try:
        pos = int(pos_str)
    except ValueError:
        return None

    return chrom, pos, ref, alt


# -----------------------------
# Load GTEx
# -----------------------------
eqtl_df = pd.read_csv(eqtl_path, sep="\t", compression="gzip")

# Parse variant_hg38 safely
parsed = eqtl_df["variant_hg38"].map(parse_variant_hg38)

# Remove None entries
mask_ok = parsed.notna()
eqtl_df = eqtl_df.loc[mask_ok].copy()

# Expand tuple into DataFrame
parsed_df = pd.DataFrame(
    parsed[mask_ok].tolist(), columns=["chrom", "pos", "ref", "alt"]
)

eqtl_df = pd.concat(
    [eqtl_df.reset_index(drop=True), parsed_df.reset_index(drop=True)], axis=1
)

# Clean spacing
eqtl_df["chrom"] = eqtl_df["chrom"].astype(str).str.strip()
eqtl_df["ref"] = eqtl_df["ref"].astype(str).str.strip()
eqtl_df["alt"] = eqtl_df["alt"].astype(str).str.strip()

# Build var_key
eqtl_df["var_key"] = (
    eqtl_df["chrom"]
    + ":"
    + eqtl_df["pos"].astype(int).astype(str)
    + ":"
    + eqtl_df["ref"]
    + ":"
    + eqtl_df["alt"]
)

# -----------------------------
# Write output
# -----------------------------
output_path = eqtl_path.replace(".tsv.bgz", "_with_varKey.tsv.gz")
eqtl_df.to_csv(output_path, sep="\t", index=False, compression="gzip")

print("Finished building GTEx var_key:", output_path)
print("Total variants parsed:", len(eqtl_df))