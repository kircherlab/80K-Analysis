#!/usr/bin/env python3
"""
One-time EMS arc preprocessing for plot_variants_v3.R.

Reads every ems_top_<tissue>.tsv.bgz file, filters to variants present in the
MPRA library, resolves each gene to its canonical TSS (GENCODE v42 gene
feature), and writes a compressed BEDPE-style TSV in the same column schema
used by the GTEx/Metabrain preprocessed files.

Output
------
data/eqtl/ems_arcs_processed.tsv.gz
    Columns (tab-separated, v_chr/v_start in columns 1-2 for fast awk lookup):
      v_chr, v_start (0-based), v_end,
      g_chr, g_start (0-based), g_end,
      pip (= ems, calibrated causal probability), beta (= ems_normalized),
      gene_name, tissue, is_brain, is_sig

Usage
-----
  cd /path/to/variant-browser
  python src/preprocess_ems_data.py

Coordinate notes
----------------
  EMS variant IDs follow the GTEx convention ("chrN_POS_REF_ALT_b38") where POS
  is 1-based (VCF).  v_start is stored as POS - 1 (0-based), matching the MPRA
  pos0 column and the GTEx preprocessed file convention.

  Gene anchors use the GTF 'gene' feature (same as process_gtex).  The TSS is
  the GTF start for + strand genes and the GTF end for - strand genes,
  converted to 0-based: g_start = tss - 1, g_end = tss.
"""

import gzip
import re
import sys
from pathlib import Path

import pandas as pd

# Configuration
DATA_BASE = Path("./")
# NOTE: these paths need to be downloaded manually from the fine-mapped EMS release and the gencode v42 release
# https://www.finucanelab.org/data
EMS_DIR   = DATA_BASE / "80K_MPRA/eQTL_fine-mapped/EMS_public"
GTF_FILE  = DATA_BASE / "80K_MPRA/gencode_v42/gencode.v42.annotation.gtf.gz"
MPRA_FILE = DATA_BASE / (
    "../../modeling/data/2605_NGN2_variants_with_model_predictions_cadd16_cadd17_alphaGenome_encoder_chrombpnet.tsv.gz"
)

OUT_DIR  = Path("data/eqtl")
EMS_OUT  = OUT_DIR / "ems_arcs_processed.tsv.gz"

EMS_SIG_THRESHOLD = 0.001  # ems (calibrated causal probability) cutoff for is_sig


# GENCODE GTF parser (mirrors process_gtex in preprocess_eqtl_data.py)

def parse_gencode_tss(gtf_path: Path) -> dict:
    """
    Return dict: ensembl_base_id → {'chrom', 'tss', 'strand', 'gene_name'}.
    TSS = chromStart (1-based) for + strand genes, chromEnd for - strand.
    Uses 'gene' features (one entry per gene, no transcript ambiguity).
    """
    print(f"Parsing GENCODE GTF: {gtf_path}", flush=True)
    re_gene_id   = re.compile(r'gene_id "([^"]+)"')
    re_gene_name = re.compile(r'gene_name "([^"]+)"')

    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    tss: dict = {}

    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue

            chrom  = parts[0]
            start  = int(parts[3])   # 1-based GTF
            end    = int(parts[4])   # 1-based GTF
            strand = parts[6]
            attrs  = parts[8]

            m_id = re_gene_id.search(attrs)
            if not m_id:
                continue
            gene_base = m_id.group(1).split(".")[0]

            m_name    = re_gene_name.search(attrs)
            gene_name = m_name.group(1) if m_name else gene_base

            tss_pos = start if strand == "+" else end
            tss[gene_base] = {
                "chrom":     chrom,
                "tss":       tss_pos,
                "strand":    strand,
                "gene_name": gene_name,
            }

    print(f"  → TSS found for {len(tss):,} genes", flush=True)
    return tss


# MPRA variant filter

def load_mpra_positions(mpra_path: Path) -> pd.DataFrame:
    """
    Return DataFrame with columns [chrom_num, pos0] for all MPRA library variants.

    Source columns (R 1-indexed → Python 0-indexed):
      col 2 (idx 1) = chr_pos_ref_alt  e.g. "10-67064084-A-T"
      col 5 (idx 4) = pos0             0-based position
    """
    print(f"Loading MPRA variants: {mpra_path}", flush=True)

    df = pd.read_csv(mpra_path, sep="\t", header=0, compression="gzip",
                     low_memory=False)
    variant_col = df.columns[1]   # chr_pos_ref_alt
    pos_col     = df.columns[4]   # pos0

    out = pd.DataFrame({
        "chrom_num": df[variant_col].astype(str).str.split("-").str[0],
        "pos0":      pd.to_numeric(df[pos_col], errors="coerce"),
    }).dropna(subset=["pos0"])
    out["pos0"] = out["pos0"].astype(int)
    out = out.drop_duplicates()

    print(f"  → {len(out):,} unique (chrom, pos0) MPRA positions", flush=True)
    return out


#  EMS processing

def process_ems(ems_dir: Path, tss_dict: dict,
                mpra_pos: pd.DataFrame, out_path: Path) -> None:

    ems_files = sorted(ems_dir.glob("ems_top_*.tsv.bgz"))
    ems_files = [f for f in ems_files if "Zone.Identifier" not in str(f)]
    print(f"\nFound {len(ems_files)} EMS tissue files", flush=True)

    # Build TSS lookup table for fast merge
    tss_rows = [
        {"gene_id_base": gid,
         "g_chr":        info["chrom"],
         "g_start":      info["tss"] - 1,   # 0-based
         "g_end":        info["tss"],
         "gene_name":    info["gene_name"]}
        for gid, info in tss_dict.items()
    ]
    tss_df = pd.DataFrame(tss_rows)

    all_chunks: list[pd.DataFrame] = []
    total_in = 0
    total_out = 0

    for ems_path in ems_files:
        tissue   = re.sub(r"\.tsv\.bgz$", "",
                          re.sub(r"^ems_top_", "", ems_path.name))
        is_brain = tissue.startswith("Brain_")

        try:
            df = pd.read_csv(
                ems_path,
                sep="\t",
                header=0,
                names=["v", "g", "rf_score_raw", "ems", "ems_normalized"],
                dtype={"v": str, "g": str,
                       "rf_score_raw": float, "ems": float, "ems_normalized": float},
                compression="gzip",   # bgz is gzip-compatible
            )
        except Exception as exc:
            print(f"  [WARN] {ems_path.name}: {exc}", flush=True)
            continue

        if df.empty:
            continue

        total_in += len(df)

        # Parse variant ID: "chr10_100009635_T_G_b38"
        # POS in EMS IDs follows GTEx/VCF 1-based convention → convert to 0-based
        v_parts       = df["v"].str.split("_", n=4, expand=True)
        df["v_chr"]   = v_parts[0]                                     # "chr10"
        df["chrom_num"] = df["v_chr"].str.replace(r"^chr", "", regex=True)
        v_pos1        = pd.to_numeric(v_parts[1], errors="coerce")
        df            = df[v_pos1.notna()].copy()
        df["v_start"] = v_pos1[v_pos1.notna()].astype(int) - 1        # 0-based
        df["v_end"]   = df["v_start"] + 1

        #  Filter: keep only MPRA library variants
        df = df.merge(
            mpra_pos,
            left_on=["chrom_num", "v_start"],
            right_on=["chrom_num", "pos0"],
            how="inner",
        ).drop(columns=["pos0"])

        if df.empty:
            continue

        #  Join canonical TSS
        df["gene_id_base"] = df["g"].str.replace(r"\.\d+$", "", regex=True)
        df = df.merge(tss_df, on="gene_id_base", how="inner")

        if df.empty:
            continue

        df["pip"]      = df["ems"]               # calibrated causal probability
        df["beta"]     = df["ems_normalized"]    # enrichment over background prior
        df["is_sig"]   = df["ems"] > EMS_SIG_THRESHOLD
        df["is_brain"] = is_brain
        df["tissue"]   = tissue

        total_out += len(df)
        print(f"  {tissue}: {len(df):,} arcs "
              f"({df['is_sig'].sum()} sig)", flush=True)

        all_chunks.append(df[[
            "v_chr", "v_start", "v_end",
            "g_chr", "g_start", "g_end",
            "pip", "beta", "gene_name", "tissue", "is_brain", "is_sig",
        ]])

    if not all_chunks:
        print(
            "\nERROR: No EMS arcs matched MPRA variants.\n"
            "  • Verify EMS variant ID coordinate system (1-based assumed).\n"
            "  • Verify MPRA pos0 column index (0-indexed col 4 assumed).\n"
            "  • Check that EMS_DIR contains ems_top_*.tsv.bgz files.",
            file=sys.stderr,
        )
        sys.exit(1)

    out = pd.concat(all_chunks, ignore_index=True)
    out = out.sort_values(["v_chr", "v_start"]).reset_index(drop=True)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False, compression="gzip")

    print(f"\nWritten {len(out):,} arcs → {out_path}")
    print(f"  Input rows across all tissues : {total_in:,}")
    print(f"  After MPRA filter + TSS join  : {total_out:,}")
    n_var = out[["v_chr", "v_start"]].drop_duplicates().shape[0]
    n_gen = out["gene_name"].nunique()
    n_tis = out["tissue"].nunique()
    print(f"  Variants: {n_var} | Genes: {n_gen} | Tissues: {n_tis}")
    print(f"  Significant (ems > {EMS_SIG_THRESHOLD}): {out['is_sig'].sum():,}")
    print(f"  Brain-tissue arcs: {out['is_brain'].sum():,}")


#  Entry point

def main():
    for path in [EMS_DIR, GTF_FILE, MPRA_FILE]:
        if not Path(path).exists():
            print(f"ERROR: required path not found:\n  {path}", file=sys.stderr)
            sys.exit(1)

    tss_dict = parse_gencode_tss(GTF_FILE)
    mpra_pos = load_mpra_positions(MPRA_FILE)
    process_ems(EMS_DIR, tss_dict, mpra_pos, EMS_OUT)

    print("\nDone.")
    print(f"  {EMS_OUT}")


if __name__ == "__main__":
    main()
