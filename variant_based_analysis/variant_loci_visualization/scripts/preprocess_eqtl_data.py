#!/usr/bin/env python3
"""
One-time eQTL preprocessing for plot_variants_v3.R.

Outputs
-------
data/eqtl/gtex_processed.tsv.gz
    GTEx per-(variant, gene, is_brain) with TSS coords, brain flag, beta_marginal.

data/eqtl/metabrain_processed.tsv.gz
    Metabrain with alt-allele-corrected beta.
    Column 'dedup_key' marks which dedup mode a row participates in:
      "cortex_EUR"  — row comes from Tissue == cortex_EUR
      "best_pvalue" — row has lowest MetaP for that (chrom, pos, GeneSymbol) pair
      Both flags can be set simultaneously (a cortex_EUR row may also be best p-value).

Usage
-----
  cd /path/to/variant-browser
  python src/preprocess_eqtl_data.py
"""

import gzip
import re
import sys
from pathlib import Path

import pandas as pd

#  Configuration
DATA_BASE = Path("./")

GENCODE_GTF  = DATA_BASE / "../../data/gencode.v42.gtf.gz"
GTEX_RAW     = DATA_BASE / "../../data/gtex_prefiltered_by_mpra_1_based.tsv.gz"
# NOTE: This file is too large to store in the repo
META_PARQUET = DATA_BASE / "neuro_Metabrain_eQTL_overlap.parquet"

OUT_DIR  = Path("data/eqtl")
GTEX_OUT = OUT_DIR / "gtex_processed.tsv.gz"
META_OUT = OUT_DIR / "metabrain_processed.tsv.gz"

NC_TO_UCSC = {
    "NC_000001.11": "chr1",  "NC_000002.12": "chr2",  "NC_000003.12": "chr3",
    "NC_000004.12": "chr4",  "NC_000005.10": "chr5",  "NC_000006.12": "chr6",
    "NC_000007.14": "chr7",  "NC_000008.11": "chr8",  "NC_000009.12": "chr9",
    "NC_000010.11": "chr10", "NC_000011.10": "chr11", "NC_000012.12": "chr12",
    "NC_000013.11": "chr13", "NC_000014.9" : "chr14", "NC_000015.10": "chr15",
    "NC_000016.10": "chr16", "NC_000017.11": "chr17", "NC_000018.10": "chr18",
    "NC_000019.10": "chr19", "NC_000020.11": "chr20", "NC_000021.9" : "chr21",
    "NC_000022.11": "chr22", "NC_000023.11": "chrX",  "NC_000024.10": "chrY",
}


#  GENCODE GTF parser

def parse_gencode_tss(gtf_path: Path) -> dict:
    """
    Return dict: ensembl_base_id → {'chrom', 'tss', 'strand', 'gene_name'}.
    TSS = chromStart (1-based) for + strand genes, chromEnd for - strand.
    """
    print(f"Parsing GENCODE GTF: {gtf_path}")
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
            start  = int(parts[3])   # 1-based
            end    = int(parts[4])   # 1-based
            strand = parts[6]
            attrs  = parts[8]

            m_id = re_gene_id.search(attrs)
            if not m_id:
                continue
            gene_base = m_id.group(1).split(".")[0]

            m_name = re_gene_name.search(attrs)
            gene_name = m_name.group(1) if m_name else gene_base

            tss_pos = start if strand == "+" else end
            tss[gene_base] = {
                "chrom":     chrom,
                "tss":       tss_pos,
                "strand":    strand,
                "gene_name": gene_name,
            }

    print(f"  → TSS found for {len(tss):,} genes")
    return tss


#  GTEx processing

def process_gtex(raw_path: Path, tss_dict: dict, out_path: Path) -> None:
    """
    Load raw GTEx filtered file, add TSS positions and brain flag, deduplicate
    to one row per (variant, gene, is_brain) keeping max-pip row, write BEDPE-style TSV.
    """
    print(f"\nProcessing GTEx: {raw_path}")
    df = pd.read_csv(raw_path, sep="\t", compression="gzip", low_memory=False)
    print(f"  Loaded {len(df):,} rows | {df['var_key'].nunique():,} unique variants")

    # Strip Ensembl version
    df["gene_base"] = df["gene"].str.replace(r"\.\d+$", "", regex=True)

    # Brain tissue flag
    df["is_brain"] = df["tissue"].str.contains("Brain", case=False, na=False)

    # Deduplicate: per (var_key, gene_base, is_brain) keep max-pip row
    idx = df.groupby(["var_key", "gene_base", "is_brain"])["pip"].idxmax()
    df  = df.loc[idx].copy().reset_index(drop=True)
    print(f"  After dedup: {len(df):,} rows")

    # Look up TSS from GENCODE
    df["g_chr_raw"] = df["gene_base"].map(lambda g: tss_dict.get(g, {}).get("chrom"))
    df["g_tss"]     = df["gene_base"].map(lambda g: tss_dict.get(g, {}).get("tss"))
    df["gene_name"] = df["gene_base"].map(lambda g: tss_dict.get(g, {}).get("gene_name", g))

    n_before = len(df)
    df = df.dropna(subset=["g_tss", "g_chr_raw"]).copy()
    print(f"  Dropped {n_before - len(df):,} rows without TSS match")

    df["g_tss"] = df["g_tss"].astype(int)

    # Variant anchor: 0-based half-open (pos is 1-based in raw file → convert)
    df["v_chr"]   = df["chrom"]
    df["v_start"] = df["pos"].astype(int) - 1   # 0-based
    df["v_end"]   = df["pos"].astype(int)

    # TSS anchor: convert GTF 1-based to 0-based half-open
    df["g_chr"]   = df["g_chr_raw"]
    df["g_start"] = df["g_tss"] - 1             # 0-based
    df["g_end"]   = df["g_tss"]

    out = df[[
        "v_chr", "v_start", "v_end",
        "g_chr", "g_start", "g_end",
        "pip", "beta_marginal", "gene_name", "tissue", "is_brain", "var_key",
    ]].rename(columns={"beta_marginal": "beta"})

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False, compression="gzip")
    print(f"  Written {len(out):,} rows → {out_path}")
    print(f"  Brain rows: {out['is_brain'].sum():,} | Non-brain: {(~out['is_brain']).sum():,}")


#  Metabrain processing

def process_metabrain(parquet_path: Path, out_path: Path) -> None:
    """
    Load Metabrain parquet, correct MetaBeta to reflect alt-allele effect,
    add UCSC chrom, and annotate each row with which dedup modes it participates in.
    """
    print(f"\nProcessing Metabrain: {parquet_path}")
    df = pd.read_parquet(parquet_path)
    print(f"  Loaded {len(df):,} rows")

    # Parse SPDI: "NC_000010.11:372603:G:A" → nc_acc, pos(0-based), ref, alt
    spdi_parts  = df["SPDI"].str.split(":", expand=True)
    df["nc_acc"] = spdi_parts[0]
    df["pos0"]   = spdi_parts[1].astype(int)   # 0-based
    df["ref"]    = spdi_parts[2]
    df["alt"]    = spdi_parts[3]

    # Map NC accession → UCSC chromosome
    df["chrom"] = df["nc_acc"].map(NC_TO_UCSC)
    n_unknown = df["chrom"].isna().sum()
    if n_unknown:
        print(f"  Warning: {n_unknown:,} rows with unknown NC accession (dropped)")
    df = df.dropna(subset=["chrom"]).copy()

    # Correct MetaBeta: if SNPEffectAllele is the REF allele, flip sign so alt = effect allele
    flip_mask             = df["SNPEffectAllele"] == df["ref"]
    df["MetaBeta_alt"]    = df["MetaBeta"].copy()
    df.loc[flip_mask, "MetaBeta_alt"] = -df.loc[flip_mask, "MetaBeta"]
    print(f"  Flipped beta for {flip_mask.sum():,} rows (SNPEffectAllele == ref)")

    # Dedup mode 1: cortex_EUR
    df["is_cortex_EUR"] = df["Tissue"] == "cortex_EUR"

    # Dedup mode 2: best_pvalue — mark which rows have the lowest MetaP per (chrom, pos0, GeneSymbol)
    df["_rank_p"] = df.groupby(["chrom", "pos0", "GeneSymbol"])["MetaP"].rank(
        method="first", ascending=True
    )
    df["is_best_pvalue"] = df["_rank_p"] == 1
    df = df.drop(columns=["_rank_p"])

    # BEDPE-style anchors (0-based)
    df["v_start"] = df["pos0"]
    df["v_end"]   = df["pos0"] + 1
    df["g_start"] = df["GenePos"]          # already 0-based per README
    df["g_end"]   = df["GenePos"] + 1

    out_cols = [
        "chrom", "v_start", "v_end",
        "chrom", "g_start", "g_end",       # both anchors on same chrom (cis-eQTL)
        "MetaBeta_alt", "MetaP", "GeneSymbol", "Tissue",
        "sig_eQTL", "is_cortex_EUR", "is_best_pvalue",
        "pos0", "GenePos",
    ]
    # Rename for clarity in R
    out = df[[
        "chrom", "v_start", "v_end",
        "g_start", "g_end",
        "MetaBeta_alt", "MetaP", "GeneSymbol", "Tissue",
        "sig_eQTL", "is_cortex_EUR", "is_best_pvalue",
        "pos0", "GenePos",
    ]].rename(columns={
        "chrom":        "v_chr",
        "MetaBeta_alt": "beta",
        "MetaP":        "p_val",
        "GeneSymbol":   "gene_name",
        "sig_eQTL":     "is_sig",
        "pos0":         "v_pos0",
        "GenePos":      "g_pos0",
    })
    out["g_chr"] = out["v_chr"]  # cis-eQTL: same chromosome

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False, compression="gzip")
    print(f"  Written {len(out):,} rows → {out_path}")
    print(f"  cortex_EUR rows: {out['is_cortex_EUR'].sum():,}")
    print(f"  best_pvalue rows: {out['is_best_pvalue'].sum():,}")
    print(f"  significant (sig_eQTL): {out['is_sig'].sum():,}")


#  Entry point

def main():
    for path in [GENCODE_GTF, GTEX_RAW, META_PARQUET]:
        if not path.exists():
            print(f"ERROR: required file not found:\n  {path}", file=sys.stderr)
            sys.exit(1)

    tss_dict = parse_gencode_tss(GENCODE_GTF)
    process_gtex(GTEX_RAW, tss_dict, GTEX_OUT)
    process_metabrain(META_PARQUET, META_OUT)

    print("\nDone. Preprocessed files:")
    print(f"  {GTEX_OUT}")
    print(f"  {META_OUT}")


if __name__ == "__main__":
    main()
