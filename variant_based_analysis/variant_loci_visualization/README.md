# Variant loci visualization

Reproduces the per‑variant genome‑browser panels in the manuscript: for a given
rsID (or `chrom/pos/ref/alt`) it renders a stacked Gviz plot of the surrounding
locus.

## Output tracks (top → bottom)

1. **Genome axis** — coordinates, with the variant position marked by a vertical line
2. **Genes** — collapsed gene models (longest transcript per gene) from local GENCODE v42
3. **MPRA regions** — tested oligo regions (blue = contains an active sequence)
4. **MPRA variants** — tested variants (green = activating, orange = repressive, red border = variant of interest)
5. **SCREEN cCREs** — ENCODE Registry V4 cCREs (blue = brain‑active, gray = other)
6. **MetaBrain eQTL** — cis‑eQTL arcs from the variant to target‑gene TSS, brain tissue (blue)
7. **NGN2 ATAC** *(optional, `show_ngn2_atac`)* — WTC11‑NGN2 ExN ATAC‑seq peaks
8. **Brain ATAC** *(optional, `brain_atac`)* — GSE113480 bigwig signal for astrocyte / cortical / hippocampus / motor
9. **GTEx eQTL / EMS** *(optional, `show_gtex` / `show_ems`)* — disabled for the manuscript figures

## Environments

| Task | conda env |
|------|-----------|
| eQTL preprocessing (`scripts/preprocess_eqtl_data.py`) | `mpra80k_python` (pandas + pyarrow) |
| plotting (`scripts/plot_variants.R`) | `mpra80k_variant_plotting` (R 4.4, Gviz, GenomicInteractions, rtracklayer, biomaRt) |

## Input data

All commands are run **from this directory**
(`variant_based_analysis/variant_loci_visualization/`).

### Already in the repo / other analysis dirs
- `../../cCRE_based_analysis/data/DIV14_WTC11-NGN2_ExN_ATAC-seq.idr0.05.bfilt.narrowPeak`
- `../../cCRE_based_analysis/data/80K_MPRA_all_sequence_regions.bed`
- `../../cCRE_based_analysis/data/brain_cCRE_open_filtered.bed.gz`
- `../../modeling/data/2605_NGN2_variants_with_model_predictions_*.tsv.gz`
- `../data/gencode.v42.gtf.gz`
- `../data/ems_public/` (only needed if `show_ems = TRUE`)

### External downloads (not tracked in git — see `.gitignore`)

```bash
# 1. SCREEN / ENCODE Registry V4 human cCREs (GRCh38)  ~124 MB
curl -fL -o ../../cCRE_based_analysis/data/GRCh38-cCREs_screen_v4.bed \
  https://downloads.wenglab.org/V4/GRCh38-cCREs.bed

# 2. GSE113480 brain cell-type ATAC-seq bigwigs  (~4 GB total)
mkdir -p data/atac_peaks
GEO=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113480/suppl
for ct in astrocyte cortical hippocampal motor; do
  curl -fL -C - -o data/atac_peaks/GSE113480_${ct}.atac-seq.bigwig \
    $GEO/GSE113480_${ct}.atac-seq.bigwig
done
```

`scripts/download_external_data.sh` runs both steps (idempotent, resumable).

### MetaBrain eQTL table

`scripts/preprocess_eqtl_data.py` expects a per‑(variant × tissue × gene)
MetaBrain cis‑eQTL table restricted to the 80K MPRA variants, at
`../data/metabrain_eqtl_full_clean.parquet` (~240 MB, columns `SPDI`,
`SNPEffectAllele`, `GeneSymbol`, `GenePos`, `Tissue`, `MetaBeta`, `MetaP`, …).
It is produced by the MetaBrain / MPRA overlap pipeline
(`metabrain_eqtl_full_clean.parquet`). Significance is defined as `MetaP < 0.01`.

## Usage

```bash
# Step 1 — one-time preprocessing (writes data/eqtl/metabrain_processed.tsv.gz, ~36 MB, committed)
conda run -n mpra80k_python python scripts/preprocess_eqtl_data.py

# Step 2 — plot
conda run -n mpra80k_variant_plotting Rscript plot_examples.R      # all manuscript examples
```

or interactively in R (env `mpra80k_variant_plotting`):

```r
source("scripts/plot_variants.R")

plot_variant(
  rsid           = "rs17572795",
  brain_atac     = c("astrocyte", "cortical", "hippocampus", "motor"),
  show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE,
  arc_lwd        = 2.5
)
# → figures/rs17572795_all_sig_split_panel.{pdf,png}
```

Key `plot_variant()` parameters:

| Parameter | Default | Notes |
|-----------|---------|-------|
| `rsid` | — | resolved to GRCh38 via Ensembl REST (BioMart fallback); or pass `chrom/pos/ref/alt` |
| `window` | `200000` | half‑window in bp |
| `metabrain_dedup` | `"all_sig"` | `"all_sig"` \| `"cortex_EUR"` \| `"best_pvalue"` |
| `arc_mode` | `"split_panel"` | `"split_panel"` (±beta stacked) \| `"color_only"` |
| `brain_atac` | `NULL` | any subset of `c("astrocyte","cortical","hippocampus","motor")` |
| `show_ngn2_atac` / `show_ccres` | `TRUE` | |
| `show_gtex` / `show_ems` | `TRUE` | off in all manuscript example calls |
