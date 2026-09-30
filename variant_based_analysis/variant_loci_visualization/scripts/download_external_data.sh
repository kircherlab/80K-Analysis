#!/usr/bin/env bash
#
# Download the external tracks required by scripts/plot_variants.R that are too
# large to store in git. Idempotent and resumable — safe to re-run.
#
# Run from the variant_loci_visualization/ directory:
#   bash scripts/download_external_data.sh
#
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ATAC_DIR="$HERE/data/atac_peaks"
CCRE_BED="$HERE/../../cCRE_based_analysis/data/GRCh38-cCREs_screen_v4.bed"
mkdir -p "$ATAC_DIR" "$(dirname "$CCRE_BED")"

log() { echo "[$(date +%H:%M:%S)] $*"; }

dl() {  # url dest
  local url="$1" dest="$2"
  if [ -s "$dest" ]; then
    log "exists, skip: ${dest##*/} ($(du -h "$dest" | cut -f1))"
    return 0
  fi
  log "downloading ${url##*/}"
  curl -fL --retry 5 --retry-delay 5 -C - -o "$dest.part" "$url"
  mv "$dest.part" "$dest"
  log "done: ${dest##*/} ($(du -h "$dest" | cut -f1))"
}

# 1. SCREEN / ENCODE Registry V4 human cCREs (GRCh38), ~124 MB
dl "https://downloads.wenglab.org/V4/GRCh38-cCREs.bed" "$CCRE_BED"

# 2. GSE113480 brain cell-type ATAC-seq bigwigs, ~4 GB total
GEO="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113480/suppl"
for ct in astrocyte cortical hippocampal motor; do
  dl "$GEO/GSE113480_${ct}.atac-seq.bigwig" "$ATAC_DIR/GSE113480_${ct}.atac-seq.bigwig"
done

log "all external data present"
