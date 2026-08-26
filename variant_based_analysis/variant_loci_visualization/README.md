## Directory to reproduce the variant loci visualization figures in the manuscript.
- Download the ATAC peaks from GEO: GSE113480
- point the required data paths to the correct location:

### Input

A single variant defined by:

| Parameter | Description | Example |
|-----------|-------------|---------|
| `chrom` | Chromosome (UCSC style) | `chr5` |
| `pos` | 0-based position | `14408058` |
| `ref` | Reference allele | `A` |
| `alt` | Alternate allele | `G` |
| `window` | Half-window in bp (optional, default 100 kb) | `100000` |

### Output tracks (top → bottom)

1. **Genome axis** — coordinates + variant position highlighted
2. **GWAS** — –log10(p) signal track
3. **eQTL** — fine-mapped eQTL arcs (GTEx v8 + Metabrain), colored by tissue (brain highlighted)
4. **EMS** — Expression Modifier Score signal track
5. **ATAC** — chromatin accessibility peaks (NGN2 ATAC)
6. **Gene models** — from cached TxDb, transcript names shown

### Example usage:

- Step 1 — one time
`python scripts/preprocess_eqtl_data.py`
`python scripts/preprocess_ems_data.py`
`python scripts/preprocess_variant_table.py`

- Step 2 — plot
```R
source("scripts/plot_variants.R")

plot_variant(rsid = 'rs1128287',    brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)
```

- NOTE: you might need to download and link different data sources like eQTL and EMS yourself so they can be preprocessed based on the variant data