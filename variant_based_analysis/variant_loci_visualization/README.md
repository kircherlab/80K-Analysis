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
`python src/preprocess_eqtl_data.py`
`python src/preprocess_ems_data.py`
`python src/preprocess_variant_table.py`

- Step 2 — plot
```R
source("src/plot_variants_v3.R")

plot_variant(rsid = 'rs1128287',    brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)
```


#### Data
- MPRA variant effects           /home/kisa/coding/80K_MPRA/80K-Analysis/09_mpra_manuscript/data/variant_effects/ 2603_NGN2_variants_with_model_predictions_cadd16_cadd17_alphaGenome_encoder.tsv.gz
- MPRA tested regions            /home/kisa/coding/80K_MPRA/WTC11_ATAC_Ahituv/80K_MPRA_all_sequence_regions.bed
- SCREEN cCREs (all)             /home/kisa/coding/80K_MPRA/cCREs_overlap/newest_screen/GRCh38-cCREs_screen_v4.bed
- SCREEN cCREs (brain subset)    /home/kisa/coding/80K_MPRA/cCREs_overlap/SCREEN_tissue/enrichment_over_tissues/brain_cCRE_open_filtered.bed.gz
- Brain ATAC bigwigs  data/atac_peaks/GSE113480_{astrocyte,cortical,hippocampal,motor}.atac-seq.bigwig
- Metabrain eQTL source /home/kisa/coding/IAWG_IGVF/data/neuron_variant/2604_metabrain_eQTL/neuro_Metabrain_eQTL_overlap.parquet
- gencode /home/kisa/coding/80K_MPRA/gencode_v42/gencode.v42.annotation.gtf.gz
- GTEx /home/kisa/coding/80K_MPRA/80K-Analysis/09_mpra_manuscript/results/variant_info/gtex_filtered/gtex_prefiltered_by_mpra_1_based.tsv.gz

#### Scripts
- plot_examples.R - plots example variants including the variant visualized in the paper:
- src/plot_variants_v3.R - plotting variant loci
- src/preprocess_eqtl_data.py - preprocesss the eQTL data
- src/preprocess_ems_data.py - preprocesss the EMS data
- src/preprocess_variant_table.py - preprocess the variant and finemapped eQTL overlap
- src/ucsc_atac_explorer.R - used if ATAC data from UCSC is desired
