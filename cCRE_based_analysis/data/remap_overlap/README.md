### Overlapping tested CREs with ReMap2022
- The remap bed file is too big but can be downloaded.
- And processed like this:

```bash
# get the target column:
awk 'BEGIN{OFS="\t"} {$NF=$NF; split($4,a,":"); print $0, a[1]}' <(zcat remap2022_nr_macs2_hg38_v1_0.bed.gz) | gzip > remap2022_nr_macs2_hg38_v1_0_target_column.bed.gz

sort the remap file:
sort -k1,1 -k2,2n <( zcat remap2022_nr_macs2_hg38_v1_0_target_column.bed.gz ) | gzip > remap2022_nr_macs2_hg38_v1_0_target_column_sorted.bed.gz
```
- using bedtools intersect to get the overlap with tested CREs:

`bedtools intersect -a ./all_tested_CRE_sorted.bed.gz -b ./remap2022_nr_macs2_hg38_v1_0_target_column_sorted.bed.gz -wo > tested_CREs_remap2022.bed.gz`

- to filter the overlaps based on the summit, use the provided script like this:

`python filter_overlap_peak_summit_overlap_tested_CREs.py`

=> `tested_80k_MPRA_CREs_remap2022_filtered_summit_overlap_no_border_overlap_min_10bp.bed.gz`

- Get the remap2022 overlap counts for the active and inactive elements:
```bash
zcat tested_80k_MPRA_CREs_remap2022_filtered_summit_overlap_no_border_overlap_min_10bp.bed.gz | \
grep -Ff 2604_NGN2_scramble_ctrl_active_elements_name.tsv | \
awk '{print $4 "\t" $8}' | \
sort | \
uniq -c > 2604_active_elements_counts.tsv
```

```bash
zcat tested_80k_MPRA_CREs_remap2022_filtered_summit_overlap_no_border_overlap_min_10bp.bed.gz | \
grep -Ff 2604_NGN2_scramble_ctrl_inactive_elements_name.tsv | \
awk '{print $4 "\t" $8}' | \
sort | \
uniq -c > 2604_not_significant_counts.tsv
```

- To generate the concordance figure with the elastic net model and the ridge regression enrichment analysis
```bash
python run_tf_activity_inactivity_concordance.py \
 --active_counts data/remap_overlap/2604_active_elements_counts.tsv \
 --inactive_counts data/remap_overlap/2604_not_significant_counts.tsv \
 --counts_no_header \
 --annotation data/remap_overlap/ngn2_element_activity_metadata_202511_bbmap_bcalm_normalization.tsv.gz \
 --elastic_net data/remap_overlap/2512_ngn2_activity_elastic_net_model_202511_bbmap_bcalm_normalization_coefficients_all_regions_removed_na_phastcons_R2_test_0.0308.tsv \
 --outdir data/remap_overlap/ \
 --n_boot 1000 --fisher_fdr 0.05 --or_enrich 2.0 --or_deplete 0.8 \
 --ridge_fdr 0.05 --n_jobs -1 --seed 42
```

