#!/bin/bash

# integrate reverse complement sequences
paste <(
    cat /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/all_missing_sequences.fa | awk '{{if ($1 ~ /^>/) {{ gsub(/[\]\[]/,"_"); print substr($1,2)"_forw"}}}}';
    cat /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/all_missing_sequences.fa | awk '{{if ($1 ~ /^>/) {{ gsub(/[\]\[]/,"_"); print substr($1,2)"_revc"}}}}';
) <(
    cat /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/all_missing_sequences.fa | awk '{{if ($1 ~ /^[^>]/) {{ seq=seq$1}}; if ($1 ~ /^>/ && NR!=1) {{print seq; seq=""}}}} END {{print seq}}';
    cat /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/all_missing_sequences.fa | awk '{{if ($1 ~ /^[^>]/) {{ seq=seq$1}}; if ($1 ~ /^>/ && NR!=1) {{print seq; seq=""}}}} END {{print seq}}' | tr ACGTacgt TGCAtgca | rev;
) > /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/missing_reference_exact.tsv


# Look up exact matches in design file
export LC_ALL=C # speed up sort
for split in {0..29}
do
    awk -v "OFS=\\t" 'NR==FNR {{a[$2] = $1; next}} {{if ($3 in a) print $2,a[$3],"270M"; }}' \
    <(
        cat /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/missing_reference_exact.tsv | awk -v "OFS=\t" '{{print $1,substr($2, 16,270)}}'
    ) \
    <(
        zcat "/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesignNoTemp/fastq/merge_split${split}.join.fastq.gz" | awk 'NR%4==2 || NR%4==1' | paste - -
    ) > "/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/match_missing_sequences/exact_match.${split}.tsv"
done
