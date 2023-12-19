#!/bin/bash

# input missing sequences (unique)
missing_seqs_reference=/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/resources/align_missing_sequences/all_missing_sequences.fa

# input reads (fastq.gz) files
read_dir=/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesignNoTemp/fastq/

# output and input missing sequences with reverse complement
missing_seqs_revcomp=/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/resources/align_missing_sequences/missing_reference_exact.tsv

# output directory for exact matches
exact_match_dir=/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/match_missing_sequences/

# integrate reverse complement sequences
paste <(
    cat $missing_seqs_reference | awk '{{if ($1 ~ /^>/) {{ gsub(/[\]\[]/,"_"); print ">"substr($1,2)"_forw"}}}}';
    cat $missing_seqs_reference | awk '{{if ($1 ~ /^>/) {{ gsub(/[\]\[]/,"_"); print ">"substr($1,2)"_revc"}}}}';
) <(
    cat $missing_seqs_reference | awk '{{if ($1 ~ /^[^>]/) {{ seq=seq$1}}; if ($1 ~ /^>/ && NR!=1) {{print seq; seq=""}}}} END {{print seq}}';
    cat $missing_seqs_reference | awk '{{if ($1 ~ /^[^>]/) {{ seq=seq$1}}; if ($1 ~ /^>/ && NR!=1) {{print seq; seq=""}}}} END {{print seq}}' | tr ACGTacgt TGCAtgca | rev;
) > $missing_seqs_revcomp


# Look up exact matches in design file
export LC_ALL=C # speed up sort
for split in {0..29}
do
    awk -F'\t' -v "OFS=\t" 'NR==FNR {{a[$2] = $1; next}} {{if ($2 in a) print $1,a[$2],$2,"270M"; }}' \
    <(
        cat $missing_seqs_revcomp | awk -v "OFS=\t" '{{print $1,substr($2, 16,270)}}'
    ) \
    <(
        zcat $read_dir"merge_split${split}.join.fastq.gz" | awk 'NR%4==2 || NR%4==1' | paste - - 
    ) > $exact_match_dir"exact_match.${split}.tsv"
done

# Merge the result in one file
exact_match_dir=/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/match_missing_sequences
cd $exact_match_dir
output_file=all_exact_match.tsv
if [ -e "$output_file" ]; then
    echo "File $output_file already exists. Are you sure you want to do that? Please remove it first."
else
    cat exact_match.*.tsv > $output_file
fi