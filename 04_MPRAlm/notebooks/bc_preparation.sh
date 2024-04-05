#!/bin/bash

config=standard
output_name=standard_with_controls


# output files:
counts_per_bc_gz_name=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/${output_name}_NGN2_counts_per_bc_sorted.tsv.gz
counts_per_bc_name=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/${output_name}_NGN2_counts_per_bc_sorted.tsv

# variants ref and alt 
sequences_to_oligos=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/${output_name}_NGN2_seqs_to_oligos_sorted.tsv.gz
joined_counts_and_seqs=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/${output_name}_NGN2_counts_sequences.tsv.gz
joined_counts_and_seqs_unzipped=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/${output_name}_NGN2_counts_sequences.tsv
single_seqs=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/${output_name}_NGN2_single_seqs.tsv
temp=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/${output_name}_NGN2_tmp.tsv
final_output=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/${output_name}_NGN2_filtered_counts_sequences.tsv.gz # 3649963


# input files:

# # standard config
# variants=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/80K_MPRA/design/variant_region_map_deduplicated.tsv.gz
# assignment=/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/experiments/${config}_bwa/assignment/assignmentFixDuplicates.tsv.gz
# count1=/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/experiments/${config}_bwa/assigned_counts/assignmentFixDuplicates/NGN2_1.merged.config.standardConfig.tsv.gz
# count2=/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/experiments/${config}_bwa/assigned_counts/assignmentFixDuplicates/NGN2_2.merged.config.standardConfig.tsv.gz
# count3=/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/experiments/${config}_bwa/assigned_counts/assignmentFixDuplicates/NGN2_3.merged.config.standardConfig.tsv.gz

# standard with variant controls
variants=/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/05_variant_region_list/resources/controls_variant_region_map.tsv
assignment=/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/experiments/${config}_bwa/assignment/assignmentFixDuplicates.tsv.gz
count1=/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/experiments/${config}_bwa/assigned_counts/assignmentFixDuplicates/NGN2_1.merged.config.standardConfig.tsv.gz
count2=/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/experiments/${config}_bwa/assigned_counts/assignmentFixDuplicates/NGN2_2.merged.config.standardConfig.tsv.gz
count3=/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/experiments/${config}_bwa/assigned_counts/assignmentFixDuplicates/NGN2_3.merged.config.standardConfig.tsv.gz


# # low config
# variants=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/80K_MPRA/design/variant_region_map_deduplicated.tsv.gz
# assignment=/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/experiments/${config}_bwa/assignment/assignmentFixDuplicates.tsv.gz

# count1=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/experiment/standard_results/results/experiments/${config}_bwa/assigned_counts/assignmentFixDuplicates/NGN2_3.merged.config.lowConfig.tsv.gz 
# count2=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/experiment/standard_results/results/experiments/${config}_bwa/assigned_counts/assignmentFixDuplicates/NGN2_2.merged.config.lowConfig.tsv.gz 
# count3=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/experiment/standard_results/results/experiments/${config}_bwa/assigned_counts/assignmentFixDuplicates/NGN2_1.merged.config.lowConfig.tsv.gz

# assign counts to barcodes and oligo ids
join <(zcat $assignment) <(zcat $count1) -t $'\t' | join - <(zcat $count2) -t $'\t' | join - <(zcat $count3) -t $'\t' | cut -f1,2,5-11 | sort -k 2 > $counts_per_bc_name 
gzip $counts_per_bc_name

# assign oligo ids to seqs and add "ref", "alt" ($1 ID, $2 ref, $3 alt)
zcat $variants | awk 'NR > 1 {print $1 "\t" $3 "\tref\n" $1 "\t" $4 "\talt"}' | sort -k 2 | gzip -c > $sequences_to_oligos


# join counts with sequence oligos file
join -1 2 -2 2  <(zcat $sequences_to_oligos) <(zcat $counts_per_bc_gz_name) -t $'\t' | cut -f 2- | sort -k 1 > $joined_counts_and_seqs_unzipped

gzip $joined_counts_and_seqs_unzipped

# storing all sequences once
zcat $joined_counts_and_seqs | awk 'NR >1{print $1","$2}' | uniq | awk -v FS="," '{print $1 "\t" $2}' | cut -f 1 |sort | uniq -c |  awk '$1 != 2' |  awk '{print $2}' > $single_seqs

# identifying sequences that have ref and alt
join -v 1 <(zcat $joined_counts_and_seqs) $single_seqs -t $'\t' | gzip -c > $temp

mv $temp $final_output

# attach header
header="variant_id\tallele\tbarcode\tDNA1\tRNA1\tDNA2\tRNA2\tDNA3\tRNA3"
echo -e $header | gzip -c | cat - $final_output > temp
mv temp $final_output


###### Files for old assignment:
# # output files:
# counts_per_bc_gz_name=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/NGN2_counts_per_bc_sorted.tsv.gz
# counts_per_bc_name=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/NGN2_counts_per_bc_sorted.tsv

# # variants ref and alt 
# sequences_to_oligos=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/NGN2_seqs_to_oligos_sorted.tsv.gz
# joined_counts_and_seqs=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/NGN2_counts_sequences.tsv.gz
# joined_counts_and_seqs_unzipped=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/NGN2_counts_sequences.tsv
# single_seqs=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/NGN2_single_seqs.tsv
# temp=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/NGN2_tmp.tsv
# final_output=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/NGN2_filtered_counts_sequences.tsv.gz

# # input files:
# variants=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/design/final_design/results/final_design/cardiac_neuro_cava_random/variant_region_map.tsv.gz
# assignment=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/experiment/results/experiments/run1counts_run2Assignment_NoDupAss/assignment/assignmentFixDuplicates.tsv.gz
# count1=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/experiment/results/experiments/run1counts_run2Assignment_NoDupAss/assigned_counts/assignmentFixDuplicates/NGN2_1.merged.config.standardConfig.tsv.gz
# count2=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/experiment/results/experiments/run1counts_run2Assignment_NoDupAss/assigned_counts/assignmentFixDuplicates/NGN2_2.merged.config.standardConfig.tsv.gz
# count3=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/experiment/results/experiments/run1counts_run2Assignment_NoDupAss/assigned_counts/assignmentFixDuplicates/NGN2_3.merged.config.standardConfig.tsv.gz
