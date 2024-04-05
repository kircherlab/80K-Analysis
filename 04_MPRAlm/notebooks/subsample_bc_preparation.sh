#!/bin/bash

# output files:
counts_per_bc_gz_name=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/test_data/NGN2_counts_per_bc.tsv.gz
sorted_counts_per_bc_name=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/test_data/NGN2_counts_per_bc.sorted.tsv

# input files:
assignment=/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/results/experiments/run1counts_run2Assignment_NoDupAss/assignment/assignmentFixDuplicates.tsv.gz
count1=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/test_data/sub_NGN2_1.merged.config.standardConfig.tsv.gz
count2=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/test_data/sub_NGN2_2.merged.config.standardConfig.tsv.gz
count3=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff/test_data/sub_NGN2_3.merged.config.standardConfig.tsv.gz

# assign counts to barcodes and oligo ids
join <(zcat $assignment) <(zcat $count1) -t $'\t' | join - <(zcat $count2) -t $'\t' | join - <(zcat $count3) -t $'\t' | cut -f1,2,5-11 | bgzip -c > $counts_per_bc_gz_name

# sort 
zcat $counts_per_bc_gz_name | sort -k 2 -o $sorted_counts_per_bc_name

# zipp 
bgzip $sorted_counts_per_bc_name
