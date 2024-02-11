#!/bin/bash

#SBATCH --job-name=use_bwa_mem2
#SBATCH --partition=medium
#SBATCH --account=hpc-ag-kircher # the SLURM account to charge
#SBATCH --nodes=1               # number of nodes to allocate
#SBATCH --ntasks=1              # number of processes (tasks) the job will start
#SBATCH --cpus-per-task=30      # number of CPUs (cores) each process (task) requires
#SBATCH --mem=128G
#SBATCH --time=00-10:00         # how long the job is permitted to run, here 2 hours

# input_ref=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/reference/reference_bwa-mem2.fa
# input_fastq=/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/match_missing_sequences/corr_matched_missing_reads.fastq # no "/" at the end
# threads=20
# output_path=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/experiment/standard_results/results/bwa2/assignment/standardAssignIGVFDesignNoTemp/bam # no "/" at the end
# log_file=/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/match_missing_sequences/bwa-mem2_test.log


# mkdir -p $output_path
# cd /data/gpfs-1/users/kisa11_c/work/coding/bwa-mem2



# ./bwa-mem2 mem -t $threads -L 80 -M -C $input_ref $input_fastq 2> $log_file
# ./bwa-mem2 mem -t $threads -L 80 -M -C $input_ref $input_fastq | samtools sort -l 0 -@ $threads > $output_path/test_bwa-mem2.bam 2> $log_file

row_number=148800

# # working tail: ref_tail_11226_server_test_1101_index
# # not working: ref_tail_11228_server_test_1101_index
# tail -n $row_number /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/resources/align_missing_sequences/reference/cor_design_no_duplicates_sequence_and_header.fa > test_data/test_ref_tail_${row_number}.fa

# ./bwa-mem2 index -p test_data/ref_tail_${row_number}_server_test_1101_index test_data/test_ref_tail_${row_number}.fa

# ./bwa-mem2 mem -t 1 test_data/ref_tail_${row_number}_server_test_1101_index /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/match_missing_sequences/corr_matched_missing_reads.fastq > test_data/ref_tail_${row_number}_server_test_1101_alignment.sam

# working head: ref_head_148798_server_test_1101_index
# not working: ref_head_148800_server_test_1101_index
head -n $row_number /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/resources/align_missing_sequences/reference/cor_design_no_duplicates_sequence_and_header.fa > test_data/test_ref_head_${row_number}.fa

./bwa-mem2 index -p test_data/ref_head_${row_number}_server_test_1101_index test_data/test_ref_head_${row_number}.fa

./bwa-mem2 mem -t 1 test_data/ref_head_${row_number}_server_test_1101_index /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/match_missing_sequences/corr_matched_missing_reads.fastq > test_data/ref_head_${row_number}_server_test_1101_alignment.sam
