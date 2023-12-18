#!/bin/bash

## [Bowtie](https://bowtie-bio.sourceforge.net/manual.shtml)
# # index reference
# bowtie-build input_reference.fasta index_prefix
# # align the reads
# bowtie [-q|-f|-r|-c] index_prefix [-1 input_reads_pair_1.[fasta|fastq] -2 input_reads_pair_2.[fasta|fastq] | input_reads.[fasta|fastq]] [options]
working_dir="/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/new_alignments/"
split=0
threads=30
cd $working_dir

bowtie-build  /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/all_missing_seqs_270_uniq.fa bowtie_index0811_uniq

bowtie -q bowtie_index0811_uniq <(gzip -dc /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesignNoTemp/fastq/merge_split${split}.join.fastq.gz ) -k 4 -v 0 -p ${threads} > /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/new_alignments/bowtie_alignments_k4_v0_split${split}.sam

## BWA
cat /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/all_missing_seqs_270_uniq.fa | awk '{{gsub(/[\]\[]/,"_")}}$0' > /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/new_alignments/bwa/all_missing_seqs_270_bwa_reference.fa;
bwa index -a bwtsw /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/new_alignments/bwa/all_missing_seqs_270_bwa_reference.fa;
samtools faidx /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/new_alignments/bwa/all_missing_seqs_270_bwa_reference.fa

# bwa mem -t ${threads} -L 80 -M -C /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/new_alignments/bwa/all_missing_seqs_270_bwa_reference.fa <(
#             gzip -dc /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesignNoTemp/fastq/merge_split${split}.join.fastq.gz
#         ) > /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/new_alignments/bwa/merge_split${split}.mapped.sam

# old, with samtools sort
bwa mem -t ${threads} -L 80 -M -C /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/new_alignments/bwa/all_missing_seqs_270_bwa_reference.fa <(
            gzip -dc /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesignNoTemp/fastq/merge_split${split}.join.fastq.gz
        )  | samtools sort -l 0 -@ ${threads} > /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/new_alignments/bam/merge_split${split}.mapped.bam

## [BBmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/) (global)
# To index and map at the same time:
# bbmap.sh in=reads.fq out=mapped.sam ref=ref.fa

#! Plan: 1. fw and rv hinzufügen
#! Plan BBmap, Minimap2, BWA mem2, Bowtie
#! Plan 1. bowtie + -k 4 + bowtie + norev (siehe max nachricht)