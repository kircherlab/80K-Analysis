#!/bin/bash
## Restart at 19.12.2023: 
# non duplicated reads config["files"]["unique_matching_reads"] = /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/match_missing_sequences/unique_missing_reads.fastq
no_duplicated_reads=/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/match_missing_sequences/unique_missing_reads.fastq
# reference: "/fast/groups/ag_kircher/MPRA/IGVF_Y1_design/resources/association_data/design_no_duplicates_sequence_and_header.fa"
uncorrected_reference=/fast/groups/ag_kircher/MPRA/IGVF_Y1_design/resources/association_data/design_no_duplicates_sequence_and_header.fa
reference=/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/resources/align_missing_sequences/reference/cor_design_no_duplicates_sequence_and_header.fa
## [Bowtie](https://bowtie-bio.sourceforge.net/manual.shtml)
# # index reference
# bowtie-build input_reference.fasta index_prefix
# # align the reads
# bowtie [-q|-f|-r|-c] index_prefix [-1 input_reads_pair_1.[fasta|fastq] -2 input_reads_pair_2.[fasta|fastq] | input_reads.[fasta|fastq]] [options]
working_dir="/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/new_alignments/bowtie/"
split=0
threads=1
cd $working_dir
# ## reference correction (removing brackets)
# cat $uncorrected_reference | awk '{{gsub(/[\]\[]/,"_")}}$0' > $reference;
## bowtie: 
bowtie_output=/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/new_alignments/bowtie/bowtie_alignment_S_k4_v0.sam
bowtie-build $reference bowtie_index1912_uniq

bowtie -q bowtie_index1912_uniq $no_duplicated_reads -S -k 4 -v 0 -p ${threads} > $bowtie_output

## old: 
# bowtie-build  /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/all_missing_seqs_270_uniq.fa bowtie_index0811_uniq

# bowtie -q bowtie_index0811_uniq <(gzip -dc /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesignNoTemp/fastq/merge_split${split}.join.fastq.gz ) -k 4 -v 0 -p ${threads} > /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/new_alignments/bowtie_alignments_k4_v0_split${split}.sam

## BWA
# bwa_output=/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/new_alignments/bwa/bwa_alignment_L80_M_C.sam
# bwa index -a bwtsw $reference
# samtools faidx $reference

# bwa mem -t ${threads} -L 80 -M -C $reference $no_duplicated_reads > $bwa_output
# echo "written to $bwa_output"
# old:
# cat /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/all_missing_seqs_270_uniq.fa | awk '{{gsub(/[\]\[]/,"_")}}$0' > /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/new_alignments/bwa/all_missing_seqs_270_bwa_reference.fa;
# bwa index -a bwtsw /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/new_alignments/bwa/all_missing_seqs_270_bwa_reference.fa;
# samtools faidx /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/new_alignments/bwa/all_missing_seqs_270_bwa_reference.fa

# bwa mem -t ${threads} -L 80 -M -C /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/new_alignments/bwa/all_missing_seqs_270_bwa_reference.fa <(
#             gzip -dc /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesignNoTemp/fastq/merge_split${split}.join.fastq.gz
#         ) > /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/new_alignments/bwa/merge_split${split}.mapped.sam

# old, with samtools sort
# bwa mem -t ${threads} -L 80 -M -C /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/new_alignments/bwa/all_missing_seqs_270_bwa_reference.fa <(
#             gzip -dc /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesignNoTemp/fastq/merge_split${split}.join.fastq.gz
#         )  | samtools sort -l 0 -@ ${threads} > /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/new_alignments/bam/merge_split${split}.mapped.bam


## BWA mem2 (installed bwa-mem2, also installed it from github) ! not possible
# bwa-mem2 index ref.fa
# bwa-mem2 mem ref.fa read1.fq read2.fq > out.sam
# bwa_work_dir=/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/new_alignments/bwa_mem2
# bwa_sam_output=$bwa_work_dir/bwa_mem2_alignment_no_flag.sam
# cd $bwa_work_dir
# # bwa-mem2 index -p $bwa_work_dir/index $reference
# # bwa-mem2 mem -t $threads $bwa_work_dir/index $no_duplicated_reads > $bwa_sam_output
/data/gpfs-1/users/kisa11_c/work/coding/bwa-mem2/bwa-mem2 index -p $bwa_work_dir/index $reference
/data/gpfs-1/users/kisa11_c/work/coding/bwa-mem2/bwa-mem2 mem -t $threads $bwa_work_dir/index $no_duplicated_reads > $bwa_sam_output


## [BBmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/) (global)
# To index and map at the same time:
# bbmap.sh in=reads.fq out=mapped.sam ref=ref.fa
# bbmap.sh in=$no_duplicated_reads out=$bwa_sam_output ref=$reference

#! Plan: 1. fw and rv hinzufügen
#! Plan BBmap, Minimap2, BWA mem2, Bowtie
#! Plan 1. bowtie + -k 4 + bowtie + norev (siehe max nachricht)