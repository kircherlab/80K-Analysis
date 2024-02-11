#!/bin/bash

#SBATCH --job-name=use_bwa_mem2
#SBATCH --partition=medium
#SBATCH --account=hpc-ag-kircher # the SLURM account to charge
#SBATCH --nodes=1               # number of nodes to allocate
#SBATCH --ntasks=30              # number of processes (tasks) the job will start
#SBATCH --cpus-per-task=2      # number of CPUs (cores) each process (task) requires
#SBATCH --time=00-02:00         # how long the job is permitted to run, here 2 hours
#SBATCH --array=0-29
#SBATCH --mem=72G

input_ref=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/reference/reference_bwa-mem2.fa
fastq_dir=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/fastq # no "/" at the end
threads=2
output_path=/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/experiment/standard_results/results/bwa2/assignment/standardAssignIGVFDesignNoTemp/bam # no "/" at the end

samples=(merge_split0.join.fastq.gz merge_split18.join.fastq.gz merge_split27.join.fastq.gz merge_split1.join.fastq.gz merge_split19.join.fastq.gz merge_split28.join.fastq.gz merge_split10.join.fastq.gz merge_split2.join.fastq.gz merge_split29.join.fastq.gz merge_split11.join.fastq.gz merge_split20.join.fastq.gz merge_split3.join.fastq.gz merge_split12.join.fastq.gz merge_split21.join.fastq.gz merge_split4.join.fastq.gz merge_split13.join.fastq.gz merge_split22.join.fastq.gz merge_split5.join.fastq.gz merge_split14.join.fastq.gz merge_split23.join.fastq.gz merge_split6.join.fastq.gz merge_split15.join.fastq.gz merge_split24.join.fastq.gz merge_split7.join.fastq.gz merge_split16.join.fastq.gz merge_split25.join.fastq.gz merge_split8.join.fastq.gz merge_split17.join.fastq.gz merge_split26.join.fastq.gz merge_split9.join.fastq.gz)

mkdir -p $output_path
cd /data/gpfs-1/users/kisa11_c/work/coding/bwa-mem2
./bwa-mem2 mem -t $threads -L 80 -M -C $input_ref <(
    gzip -dc $fastq_dir/${samples[$SLURM_ARRAY_TASK_ID]}
)  | samtools sort -l 0 -@ $threads > $output_path/merge_split${SLURM_ARRAY_TASK_ID}.bwa-mem2.bam