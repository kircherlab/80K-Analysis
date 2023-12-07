#!/bin/bash
#SBATCH --job-name=match_missing_sequences
#SBATCH --partition=short
#SBATCH --account=hpc-ag-kircher # the SLURM account to charge
#SBATCH --nodes=1               # number of nodes to allocate
#SBATCH --ntasks=30              # number of processes (tasks) the job will start
#SBATCH --cpus-per-task=1      # number of CPUs (cores) each process (task) requires
#SBATCH --time=00-02:00         # how long the job is permitted to run, here 2 hours
#SBATCH --array=0-29
#SBATCH --mem-per-cpu=5G


samples=(merge_split0.join.fastq.gz merge_split18.join.fastq.gz merge_split27.join.fastq.gz merge_split1.join.fastq.gz merge_split19.join.fastq.gz merge_split28.join.fastq.gz merge_split10.join.fastq.gz merge_split2.join.fastq.gz merge_split29.join.fastq.gz merge_split11.join.fastq.gz merge_split20.join.fastq.gz merge_split3.join.fastq.gz merge_split12.join.fastq.gz merge_split21.join.fastq.gz merge_split4.join.fastq.gz merge_split13.join.fastq.gz merge_split22.join.fastq.gz merge_split5.join.fastq.gz merge_split14.join.fastq.gz merge_split23.join.fastq.gz merge_split6.join.fastq.gz merge_split15.join.fastq.gz merge_split24.join.fastq.gz merge_split7.join.fastq.gz merge_split16.join.fastq.gz merge_split25.join.fastq.gz merge_split8.join.fastq.gz merge_split17.join.fastq.gz merge_split26.join.fastq.gz merge_split9.join.fastq.gz)

awk -v "OFS=\t" 'NR==FNR {{a[$2] = $1; next}} {{if ($3 in a) print $1,a[$3],"270M";}}' \
<(
    cat /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/missing_reference_exact_uniq.tsv | awk -v "OFS=\t" '{{print $1,substr($2, 16,270)}}'
) \
<( zcat /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesignNoTemp/fastq/${samples[$SLURM_ARRAY_TASK_ID]} | awk 'NR%4==2 || NR%4==1' | paste - -
) > /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/match_missing_sequences/exact_match.${SLURM_ARRAY_TASK_ID}.tsv


