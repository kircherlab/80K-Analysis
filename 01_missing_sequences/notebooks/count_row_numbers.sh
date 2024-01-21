# bash script which calls samtools view and counts the numbers of rows in each bam file

missing_sequences_bam_list=('/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/identified_missing_sequences/missing_sequences_min_quality.bam' '/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/identified_missing_sequences/missing_sequences_alignment_start.bam' '/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/identified_missing_sequences/missing_sequences_alignment_end.bam' '/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/identified_missing_sequences/missing_sequences_sequence_length_min.bam' '/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/identified_missing_sequences/missing_sequences_sequence_length_max.bam')

for bam in "${missing_sequences_bam_list[@]}"; do
    echo $bam
    samtools view -c $bam
    echo "-------------------"
done