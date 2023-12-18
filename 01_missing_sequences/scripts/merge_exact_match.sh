exact_match_dir=/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/match_missing_sequences
cd $exact_match_dir
output_file=all_exact_match.tsv
if [ -e "$output_file" ]; then
    echo "File $output_file already exists. Are you sure you want to do that? Please remove it first."
else
    cat exact_match.*.tsv > $output_file
fi