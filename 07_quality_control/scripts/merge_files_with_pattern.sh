# This script merges all files with a defined pattern in the directory and writes the header only once

directory=/home/kisa/coding/80K_MPRA/enformer_data/enformer_matched_with_variants_enformer_class
output_file=merged_enformer_variant_all_max_no_filter.tsv
cd $directory
first=1

for file in *prioritized_variants_enformer_class.tsv; do
    if [ "$first" ]
    then
        cat "$file"
        first=
    else
        cat "$file" | tail -n +2
    fi
done > $output_file
