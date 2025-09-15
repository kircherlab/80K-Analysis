reference=/home/kisa/coding/80K_MPRA/reference_genome_match/hg38.fa
for i in /home/kisa/coding/80K_MPRA/80K-Analysis/07_quality_control/notebooks/control_metadata/80k_hg38_dnase_controls_no_shuffled.bed; do
    name=`basename $i .bed`
    bedtools getfasta -nameOnly -fi $reference -bed $i | \
    awk '{if ($0 ~ /^>/){print $0} else { print toupper($0)}}' > $name.tmp

    fasta-dinucleotide-shuffle -f $name.tmp -s 42 -c 200 | \
    python /home/kisa/coding/80K_MPRA/MPRAOligoDesign/workflow/scripts/kMerFilter.py -i $name.tmp -k 6 > $name.diShuffled.fa
    rm $name.tmp
done;