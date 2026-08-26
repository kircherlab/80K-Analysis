Using H13 clustered and Fimo to get the overlapping TFs for the sequences of interest using the following command and place the fimo.tsv here:

```bash

# cluster H13CORE meme format
awk 'NR==FNR{a[$1]++; next}FNR < 10 {print $0}FNR>=10 {if (a[$2]) toprint=1; if ($0=="" && toprint) {print ""; toprint=0} if(toprint) print $0}' resources/motif_analysis/cluster_list.tsv resources/motif_analysis/H13CORE_meme_format.meme > resources/motif_analysis/H13CORE_meme_format_clustered.meme


conda create -n fimo_env bioconda::meme
conda activate fimo_env

# use fimo: https://meme-suite.org/meme/doc/fimo.html
fimo --max-stored-scores 1000000 --o ./all_neuro_ctrls_scrambled_tested_h13  ./H13CORE_meme_format_clustered.meme <( zcat 80K_tested_neuro_ctrls_scrambled_fimo_input.fa.gz )
```