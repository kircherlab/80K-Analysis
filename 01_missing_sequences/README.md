File explanation:
all_missing_sequences.fa: result from missing_sequences_per_label.ipynb (are all the missing sequences in the 80K dataset?!)
bowtie*: output from bowtie
all_missing_seqs*: Preparation for the mappers

*.ipynb: analysis steps, work in progress
*.py: finished analysis steps (usable for workflows with small changes (adding usage of snakemake variable))

`run_seq_match_sbatch.sh`: generates the reverse complement of all missing sequenes for exact matching with `match_missing_seqs_sbatch.sh` -> result missing_reference_exact.tsv

missing_sequence_lookup: slow python version of exact match

/match_missing_sequences/*: result of fasta bash script for exact matching reads to missing sequences

/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/match_missing_seqs_sbatch.sh: also computing exact matches (fast)

different_aligners.sh: aligning the sequences with different aligners

`/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/missing_reference_exact_uniq.tsv`: includes fw and rv in the end for forward and reverse complement sequences -> all are unique

example `cardiac_neuro_cava_random:REF_TRIO|ENSG00000038382.23|EH38E2358113_fwd_tile1-1_fw`: can be found in missing_reference_exact_uniq.tsv and exact_match.0.tsv
