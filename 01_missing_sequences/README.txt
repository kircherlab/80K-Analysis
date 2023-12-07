File explanation:

bowtie*: output from bowtie
all_missing_seqs*: Preparation for the mappers

*.ipynb: analysis steps, work in progress
*.py: finished analysis steps (usable for workflows with small changes (adding usage of snakemake variable))

missing_sequence_lookup: slow python version of exact match

/match_missing_sequences/*: result of fasta bash script for exact matching reads to missing sequences

/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/match_missing_seqs_sbatch.sh: also computing exact matches 