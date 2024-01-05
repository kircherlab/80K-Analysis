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

tmp/ folder for missing_reference_exact_uniq.tsv / missing_seq_dict.json / old_missing_reference_exact.tsv (it can be interesting for downstream analysis but it is not needed for the workflow)

Find all reads matched to the missing sequences
- found problem with names in reads: adding "_" in reads name: e.g. `zcat $read_dir"merge_split0.join.fastq.gz" | head -n 30 | awk 'NR%4==2 || NR%4==1' | paste - - | sed -e 's/ /_/g' > $exact_match_dir"exact_match.0.tsv`

align these reads to the reference genome

Check for each aligner: how many reads are aligned to the missing sequences and how many reads are matched to the other sequences?

Found a lot of different reads with the same sequence: 1812292 (first sequence in `/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/match_missing_sequences/matched_missing_reads.fa` is present ~730 times)

non duplicated reads are in config["files"]["unique_matched_reads"] (2800) -> is fine, found 2800 matching sequences

## Woring on the alignment of the reads to the reference genome
### Bowtie:
- Output: 3.2M
    - Total time for backward call to driver() for mirror index: 00:00:21
    - Setting the index via positional argument will be deprecated in a future release. Please use -x option instead.
        - reads processed: 2800
        - reads with at least one alignment: 2800 (100.00%)
        - reads that failed to align: 0 (0.00%)
    - Reported 4425 alignments => `/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/new_alignments/bowtie/bowtie_alignment_k4_v0.sam`

### BWA:
- Output: 12.89 MB
    - [main] Real time: 23.935 sec; CPU: 20.033 sec
    - [M::bwa_idx_load_from_disk] read 0 ALT contigs
    - [M::process] read 2800 sequences (756000 bp)...
    - [M::mem_process_seqs] Processed 2800 reads in 5.543 CPU sec, 5.446 real sec
    - [main] Version: 0.7.17-r1188
    - [main] CMD: `bwa mem -t 10 -L 80 -M -C /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/resources/align_missing_sequences/reference/cor_design_no_duplicates_sequence_and_header.fa /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/match_missing_sequences/unique_missing_reads.fastq`
    - [main] Real time: 5.712 sec; CPU: 5.669 sec

### BWA mem 2: (index works but mem does not)
- Output (mem - failing) (Environment: 280 GB mem and 12 Threads)
```bash 
*** buffer overflow detected ***: bwa-mem2 terminated different_aligners.sh: line 59: 1493544 Aborted (core dumped) bwa-mem2 mem -t $threads $bwa_work_dir/index $no_duplicated_reads > $bwa_sam_output 
``` 

### BBMap: 
- Output:
    - low-Q discards:           0.0000%               0         0.0000%                  0
    - perfect best site:      100.0000%            2800       100.0000%             756000
    - semiperfect site:       100.0000%            2800       100.0000%             756000
    - Match Rate:                   NA               NA       100.0000%             756000
    - Error Rate:               0.0000%               0         0.0000%                  0
    - Sub Rate:                 0.0000%               0         0.0000%                  0
    - Del Rate:                 0.0000%               0         0.0000%                  0
    - Ins Rate:                 0.0000%               0         0.0000%                  0
    - N Rate:                   0.0000%               0         0.0000%                  0 
    - Total time:             11.919 seconds.

## Running bowtie in the MPRAsnakeflow
- remove temp() of the rules: assignment_mapping, assignment_merge (fastq/merge... bam/merge_split...)
- run MPRAsnakeflow with standardConfig
## Comparing low config and number of variants with standard config results of #variants
- check reference:
    - Number of sequences: 80215
    - list of all variants: 46458 `cat /fast/groups/ag_kircher/MPRA/IGVF_Y1_design/resources/association_data/design_no_duplicates_sequence_and_header.fa | grep ">" | grep -E "\-[A-Z]+-[A-Z]+" | wc -l`
      - write list to tsv file
    - Number of SLEA in header: 200