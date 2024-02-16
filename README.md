# File explaining the content of this analysis folder

## Summary
- Code I used to analyze the 80K data
  - work in progress is done in .ipynb
  - GOAL: simplified scripts (.py) for the final analysis pipeline
- Directories:
  - `01_missing_sequences`: analysis of the missing sequences (mapping problem)
  - `02_controls`: plotting the controlls and how well they perform
  - `03_config_check`: try different configurations for the analysis
  - `04_MPRAlm`: Analyze 80K with MPRAlm (e.g. compare to bc_MPRAlm; plots for tiger team)
  - `05_variant_region_list`: get a table which has for each sequence tested the meta data and the additional information (e.g. enformer class, cell-type)
- Resulting data is stored in the group folder `/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/project/<folder_name>/`

### Code
- missing_sequence_lookup.py: script for finding the missing sequences in the read files by exact matches
  - input: one merged read file (fastq) and the reference set as well as the assignement
  - output: hits of the missing sequences in the read file
- `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesign/missing_sequences_per_label.ipynb` preperation of the missing_sequence_lookup script

### Data
- sorted_counts.tsv
  - counts of the BC and oligo assignment (10451796 lines)
- `missing_sequence_dict.json`: dictionary of the missing sequences and their header
- `match_missing_sequences/exact_match.{split}.tsv`: matched missing sequences to the reference (for used code see `labintern/notes/kilian/2023/MPRA_80K/investigate_files.md` (`run_seq_match_sbatch.sh` similar `80K_analysis/01_missing_sequences/match_missing_seqs_sbatch.sh`))


