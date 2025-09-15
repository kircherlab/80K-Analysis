### Snakemake workflow for the analysis of TFBS enrichment between groups of sequences using a single master bed of high-confidence TFBSs
#### Step 1: Prepare the input files
1. Single master bed file of high-confidence TFBSs (e.g. from ReMap 2022 non-redundant peaks (using the peak summits))
2. Two bed files of sequences to compare

#### Workflow steps summary:
1. Generate the overlap files using bedtools intersect
2. Process the overlap files for each TFBS to create a contingency table => perform Fisher's exact test => FDR correction
3. Plot the resulting data for a dot plot of the significant TFBSs using the odds ratio on the x-axis, size (number of findings), shape (significants) in a single seaborn, matplotlib plot.

#### Configuration file (config.yaml) example:
```yaml
# Input files
master_TFBS_bed: "path/to/master_TFBS.bed"  #
sequence_set_A_bed: "path/to/sequence_set_A.bed"
sequence_set_B_bed: "path/to/sequence_set_B.bed"

# Names for plotting and output files
active_label: "Active"
inactive_label: "Inactive"
output_prefix: "TFBS_enrichment"
# Statistical parameters
fdr_threshold: 0.05
min_findings: 5  # Minimum number of findings to consider a TFBS for plotting
# Plotting parameters
plot_title: "TFBS Enrichment Analysis"
plot_xlabel: "Odds Ratio"
plot_ylabel: "Transcription Factors"
plot_size: [10, 8]  # Width, Height in inches
plot_dpi: 300
plot_format: "png"  # e.g., png, pdf, svg
```

#### Step 2: Compute counts for each TFBS in both sets of sequences and perform statistical tests

- rules
    - compute_TFBS_counts
    - tfbs_fisher_exact_test

- detailed: TFBS Enrichment in Sequence set A vs. Sequence set B
    goal: Efficiently calculate TFBS enrichment using a single master file and prepare a consolidated data table for plotting. Data table output structure: | TFBS_name | odds_ratio | adjusted_pvalue | active_overlaps | comparison_group |

    Tools: Python with pandas and scipy.stats, bedtools intersect.

    Procedure:

        Generate Overlap File:

            Use bedtools intersect to find overlaps between your active and inactive sequences and the master TFBS BED file. The -wao (write all overlaps) flag is critical here as it keeps all information from both files, including the TF names.

            bedtools intersect -a active_sequences.bed -b master_TFBS.bed -wao > active_overlaps.txt

            bedtools intersect -a inactive_sequences.bed -b master_TFBS.bed -wao > inactive_overlaps.txt

        Process Overlaps and Build Contingency Tables:

            Write a Python script to read active_overlaps.txt and inactive_overlaps.txt.

            For each TFBS, count the number of sequences it overlaps in both the active and inactive sets.

            The script should create a pandas DataFrame to store the counts for each TFBS. The columns of this DataFrame should be: TFBS_name, active_overlaps, inactive_overlaps, active_non_overlaps, inactive_non_overlaps.

            The non_overlaps counts can be derived from the total number of sequences in each group. For example, active_non_overlaps would be (Total Active Sequences) - (active_overlaps).

        Perform Fisher's Exact Test and FDR Correction:

            Iterate through the DataFrame from the previous step. For each TFBS, use scipy.stats.fisher_exact to perform the test on the 2×2 contingency table. Store the resulting p-value and odds ratio.

            After testing all TFBS, apply a multiple comparison correction (e.g., Benjamini-Hochberg) to the p-values.

Assumptions of Fisher's Exact Test

Fisher's Exact Test is a non-parametric test and does not assume a normal distribution. It is particularly well-suited for this use case. The main assumptions are:

    Categorical Data: The data must be in the form of counts in a 2×2 contingency table, which is exactly what we have (TFBS Present vs. TFBS Absent, and Active vs. Inactive).

    Independence: The observations must be independent. In your case, each sequence is an independent observation, and the presence or absence of a TFBS within it is a single event. Therefore, the assumptions are met.

Step 3: Addressing Question 2: Enrichment in Open-Chromatin vs. Non-Open-Chromatin Sequences (Revised)

Goal: Use the same efficient approach to compare enrichment within the active sequence subgroups.

    Tools: Python with pandas and scipy.stats, bedtools intersect.

    Procedure:

        Generate Overlap Files for Subgroups:

            You'll need active_oc.bed and active_non_oc.bed as created in the previous plan.

            Run bedtools intersect for each of these with your master_TFBS.bed file:

            bedtools intersect -a active_oc.bed -b master_TFBS.bed -wao > active_oc_overlaps.txt

            bedtools intersect -a active_non_oc.bed -b master_TFBS.bed -wao > active_non_oc_overlaps.txt

        Process and Build Contingency Tables (Replicated Logic):

            Use the same Python script logic from Step 2, but this time process the _oc_overlaps.txt and _non_oc_overlaps.txt files.

            The new DataFrame will have columns for TFBS_name, active_oc_overlaps, active_non_oc_overlaps, and their respective non_overlaps counts.

        Perform Fisher's Exact Test and FDR Correction:

            Repeat the statistical testing and FDR correction steps on this new DataFrame.

#### Step 3: Plot the significant TFBSs and their enrichment
- input
    - resulting data table
    - label A and label B from config.yaml
    - ouptput directory, output prefix, fdr threshold, min findings, plot size from config.yaml
```
example:
| TFBS_name | odds_ratio | adjusted_pvalue | active_overlaps | comparison_group |
|---|---|---|---|---|
| TF_A | 5.2 | 0.0001 | 150 | active_vs_inactive |
| TF_B | 3.1 | 0.002 | 85 | active_vs_inactive |
```
- rules
    - plot_TFBS_enrichment (scripts/plot_enriched_tfbs.py)
    - input