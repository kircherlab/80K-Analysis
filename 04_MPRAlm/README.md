# Directory for performing MPRAlm and bc_MPRAlm

## Outline
- Preparing data
- results and preparation script: `/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/bc_tradeoff`

### Preparing data: input for (bc_)MPRAlm
- File from Pia: `/data/gpfs-1/work/groups/ag_kircher/MPRA/ProximalPromoter/projects/tt_analysis/README.md`
# Malacoda quantification of proximal promoter  

Responsible: Pia

As part of the tiger team MPRA analysis.

Malacoda requires the counts per barcode, not per oligo. 
This folder only contains files to create the correct input for malacoda. All actual code and calculations are done on the Omics cluster in Lübeck.

/data/gpfs-1/work/groups/ag_kircher/MPRA/ProximalPromoter/projects/malacoda

## Creating input files

```sh
variants=/data/gpfs-1/work/groups/ag_kircher/MPRA/ProximalPromoter/experiments/second_approach/mprasnakeflow/results.first_results/experiments/run1+2+3HepG2/variants/standardRun1+2+3/standard/HepG2_variantTable.tsv.gz
assignment=/data/gpfs-1/work/groups/ag_kircher/MPRA/ProximalPromoter/experiments/second_approach/mprasnakeflow/results.first_results/experiments/run1+2+3HepG2/assignment/standardRun1+2+3.tsv.gz
count1=/data/gpfs-1/work/groups/ag_kircher/MPRA/ProximalPromoter/experiments/second_approach/mprasnakeflow/results.first_results/experiments/run1+2+3HepG2/assigned_counts/standardRun1+2+3/HepG2_1.merged.config.standard.tsv.gz
count2=/data/gpfs-1/work/groups/ag_kircher/MPRA/ProximalPromoter/experiments/second_approach/mprasnakeflow/results.first_results/experiments/run1+2+3HepG2/assigned_counts/standardRun1+2+3/HepG2_2.merged.config.standard.tsv.gz
count3=/data/gpfs-1/work/groups/ag_kircher/MPRA/ProximalPromoter/experiments/second_approach/mprasnakeflow/results.first_results/experiments/run1+2+3HepG2/assigned_counts/standardRun1+2+3/HepG2_3.merged.config.standard.tsv.gz

# assign counts to barcodes and oligo ids
join <(zcat $assignment) <(zcat $count1) -t $'\t' | join - <(zcat $count2) -t $'\t' | join - <(zcat $count3) -t $'\t' | cut -f1,2,5-11 | bgzip -c > HepG2_counts_per_bc.tsv.gz

# assign oligo ids to seqs and add "ref", "alt"
zcat $variants | awk 'NR > 1 && $10 > 0 && $17 > 0 {print $1 "\t" $2 "\tref\n" $1 "\t" $3 "\talt"}'  | bgzip -c > HepG2_seqs_to_oligos.tsv.gz

# sort before joining
zcat HepG2_seqs_to_oligos.tsv.gz | sort -k 2 -o HepG2_seqs_to_oligos.tsv
zcat HepG2_counts_per_bc.tsv.gz | sort -k 2 -o HepG2_counts_per_bc.tsv
bgzip HepG2_seqs_to_oligos.tsv
bgzip HepG2_counts_per_bc.tsv

# join counts with sequence oligos file
join -1 2 -2 2  <(zcat HepG2_seqs_to_oligos.tsv.gz) <(zcat HepG2_counts_per_bc.tsv.gz) -t $'\t' | cut -f 2- | bgzip -c > malacoda_input_HepG2.tsv.gz

# remove sequences that miss either ref or alt
sort -k 1 <(zcat malacoda_input_HepG2.tsv.gz) > malacoda_input_HepG2.tsv
bgzip malacoda_input_HepG2.tsv
zcat malacoda_input_HepG2.tsv.gz | awk 'NR >1{print $1","$2}' | uniq | awk -v FS="," '{print $1 "\t" $2}' | cut -f 1 |sort | uniq -c |  awk '$1 != 2' |  awk '{print $2}' > single_seqs.temp.tsv
join -v 1 <(zcat malacoda_input_HepG2.tsv.gz) single_seqs.temp.tsv -t $'\t' | bgzip -c > temp
mv temp malacoda_input_HepG2.tsv.gz

# attach header
header="variant_id\tallele\tbarcode\tDNA1\tRNA1\tDNA2\tRNA2\tDNA3\tRNA3"
echo -e $header | bgzip -c | cat - malacoda_input_HepG2.tsv.gz > temp
mv temp malacoda_input_HepG2.tsv.gz
```