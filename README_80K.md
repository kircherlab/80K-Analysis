## Notes about analyzing 80K
- Max did the fist step of the analysis without a labeling file
- I did the second step: 
  - used the following command for the computing the labeling file and added it then to the config (located here)
```bash
#!/bin/bash

# Set the input and output file names
input_file="/fast/groups/ag_kircher/work/MPRA/IGVF_Y1_design/resources/association_data/design_no_duplicates.fa"
output_file="/fast/groups/ag_kircher/work/MPRA/IGVF_Y1_design/resources/association_data/design_no_duplicates_label.tsv"

# grep all headers and print them as tsv with the fist part as label
grep "^>" $input_file | sed 's/^>//g' | awk -F ":" '{print $1":"$2 "\t" $1}' > $output_file
```
  - used the following command for testing the snakemake workflow with the labeling file
``` bash
snakemake --use-conda  --configfile config.yaml --snakefile /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/workflow/Snakefile --conda-prefix ../../MPRAsnakeflow_projects/conda --keep-going --cluster-config /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/config/sbatch.yml --cluster-status status.py --cluster "sbatch --parsable --nodes=1 --ntasks-per-node={cluster.threads} --mem {cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output} -e {cluster.error}"  --jobs 10 --cluster-cancel scancel -n --quiet
```
```bash
snakemake --use-conda  --configfile low_config.yaml --snakefile /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/workflow/Snakefile --conda-prefix ../../MPRAsnakeflow_projects/conda --keep-going --cluster-config /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/config/sbatch.yml --cluster-status status.py --cluster "sbatch --parsable --nodes=1 --ntasks-per-node={cluster.threads} --mem {cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output} -e {cluster.error}"  --jobs 40 --cluster-cancel scancel -n --quiet```

Check if lower bc levels have an influence: all\_config\_test.yaml (lowConfig:
for low bc counts)

- run MPRAsnakeflow with bowtie (only bowtieRun config (standard config but added bowtie call to the MPRAsnakeflow))
- Bowtie call `bowtie -x <index-base> -m <allow-multiple-alignments> --best (only the best alignement) --strata (only in the best strata)`
  - -x is not there; lets use -n alignment option
- (TODO:) Question: Am I interested in all alignments or only in the best alignment?
```bash 
snakemake --use-conda  --configfile bowtie_config.yaml --snakefile /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/workflow/Snakefile --conda-prefix ../../MPRAsnakeflow_projects/conda --keep-going --cluster-config /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/config/sbatch.yml --cluster-status status.py --cluster "sbatch --parsable --nodes=1 --ntasks-per-node={cluster.threads} --mem {cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output} -e {cluster.error}"  --jobs 40 --cluster-cancel scancel -n --quiet
```
- run again standard config to have a flawless set of data (checked out MPRAsnakeflow at standardMPRAsnakeflow, removed temp in rules (assignment_mapping, assignment_merge))
```bash
snakemake --use-conda  --configfile standard_config.yaml --snakefile /data/gpfs-1/users/kisa11_c/work/coding/standardMPRAsnakeflow/workflow/Snakefile --conda-prefix ../../MPRAsnakeflow_projects/conda --keep-going --cluster-config /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/config/sbatch.yml --cluster-status status.py --cluster "sbatch --parsable --nodes=1 --ntasks-per-node={cluster.threads} --mem {cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output} -e {cluster.error}"  --jobs 40 --cluster-cancel scancel -n --quiet
```

- testing bowtie code for split 17: `bowtie --best -x results/assignment/assignIGVFDesignNoTemp/reference/bowtie -q <(gzip -dc results/assignment/assignIGVFDesignNoTemp/fastq/merge_split17.join.fastq.gz) -S -p 4 | samtools sort -l 0 -@ 4 > bowtie_test_split_17.bam`
- # reads processed: 3803775
- # reads with at least one alignment: 3463631 (91.06%)
- # reads that failed to align: 340144 (8.94%)
- Reported 3463631 alignments
- [bam_sort_core] merging from 0 files and 4 in-memory blocks...

### Bowtie + bowtie2 run:
- Error in rule assigned_counts_filterAssignment:
    - jobid: 232
    - input: results/assignment/assignIGVFDesignNoTempBowtie/assignment_barcodes.standardConfig.sorted.tsv.gz, /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/workflow/rules/../scripts/count/samplerer_assignment.py
    - output: results/experiments/run1counts_run2Assignment_NoDupAss_bowtie/assignment/assignmentFixDuplicates.tsv.gz
    - log: results/logs/assigned_counts/filterAssignment.run1counts_run2Assignment_NoDupAss_bowtie.assignmentFixDuplicates.log (check log file(s) for error details)
    - conda-env: /data/gpfs-1/work/groups/ag_kircher/MPRA/MPRAsnakeflow_projects/conda/ada89664f6498dde80c90c8465b7875c_
    - shell: python /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/workflow/rules/../scripts/count/samplerer_assignment.py         --input results/assignment/assignIGVFDesignNoTempBowtie/assignment_barcodes.standardConfig.sorted.tsv.gz                                    --output results/experiments/run1counts_run2Assignment_NoDupAss_bowtie/assignment/assignmentFixDuplicates.tsv.gz &> results/logs/assigned_counts/filterAssignment.run1counts_run2Assignment_NoDupAss_bowtie.assignmentFixDuplicates.log
    - (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)
    - cluster_jobid: 2988001
    - log said: empty input data 
    - `results/assignment/assignIGVFDesignNoTempBowtie/assignment_barcodes.standardConfig.sorted.tsv.gz` is empty
    - previous file: from rule has rows 556349 (`0	other	NA` but also other `1	C_SLEA:SLEA_hg18:chr2:210861483-210861650|103:V_HNF4_Q6:AAGGTCCAG;155:V_HNF4_Q6:AAGGTCCAG	16;268M;XA:i:0;MD:Z:266A0A0;255`)
	- compare it to the standard results: current state not working ()
	- run it locally in order to understand the error (on MPRA_files bottom)
		- results/assignment/assignIGVFDesignNoTempBowtie/barcodes_incl_other.sorted.tsv.gz (number of rows: 114104369)
  - `/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesignNoTempBowtie/bam/`
    - `merge_split13.mapped.bam` (0-29) all 1.9GB 
    - NOTE: perform snakemake workflow using the output of one rule: `snakemake ... -R <rule_name>`
    - call for summary of mappings on the reference sequences `samtools idxstats merge_split0.mapped.bam`
    - call view bam: `samtools view merge_split0.mapped.bam | less`
    - combine the bam files of the splits `cd results/assignment/assignIGVFDesignNoTempBowtie/bam/ && samtools merge -@ 38 -o merged.bam merge_split0.mapped.bam merge_split12.mapped.bam merge_split16.mapped.bam merge_split2.mapped.bam merge_split23.mapped.bam merge_split27.mapped.bam merge_split4.mapped.bam merge_split8.mapped.bam merge_split1.mapped.bam merge_split13.mapped.bam merge_split17.mapped.bam merge_split20.mapped.bam merge_split24.mapped.bam merge_split28.mapped.bam merge_split5.mapped.bam merge_split9.mapped.bam merge_split10.mapped.bam merge_split14.mapped.bam merge_split18.mapped.bam merge_split21.mapped.bam merge_split25.mapped.bam merge_split29.mapped.bam merge_split6.mapped.bam merge_split11.mapped.bam merge_split15.mapped.bam merge_split19.mapped.bam merge_split22.mapped.bam merge_split26.mapped.bam merge_split3.mapped.bam merge_split7.mapped.bam`
    - create idx: `cd results/assignment/assignIGVFDesignNoTempBowtie/bam/ && samtools index -b -@ 38 merged.bam`
    - use samtools idxstats and investigate how many reads are found `cd results/assignment/assignIGVFDesignNoTempBowtie/bam/ && samtools idxstats merged.bam > idxstats_bowtie.tsv`
    - number of idxstat_bowtie.tsv: 80216 (`cat idxstats_bowtie.tsv | wc -l`) 
    - number of idxstat_bowtie.tsv: $3 (number of hits >10: 77527) (`awk '$3>10 {print $0}' idxstats_bowtie.tsv | wc -l`)
    - number of idxstat_bowtie.tsv: $3 (number of hits >0: 78613) (`awk '$3>0 {print $0}' idxstats_bowtie.tsv | wc -l`)
    - compare to bwa results in standard_results/: 
      - get all files to file: `find ./  -printf "%f\n" | grep "mapped.bam" > list_of_bams.tsv`
      - `samtools merge -@ 38 -o merged_bwa.bam -b list_of_bams.tsv` => 11G bam file
      - `ls -lah`
      - `samtools index -b -@ 38 merged_bwa.bam`
      - `samtools idxstats merged_bwa.bam > idxstats_bwa.tsv`
      - `awk '$3>0 {print $0}' idxstats_bwa.tsv | wc -l`
- run bowtie2 rule:
  - generate: `find ./  -printf "%f\n" | grep "mapped.bowtie2.bam" > list_of_bowtie2_bams.tsv`
  - merge bams: `samtools merge -@ 38 -o merged_bowtie2.bam -b list_of_bowtie2_bams.tsv`
  - index bam: `samtools index -b -@ 38 merged_bowtie2.bam`
  - idxstat: `samtools idxstats merged_bowtie2.bam > idxstats_bowtie2.tsv`
  - `awk '$3>0 {print $0}' idxstats_bowtie2.tsv | wc -l`
### Configs and standard_results/
- standard_config
- first_bowtie_config: normal config in prior snakemake results (with bowtie as mapper)
- bowtie_config: new experiment and assignment with bowtie
- Start snakemake with standardMPRAsnakeflow directory: `snakemake --use-conda  --configfile low_config.yaml --snakefile /data/gpfs-1/users/kisa11_c/work/coding/standardMPRAsnakeflow/workflow/Snake
file --conda-prefix ./conda --keep-going --cluster-config /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/config/sbatch.yml --cluster-status status.py --cluster "sbatch --parsable --
nodes=1 --ntasks-per-node={cluster.threads} --mem {cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output} -e {cluster.error}"  --jobs 40 --cluster-cancel scancel -n`

- next problem: counts_umi_raw_counts test it locally
Error in rule counts_umi_create_BAM:                                                                                            ```jobid: 172                                                                                                                  input: /fast/groups/ag_kircher/MPRA/IGVF_Y1_design/data/counts/Ngn2-DNA-1_S1_R1_001.fastq.gz, /fast/groups/ag_kircher/M$
RA/IGVF_Y1_design/data/counts/Ngn2-DNA-1_S1_R3_001.fastq.gz, /fast/groups/ag_kircher/MPRA/IGVF_Y1_design/data/counts/Ngn2-D$
A-1_S1_R2_001.fastq.gz, /data/gpfs-1/users/kisa11_c/work/coding/standardMPRAsnakeflow/workflow/rules/../scripts/count/FastQ$
doubleIndexBAM.py, /data/gpfs-1/users/kisa11_c/work/coding/standardMPRAsnakeflow/workflow/rules/../scripts/count/MergeTrimR$adsBAM.py
    output: results/experiments/run1counts_run2Assignment_NoDupAss/counts/useUMI.NGN2_1_DNA.bam
    log: results/logs/counts/umi/create_BAM.run1counts_run2Assignment_NoDupAss.NGN2.1.DNA.log (check log file(s) for error $etails)
    conda-env: /data/gpfs-1/work/groups/ag_kircher/MPRA/IGVF_Y1_design/experiment/standard_results/conda/a5afef128c1a485ea0$d3d7e3c7a203a_
    shell:                                                                                                                               set +o pipefail;

        fwd_length=`zcat /fast/groups/ag_kircher/MPRA/IGVF_Y1_design/data/counts/Ngn2-DNA-1_S1_R1_001.fastq.gz | head -2 | tail -1 | wc -c`;
        fwd_length=$(expr $(($fwd_length-1)));

        rev_start=$(expr $(($fwd_length+1)));

        minoverlap=`echo ${fwd_length} ${fwd_length} 15 | awk '{print ($1+$2-$3-1 < 11) ? $1+$2-$3-1 : 11}'`;

        echo $rev_start >> results/logs/counts/umi/create_BAM.run1counts_run2Assignment_NoDupAss.NGN2.1.DNA.log
        echo $minoverlap >> results/logs/counts/umi/create_BAM.run1counts_run2Assignment_NoDupAss.NGN2.1.DNA.log

        paste <( zcat /fast/groups/ag_kircher/MPRA/IGVF_Y1_design/data/counts/Ngn2-DNA-1_S1_R1_001.fastq.gz ) <( zcat /fast/groups/ag_kircher/MPRA/IGVF_Y1_design/data/counts/Ngn2-DNA-1_S1_R3_001.fastq.gz  ) <( zcat /fast/groups/ag_kircher/MPRA/IGVF_Y1_design/data/counts/Ngn2-DNA-1_S1_R2_001.fastq.gz ) |         awk '{if (NR % 4 == 2 || NR % 4 == 0) {print $1$2$3} else {print $1}}' |         python /data/gpfs-1/users/kisa11_c/work/coding/standardMPRAsnakeflow/workflow/rules/../scripts/count/FastQ2doubleIndexBAM.py -p -s $rev_start -l 0 -m 17 --RG NGN2_1_DNA |         python /data/gpfs-1/users/kisa11_c/work/coding/standardMPRAsnakeflow/workflow/rules/../scripts/count/MergeTrimReadsBAM.py --FirstReadChimeraFilter '' --adapterFirstRead '' --adapterSecondRead '' -p --mergeoverlap --minoverlap $minoverlap > results/experiments/run1counts_run2Assignment_NoDupAss/counts/useUMI.NGN2_1_DNA.bam 2>> results/logs/counts/umi/create_BAM.run1counts_run2Assignment_NoDupAss.NGN2.1.DNA.log

        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)
    cluster_jobid: 2987694
```
- added mkdir -p results/experiments/{params.project}/counts/ to the shell command and added project to params with lambda
- error in the end because concatenation of the files is not working (only one file is created)
- investigate the results:
  - use missing_sequenes_per_label.ipynb: create again header dataframe: `cat results/assignment/{assignment}/reference/reference.fa | grep ">" | awk '{print substr($0,2)}' > results/assignment/{assignment}/reference/all_headers.tsv` (`cat reference.fa | grep ">" | awk '{print substr($0,2)}' > all_headers.tsv`)
  - problem bowtie does not include the bc information in the sam file => there is a flag for bowtie2 

### Finding: BWA new has less not found sequences then the initial run
- try to regenerate the 5084 sequences from the only run
  - directory: `/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesign/bam`
  - merge the bam (`samtools merge -@ 3 -o merged_bwa.bam merge_split*.mapped.bam`)
  - index the bam (`samtools index -b -@ 30 merged_bwa.bam`)
  - idxstats the bam (`samtools idxstats merged_bwa.bam > idxstats_bwa_old.tsv`)
  - `awk '$3==0 {print $0}' idxstats_bwa_old.tsv | wc -l` => 1603 (probably old bowtie results)
  - run mit mehreren alignment: results/all_alignments/standardAssignIGVFDesignNoTemp/bam/aligned_merged_reads.bam => 387 unmatched sequences

#### Try to redo the missing sequence finding (find_missing_sequences.py)
- bwa: `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesignNoTemp/assignment_barcodes.standardConfig.sorted.tsv.gz` 5084 sequences missing
  - `zcat /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesignNoTemp/barcodes_incl_other.sorted.tsv.gz | wc -l` => 
- bowtie: `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesignNoTempBowtie/assignment_barcodes.standardConfig.sorted.tsv.gz` 2325 sequences missing
  - barcodes_incl_other.sorted.tsv.gz: 114104369 rows
- Problem seem to be the way from bam to assignment_barcodes (intermediate barcodes_incl_other): currently: `samtools view -F 1792 {input}` -F <FLAG> (exclude if a bit is correct for the element at hand) 1792 (not primary alignment, read fails platform/quality check/read is PCR or optical duplicate) `https://broadinstitute.github.io/picard/explain-flags.html`

### check the barcode_others step for bwa and bowtie
- where is the (merged) bam file?
- `samtools merge -@ 8 -b /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/list_of_bams.tsv -o /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/bwa_merged.bam`
  - `samtools index -b -@ 8 /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/bwa_merged.bam`
  - samtools idxstats /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/bwa_merged.bam > /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/idxstats_bwa_1201.tsv

- bam: bowtie: `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesignNoTempBowtie/bam/mapped_bowtie.bam`
- bam: bowtie2: `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesignNoTempBowtie/bam/merged_bowtie2.bam`

- bam: 11G bwa-mem2: `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/bwa2/assignment/standardAssignIGVFDesignNoTemp/bam/aligned_merged_reads_bwa-mem2.bam` 
  - index: samtools index -b -@ 8 /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/bwa2/assignment/standardAssignIGVFDesignNoTemp/bam/aligned_merged_reads_bwa-mem2.bam
  - idxstats: samtools idxstats /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/bwa2/assignment/standardAssignIGVFDesignNoTemp/bam/aligned_merged_reads_bwa-mem2.bam > /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/bwa2/assignment/standardAssignIGVFDesignNoTemp/bam/idxstats_bwa-mem2.tsv
  - unmapped: awk '$3==0 {print $0}' /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/bwa2/assignment/standardAssignIGVFDesignNoTemp/bam/idxstats_bwa-mem2.tsv | wc -l # (925)
  - mapped: awk '$3>0 {print $0}' /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/bwa2/assignment/standardAssignIGVFDesignNoTemp/bam/idxstats_bwa-mem2.tsv | wc -l # 79291
  # => Result of bwa-mem2 is same as bwa 

# normal bwa mem + samtools view -F 1792
- samtools view -F 1792 -b /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/bwa_merged.bam > /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/bwa_merged_view_1792_output.bam
  - index the resulting bam and idxstats it:
    - samtools index -b -@ 8 /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/bwa_merged_view_1792_output.bam
    - samtools idxstats /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/bwa_merged_view_1792_output.bam > /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/idxstats_bwa_merged_view_1792_output.tsv

### Investigating difference of barcode assignment
- only filtering bwa result with samtools view 1792: `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/barcodes_incl_other.bwa_view_1792.tsv`
  - ~16GB (114104369 lines)
- Filtering with hard limits: (see rule `assignment_getBCs`) `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/barcodes_incl_other.sorted.tsv.gz`
  - ~600MB (2665005 lines) 
#### Find good examples of missing sequences with exact matches
- at the moment there is no cardiac_neuro_cava_random sequence in the barcode_bwa (which is only filtered by 1792) => this filter might be the problem? How can we test that? 
  - Check the bam file if it uncludes the exact matching reads (do they have a valid alignment? to which sequence? Are there sequences which are duplicated?)
  - Check if `/home/kisa/coding/80K_MPRA/80K-Analysis/05_variant_region_list/resources/design_no_duplicates_sequence_and_header.fa` has duplicated sequences
    - rows starting with ">" == rows not starting with ">" => header and sequence are alternating (Yes all sequences are unique)
#### Find them in the bam file and check which quality measure they fail
- TODO!!
#### Generate the table of all variants and regions
- at 05_variant_region_list: variant_region_list.ipynb
  - writing variants, regions and references to files (smaller subsets for generating one table with meta data)
  - Numbers are 8900 regions, 18582 references and 46456 variants (more variants than expected, less regions than expected, no expectation)

