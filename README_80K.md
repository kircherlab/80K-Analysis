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
    - config (aka StandardConfig (from Max)):
``` bash
snakemake --use-conda  --configfile config.yaml --snakefile /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/workflow/Snakefile --conda-prefix ../../MPRAsnakeflow_projects/conda --keep-going --cluster-config /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/config/sbatch.yml --cluster-status status.py --cluster "sbatch --parsable --nodes=1 --ntasks-per-node={cluster.threads} --mem {cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output} -e {cluster.error}"  --jobs 10 --cluster-cancel scancel -n --quiet
```
- config with lower thresholds on Barcode counts
```bash
snakemake --use-conda  --configfile low_config.yaml --snakefile /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/workflow/Snakefile --conda-prefix ../../MPRAsnakeflow_projects/conda --keep-going --cluster-config /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/config/sbatch.yml --cluster-status status.py --cluster "sbatch --parsable --nodes=1 --ntasks-per-node={cluster.threads} --mem {cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output} -e {cluster.error}"  --jobs 40 --cluster-cancel scancel -n --quiet
```

### Check if lower bc levels have an influence: all\_config\_test.yaml (lowConfig: for low bc counts)
- run MPRAsnakeflow with bowtie (only bowtieRun config (standard config but added bowtie call to the MPRAsnakeflow))
- Bowtie call `bowtie -x <index-base> -m <allow-multiple-alignments> --best (only the best alignement) --strata (only in the best strata)`
  - -x is not there; lets use -n alignment option
- (TODO:) Question: Am I interested in all alignments or only in the best alignment? -> best but (spoiler) the output was not usable from bowtie because barcodes were not streamlined
```bash 
snakemake --use-conda  --configfile bowtie_config.yaml --snakefile /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/workflow/Snakefile --conda-prefix ../../MPRAsnakeflow_projects/conda --keep-going --cluster-config /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/config/sbatch.yml --cluster-status status.py --cluster "sbatch --parsable --nodes=1 --ntasks-per-node={cluster.threads} --mem {cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output} -e {cluster.error}"  --jobs 40 --cluster-cancel scancel -n --quiet
```
- run again standard config to have a flawless set of data (checked out MPRAsnakeflow at standardMPRAsnakeflow, removed temp in rules (assignment_mapping, assignment_merge))
```bash
snakemake --use-conda  --configfile standard_config.yaml --snakefile /data/gpfs-1/users/kisa11_c/work/coding/standardMPRAsnakeflow/workflow/Snakefile --conda-prefix ../../MPRAsnakeflow_projects/conda --keep-going --cluster-config /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/config/sbatch.yml --cluster-status status.py --cluster "sbatch --parsable --nodes=1 --ntasks-per-node={cluster.threads} --mem {cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output} -e {cluster.error}"  --jobs 40 --cluster-cancel scancel -n --quiet
```

- testing bowtie code for split 17: `bowtie --best -x results/assignment/assignIGVFDesignNoTemp/reference/bowtie -q <(gzip -dc results/assignment/assignIGVFDesignNoTemp/fastq/merge_split17.join.fastq.gz) -S -p 4 | samtools sort -l 0 -@ 4 > bowtie_test_split_17.bam`
- Output of bowtie:
  - reads processed: 3803775
  - reads with at least one alignment: 3463631 (91.06%)
  - reads that failed to align: 340144 (8.94%)
  - Reported 3463631 alignments
  - [bam_sort_core] merging from 0 files and 4 in-memory blocks...
  - ...

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
- Start snakemake with standardMPRAsnakeflow directory:
```bash 
snakemake --use-conda  --configfile low_config.yaml --snakefile /data/gpfs-1/users/kisa11_c/work/coding/standardMPRAsnakeflow/workflow/Snake
file --conda-prefix ./conda --keep-going --cluster-config /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/config/sbatch.yml --cluster-status status.py --cluster "sbatch --parsable --
nodes=1 --ntasks-per-node={cluster.threads} --mem {cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output} -e {cluster.error}"  --jobs 40 --cluster-cancel scancel -n
```
- next problem: counts_umi_raw_counts test it locally Error in rule counts_umi_create_BAM: (solved by creating result directory)                                  
```bash 
jobid: 172
input: /fast/groups/ag_kircher/MPRA/IGVF_Y1_design/data/counts/Ngn2-DNA-1_S1_R1_001.fastq.gz, /fast/groups/ag_kircher/M$
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

### Finding: BWA new has less not found sequences then the initial run (spoiler: this was because I took the bam as comparison and didn't look at the subsequent filtering step)
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
- merge list of bams: `samtools merge -@ 8 -b /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/list_of_bams.tsv -o /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/bwa_merged.bam`
- index merged bam: `samtools index -b -@ 8 /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/bwa_merged.bam`
- idxstats of bam: `samtools idxstats /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/bwa_merged.bam > /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/idxstats_bwa_1201.tsv`
- bam: bowtie: `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesignNoTempBowtie/bam/mapped_bowtie.bam`
- bam: bowtie2: `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesignNoTempBowtie/bam/merged_bowtie2.bam`

- bam: 11G bwa-mem2: `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/bwa2/assignment/standardAssignIGVFDesignNoTemp/bam/aligned_merged_reads_bwa-mem2.bam` 
  - index: `samtools index -b -@ 8 /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/bwa2/assignment/standardAssignIGVFDesignNoTemp/bam/aligned_merged_reads_bwa-mem2.bam`
  - idxstats: `samtools idxstats /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/bwa2/assignment/standardAssignIGVFDesignNoTemp/bam/aligned_merged_reads_bwa-mem2.bam > /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/bwa2/assignment/standardAssignIGVFDesignNoTemp/bam/idxstats_bwa-mem2.tsv`
  - unmapped: `awk '$3==0 {print $0}' /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/bwa2/assignment/standardAssignIGVFDesignNoTemp/bam/idxstats_bwa-mem2.tsv | wc -l` # (925)
  - mapped: `awk '$3>0 {print $0}' /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/bwa2/assignment/standardAssignIGVFDesignNoTemp/bam/idxstats_bwa-mem2.tsv | wc -l` # 79291
  - => Result of bwa-mem2 is same as bwa 

### normal bwa mem + samtools view -F 1792
- filtering samtools 1792: `samtools view -F 1792 -b /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/bwa_merged.bam > /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/bwa_merged_view_1792_output.bam`
  - index the resulting bam and idxstats it:
    - `samtools index -b -@ 8 /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/bwa_merged_view_1792_output.bam`
    - `samtools idxstats /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/bwa_merged_view_1792_output.bam > /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/idxstats_bwa_merged_view_1792_output.tsv`

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
- check if all sequences are upper case: ` cat /fast/groups/ag_kircher/MPRA/IGVF_Y1_design/resources/association_data/design_no_duplicates_sequence_and_header.fa | grep -v ">" | awk '$0 !~ /[A-Z]/  {print $0}'`
- Run `80K_analysis/01_missing_sequences/notebooks/missing_sequences_in_bam.py` => alignments of missing (results/identified_missing_sequences/missing_sequences.bam) => test this with the quality control rule of MPRASnakeflow => the output should be empty
  - Output of the script: (command line)
    - found 1579132 reads with missing sequences
    - found 2639 different oligo_names with missing sequences: max: 2639
- currently: preparing a python script doing the same as the rule 
  - run rule call (can be found in `redo_assignment_getBC.py`) => No output was prepared 
  - run the samtools command in the command line => nothing found so all failed (was expected)
  - Write script which tests the condition and writes which fail at which quality measure (for merged_bam.bam from BWA)
    - quality: `/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/identified_missing_sequences/missing_sequence_min_quality.bam`
      - 4227071 reads (weird looking long alignments)
    - alignment start: `/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/identified_missing_sequences/missing_sequence_alignment_start.bam` 
      - 394 reads
    - alignment end: `/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/identified_missing_sequences/missing_sequences_alignment_end.bam`, 
      - 206245 reads
    - length min: `/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/identified_missing_sequences/missing_sequences_sequence_length_min.bam`
      - 2403378
    - length maximum: `/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/identified_missing_sequences/missing_sequences_sequence_length_max.bam`:
      - 0 reads
#### Parse Bam file because bwa returns mapping quality 0 sequences
- examples for multi mapping reads NOS3 in e.g. `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/bam/bwa_merged.bam`
- Plan: Postprocess the bam file (found duplicates)
  - check if unmapped: continue
  - check identity (if ~95% of sequence should be matched) with cigar and md
    - if mapping quality ($5) >= 1 
      - optional: expected length of the read (configuration) ~> 265
    - if mapping quality ($5) < 0:
      - if mapping on forward ($2 == 0)
        - check if AS > XS (check how many cases have difference of 5) -> keep and set mapping quality ($5) = 1 (is the best alignment)
      - if mapping with reversed read ($2 == 16) and strandedness (configuration) is set:
        - try recovering / finding the correct match
          - if AS == XS (the line should have a XA)
            - it can be recovered if only one of the alternatives at XA (separated by ";") can be matched on +strand (XA:Z:<oligo_name>,+/-<pos>,<cigar>,<number_differences>;)
            - if more than one alternative can be matched on + than we count and throw a warning of the overall number
  - Write all of this in a file with barcode \t reference_name \t position;cigarstring;NM:i:<number_of_differences>;MD:Z:<MD>;<number I do not understand>
- example code: `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/bin`
  - `removeSequenceErrors.py` - cigar found
- understood the MD parsing and cigar parsing
- wrote identity check function (computing number of matches from the cigar)
  - how to check the number of mismatches?
  - for one split:
    - unmapped count:  762
    - high_identity_count:  3791023
    - count:  3803677
    - proportion of high quality:  0.9966732190982568
  - weird examples: 
    - NB501960:812:HH53WAFX5:1:11309:25329:20408      0       cardiac_neuro_cava_random:REF_AHDC1|ENSG00000126705.15|EH38E1331564_rev_tile1-1 16      3       270M    *       0       0       TCTCTGGGCCTTGGTTTCCTTCCTCCTGCATAGCGAAGGAGGTTGGATTAGTGGCTCCCTAAGGCACCTTCTAGCTCTGACAGGCTCCAAGCCTGTGTTGACTGATGTGTCCTAGGAGATAGGCGCACACAGAGAACCAAGTCAGCTCCGAGAATCCTGTGAAGGTATCGCCACCCCACCCCCAGATGGCTGGAGTGCCTCCCTTCCTGAGACACACCCTTCATGGATACTGGTGGAGGTTGTGGTGGATGGAGGGGGCTTATCACCCAA  AAAAAEE/E/E<EEEEAEEE/AEEEEEEEEE/E/EEEAEEE/EEEEEAE/AEEE/E<EEE/E<E//AEEEA/AEAAA<E/EEE<EE/E///</EEEEE/EAAEE/AE<E//6<AE/AEEEAEEE:B"CHGGHHGEHHHHD<GGHFD<////<</A//////<A<AA<EA/AE/AA<</EAEAA/<E/EEE<EEEE/EA/A/A</EEEE/EE/E/AA/EEEAE6EE/EEAEEAA//EEEA/EEEEE/AEEEEEE/EEEEEE/AEEAAAAAA  NM:i:6  MD:Z:148G4T6C8C16A36G46 AS:i:240        XS:i:235        XA:Z:cardiac_neuro_cava_random:ALT_AHDC1|ENSG00000126705.15|EH38E1331564_rev_tile1-1_AHDC1|ENSG00000126705.15|EH38E1331564|1-27573227-A-G,+16,270M,7;cardiac_neuro_cava_random:ALT_AHDC1|ENSG00000126705.15|EH38E1331564_rev_tile1-1_AHDC1|ENSG00000126705.15|EH38E1331564|1-27573148-C-T,+16,270M,7;   XI:Z:ACGGGAGTCAGTATG,YI:Z://AAAAEEAEEE/EE
    - NB501960:812:HH53WAFX5:1:11105:15385:4842       0       cardiac_neuro_cava_random:NRXN3|ENSG00000021645.20|EH38E1730612_fwd_tile1-1     16      0       270M    *       0       0       CATGTGTGTGCGTGTGTAAGAGCAATAGATGGAGCAGCCAATGCAGGGAGAAGATCTGAAATTTGGAGCCTATTTAAGAGTGAAAGAACAGCTGCCTTTTAAATGGTTGCTAATTCTGACAATTAATGCTGTTCTGTGAGTGTGTGATTTTCTGTGATAGCAGTTAGGAGGGATGATGGTGTGTAATAGCTCATTTCCTAAGCTTTATGGTAATGTAGCATAGTCAAGACCATAGAACTGAAAGAGAACTGGGTCCCACATGCTCTGCAG  AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEAEEAEEEEEEEEEEEAEEEEEEEEEEEAEEEEEEEEEEEAEEEAEEEEEE<EEEEEGHHGGGGHHGHHHHGHHHGHHHEEEEAEEAEAAEAEAEA<EAAAA6AAA<6AAAEA<EAEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA  NM:i:0  MD:Z:270        AS:i:270        XS:i:268        XA:Z:GC_Vista:fb;mb__vistaElementControl|chr14:78308172-78308421,+18,270M,2;    XI:Z:CTGAACTATAATCAG,YI:Z:AAAAAEEEEEEEEEE
    - NB501960:812:HH53WAFX5:3:11405:4093:19168  16  cardiac_neuro_cava_random:REF_DRD4|ENSG00000069696.7|EH38E2937745_fwd_tile1-1  16  0  270M  *  0  0  TCTGGGAAGGGGGCCCCACAGGCAGCACCTACCGCAGGAGCTCAGTGTGAGCCACTGTCGGCCTGTGGGTGTGTGTCGTGTGTCGCAGCCAGCTCAGTTGCTGTCAGGATTCCACAGGCTGGGCAGCGTAAACCGCAGGTCTTCCTTTTCTTGGTTTTGGAAACTAGACATCTGAGATACCAGCAGGCCTGGTTCCTGGGGAGGCCCCCTTCCTGACTTACAGACGGCCGCCTCCCCGCTGTGACCTCACCTGGCCTTCCCTCATCTAAG  array('B', [32, ..., 32])  [('NM', 1), ('MD', '110A159'), ('AS', 265), ('XS', 265), ('XI', 'CTGAGAGTGACGTTG,YI:Z:AAAAAEEEEEEEEEE')]
    - I looked into examples of 
#### Generate the table of all variants and regions
- at 05_variant_region_list: variant_region_list.ipynb
  - writing variants, regions and references to files (smaller subsets for generating one table with meta data)
  - Numbers are 8900 regions, 18582 references and 46456 variants (more variants than expected, less regions than expected, no expectation)

- for checking if subsequent steps work with the new `assignment_barcodes_incl_other.sorted.tsv.gz` file: 
```bash
snakemake --use-conda  --configfile standard_config.yaml --snakefile /data/gpfs-1/users/kisa11_c/work/coding/standardMPRAsnakeflow/workflow/Snakefile --conda-prefix ../../MPRAsnakeflow_projects/conda --keep-going --cluster-config /data/gpfs-1/users/kisa11_c/work/coding/standardMPRAsnakeflow/config/sbatch.yml --cluster-status status.py --cluster "sbatch --parsable --nodes=1 --ntasks-per-node={cluster.threads} --mem {cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output} -e {cluster.error}"  --jobs 40 --cluster-cancel scancel -n --quiet
```
- it installs python3 packages

#### Run the new assignment process again (testing with low config) and check how many sequences get through the filter (if we trust the aligner)
- current results: we want to have high identity alignments
- now: we trust him but count how many times we would have thrown this read away
- run MPRAsnakeflow with low config (2 bc 1 dna 1 rna) and the new assignment file
```bash
snakemake --use-conda  --configfile low_config.yaml --snakefile /data/gpfs-1/users/kisa11_c/work/coding/standardMPRAsnakeflow/workflow/Snakefile --conda-prefix ../../MPRAsnakeflow_projects/conda --keep-going --cluster-config /data/gpfs-1/users/kisa11_c/work/coding/standardMPRAsnakeflow/config/sbatch.yml --cluster-status status.py --cluster "sbatch --parsable --nodes=1 --ntasks-per-node={cluster.threads} --mem {cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output} -e {cluster.error}"  --jobs 50 --cluster-cancel scancel -n
```
- 135 jobs
- check the log files
- got killed at 1 by slurm (thought it would take only 5h)
- identified
```bash
snakemake --use-conda  --configfile low_config.yaml --snakefile /data/gpfs-1/users/kisa11_c/work/coding/standardMPRAsnakeflow/workflow/Snakefile --conda-prefix ../../MPRAsnakeflow_projects/conda --keep-going --cluster-config /data/gpfs-1/users/kisa11_c/work/coding/standardMPRAsnakeflow/config/sbatch.yml --cluster-status status.py --cluster "sbatch --parsable --nodes=1 --ntasks-per-node={cluster.threads} --mem {cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output} -e {cluster.error}"  --jobs 50 --cluster-cancel scancel --rerun-incomplete -n
```
- 107 jobs
- `barcodes_incl_other.sorted.tsv.gz` with new assignment script: `93159939` (not matched barcodes are missing) => remove the BCs again and redo the workflow
- fixed script and added the barcode with "other" and "NA" if it cannot be assigned 
- How many reads do not have an XI tag? # started job with job id 4010688 slurm log in ~/slurm*4010688* => found error in debug script => restart: 4016253

### Work on variant region map (match manually with all the sequences)
- result can be found in `/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/04_MPRAlm/notebooks/investigate_80K_data.ipynb` 
  - variant_region_map_var_merged["sequence"].isna().sum() # 46374  (left_on="Variant") between header of design fasta and Variant column is no overlap
  - variant_region_map_region_merged["sequence"].isna().sum() # 46374  (left_on="Region") between header of design fasta and Region column is no overlap
  - variant_region_map_ref_merged["sequence"].isna().sum() # 778  (left_on="REF_ID") between header of design fasta and REF_ID column is overlap (not for 778 sequences)
  - variant_region_alt_merged["sequence"].isna().sum() # 1374  (left_on="ALT_ID") between header of design fasta and ALT_ID column is overlap (not for 1374 sequences)
- Investigate how it was done for 80K: MPRAOligodesign workflow / github repository 
  - rule producing: variant_region_map: https://github.com/kircherlab/MPRAOligoDesign/blob/master/workflow/rules/oligo_design.smk
  - script producing: variant_region_map: https://github.com/kircherlab/MPRAOligoDesign/blob/master/workflow/scripts/oligo_design/getSequencesInclVariants.py
  - Input of script:
    - regions: bed file 
    - variants: vcf file
    - reference/design: fasta file
    - outputs for: variants, variants-removed, regions, regions removed, design, design variant map
    - variant-edge-exclusion: What was given here?
  - What does this script do? 
  0. find config for MPRAOligoDesign
    - 5K_config: `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/final_design/design_5K.conf.yml`
  1. find vcf and bed file for 80K data
    - vcf: `/data/gpfs-1/groups/ag_kircher/MPRA/IGVF_Y1_design/design/final_design/input/variants_5K.vcf`
    - bed: `/data/gpfs-1/groups/ag_kircher/MPRA/IGVF_Y1_design/design/final_design/input/regions_5K.bed`
 
- Table for mapping headers: sanity checking with Mohan: `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/scripts/check_sanity_oligomap.ipynb`
- Changed the IDs in the variant_region_map to the headers of the new design file
  - script `/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/scripts/modify_variant_table.py`
- check now the results of (bc)MPRAlm
  - lowConfig results: `results/experiments/lowConfig_bwa/statistic/barcode/assigned_counts/assignmentFixDuplicates/NGN2_lowConfig_barcode_correlation.tsv`
  - generate new MPRAlm results with the new `variant_region_map` (`/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/80K_MPRA/design/variant_region_map_deduplicated.tsv.gz`)
  1. Current results new assignment (found error, old assignment was used)
  2. New variant region map (whats the influence on the variant numbers)
  3. New variant region map + lowConfig during MPRAsnakeflow (whats the influence on the region numbers)
    - run bc_preparation.sh (7min) (job id 4365763) (`sbatch --mem=20G --time=0-09:00 -c 2 --wrap="bash /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/bc_preparation.sh"`) (820 single sequences)
    `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/lowConfig_NGN2_filtered_counts_sequences.tsv.gz` (5827836 rows)
    - run preprocess: `/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/04_MPRAlm/bc_MPRAlm/scripts/preprocess.py` => e.g. `bc_tradeoff/results/preprocess/lowConfig_mpralm_input_NGN2.tsv` (1628686 rows)
    - run MPRAlm `/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/04_MPRAlm/MPRAlm/scripts/run_mpralm.r` => `MPRAlm/weights_lowConfig_NGN2.feather` (~25000 variants can be tested)
    - understood what the preprocessing does
    - results in `25411` variant ids
- Start the standard workflow again because wrong filtering was used:
```bash
snakemake --use-conda  --configfile standard_config.yaml --snakefile /data/gpfs-1/users/kisa11_c/work/coding/standardMPRAsnakeflow/workflow/Snakefile --conda-prefix ../../MPRAsnakeflow_projects/conda --keep-going --cluster-config /data/gpfs-1/users/kisa11_c/work/coding/standardMPRAsnakeflow/config/sbatch.yml --cluster-status status.py --cluster "sbatch --parsable --nodes=1 --ntasks-per-node={cluster.threads} --mem {cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output} -e {cluster.error}"  --jobs 40 --cluster-cancel scancel -n --quiet
```
- Run bc_preparation again: `sbatch --mem=20G --time=0-09:00 -c 2 --wrap="bash /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/bc_preparation.sh"`
  - job id: 4382367 (10 min) 1058 single sequences (not found ref and alt)
  - results in `25395` variant ids
- Checking the MPRAsnakeflow results: there is combined output (less preprocessing needed)
  - low config:
    - `/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/experiment/standard_results/results/experiments/lowConfig_bwa/assigned_counts/assignmentFixDuplicates/lowConfig/NGN2_2_barcode_assigned_counts.tsv.gz` (14518686 rows) (BC, assigned oligo, dna_count, rna_count)
    - `/data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/experiment/standard_results/results/experiments/lowConfig_bwa/assigned_counts/assignmentFixDuplicates/lowConfig/NGN2_allreps_merged.combined.tsv.gz` (14518686 rows) (BC, assigned oligo, dna_count, rna_count)
      - Question at hand: does this contain all information I need for variants and how many are included? (join with variant_region_map)
        - I have xx variants see presentation
- Running preprocessing: standard: `41284` variants with > 2 barcodes
- sanity check because shitty understanding of custom script: 
  - job id: `4538047` for lowConfig (bc_preparation.sh): 42598 
  - job id: `4538055` for standard (bc_preparation.sh): 42129
  - Preprocessing: standard: `4538254` => in table: `41284`
  - Preprocessing: lowConfig: `4538258` => in table: `42137`
  - generate new task: 06_variant_analysis
  - Questions: try to get the numbers of positive and negative effects
    - Add enformer classes (problems?)
#### metadata format in variant_region_list.ipynb
- output in `MPRA/IGVF_Y1_design/80K_MPRA/design`
- variant and sequence position is missing
### Check quality control of MPRAsnakeflow
- DNA plot less than RNA plot (`scripts/plot_perBCCounts_correlation.R`)
  - problems with on demand portal: tidyverse not installable and session crashes alot
  - Testing only with 2 replicates (all barcodes in merged files: e.g. `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/experiments/standard_bwa/assigned_counts/assignmentFixDuplicates/NGN2_1.merged.config.standardConfig.tsv.gz`)
```R
% DNA
correlation_plots <- cowplot::plot_grid(plotlist = plots_correlations_dna, ncol = 1)
ggsave("test_plot_corrlation_dna.png", correlation_plots, width = 15, height = 10 * length(plots_correlations_dna), dpi=96, type="cairo")
% rna:
correlation_plots_rna <- cowplot::plot_grid(plotlist = plots_correlations_rna, ncol = 1)
ggsave("test_plot_corrlation_rna.png", correlation_plots_rna, width = 15, height = 10 * length(plots_correlations_rna), dpi=96, type="cairo")
```
  - Checking for reasons why less dots are in the DNA barcode plot than in the RNA plot (data which is plotted in the script (`res_plot`) is 100000 rows long and has for each value (`DNA_normalized_log2.x`, `DNA_normalized_log2.y`) and (`RNA_normalized_log2.x`, `RNA_normalized_log2.y`) entries 
    - Possibility: more duplicates in DNA_normalized_log2 x and y values than in RNA
      - `sum(duplicated(res_plot[c("DNA_normalized_log2.x", "DNA_normalized_log2.y")]))` => 99622
      - `sum(duplicated(res_plot[c("RNA_normalized_log2.x", "RNA_normalized_log2.y")]))` => 98417
        - We plot 100000 values the DNA and RNA values differ in percent of duplicates: DNA: 99.622% and RNA: 98.417%
        - 1 - 0.99622 is the proportion of dots in the DNA plot
        - check the number of unique values `DNA_normalized_log2` and `DNA_normalized_log2`
          - `nrow(unique(res_plot[c("DNA_normalized_log2.x", "DNA_normalized_log2.y")]))` => 378 (for all data: 876)
          - `nrow(unique(res_plot[c("RNA_normalized_log2.x", "RNA_normalized_log2.y")]))` => 1583 (for all data: 5040)
  - in RNA_plot per barcode found error (?) (it fixes the error)
  ```R  
  min <- min(data$`RNA_normalized.x_log2`) # should be min(data$`RNA_normalized_log2.x`)
  max <- max(data$`RNA_normalized.y_log2`) # should be max(data$`RNA_normalized_log2.y`)
  ``` 
  - check the unique values of log DNA and RNA normalized log2 values for this data
```R
length(unique(plots_correlations_dna_list[[1]]$data$dna_normalized_log2.x)) # => 12922
length(unique(plots_correlations_dna_list[[1]]$data$dna_normalized_log2.y)) # => 12056
length(unique(plots_correlations_dna_list[[1]]$data$rna_normalized_log2.y)) # => 25786
length(unique(plots_correlations_dna_list[[1]]$data$rna_normalized_log2.x)) # => 25022
```
- have the warning for RNA plot: (it is because ggplot removes some values to have a nicer looking plot but all data is used to compute (summary) statistics according to [stack](https://stackoverflow.com/questions/32505298/explain-ggplot2-warning-removed-k-rows-containing-missing-values))
```
Warnmeldungen:
1: Removed 1533 rows containing missing values or values outside the scale range
(`geom_point()`). 
2: Removed 697 rows containing missing values or values outside the scale range
(`geom_point()`).
```
- tested all plots with different number of labels and looks better with lables underneath the plot
- resulting plots locally in documents folder (both plotting objects have 100000 elements but the DNA plot is not showing the same amount of barcodes)
- label problem for inserts
  - rule: statistic_correlation_calculate: `workflow/scripts/count/plot_perInsertCounts_correlation.R`
    - input: `results/experiments/{{project}}/assigned_counts/{{assignment}}/{{config}}/{{condition}}_{replicate}_merged_assigned_counts.tsv.gz` for replicate 1: 74691 rows (header: `name    dna_counts      rna_counts      dna_normalized  rna_normalized  ratio   log2    n_obs_bc`)
    - with legend position bottom plot looks better
    - increased size of plot: looks nice (stored locally) but increase of file size
    - For inserts the number is not different because of zero counts (Pia idea)

### Compute plots for 80K data
- run MPRAsnakeflow + Report 
- get everything together in excel file ([in sharepoint](https://charitede-my.sharepoint.com/:x:/r/personal/kilian_salomon_bih-charite_de/_layouts/15/Doc.aspx?sourcedoc=%7BD586D4C6-D692-4638-A51C-943F1BAECE8B%7D&file=sanity_check_NGN2.xlsx&action=default&mobileredirect=true))

### Investigating Controls
- Found directory with files for controls: `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/resources/controls`
  - README points you to the files in the folder 
    - Seem to be the original files
    - e.g. fasta folder: 
      - andrew_cardio_positive-controls.fasta => all spaces (" ") are replaced by ":" in the current version
      - went through all files and matched a subset to the original header => only andrew_cardio_positive-controls.fasta: spaces replaced by ":"
    - e.g. meta folder: 
      - bed: `in-house.bed` => no ID column (Not matchable?)
- SLEA controls: n=200 100 positive an 100 negative (tested in HepG2)
- control groups with variants: `grep "variant" '/home/kisa/coding/80K_MPRA/80K-Analysis/05_variant_region_list/metadata_example_1903.tsv' | grep "control" | cut -f 14 | sort | uniq | wc -l` => 8
  - C_positive_heart_CAD
  - C_positive_neuron_CD => has references without matching alternative and alternatives without matching references: Mohan said these are element variants as well
  - GC_Atrial_fib
  - GC_Kircher
  - GC_Liang
  - GC_Mendelian_variants
  - GC_Mohlke
  - GC_Selvarajan (199 ALT)
- C_positive_neuron_CD: According to Mohan: conversation gitter 20.03.2024 C_positive_neuron_CD: has no variant controls (we got 100 positive controls with highest effect in NGN2 from Chengyu)
  - Reference without alternative: 'C_positive_neuron_CD:n2_rs11876_C_T_ref_50::chr16:2038176-2038446-mean_ratio1.93' (looked for matching rsID)
  - Alternative without reference: 'C_positive_neuron_CD:p1_rs55985730_T_G_alt_50_T::chr7:128776855-128777125-mean_ratio2.34' (looked for matching rsID)
- haveing 6275 controls on the experiment
  - element inactive control    5019
  - variant negative control     967
  - element active control      198
  - variant positive control     91

#### Questions about controls:
- GC_Mendelian_variants: 
  - e.g.:
    - REF1: GC_Mendelian_variants:REF_chr10:23219434A*G|PTF1A, REF2: GC_Mendelian_variants:REF_chr10:23219376A*C|PTF1A 
    - ALT: GC_Mendelian_variants:ALT_chr10:23219434A*G|PTF1A_chr10:23219376A*C|PTF1A
