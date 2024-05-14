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

### Combine enformer class files (information from the prioritized variants)
- enformer files: `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k`
- `merged_enformer_all_max_value.tsv`
- following code leads to `/home/kisa/coding/80K_MPRA/server_results/enformer_results/enformer_class/merged_enformer_class.tsv`
- first remove unused columns of the enformer file: results in: `{gene_type}.{variant_type}_680_columns.tsv`
- merge max_value enformer: 
```bash
first=1
for f in /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/cardiac/enformer/cardiac.ultra-rare_max_values.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/cardiac/enformer/cardiac.singleton_max_values.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/cava/enformer/cava.ultra-rare_max_values.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/cava/enformer/cava.singleton_max_values.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/neuro/enformer/neuro.ultra-rare_max_values.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/neuro/enformer/neuro.singleton_max_values.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/random/enformer/random.ultra-rare_max_values.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/random/enformer/random.singleton_max_values.tsv

do
    if [ "$first" ]
    then
        zcat "$f"
        first=
    else
        zcat "$f" | tail -n +2
    fi
done > /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/merging_enformer_info/merged_enformer_all_max_value.tsv
```
- this was the old version locally (did the sampling again myself but no seed was given => could not reproduce the sampling like mohan did)
```bash
first=1
for f in /home/kisa/coding/80K_MPRA/server_results/enformer_results/enformer_class/cava.singleton_enformer_class.tsv /home/kisa/coding/80K_MPRA/server_results/enformer_results/enformer_class/random.ultra-rare_enformer_class.tsv /home/kisa/coding/80K_MPRA/server_results/enformer_results/enformer_class/cardiac.singleton_enformer_class.tsv /home/kisa/coding/80K_MPRA/server_results/enformer_results/enformer_class/neuro.ultra-rare_enformer_class.tsv /home/kisa/coding/80K_MPRA/server_results/enformer_results/enformer_class/neuro.singleton_enformer_class.tsv /home/kisa/coding/80K_MPRA/server_results/enformer_results/enformer_class/cava.ultra-rare_enformer_class.tsv /home/kisa/coding/80K_MPRA/server_results/enformer_results/enformer_class/cardiac.ultra-rare_enformer_class.tsv
do
    if [ "$first" ]
    then
        cat "$f"
        first=
    else
        cat "$f" | tail -n +2
    fi
done > /home/kisa/coding/80K_MPRA/server_results/enformer_results/enformer_class/merged_enformer_class.tsv
```

- combine all variants 

```bash
first=1
for f in /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/merging_variants/cardiac.common_variants.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/merging_variants/cardiac.rare_variants.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/merging_variants/cardiac.ultra-rare_variants.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/merging_variants/cardiac.singleton_variants.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/merging_variants/cava.common_variants.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/merging_variants/cava.rare_variants.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/merging_variants/cava.ultra-rare_variants.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/merging_variants/cava.singleton_variants.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/merging_variants/neuro.common_variants.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/merging_variants/neuro.rare_variants.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/merging_variants/neuro.ultra-rare_variants.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/merging_variants/neuro.singleton_variants.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/merging_variants/random.common_variants.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/merging_variants/random.rare_variants.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/merging_variants/random.ultra-rare_variants.tsv /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/merging_variants/random.singleton_variants.tsv
do
    if [ "$first" ]
    then
        cat "$f"
        first=
    else
        cat "$f" | tail -n +2
    fi
done > /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/merging_variants/merged_prioritized_variants.tsv
```
- number of variants: 71410 (what to do with elements which are in two groups e.g. neuro and cava)? Two entries at the moment
- 13.05: new way to combine the variants of common and enformer prioritized into one file with unique ids: merging of enformer variant information into a list of tuples (gene_set, variant_type, enformer_class)
- Idea merge to the variant region map by cutting the group

##### Investigate enformer merged file and how to connect the ID to the oligos
- found duplicates in the variant region list: Max: 06.05.2024 said that this is because variants are connected to different genes
- For the workflow remove one pair of the ref and alt from the fasta file


#### Get sequences of variants (1KB around the variant) for enformer predictions
- how many variants do I expect?
- cat /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/resources/association_data/design_no_duplicates_sequence_and_header.fa | grep -E 'ALT_|REF_' | grep "~" | wc -l
- wanted to add region information for the variant sequences: for the tested sequences this was simple and for the variant controls (`GC_Mendelian_variants` and `C_positive_heart_CAD`) this is not possible with the region information in `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/final_design/results/final_design` because no region information for these control groups exists
- solution: create region file first with blat (see "Regions of controls" below)

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
  - info                   class                     number of variants   after assignment (with elements)
  - C_positive_heart_CAD   variant negative control     97                    92
  - C_positive_neuron_CD   variant positive control     91                    88 (has references without matching alternative and alternatives without matching references: Mohan said these are element variants conversation: 20.03)
  - GC_Atrial_fib          variant negative control     44                    45
  - GC_Kircher             variant negative control    203                    185
  - GC_Liang               variant negative control     16                    16
  - GC_Mendelian_variants  variant negative control    209                    185
  - GC_Mohlke              variant negative control     34                    31
  - GC_Selvarajan          variant negative control    364                    356 (364 in design and 356 after assignment with elements (199 ALT in design and 194 in assignment (standard_bwa))
#### Additional notes about controls
 C_positive_neuron_CD: According to Mohan: conversation gitter 20.03.2024 C_positive_neuron_CD: has no variant controls (we got 100 positive controls with highest effect in NGN2 from Chengyu)
  - Reference without alternative: 'C_positive_neuron_CD:n2_rs11876_C_T_ref_50::chr16:2038176-2038446-mean_ratio1.93' (looked for matching rsID)
  - Alternative without reference: 'C_positive_neuron_CD:p1_rs55985730_T_G_alt_50_T::chr7:128776855-128777125-mean_ratio2.34' (looked for matching rsID)
- haveing 6275 controls on the experiment
  - element inactive control    5019
  - variant negative control     967
  - element active control      198
  - variant positive control     91
- many variant control get lost on the way to MPRAlm (filtering for barcodes per allele)
  - before filtering for bc_MPRAlm:
    - C_positive_heart_CAD	49
    - GC_Atrial_fib	23
    - GC_Kircher	198
    - GC_Liang	8
    - GC_Mendelian_variants	174
    - GC_Mohlke	20
	  - GC_Selvarajan	198
  - after filtering:
    - C_positive_heart_CAD     43
    - GC_Atrial_fib            17
    - GC_Liang                  7
    - GC_Mohlke                11
    - GC_Selvarajan           112

- GC_Liang: eQTLs from neural 8 variants: 7 matching the ref/alt >= 2 barcodes
  - all non-significant
- GC_Mohlke: Hep and NonHep controls
- GC_Selvarajan: HepG2 controls
- GC_Atrial_fib: Atrial fibroblast controls
- during adding variant controls to bc_MPRAlm matching problem with headers again => check and fix it => mendelian variants and gc kircher variants in bc_MPRAlm
#### Questions about controls:
- _headerDuplicate: 2 from Mendelian_variants and 1 Mohlke (what is the reason?)
- looking for the reason of variant_ids: rs1218584|KCNN3|STARR-seq-AF~rs76749863|KCNN3|STARR-seq-AF_fwd_tile1-1
  - I find alt with multiple rsIDs (what does this mean?) GC_Atrial_fib:ALT_rs1218584|KCNN3|STARR-seq-AF~rs76749863|KCNN3|STARR-seq-AF_fwd_tile1-1_rs76749863
- GC_Mendelian_variants: 
  - e.g.:
    - REF1: GC_Mendelian_variants:REF_chr10:23219434A*G|PTF1A, REF2: GC_Mendelian_variants:REF_chr10:23219376A*C|PTF1A 
    - ALT: GC_Mendelian_variants:ALT_chr10:23219434A*G|PTF1A_chr10:23219376A*C|PTF1A
  - in current design file: changed A>C to A*C
  - problems with: GC_Mendelian_variants:REF_chr7:156791472C*T|SHH: seems like only ALTs are there (no REFs)
- GC_Selvarajan: has one variant which does not match a reference
  - `GC_Selvarajan:ALT_rs2297787|STARR-seq-HepG2_fwd_tile1-1_rs2297787`
- Is this a reference? Which variant is associated? 
  - reference: `>C_negative_neuron_MK:rdhs_475428_chr3_148289221_148289490_reference__0.420408294159496`
  - alt: `>C_negative_neuron_MK:tile_22397_chr2_103176684_103176953_A_T_257__0.418330976255658`
- `GC_Mohlke:REF_NC000001.11|230158967|C|A|MohlkeHepControls,NC000001.11|230159168|C|T|MohlkeHepControls,NC000001.11|230159329|CTTAAAGTGTTCAGCACTCCCCT|CT|MohlkeHepControls_fwd_tile3-3` is a duplicated header in the variant region map (REF_ID) => cannot be matched

- check matchable headers with modified headers
  - ref: 
    - "GC_Mohlke:REF_NC000001.11|230158967|C|A|MohlkeHepControls,NC000001.11|230159168|C|T|MohlkeHepControls,NC000001.11|230159329|CTTAAAGTGTTCAGCACTCCCCT|CT|MohlkeHepControls_fwd_tile1-3" - tiling what does this mean?
    - GC_Mendelian_variants:ALT_chr7:156791472C>T|SHH_chr7:156791413A>C|SHH (no GC_Mendelian_variants:ALT_chr7:156791472C>T)
    - GC_Mohlke:ALT_NC000010.11|100315721|G|A|MohlkeNonHepControl_fwd_tile1-1_NC000010_11_100315721_G_A
  - Groups
    - From variant_maps
       - label
        - GC_Selvarajan            198
        - GC_Kircher               198
        - GC_Mendelian_variants    174
        - C_positive_heart_CAD      49
        - GC_Atrial_fib             23
        - GC_Mohlke                 20
        - GC_Liang                   8
    - matchable
      - ref 
        - label
          - GC_Selvarajan            198
          - GC_Kircher               198
          - GC_Mendelian_variants    163
          - C_positive_heart_CAD      48
          - GC_Atrial_fib             23
          - GC_Mohlke                 18
          - GC_Liang                   8
      - alt 
        - label
          - GC_Selvarajan            198
          - GC_Kircher               198
          - GC_Mendelian_variants    161
          - C_positive_heart_CAD      48
          - GC_Atrial_fib             23
          - GC_Mohlke                 17
          - GC_Liang                   8
    - final numbers from current design 16042024
      - label
        - cardiac_neuro_cava_random    46374
        - GC_Selvarajan                  198
        - GC_Kircher                     198
        - GC_Mendelian_variants          161
        - C_positive_heart_CAD            48
        - GC_Atrial_fib                   23
        - GC_Mohlke                       17
        - GC_Liang                         8
  - identifying discrepance between number of headers before matching in the old variant map and after matching with the design file (17 not matching) 4 references
    - expectation: filtering of these sequences resulted in removing many sequences
    - manually checked examples
      - not in design anymore
        - GC_Mohlke:ALT_NC000001.11|230158967|C|A|MohlkeHepControls,NC000001.11|230159168|C|T|MohlkeHepControls,NC000001.11|230159329|CTTAAAGTGTTCAGCACTCCCCT|CT|MohlkeHepControls_fwd_tile1-3_NC000001_11_230158967_C_A
        - 'C_positive_heart_CAD:ALT_rs12721051_rs12721051' not there anymore but same rsid is in another group (Selvarajan)
        - GC_Mohlke:ALT_NC000010.11|100315721|G|A|MohlkeNonHepControl_fwd_tile1-1_NC000010_11_100315721_G_A (gone because ref is gone)
      - 'GC_Mendelian_variants:ALT_chr7:156791472C>T not in the design anymore (ref gone)
        - 'GC_Mendelian_variants:ALT_chr7:156791472C>T|SHH_chr7:156791413A>C|SHH',
        - 'GC_Mendelian_variants:ALT_chr7:156791472C>T|SHH_chr7:156791459T>C|SHH',
        - 'GC_Mendelian_variants:ALT_chr7:156791472C>T|SHH_chr7:156791472C>G|SHH',
        - 'GC_Mendelian_variants:ALT_chr7:156791472C>T|SHH_chr7:156791472C>T|SHH',
        - 'GC_Mendelian_variants:ALT_chr7:156791472C>T|SHH_chr7:156791474G>A|SHH',
        - 'GC_Mendelian_variants:ALT_chr7:156791472C>T|SHH_chr7:156791480G>A|SHH',
        - 'GC_Mendelian_variants:ALT_chr7:156791472C>T|SHH_chr7:156791542A>C|SHH',
        - 'GC_Mendelian_variants:ALT_chr7:156791472C>T|SHH_chr7:156791547A>G|SHH',
        - 'GC_Mendelian_variants:ALT_chr7:156791472C>T|SHH_chr7:156791571T>A|SHH',
        - 'GC_Mendelian_variants:ALT_chr7:156791472C>T|SHH_chr7:156791579C>T|SHH',
        - 'GC_Mendelian_variants:ALT_chr7:156791472C>T|SHH_chr7:156791581A>G|SHH',
      - Kombination of two alts and one is not a SNV?
        - GC_Mendelian_variants:ALT_chr8:11703860G>T|GATA4_chr8:11703890AG>A|GATA4 (Was soll das heißen??? Alt kombiniert mit anderem alt?)
        - GC_Mendelian_variants:ALT_chr8:11703890AG>A|GATA4_chr8:11703890AG>A|GATA4
  - in new variant region map of controls: 653 control variants
    - label
      - GC_Selvarajan            198
      - GC_Kircher               198
      - GC_Mendelian_variants    161
      - C_positive_heart_CAD      48
      - GC_Atrial_fib             23
      - GC_Mohlke                 17
      - GC_Liang                   8
  - controls have non unique variant ids: (Problem I dont't know how the variant Ids are computed)
    - GC_Selvarajan Group has non unique variant ids: (all, unique) 198 170 (because of tiling I think)
    - GC_Kircher  has unique variant ids
    - GC_Mendelian_variants Group has non unique variant ids: (all, unique) 161 50 (probably also tiling)
    - C_positive_heart_CAD  has unique variant ids
    - GC_Atrial_fib  has unique variant ids
    - GC_Mohlke Group has non unique variant ids: (all, unique) 17 16
    - GC_Liang  has unique variant ids

#### Using control match table to get them into barcode MPRAlm
- first attempt: no output in final table (filtering too harsh?)
- inverstigating step by step: 
  - counts_per_bc_sorted (basically join) => worked
  - found out that controls variant map does not include a region column but script assumed one => fixed 
- we have 181 alternative sequences in GC_Selvarajan sequences => 
- merge mpralm_dna and mpralm_rna tables for test variants and controls:
  - variant table: /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/bc_preparation/controls/GC_Selvarajan_standard_NGN2_filtered_counts_sequences.tsv.gz 
```bash
first=1
for f in /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/standard_NGN2_filtered_counts_sequences.tsv.gz /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/bc_preparation/controls/GC_Selvarajan_standard_NGN2_filtered_counts_sequences.tsv.gz
do
    if [ "$first" ]
    then
        zcat "$f"
        first=
    else
        zcat "$f" | tail -n +2
    fi
done > /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/bc_preparation/test_controls_concat/GC_Selvarajan_standard_NGN2_filtered_counts_sequences.tsv
```
- concatinated variant table has 42305 unique alternative sequences
- add chrom, pos, ref alt to the table for controls and tested oligos
  - bed file for tested sequences `(/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/final_design/input/regions_5K.bed)` (not for controls!) (after replacing ~ with , in the names matching was possible)
  - where is the bed for the controls?
  - do same for elements and find control bed file 
  - ``found`` bed files for 
    - GC_Atrial_fib
    - GC_Liang
    - GC_Selvarajan
    - GC_Mohlke
    - GC_Kircher       
    - GC_Cort_Chengyu
    - GC_GABA_Chengyu
    - GC_Glut_Chengyu
    - GC_Hon
    - GC_Mendelian_variants
    - C_positive_heart_CAD
    - GC_Vista
    - GC_DNase_positive
    - GC_DNase_negative_brain
    - GC_DNase_negative_blood      
    - MK
- bed files ``missing`` for:
  - C_negative_heart_MK
  - C_negative_neuron_MK
  - C_negative_neuron_NP
  - C_positive_heart_MK
  - C_positive_neuron_CD
  - C_positive_neuron_MK
  - C_positive_neuron_NP
  - C_positive_heart_AB
  - C_SLEA
  - GC_DNase_positive_shuffeled
  - GC_DNase_negative_brain_shuffeled
  - GC_DNase_negative_blood_shuffeled
- All variant related barcode sequences: generated from (after filtering: 16264 barcodes) (path: `/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/standard_with_controls_NGN2_filtered_counts_sequences.tsv.gz`) 
  - filtered out because single sequences 5; => 194 different variant ids (different variants (considering ref and alt as different: 387 sequences))
    - combine: /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/standard_with_controls_NGN2_filtered_counts_sequences.tsv.gz and /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/standard_NGN2_filtered_counts_sequences.tsv.gz
```bash
first=1
for f in /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/standard_NGN2_filtered_counts_sequences.tsv.gz /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/standard_with_controls_NGN2_filtered_counts_sequences.tsv.gz
do
    if [ "$first" ]
    then
        zcat "$f"
        first=
    else
        zcat "$f" | tail -n +2
    fi
done > /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/standard_with_variant_controls_NGN2_filtered_counts_sequences.tsv
gzip /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/standard_with_variant_controls_NGN2_filtered_counts_sequences.tsv
```
  - Summary: `standard_NGN2_filtered_counts_sequences.tsv.gz`:
    - different variant ids: 42130 (ref and alt different: 84242) (41284 variants end in mpralm)
    - 3531738 different barcodes
  - filtered_counts_sequences based on controls: now combine with the tested sequence result and go to mpralm (preprocessing)
  - use python script: `/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/04_MPRAlm/bc_MPRAlm/scripts/preprocess.py`
    - generate input for mpralm from filtered counts (barcodes per variant sequences)
    - input number: 42322: filtering: at least 2 barcodes for variant and ref and alt (=> 41474 still valid)
    - preprocessing: 41474 variants with >= 2 barcodes in alt and ref
    - after bc_mpralm: with all variant controls (filtering for 2 barcodes per ref and alt): 41474
    - group distribution with variant controls:
      - cardiac_neuro_cava_random    41284
      - GC_Selvarajan                  112
      - C_positive_heart_CAD            43
      - GC_Atrial_fib                   17
      - GC_Mohlke                       11
      - GC_Liang                         7

- using variant region map with variant controls (all sets)
  - filtered_counts_sequences (barcode to ref and alt association and counts for each replicate): 3548001
  - unique varaint_id + allele: 84628
  - after MPRAsnakeflow and variants with ref and alt: 42611
  - filtering for 10 barcodes for ref and alt: 34933
  - for a table of these controls see share point progress from Kilian (https://charitede-my.sharepoint.com/:p:/r/personal/kilian_salomon_bih-charite_de/_layouts/15/Doc.aspx?sourcedoc=%7BC8A28A99-FF4A-4A83-AABD-8BB0C7606192%7D&file=20240301_0601progress.pptx&action=edit&mobileredirect=true)

#### Regions of controls (solvable with blat)
- current questions: missing controls have chrom, start end in header
- Beispiel header von `MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/controls/negative_250_heart_MK.fa`
```
>tile_6903_chr11_9614045_9614314_reference__0.958461950470297
>tile_35543_chr5_173245685_173245954_reference__0.556701926768965
>tile_36416_chr6_35490266_35490535_reference__0.494964379737065
>tile_16405_chr16_58119639_58119908_reference__0.479510907361922
>tile_27236_chr22_37986181_37986450_reference__0.476821379685758
```
- Beispiel header von `MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/controls/negative_250_neuron_MK.fa`
```
>tile_14444_chr15_67066278_67066547_reference__1.1385203581298
>tile_9599_chr12_65742272_65742541_reference__1.12348297546129
>tile_9256_chr12_26114182_26114451_reference__1.06883601881989
>tile_44145_chr8_80578146_80578415_reference__1.04321190480237
>tile_1999_chr1_88462556_88462825_G_C_19__0.955831785228529
```
- Beispiel von `MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/controls/positive_100_neuron_CD.fa`
```
>p1_rs9975055_T_G_alt_50_G::chr21:44929969-44930239-mean_ratio2.18
>n1_rs2279982_G_A_alt_50::chr2:164841902-164842172-mean_ratio2.16
>p1_rs34761481_G_A_alt_50_G::chr7:129161739-129162009-mean_ratio2.14
>p1_rs7115714_G_A_alt_50_A::chr11:120424017-120424287-mean_ratio2.13
```
- Beispiel von `MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/controls/positive_100_neuron_NP.fa`
```
>GW18_PFC_ABC_chr15_89400286_89400556_0.830617698776558
>NGN2_iPSC_ABC_chr4_112626228_112626498_1.11594727946851
>Midfetal_Cortex_Trevino_chr20_63366571_63366841_1.30229686520284
>NGN2_iPSC_ABC_chr2_2759136_2759406_1.42242758369854
```
- Use blat to get control regions: (first approach was not so good, because not filtered for reference sequences: I will show here only improved process)
  - all control sequences (meant to be variant controls (problems in variant map therefore not entirely sure)) with REF in name from design: `/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/05_variant_region_list/results/control_sequences/control_variant_reference_sequence.fasta` # 312
  - all sequences are unique and all sequences have a length of 300 with 15bp adapter
  - needed to remove the adapter and filter for the ":REF_" pattern
  - use blat to find perfect matches in hg38
    - `blat /data/cephfs-2/unmirrored/groups/ag-kircher/MPRA/kilian_projects/hg38/ncbi_dataset/data/GCF_000001405.26/GCF_000001405.26_GRCh38_genomic.fna /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/05_variant_region_list/results/control_sequences/control_variant_reference_sequence.fasta blat_variant_control_sequences_matched_file.psl`
    - generate table: `tail -n +6 blat_variant_control_sequences_matched_file.psl > blat_variant_control_sequence_result_without_header.tsv`
    - header list: `['match', 'mismatch', 'rep_match', 'Ns', 'Q_gap_count', 'Q_gap_bases', 'T_gap_count', 'T_gap_bases', 'strand', 'Q_name', 'Q_size', 'Q_start', 'Q_end', 'T_name', 'T_size', 'T_start', 'T_end', 'block_count', 'blockSizes', 'qStarts', 'tStarts']`
  - Next step: 08.04: generate tsv and only get the 270 matchings (check if only 1 per sequence and check the strand information)
  - 'GC_Kircher' make still problems: header too long for blat and in `/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/05_variant_region_list/resources/controls_variant_region_map.tsv.gz` different variant ids than in bed and vcf files (I think I can merge them directly to vcf)
  - merged with vcf IDs and resulting table has 203 unique headers but only 100 unique IDs which is odd
    - I cannot take these as IDs
#### Using blat to get genomic coordinates of the sequences for the metadata file (used: GCF_000001405.26_GRCh38_genomic.fna)
- code: `/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/05_variant_region_list/notebooks/blat_for_design.ipynb`
- made matching correct for elements and variants (found that bed file is not what I assumed it is)
- Found region file is not helpful: `MPRA/IGVF_Y1_design/design/final_design/input/regions_5K.bed` (not 270 bp but from 150 - 370 long regions)
  - create Sam file from matches of sequences without adapter yourself (use blat locally) `blat /data/cephfs-2/unmirrored/groups/ag-kircher/MPRA/kilian_projects/hg38/ncbi_dataset/data/GCF_000001405.26/GCF_000001405.26_GRCh38_genomic.fna /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/05_variant_region_list/tested_ref_and_elements.fasta /data/cephfs-2/unmirrored/groups/ag-kircher/MPRA/kilian_projects/hg38/ncbi_dataset/data/GCF_000001405.26/blat_matched_file.psl` => (generate table: `tail -n +6 blat_matched_file.psl > blat_result_without_header.tsv`)
  - header list: `['match', 'mismatch', 'rep_match', 'Ns', 'Q_gap_count', 'Q_gap_bases', 'T_gap_count', 'T_gap_bases', 'strand', 'Q_name', 'Q_size', 'Q_start', 'Q_end', 'T_name', 'T_size', 'T_start', 'T_end', 'block_count', 'blockSizes', 'qStarts', 'tStarts']`
  - start positions have small differences (e.g. 39 bases)
    - 39    1955
    - 38    1508
    - 35    1380
    - 40    1364
    - 37    1231
  - sanity check:
    - are all sequences of blat inbetween the regions of region_5k.bed (old file):
      - No: 11538 sequences are not inbetween because they are not as large as the 270bps (their length is <= 271)
    - I checked in UCSC: 
      - `cardiac_neuro_cava_random:PRDM16|ENSG00000142611.17|EH38E2779779_fwd_tile1-1` one example region where the ucsc does not show any information about a cCRE from encode: 
        - blat regions: NC_000001.11 (chr1) 3210194 3210464
        - region_5k regions: chr1 3210229 3210430
        - [UCSC with blat regions](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A3210194%2D3210464&hgsid=2082689898_VAczbUehrymUlkaWiiYiCorrQZsL)
- 22.04: still problems with the headers and afraid of having matching problems if done too fast (deadline for bed to Thorben)
  - create fasta of sequences for "interesting" sequences (are in NP, MK or CD neuro positive / negative control group)
    - `/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/05_variant_region_list/results/control_sequences/blat_interesting_neuro_controls.fa`
    - `blat /data/cephfs-2/unmirrored/groups/ag-kircher/MPRA/kilian_projects/hg38/ncbi_dataset/data/GCF_000001405.26/GCF_000001405.26_GRCh38_genomic.fna /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/05_variant_region_list/results/control_sequences/blat_interesting_neuro_controls.fa /data/cephfs-2/unmirrored/groups/ag-kircher/MPRA/kilian_projects/hg38/ncbi_dataset/data/GCF_000001405.26/blat_interesting_neuro_controls_matched.psl`
    - `tail -n +6 /data/cephfs-2/unmirrored/groups/ag-kircher/MPRA/kilian_projects/hg38/ncbi_dataset/data/GCF_000001405.26/blat_interesting_neuro_controls_matched.psl > /data/cephfs-2/unmirrored/groups/ag-kircher/MPRA/kilian_projects/hg38/ncbi_dataset/data/GCF_000001405.26/blat_interesting_neuro_controls_matched.tsv`
  - investigate interesting control headers and create parsing strategy for each group (for neuro use _chr as split)
    - `MK`: (coordinates not checked yet)
      - `>MK:SNRNP70|chr19-49085144+49085413|reference`
    - According to blat: coordinates are the same for MK but for NP: header_start-1 = blat_start
    - `C_positive_neuron_NP`: _chr_start_end (start, end]
      - `>C_positive_neuron_NP:GW18_PFC_ABC_Midfetal_Cortex_Trevino_Midfetal_Cortex_Ziffra_ABC_chr10_79221115_79221385_5.21059445249138`
    - `C_positive_neuron_MK`: chr_ start _ end (both within the range [start,end])
      - `>C_positive_neuron_MK:tile_35742_chr6_3247831_3248100_reference_0.892141141777512`
    - `C_negative_neuron_MK`: chr_ start _ end (both within the range [start,end])
    - `C_positive_neuron_CD`: ::chr: start - end (start, end]
      - `>C_positive_neuron_CD:p1_rs55985730_T_G_alt_50_T::chr7:128776855-128777125-mean_ratio2.34`
    - excluded cases: 
      - `C_positive_neuron_CD`:
        - `>C_positive_neuron_CD:c1_NA_NA_NA_NA::72hr_top_94-mean_ratio2.31`
#### Found gene names per gene set: 
- Github repo [MPRA_design](https://github.com/kircherlab/MPRA_design/tree/main) resources: `MPRA_design/resources/gene_lists`

- Todo: metadata table
  - investigated variant information (vcf like file) and identifier from there are without "," or "~" so remember to split the ids by "~" later in order to match the variant information
  - check control region (length check first)
    - Not all controls have region files but for these who have: region files are 270 bp long
  - verify if regions are fine (after mapping to the reference)
  - add regions to the alternative sequences with alt_ref_map_table
    - go over variant_region_map and make dict: ref: [corresponding alt]
  - add SPDI (for variants)
  - add for each reference the SPDI of each variant it is responsible for
- Todo: variant analysis
  - Check numbers of controls
  - Add all controls to bcMPRAlm
    - Generated alt_ref_maps for analysis: `80K-Analysis/05_variant_region_list/results/variant_control_map`
    - combine these with variant region map
  - make workflow for bcMPRAlm data preparation
  - make elementwise bcMPRAlm (notes on how to do it in perform bc R script)
- check if control bed files are what they look like
  - GC_Cort_Chengyu has only 270bp regions
  - GC_GABA_Chengyu has only 270bp regions
  - GC_Glut_Chengyu has only 270bp regions
  - GC_Hon has only 270bp regions
  - GC_Vista has only 270bp regions
  - GC_DNase_positive has only 270bp regions
  - GC_DNase_negative_brain has only 270bp regions
  - GC_DNase_negative_blood has only 270bp regions

#### Significant results 
- neuro: 
  - high positive result: e.g. 'KCNT2|ENSG00000162687.19|EH38E2855757|1-196494459-G-A'
    - in R could only find one high logfc alt, Na at ref (2 vs NA) and a bit higher alt vs ref: 1.7 vs 1.41

#### 09.04.2024: Meeting with Mohan about matching and bc_MPRAlm
- control: ref alt id match: check if ref and alt id are not there than they did not make it in the design
- check the header modifications
- enformer results:
  - he will prepare enformer predictions from fasta like tsv files
  - Process for variant effects (delta)
    - default by enformer with variants: takes the whole context (>190KB) from vcf
    - NOT 1kb and not only sequence
  - enformer is doing a alt - ref => positive if alt is higher than ref and negative if ref is higher than alt

### 11.04.2024: investigate bc_MPRAlm results
- all bc_MPRAlm sequences have 10 barcodes: no, 4300 don't have
  - why? looking into mprasnakeflow: rule assigned_counts_make_master_tables (`/workflow/rules/assigned_counts.smk`) uses the threshold and produces
    - Input for bc_MPRAlm_preparation: count files from mprasnakeflow (unfiltered!!) (after producing these the workflow filters for the specified thresholds in the config file)
    - `bc_merged="results/experiments/{project}/assigned_counts/{assignment}/{config}/{condition}_allreps_merged_barcode_assigned_counts.tsv.gz",` is the filtered and merged file
    - current filter function was only checking for at least 2 barcodes (42322 is number of sequences with at least 2 barcodes (with part of control))
    - if filtered for 10 barcodes only 34683 variants can be found 
    - if filtered for 5 barcodes only 38858 variants can be found


### 18.04.2024: 

#### download re-sequencing
- infor from GNU ftp wget manual (https://ftp.gnu.org/old-gnu/Manuals/wget-1.8.1/html_mono/wget.html)
- ftp://user:password@host:port it needs pretty long to connect (what means long?) over 5 minutes
- checksums: (found all manually)
- md5sum Ngn2-DNA-1_S1_R1_001.fastq.gz
  - cb19ab721a996189e721641641c7f64a  Ngn2-DNA-1_S1_R1_001.fastq.gz
  - 1745e7c6ae3a9e2a6b52668fc2b8f2d4  Ngn2-DNA-1_S1_R2_001.fastq.gz


#### investigate why only 41 Mendelian variant controls in the final design
- current analysis: check how many will pass the quality standard of 75% of the barcodes need to be matched to one sequence
- check from assignment file how many sequences would have at least 10 barcodes => 197
- Why do we only find 41 variants: 
```bash
awk -F'\t' '{
    split($NF, arr, "/");
    ratio = arr[2] / arr[1];
    if (ratio > 0.75) print $0
}' <( zcat /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/experiments/standard_bwa/assignment/assignmentFixDuplicates.tsv.gz | grep "Mendelian" ) | cut -f2 | sort | awk '{ count[$0]++ } END { for (entry in count) if (count[entry] >= 10) print entry }' | less -S
```
- check counts + assignment:
```bash
zcat /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/bc_preparation/standard_with_all_variant_controls_NGN2_counts_per_bc_sorted.tsv.gz | grep "Mendelian" | cut -f2 | sort | uniq | wc -l
# 204 / 209 design
zcat /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/bc_preparation/standard_with_all_variant_controls_NGN2_counts_per_bc_sorted.tsv.gz | grep "Mendelian" | cut -f2 | sort | uniq | grep "ALT" | wc -l
# 157 / 161 design
zcat /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/bc_preparation/standard_with_all_variant_controls_NGN2_counts_per_bc_sorted.tsv.gz | grep "Mendelian" | cut -f2 | sort | uniq | grep "REF" | wc -l
# 47 / 49 design
```
- counts_sequences 
```bash
zcat /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/bc_preparation/standard_with_all_variant_controls_NGN2_counts_sequences.tsv.gz | grep "Mendelian" | cut -f 1 | sort | uniq | wc -l
# 49
```
- new variant controls
```bash
zcat /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/projects/bc_tradeoff/results/bc_preparation/standard_with_all_variant_controls_mendelian_NGN2_counts_sequences.tsv.gz | grep "Mendelian" | cut -f 1 | sort |uniq | wc -l
# 160
```

### 23.04.2024
- found output workflow in `final_design/results/final_design/region.bed`
- label
  - cardiac_neuro_cava_random    27556
  - GC_Vista                       256
  - GC_Cort_Chengyu                185
  - GC_Selvarajan                  166
  - GC_GABA_Chengyu                 85
  - GC_Glut_Chengyu                 83
  - GC_Atrial_fib                   22
  - GC_Mohlke                       18
  - GC_Liang                         8
  - GC_Hon                           6
  - GC_Kircher                       5

- found non snv in final region.bed 
  - 'C_positive_heart_CAD:rs34091558'
    - 'rs34091558'
  - 'GC_Mendelian_variants:chr1:209816133C>CA|IRF6'
    - 'chr1:209816133C>CA|IRF6'
  - 'GC_Mohlke:NC000001_11_230159329_CTTAAAGTGTTCAGCACTCCCCT_CT'
    - 'NC000001.11|230158967|C|A|MohlkeHepControls,NC000001.11|230159168|C|T|MohlkeHepControls,NC000001.11|230159329|CTTAAAGTGTTCAGCACTCCCCT|CT|MohlkeHepControls_fwd_tile3-3'
  - 'GC_Mendelian_variants:chr7:156791274T>TTAAGGAAGTGATT|SHH'
    - 'chr7:156791255G>C|SHH', 'chr7:156791257G>A|SHH', 'chr7:156791274T>TTAAGGAAGTGATT|SHH'
  - 'GC_Mendelian_variants:chr8:11703890AG>A|GATA4'
    - 'chr8:11703860G>T|GATA4', 'chr8:11703890AG>A|GATA4'

- found that regions in vcf file are not the same as regions id in bed file
- changing bed file (adding new column without group information ( split ":" ))
- multiple regions can be found (641 in cardiac_neuro_cava_random)
- both regions do not need to be the same
- I tried to map one variant position to the region of interest but over 700 variants have multiple regions and some of the listed region ids from the vcf file do not exist in the bed file (=> this in not unique)
  - Example: `cardiac_neuro_cava_random:ARV1|ENSG00000173409.14|EH38E2873164|1-231021494-T-C`
- Procedure I will only take the first region and hope the best
- Found 136 regions in the vcf file which are not in the region.bed they remain because of a bug in the MPRAOligoDesign
- region cannot be found for some sequences: 
  - For all sequences: (example: GC_Mendelian_variants:chr1:21564170G>A|ALPL)
    - Empty region_list:  0
    - Region not findable:  65
    - Skipped, because multiple regions 789 
  - Only for tested:
    - Region not findable:  0
    - Skipped, because multiple regions 641)

- Found 94 ref and alt sequences from the old designed headers / I havent checked all but they apear in the new design file as well but neither in the variant map nor in the variant.vcf are not in the variant region map of the new design
- Check in the MPRAoligo design code: check previous unfiltered results and where these sequences get lost. (Issue in the code?)
  - variants which are filtered out remain in the design.fasta (e.g. ALT_LMNA|ENSG00000160789.24|EH38E2840374_fwd_tile1-1_LMNA|ENSG00000160789.24|EH38E2840374|1-156033864-A-T)
- Filtered best 100, low 100 and 200 random for >50 barcodes
  - all these we have regions for
- Number of tested sequences we have regions for: 73846
- Overall number of sequences in the deduplicated design: 73940
- Answer Chengyu
Hi Chengyu, Hi Nadav,
Thank you for the additional sequencing data, after some issues with our firewall, I was now able to download it. I will keep you posted with the results after running MPRAsnakeflow on the combined data. 
Btw, apparently, I have made some incorrect assumptions about your control set (aka C_positive_neuron_CD). I had thought that you had given us 100 positive and negative sequences during my presentation at the call two weeks ago. That was wrong, I confused it with the SLEA control group where we have 100 positive and 100 negative sequences. Turns out that these were supposed to be only positive controls. 
I have looked at the correlation between the values in the headers and our mean ratios of the 3 replicates. Unfortunately, the correlation between these two experiments does not look so good. However, if I increase the number of required barcodes per insert (e.g. to 70), then the correlation increases (see plots attached for a min. of 15, 30 and a min. of 70 barcodes per insert). However, I am wondering whether I am using the correct values for the previous read-out of these elements. I assumed that the numerical value in the header is the value from your experiment. For example:
>C_positive_neuron_CD:p1_rs6813360_A_C_ref_50_A::chr4:155359187-155359457-mean_ratio2.42
Is this correct?
Best,
Kilian
>C_positive_neuron_CD:p1_rs6813360_A_C_ref_50_A::chr4:155359187-155359457-mean_ratio2.42
>C_positive_neuron_CD:p1_rs55985730_T_G_alt_50_T::chr7:128776855-128777125-mean_ratio2.34
- investigating neuro insteresting
  - C_negative_neuron_MK:rdhs_441080_chr3_13042165_13042434_reference__0.47708598278664
    - negative strand hit coordinates not (start, end]

### Resequencing data seems like only 15bp reads are there

- Ngn2-DNA-1_S1_R1_001.fastq.gz
- @NB501960:958:H3FCMBGXW:1:11101:19983:1040 1:N:0:TTGAACTNAG
- TGGGTCCATGCCTGA
- +
- AAAAAEEEEEEEEEE

- Ngn2-DNA-1_S1_R2_001.fastq.gz
- @NB501960:958:H3FCMBGXW:1:11101:19983:1040 2:N:0:TTGAACTNAG
- CCAGGGCACNGGCTAC
- +
- AAAAAEEEE#EEEEEE

- Ngn2-DNA-1_S1_R3_001.fastq.gz
- @NB501960:958:H3FCMBGXW:1:11101:19983:1040 3:N:0:TTGAACTNAG
- TCAGGCNTGNACCCN
- +
- AAAAAE#EE#EEEE#


- Ngn2-DNA-2_S2_R1_001.fastq.gz
- @NB501960:958:H3FCMBGXW:1:11101:18158:1040 1:N:0:TATAATTNAA
- TTTGACAAACTTGCG
- +
- AAAAAEEEEEEEEEE

- Ngn2-DNA-3_S3_R1_001.fastq.gz
- @NB501960:958:H3FCMBGXW:1:11101:25544:1042 1:N:0:GAAGCTGNAG
- CGCCCTGGTTATAAT
- +
- AAAAAEAEEEEAEEA

- old 80K data: 146 basepairs long
- FW: /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/data/association/IGVF_neuro_S1_R1_001.fastq.gz
- @NB501960:812:HH53WAFX5:1:11101:26148:1108 1:N:0:NCGTAGACCA
GGCTTCTGATAAGCCGCCAATTCATAGTGTGGGTTTGGAGAGCTGGAGACGGGGTGAGAAAGCTGAGGCCTTTGCAAAGTCTATTTACACAGTGGCCAGAGTCCCTTTCACCCTCTCCGGCAAACAGGCCCTGGAGCACAGGCATT
+
- BC: /data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/data/association/IGVF_neuro_S1_R2_001.fastq.gz
- check length of reads with unix command 
```bash
file_list="Ngn2-DNA-1_S1_R1_001.fastq.gz Ngn2-DNA-1_S1_R2_001.fastq.gz Ngn2-DNA-1_S1_R3_001.fastq.gz Ngn2-DNA-2_S2_R1_001.fastq.gz Ngn2-DNA-2_S2_R2_001.fastq.gz Ngn2-DNA-2_S2_R3_001.fastq.gz Ngn2-DNA-3_S3_R1_001.fastq.gz Ngn2-DNA-3_S3_R2_001.fastq.gz Ngn2-DNA-3_S3_R3_001.fastq.gz"
for file in $file_list; do
  echo "$file"
  awk 'NR%4 == 2 {if (length($0) > 15) print length($0)}' <(zcat $file)
done
```
- Ngn2-DNA-1_S1_R1_001.fastq.gz: 15 length
- Ngn2-DNA-1_S1_R2_001.fastq.gz: 16 length
- non of the resequences DNA files contains a read above 100
- check RNA:
```bash
file_rna_list="Ngn2-RNA-1_S4_R1_001.fastq.gz Ngn2-RNA-1_S4_R2_001.fastq.gz Ngn2-RNA-1_S4_R3_001.fastq.gz Ngn2-RNA-2_S5_R1_001.fastq.gz Ngn2-RNA-2_S5_R2_001.fastq.gz Ngn2-RNA-2_S5_R3_001.fastq.gz Ngn2-RNA-3_S6_R1_001.fastq.gz Ngn2-RNA-3_S6_R2_001.fastq.gz Ngn2-RNA-3_S6_R3_001.fastq.gz"
for file in $file_rna_list; do
  echo "$file"
  awk 'NR%4 == 2 {if (length($0) > 100) print length($0)}' <(zcat $file)
done
```

### Modify current design file
```bash
#!/bin/bash

# Set the input and output file names
input_file="/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/unifying_headers_kilian_042024/design_no_duplicates_old_naming.fa"
output_file="/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/unifying_headers_kilian_042024/GC_Kircher_shortened_headers_match_table.tsv"

# grep all headers and print them as tsv with the fist part as label
grep "^>" $input_file | sed 's/^>//g' | awk -F ":" '{print $1":"$2 "\t" $1}' > $output_file
```
```bash
snakemake --use-conda  --configfile standard_config.yaml --snakefile /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/workflow/Snakefile --conda-prefix /data/gpfs-1/users/kisa11_c/work/coding/MPRA/MPRAsnakeflow_projects/conda --keep-going --cluster-config /data/gpfs-1/users/kisa11_c/work/coding/MPRAsnakeflow/config/sbatch.yml --cluster-status status.py --cluster "sbatch --parsable --nodes=1 --ntasks-per-node={cluster.threads} --mem {cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output} -e {cluster.error}"  --jobs 40 --cluster-cancel scancel -n --quiet
```

### Sanity checking MPRAOligoDesign
- location of the data to be checked: `MPRA/IGVF_Y1_design/design/final_design/results/oligo_design/cardiac_neuro_cava_random`
- Both examples are in the design file
- e.g. cardiac_neuro_cava_random:ALT_MTHFR|ENSG00000177000.13|EH38E1319148_rev_tile1-1_MTHFR|ENSG00000177000.13|EH38E1319148|1-11841412-C-T (is in design and could be matched) => is in the design and the variant region map and the variant file
- e.g. cardiac_neuro_cava_random:ALT_DOCK7|ENSG00000116641.18|EH38E1354059_rev_tile1-1_DOCK7|ENSG00000116641.18|EH38E1354059|1-62589630-A-C (is in not fitting) => is in the design file but not in the variant region map or variant file
  - is still in `design_variants.variant_region_map.tsv.gz` but not in filtered `design_variants_filtered.variant_region_map.tsv.gz` 


### Investigate the ultraconserved variants:
- found 69 variants in the 80K design which are within the ultraconserved regions
- sanity check: found all of these in the variant region map
- look for significant variants: 1 is significant but no paper over the variant found

### 