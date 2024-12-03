# create design_no_duplicates_sequence_and_header.fa

we have duplicate headers:

```bash
cat /fast/groups/ag_kircher/MPRA/IGVF_Y1_design/resources/association_data/design_no_duplicates.fa | \
grep ">" | sort | uniq -c | sort -nr | head -n 10
```

And we have spaces which do not belong to the header:

```bash
cat /fast/groups/ag_kircher/MPRA/IGVF_Y1_design/resources/association_data/design_no_duplicates.fa | \
grep ">" | grep " " | head -n 10
```

```bash
cat /fast/groups/ag_kircher/MPRA/IGVF_Y1_design/resources/association_data/design_no_duplicates.fa | \
sed 's/ /:/g' | \
sed '0,/>GC_Mohlke:REF_NC000001.11|230158967|C|A|MohlkeHepControls~NC000001.11|230159168|C|T|MohlkeHepControls~NC000001.11|230159329|CTTAAAGTGTTCAGCACTCCCCT|CT|MohlkeHepControls_fwd_tile3-3/s//>GC_Mohlke:REF_NC000001.11|230158967|C|A|MohlkeHepControls~NC000001.11|230159168|C|T|MohlkeHepControls~NC000001.11|230159329|CTTAAAGTGTTCAGCACTCCCCT|CT|MohlkeHepControls_fwd_tile3-3_headerDuplicate1_2/' | \
sed '0,/>GC_Mohlke:REF_NC000001.11|230158967|C|A|MohlkeHepControls~NC000001.11|230159168|C|T|MohlkeHepControls~NC000001.11|230159329|CTTAAAGTGTTCAGCACTCCCCT|CT|MohlkeHepControls_fwd_tile3-3/s//>GC_Mohlke:REF_NC000001.11|230158967|C|A|MohlkeHepControls~NC000001.11|230159168|C|T|MohlkeHepControls~NC000001.11|230159329|CTTAAAGTGTTCAGCACTCCCCT|CT|MohlkeHepControls_fwd_tile3-3_headerDuplicate2_2/' | \
sed '0,/>GC_Mendelian_variants:REF_chr8:11703890AG\*A|GATA4/s//>GC_Mendelian_variants:REF_chr8:11703890AG\*A|GATA4_headerDuplicate1_2/' | \
sed '0,/>GC_Mendelian_variants:REF_chr8:11703890AG\*A|GATA4/s//>GC_Mendelian_variants:REF_chr8:11703890AG\*A|GATA4_headerDuplicate2_2/' | \
sed '0,/>GC_Mendelian_variants:REF_chr8:11703860G\*T|GATA4/s//>GC_Mendelian_variants:REF_chr8:11703860G\*T|GATA4_headerDuplicate1_2/' | \
sed '0,/>GC_Mendelian_variants:REF_chr8:11703860G\*T|GATA4/s//>GC_Mendelian_variants:REF_chr8:11703860G\*T|GATA4_headerDuplicate2_2/' \
> /fast/groups/ag_kircher/MPRA/IGVF_Y1_design/resources/association_data/design_no_duplicates_sequence_and_header.fa
```

### create design_no_duplicates_sequence_and_header.fa
- still duplicates in sequences
- downstream problem: in originial design file headers had "," in them. Mohan said they were replaced with "~" because they made problems with a MPRAsnakeflow alignment tool (see gitter from 4.4.2024) but they were not changed for the enformer analysis

### create removed_brackets_design_no_duplicates_sequence_and_header.fa
```bash
 sed 's/\[.\]/_./g' /data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/resources/association_data/design_no_duplicates_sequence_and_header.fa > /data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/resources/association_data/removed_brackets_design_no_duplicates_sequence_and_header.fa
```
- make label file fitting removed brackets design file:
```bash
( echo -e "name\tlabel"; grep "^>"
removed_brackets_design_no_duplicates_sequence_and_header.fa | sed 's/^>//g'
| awk -F ":" '{print $0 "\t" $1}'; )
> removed_brackets_design_no_duplicates_sequence_and_header_label.fa
```
- removed adapters because of forward collisions:
  `design_no_duplicates_sequence_and_header_no_adapter_no_brackets.fa`
- added random adapters (performs better in bwa?)
  - `scp /home/kisa/coding/80K_MPRA/design_data/design_no_duplicates_sequence_and_header_random_adapter_no_brackets.fa kisa11_c@172.16.35.182:/data/cephfs-2/unmirrored/groups/kircher/users/kisa11_c/projects/MPRA/IGVF_Y1_design/resources/association_data/`

### handle collisions in designed sequences
- collisions from check_design_file rule of MPRAsnakeflow (29.07.2024) were
  used and following headers were handled seperately:
```python
removable_header_bc_of_collisions = ['MK:tile_19098|chr18-25465452+25465722|G-A-1', 'MK:tile_19098|chr18-25465452+25465722|A-G-268', 'MK:tile_37449|chr6-98703408+98703677|reference', 'GC_Mendelian_variants:ALT_chr7:156791255G*C|SHH_chr7:156791274T*TTAAGGAAGTGATT|SHH', 'GC_Mendelian_variants:ALT_chr7:156791581A*G|SHH_chr7:156791579C*T|SHH',
 'GC_Mendelian_variants:ALT_chr7:156791579C*T|SHH_chr7:156791542A*C|SHH', 'MK:tile_47638|chr16-52435608+52435878|A-T-1', 'GC_Mendelian_variants:ALT_chr7:156791579C*T|SHH_chr7:156791581A*G|SHH', 'MK:tile_47629|chr14-103542886+103543156|C-G-33', 'MK:tile_47742|chr5-87944920+87945190|G-C-165', 'GC_Mendelian_variants:REF_chr7:156791472C*G|SHH',
 'GC_Mendelian_variants:REF_chr7:156791579C*T|SHH', 'MK:tile_47629|chr14-103542886+103543156|reference', 'GC_Vista:fb;fm_mm1912_vistaElementControl|chr7:42140681-42140930', 'MK:tile_47627|chr14-103542806+103543076|T-A-248', 'GC_Mendelian_variants:REF_chr10:23219434A*G|PTF1A', 'MK:tile_47626|chr14-103542805+103543075|reference',
 'MK:tile_37403|chr6-98416848+98417117|reference', 'GC_Mendelian_variants:ALT_chr7:156791472C*G|SHH_chr7:156791459T*C|SHH', 'GC_Mendelian_variants:ALT_chr7:156791255G*C|SHH_chr7:156791255G*C|SHH', 'MK:tile_47744|chr5-87945000+87945270|T-A-62', 'MK:tile_47743|chr5-87944921+87945191|T-A-141', 'MK:tile_47626|chr14-103542805+103543075|C-G-114',
 'MK:tile_47745|chr5-87945001+87945271|reference', 'GC_Selvarajan:ALT_rs499966|STARR-seq-HepG2_fwd_tile1-1_rs499966', 'GC_Mendelian_variants:ALT_chr7:156791474G*A|SHH_chr7:156791571T*A|SHH', 'GC_Mendelian_variants:ALT_chr10:23219436A*G|PTF1A_chr10:23219434A*G|PTF1A', 'C_positive_neuron_MK:tile_34824_chr5_141041068_141041337_A_T_1_0.562204129450572', 'C_positive_neuron_MK:tile_34824_chr5_141041068_141041337_G_C_269_0.550310056665328', 'MK:tile_47785|chr6-164344756+164345026|G-C-269', 'MK:tile_47725|chr4-125519552+125519822|T-A-269', 'GC_Mendelian_variants:ALT_chr10:23219434A*G|PTF1A_chr10:23219517A*C|PTF1A',
 'GC_Mendelian_variants:ALT_chr7:156791579C*T|SHH_chr7:156791547A*G|SHH', 'MK:tile_47629|chr14-103542886+103543156|T-A-168', 'GC_Mendelian_variants:ALT_chr7:156791472C*G|SHH_chr7:156791472C*G|SHH', 'GC_Mendelian_variants:ALT_chr10:23219434A*G|PTF1A_chr10:23219376A*C|PTF1A', 'MK:tile_47562|chr1-52663056+52663326|C-G-1', 'MK:tile_47791|chr6-170438656+170438926|G-C-1',
 'GC_Mendelian_variants:ALT_chr7:156791472C*G|SHH_chr7:156791480G*A|SHH', 'GC_Mendelian_variants:ALT_chr7:156791579C*T|SHH_chr7:156791472C*T|SHH', 'GC_Mendelian_variants:ALT_chr7:156791474G*A|SHH_chr7:156791542A*C|SHH', 'GC_Vista:fb_hs262_vistaElementControl|chr5:77645014-77645263', 'GC_Mendelian_variants:ALT_chr7:156791474G*A|SHH_chr7:156791581A*G|SHH',
 'GC_Mendelian_variants:ALT_chr7:156791255G*C|SHH_chr7:156791257G*A|SHH', 'MK:tile_47744|chr5-87945000+87945270|G-C-265', 'MK:tile_47747|chr5-87945240+87945510|T-A-1', 'MK:tile_47627|chr14-103542806+103543076|C-G-129', 'GC_Mendelian_variants:ALT_chr7:156791579C*T|SHH_chr7:156791474G*A|SHH', 'MK:tile_47598|chr12-102961789+102962059|A-T-1',
 'GC_Mendelian_variants:ALT_chr7:156791472C*G|SHH_chr7:156791472C*T|SHH', 'GC_Mendelian_variants:ALT_chr7:156791472C*G|SHH_chr7:156791579C*T|SHH', 'MK:rdhs_334645|chr2-6772346+6772615|reference', 'MK:tile_47744|chr5-87945000+87945270|T-A-241', 'GC_Selvarajan:REF_rs499966|STARR-seq-HepG2_fwd_tile1-1', 'GC_Mendelian_variants:ALT_chr7:156791474G*A|SHH_chr7:156791413A*C|SHH',
 'GC_Mendelian_variants:ALT_chr7:156791581A*G|SHH_chr7:156791472C*G|SHH', 'MK:tile_44725|chr8-127646130+127646399|reference', 'MK:tile_47629|chr14-103542886+103543156|C-G-49', 'MK:newcore_405894|chrX-71182006+71182275|reference', 'GC_Mendelian_variants:ALT_chr7:156791474G*A|SHH_chr7:156791547A*G|SHH', 'C_negative_neuron_NP:Fetal_Cerebrum_Cicero_chr10_11131297_11131567_2.37359878574768',
 'MK:tile_985|chr1-33363965+33364235|C-T-268', 'MK:tile_985|chr1-33363965+33364235|T-C-1', 'GC_Mendelian_variants:ALT_chr10:23219436A*G|PTF1A_chr10:23219436A*G|PTF1A', 'GC_Mendelian_variants:ALT_chr7:156791579C*T|SHH_chr7:156791571T*A|SHH', 'GC_Mendelian_variants:ALT_chr7:156791472C*G|SHH_chr7:156791474G*A|SHH', 'GC_Vista:fb;mb__vistaElementControl|chr14:78308172-78308421',
 'GC_Mendelian_variants:REF_chr7:156791255G*C|SHH', 'MK:tile_47743|chr5-87944921+87945191|reference', 'MK:tile_10173|chr12-113932852+113933122|C-T-268', 'MK:tile_10173|chr12-113932852+113933122|G-A-1', 'GC_Mendelian_variants:ALT_chr10:23219434A*G|PTF1A_chr10:23219508A*G|PTF1A', 'MK:tile_47744|chr5-87945000+87945270|G-C-85', 'GC_Mendelian_variants:ALT_chr7:156791581A*G|SHH_chr7:156791480G*A|SHH'
]
```
- with adapter:
  - `design_no_duplicates_sequence_and_header_with_adapter_no_brackets_collisions.fa`
  - `design_no_duplicates_sequence_and_header_with_adapter_no_brackets_no_collisions`
- without adapter:
  - `/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/resources/association_data/design_no_duplicates_sequence_and_header_no_adapter_no_brackets_no_collisions.fa`
  - `/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/resources/association_data/design_no_duplicates_sequence_and_header_no_adapter_no_brackets_collisions.fa`

- without collisions and with adapter but without brakets: created locally within state_of_data.ipynb 10.09.2024 `renamed_design_no_duplicates_sequence_and_header_with_adapter_no_brackets_no_collisions.fa`

- try to remove collisions with adapters added to the fastq files of the
  assignment
- currently newest fasta and label file:
  - `/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/resources/association_data/removed_brackets_design_no_duplicates_sequence_and_header.fa`
  - `/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/resources/association_data/removed_brackets_design_no_duplicates_sequence_and_header_label.tsv`

#### Renamed three sequences from GC_Selvarajan
- in notebook (80K-Analysis/07_quality_control/notebooks/state_of_data.ipynb) => wrote again to /data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/resources/association_data/removed_brackets_design_no_duplicates_sequence_and_header.fa

## NOTE: During preparing the Metadata file we noticed that the design_no_duplicates has removed the headers without merging them
- this means these sequences cannot be assessed later
  remove the duplicated sequences by concatenating the headers with '#'
- see generate_deduplicated_design_file.py

#### Deduplicated design file:
- no_duplicated_sequences_design.fa: headers with same sequences are
  concatenated with a '#' as separator (note: if you want to match it again
with the variant map or the regions.bed split the headers within a pandas
dataframe shown in the notebook)
- found duplicated sequence: header is reference but sequence is from alternative (GC_Mohlke:REF_NC000001.11|230158967|C|A|MohlkeHepControls,NC000001.11|230159168|C|T|MohlkeHepControls,NC000001.11|230159329|CTTAAAGTGTTCAGCACTCCCCT|CT|MohlkeHepControls_fwd_tile3-3)
- use `zcat /data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/final_design/results/final_design/design.fa.gz | sed 's/ /:/g' > design_removed_spaces.fa`

