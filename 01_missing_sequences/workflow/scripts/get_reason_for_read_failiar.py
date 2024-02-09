# Needs to have a bam file with an index
import os
import pandas as pd
import pysam
import yaml

# get the input file (bam file)
config_path = '/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/config/config.yml'
config_path = '/home/kisa/coding/80K_MPRA/80K-Analysis/01_missing_sequences/config/config.yml' # local

with open(config_path) as conf:
    config = yaml.load(conf, Loader=yaml.FullLoader)
    conf.close()

# merged_bam = config['files']['merged_bam']

missing_sequence_bam = config['files']['missing_sequence_bam']

# quality measures from config (/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/standard_config.yaml):
mapping_quality_min = config['quality']['min_mapping_quality']
alignment_start_min = config['quality']['alignment_start']['min']
alignment_start_max = config['quality']['alignment_start']['max']
sequence_length_min = config['quality']['sequence_length']['min']
sequence_length_max = config['quality']['sequence_length']['max']
missing_sequences_min_quality = os.path.join(config['general']['output_dir'], "identified_missing_sequences", "missing_sequences_min_quality.bam")        
missing_sequences_alignment_start = os.path.join(config['general']['output_dir'], "identified_missing_sequences", "missing_sequences_alignment_start.bam")
missing_sequences_alignment_end = os.path.join(config['general']['output_dir'], "identified_missing_sequences", "missing_sequences_alignment_end.bam")
missing_sequences_sequence_length_min = os.path.join(config['general']['output_dir'], "identified_missing_sequences", "missing_sequences_sequence_length_min.bam")
missing_sequences_sequence_length_max = os.path.join(config['general']['output_dir'], "identified_missing_sequences", "missing_sequences_sequence_length_max.bam")

samfile = pysam.AlignmentFile(missing_sequence_bam, "rb")
min_quality_bam = pysam.AlignmentFile(missing_sequences_min_quality, "wb", template=samfile)
alignment_start_bam = pysam.AlignmentFile(missing_sequences_alignment_start, "wb", template=samfile)
alignment_end_bam = pysam.AlignmentFile(missing_sequences_alignment_end, "wb", template=samfile)
sequence_length_min_bam = pysam.AlignmentFile(missing_sequences_sequence_length_min, "wb", template=samfile)
sequence_length_max_bam = pysam.AlignmentFile(missing_sequences_sequence_length_max, "wb", template=samfile)

for read in samfile.fetch():
    # # implemented logic
    # if read < mapping_quality_min => write
    # if read ($4) < alignment_start_min => write
    # if read ($4) > alignment_start_max => write
    # print(read.reference_end)
    # if sequence length ($10) < sequence_length_min => write
    # if sequence length ($10) > sequence_length_max => write

    if read.mapping_quality < mapping_quality_min:
        prepared_read_name = f'@{read.query_name} {read.tags[-1][0]}:Z:{read.tags[-1][1]}'
        # write to min_quality_bam
        min_quality_bam.write(read)

    if read.reference_start < alignment_start_min:
        prepared_read_name = f'@{read.query_name} {read.tags[-1][0]}:Z:{read.tags[-1][1]}'
        # write to alignment_start
        alignment_start_bam.write(read)

    if read.reference_start > alignment_start_max:
        prepared_read_name = f'@{read.query_name} {read.tags[-1][0]}:Z:{read.tags[-1][1]}'
        # write to alignment_end
        alignment_end_bam.write(read)

    if len(read.query_sequence) < sequence_length_min:
        prepared_read_name = f'@{read.query_name} {read.tags[-1][0]}:Z:{read.tags[-1][1]}'
        # write to sequence_length_min
        sequence_length_min_bam.write(read)

    if len(read.query_sequence) > sequence_length_max:
        prepared_read_name = f'@{read.query_name} {read.tags[-1][0]}:Z:{read.tags[-1][1]}'
        # write to sequence_length_max
        sequence_length_max_bam.write(read)

samfile.close()
min_quality_bam.close()
alignment_start_bam.close()
alignment_end_bam.close()
sequence_length_min_bam.close()
sequence_length_max_bam.close()