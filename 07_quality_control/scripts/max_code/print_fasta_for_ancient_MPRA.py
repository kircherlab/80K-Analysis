"""

Add headers to the table of interest in order to write a fasta file

:Author: Kilian Salomon
:Contact: kilian.salomon@bih-charite.de
:Date: *25.04.2024
:Type: tool
:Input: Dataframe with matchable headers, design fasta file
:Output: Fasta file of the subset of sequences
"""

# import 
import pandas as pd
import yaml

# load helpful functions
import sys
sys.path.append('../00_helpful_functions')
import helpful_functions as hf


# config_path = "/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/global80K_config.yaml"
config_path = "/home/kisa/coding/80K_MPRA/80K-Analysis/global80K_config.yaml"
# load config file
with open(config_path, "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
# data

design_df = config['files']['final_design']['design_fasta']

df_of_interest = "/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/07_quality_control/results/ancient_mpra/filtered_tested_sequences.tsv"

# output
output_fasta = "/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/07_quality_control/results/ancient_mpra/filtered_tested_sequences.fa"


# helpful functions
def write_fasta(sequence_df, output_path, header=['header', 'sequence']):
    """
    Write the fasta file with the header and sequence in one row
    """
    with open(output_path, 'w') as f:
        for index, row in sequence_df.iterrows():
            f.write('>' + row[header[0]] + '\n' + row[header[1]] + '\n')
    return True






# write fasta file 
write_fasta(df_of_interest_merged, output_fasta, ['name_mean_log2', 'sequence'])