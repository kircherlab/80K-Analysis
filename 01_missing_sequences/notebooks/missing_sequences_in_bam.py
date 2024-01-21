# script to find missing sequences in bam files
# Date: 18-01-2024


import pandas as pd
from Bio import SeqIO
import yaml
import pysam

def fasta2pandasDF(fasta_file):
    """Read a fasta file but the sequence is in one line"""
    # read the fasta file with the sequences and prepare a tsv with header and sequence using biopython
    records = list(SeqIO.parse(fasta_file, "fasta"))
    design_df = pd.DataFrame(columns=['header', 'sequence'])
    header = [] 
    sequence = []
    for record in records:
        header.append(record.id)
        sequence.append(str(record.seq))

    design_df['header'] = header
    design_df['sequence'] = sequence
    return design_df

def write_fasta(sequence_df, output_path):
    """
    Write the fasta file with the header and sequence
    """
    with open(output_path, 'w') as f:
        for index, row in sequence_df.iterrows():
            f.write('>' + row['header'] + '\n' + row['sequence'] + '\n')
    return True
    
def read_raw_oligo_read_table(file_path=None):
    """
    Reads the read_oligo table 
    if file path == None: read example
    """
    if not file_path:
        reads = ["@NB501960:812:HH53WAFX5:1:11101:21622:1110 XI:", "@NB501960:812:HH53WAFX5:1:11101:4019:1255 XI:Z", "@NB501960:812:HH53WAFX5:1:11101:3848:1256 XI:Z",
        "@NB501960:812:HH53WAFX5:1:11101:3848:1256 XI:Z", "@NB501960:812:HH53WAFX5:1:11101:5659:1395 XI:Z"]

        oligo = [">cardiac_neuro_cava_random:REF_PPP2R2A|ENSG000", ">cardiac_neuro_cava_random:REF_ANK2|ENSG000001", ">cardiac_neuro_cava_random:REF_LMNA|ENSG000001",
        ">cardiac_neuro_cava_random:REF_ANK2|ENSG000001",
        ">MK:newcore_110746|chr13-80235318+80235587|ref"]

        reads_oligo_dict = {
            "read_name": reads,
            "oligo_name": oligo 
        }
        raw_read_oligo_df = pd.DataFrame(reads_oligo_dict)
    else:
        raw_read_oligo_df = pd.read_csv(file_path, sep="\t", header=None)
        raw_read_oligo_df.columns = ["read", "oligo", "sequence", "match_length"]
        # the oligo is ><oligo_name>_revc/forw (remove the > and split at last _)
        raw_read_oligo_df["oligo_name"] = raw_read_oligo_df["oligo"].str.split("_").str[:-1].str.join("_").str.lstrip(">")
    return raw_read_oligo_df



config_path = '/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/config/config.yml'

with open(config_path) as conf:
    config = yaml.load(conf, Loader=yaml.FullLoader)
    conf.close()

## input:
assigned_barcodes = config['files']['assigned_barcodes']
all_headers = config['files']['all_headers']

# all assigned sequences with barcode
assigned_seq_df = pd.read_csv(assigned_barcodes, sep='\t', header=None)
assigned_seq_df.columns = ['barcode', 'oligo_name', 'quality', 'number of matches']

# all sequences / headers in the design
all_seq_df = pd.read_csv(all_headers, sep='\t', header=None)
all_seq_df.columns = ['oligo_name']
# add the label to the all_seq_df
all_seq_df['label'] = all_seq_df['oligo_name'].str.split(':').str[0]

# left join on reference by oligo_name (to check which oligos are only in reference)
all_seq_with_assignment = all_seq_df.merge(assigned_seq_df, on='oligo_name', how='left')

# get only the missing sequence names and store them in a tsv file (name and label) (barcode is na because left join did not find a match)
missing_sequences = all_seq_with_assignment[all_seq_with_assignment['barcode'].isna()][['oligo_name', 'label']]

# reading the reads and match them with the missing sequences in order to get the reads matching the missing sequences with exact matches
all_exact_match = config['files']['all_exact_match']
raw_read_oligo_df = read_raw_oligo_read_table(all_exact_match) # replace None by path to all_exact_match


# match this table to the table of missing sequences by oligo_name (left join)
merged_read_missing_oligo_df = raw_read_oligo_df.merge(missing_sequences, how="left", on="oligo_name")

if config['general']['cardiac_neuro_cava_random']:
    # # filter for results with lable == cardiac_neuro_cava_random
    merged_read_missing_oligo_df = merged_read_missing_oligo_df[merged_read_missing_oligo_df["label"] == "cardiac_neuro_cava_random"]
    merged_read_missing_oligo_df = merged_read_missing_oligo_df[merged_read_missing_oligo_df["label"].notna()]

# create the dict of missing sequences per read and reads per missing sequence
missing_sequence_list_per_read = merged_read_missing_oligo_df.groupby("read")["oligo_name"].apply(list) # 1579132 reads (only cardiac_neuro_cava_random)
missing_sequence_dict_per_read = missing_sequence_list_per_read.to_dict()
read_list_per_missing_sequence = merged_read_missing_oligo_df.groupby("oligo_name")["read"].apply(list) # 2639 sequences (only cardiac_neuro_cava_random)
read_dict_per_missing_sequence = read_list_per_missing_sequence.to_dict()

# read the merged bam and the output bam of all the results from the missing sequences
merged_bam = config['files']['merged_bam']

missing_sequence_bam = config['files']['missing_sequence_bam']
if config['general']['cardiac_neuro_cava_random']:
    missing_sequence_bam =  config['files']['missing_sequence_bam_cardiac']

read_count = 0
read_set = set()
samfile = pysam.AlignmentFile(merged_bam, "rb")
missing_sequence_alignments = pysam.AlignmentFile(missing_sequence_bam, "wb", template=samfile)
print("Start reading the bam file and writing the missing sequences to a new bam file...")
for read in samfile.fetch():
    # TODO: add a function which gets the read and the sequences of interest table and returns if the read is a exact matching read or not 
    # see https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment for examples of pysam functionality
    prepared_read_name = f'@{read.query_name} {read.tags[-1][0]}:Z:{read.tags[-1][1]}'
    # try to find the prepared read name in the missing_sequence_dict_per_read
    if prepared_read_name in missing_sequence_dict_per_read.keys():
        read_count += 1
        oligo_name = missing_sequence_dict_per_read[prepared_read_name][0]
        read_set.add(oligo_name)
        # print(prepared_read_name) # debug (slows computation down)
        missing_sequence_alignments.write(read)
    if len(read.tags[-1]) != 2:
        print(f"read.tags is different than 2: {read.tags}\n{read.query_name}\n{prepared_read_name}\n")

missing_sequence_alignments.close()
samfile.close()
print(f"found {read_count} reads with missing sequences")
print(f"found {len(read_set)} different oligo_names with missing sequences: max: {len(read_dict_per_missing_sequence.keys())}")


# for later: output_samfile.write(AlignedSegment_to_be_written)