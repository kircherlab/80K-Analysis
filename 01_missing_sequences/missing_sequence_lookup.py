## python script with command line input for lookup of missing sequences in reads
# expected call: missing_sequence_lookup.py <path to reads> <sequence length> <recompute (True/[False])> <output file path>


# imports
import sys
import gzip
import matplotlib.pyplot as plt
import pandas as pd
import json 
import yaml


# load reference fasta specified in config
config_path = '/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/config.yaml'
assignment_name = 'assignIGVFDesignNoTemp'

with open(config_path) as conf:
    config = yaml.load(conf, Loader=yaml.FullLoader)
    conf.close()

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# command line input
_reads_path = sys.argv[1]
_seq_length = int(sys.argv[2])
_recompute = str2bool(sys.argv[3])

_output_file_path = sys.argv[4]

_all_ref_seqs = '/fast/work/groups/ag_kircher/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesign/all_ref_sequences.tsv'
_all_assignment = '/fast/work/groups/ag_kircher/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesign/assignment_barcodes.standardConfig.sorted.tsv.gz'
_default_seq_dict_path = '/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/missing_seq_dict.json'
_ref_fasta_path = config['assignments'][assignment_name]['reference'] 


# functions
def prepare_missing_seq_dict(ref_fasta_path, missing_seqs_list, max_seq_number=145): # 5084
    with open(ref_fasta_path) as file:
    # iterate over lines in fasta file and put the fasta sequences in a dict where the key is the header and included in missing_seqs and the value is the sequence
        missing_seq_dict = {}
        header = ''
        seq = ''
        included = False
        for line in file:
            line = line.rstrip()
            if line.startswith('>'): 
                if line[1:] in missing_seqs_list:
                    header = line
                    missing_seq_dict[header] = ''
                    included = True
                else:
                    included = False
                    continue
            elif included:
                missing_seq_dict[header] += line.upper()
    return missing_seq_dict

def build_missing_seq_dict(ref_fasta_path, max_seq_number=145):
    '''Prepares a list of missing sequences from the reference fasta file'''
    # load reference (names of oligos) and the assignment file
    all_seq_df = pd.read_csv(_all_ref_seqs, header=None, sep='\t')
    all_seq_df.columns = ['oligo_name']
    all_seq_df['label'] = all_seq_df['oligo_name'].str.split(':').str[0]

    # all assigned sequences with barcode
    assigned_seq_df = pd.read_csv(_all_assignment, sep='\t', header=None)
    assigned_seq_df.columns = ['barcode', 'oligo_name', 'quality', 'number of matches']

    # left join to get the missing sequences
    all_seq_with_assignment = all_seq_df.merge(assigned_seq_df, on='oligo_name', how='left')
    missing_seqs = all_seq_with_assignment.loc[all_seq_with_assignment['barcode'].isna()]
    # drop all columns except oligo_name and label
    missing_seqs = missing_seqs.drop(['barcode', 'quality', 'number of matches'], axis=1)
    # missing sequence headers as list and prepare dict (with sequences) 
    missing_seqs_list = missing_seqs['oligo_name'].values.tolist()
    missing_seq_dict = prepare_missing_seq_dict(ref_fasta_path, missing_seqs_list, max_seq_number)
    return missing_seq_dict


def prep_query_seq(seq, length, reverse=False):
    '''
    Prepares the query sequence for the search
    @param: sequence with full length; length: disired length; reverse: bool if reverse sequence is required
    @output: query sequence
    '''
    if reverse:
        query_seq = seq[::-1]
    return seq[:length]
        
# feature assignment exact matches 
def exact_match(fastq_file_path, seq_dict, seq_length, reverse=False):
    '''Tries to find all sequences in one read sequence. Iterates each missing sequence for an exact match. The sequence length is given by the user.'''
    seperation = ''.join(['-'] * 20)
    counter = 0
    with gzip.open(fastq_file_path, 'rt') as reads_fastq:
        with open(_output_file_path, 'w') as output_file:
            for line in reads_fastq:
                if line.startswith('@'):
                    header = line
                elif line.startswith('+'):
                    pass
                else:
                    read_seq = line
                    for header, seq in seq_dict.items():
                        query_seq = prep_query_seq(seq, seq_length, reverse)
                        if query_seq in read_seq:
                            counter += 1
                            output_file.write(f'{seperation}\n\nFound match\n')
                            output_file.write(f'In header: {header}\n')
                            output_file.write(f'match of length {seq_length} of\nquery: {header}\n\nseq:{query_seq}\n\n')
                            output_file.write(f'in sequence: {read_seq}\n')
                            output_file.write(f'{seperation}\n')
                            # print(f'{seperation}\n\nFound match\n')
                            # print('In header: ',header)
                            # print(f'match of length {seq_length} of ', query_seq)
                            # print('in sequence: ', seq)
                            # print(f'{seperation}')
            output_file.write(f'Found {counter} matches in {fastq_file_path}\n')
    print(f'Found {counter} matches')
# main
def main():
    # compute missing sequence dict 
    if _recompute:
        missing_seq_dict = build_missing_seq_dict(_ref_fasta_path, 145)
    else:
        missing_seq_dict = json.load(open(_default_seq_dict_path,))

    exact_match(_reads_path, missing_seq_dict, _seq_length, False)    
    #! Start it in parallel on the cluster

if __name__ == '__main__':
    main()