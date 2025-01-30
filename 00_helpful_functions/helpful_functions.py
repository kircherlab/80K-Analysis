# """
# Author: Kilian Salomon (Kilian.Salomon@bih-charite.de)
# *: 03.2024
# """

import gzip # zipped files
import numpy as np # isna
import re # regex for header matching
import pandas as pd
from Bio import SeqIO

def get_region_info(header):
    """
    gets region id (screen enhancer identifier) from the header of tested sequences
    cardiac_neuro_cava_random:TRIO|ENSG00000038382.23|EH38E2358394|5-14408059-A-G => EH38E2358394
    """
    if 'cardiac_neuro_cava_random' not in header:
        raise ValueError('Only defined for tested headers:', header)
        return 'no_match'
    pattern = r"EH38E[\d]+" # EH38E2358394
    region_match = re.search(pattern, header)
    if region_match:
        return region_match.group()
    else:
        print('Header does not contain region info:', header)
        return 'no_match'


def get_gene_name(header, with_controls=True):
    """Returns the gene name: cardiac_neuro_cava_random:SKI|ENSG00000157933.11|EH38E2778533_fwd_tile1-1 => SKI"""
    if 'cardiac_neuro_cava_random' not in header:
        if not with_controls:
            print(header)
            raise ValueError("Function only defined for tested headers")
        return "Not given - control sequence"


    if 'ALT_' in header:
        return header.split(":ALT_")[1].split('|')[0]

    if 'REF_' in header:
        return header.split(":REF_")[1].split('|')[0]

    return header.split(':')[1].split('|')[0]


def write_fasta(sequence_df, output_path, header=['header', 'sequence']):
    """
    Write the fasta file with the header and sequence
    """
    with open(output_path, 'w') as f:
        for index, row in sequence_df.iterrows():
            f.write('>' + row[header[0]] + '\n' + row[header[1]] + '\n')
    return True


def read_zipped_fasta(fasta_file):
    """
    Read zipped fasta file with Biopython.
    """
    handle = gzip.open(fasta_file, 'rt')
    fasta_sequences = SeqIO.parse(handle,'fasta')
    return fasta_sequences


def fasta_to_dataframe(fasta_file, columns=[]):
    """
    Convert a fasta file to a pandas dataframe.
    """
    # case for ziped files:
    if fasta_file.endswith('.gz'):
        fasta_sequences = read_zipped_fasta(fasta_file)
    else:
        fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    header = []
    sequence = []
    for fasta in fasta_sequences:
        header.append(fasta.id)
        sequence.append(str(fasta.seq))
    df = pd.DataFrame({'header': header, 'sequence': sequence})
    if columns != [] and len(columns) == 2:
        df.columns = columns
    return df


def get_label(header, check_for_multi_header=False, seperator=";"):
    """
    Returns the label of a header. Expected format: <label>:<rest_of_header>
    Can be used by applying on header column of a dataframe.
    """
    if check_for_multi_header: # seperate name by seperator
        if seperator not in header:
            return header.split(":")[0]
        return [get_label(name) for name in header.split(seperator)]

    return header.split(":")[0]


def is_tested(header):
    """
    Returns for a header if it is a tested sequence.
    """
    label = get_label(header, check_for_multi_header=True, seperator=";")
    if isinstance(label, list):
        return 'cardiac_neuro_cava_random' in label

    return label == 'cardiac_neuro_cava_random'


def is_control(header):
    """
    Returns for a header if it is from a control sequence.
    """
    return not is_tested(header)

### Function to get a unique identifier for variants to merge unique variants to the variant vcf or the mpralm results
def has_variant_info_in_header(header):
    """Check if the header has the variant info (chr-pos-ref-alt) in the ending of the header"""
    pattern = r'([A-Z]|[0-9]+)-[0-9]+-[A-Z]-[A-Z]'
    matches = re.search(pattern, header)
    if matches:
        return True
    else: return False


def get_chrom_pos_ref_alt_pattern(header):
    """Extracts the chrom pos ref alt information from the end of the header"""
    if header is None:
        return "NA"
    if has_variant_info_in_header(header):
        pattern = r'([A-Z]|[0-9]+)-[0-9]+-[A-Z]-[A-Z]'
        matches = re.search(pattern, header)
        if matches:
            return matches.group()
    return "NA"


def is_reference(allele):
    """
    Since allele is list: returns if 'ref' in allele list
    - Different cases:
      - NaN (float)
      - list of either 'ref' or 'alt'
    """
    current_allele = allele

    # check case if allele is NA (float) or None
    if allele is None or allele == 'NA':
        return False
    if isinstance(allele, float):
            return False
    if type(allele) is list:
        current_allele = allele[0]

    if 'ref' in current_allele:
        return True
    return False


def is_alternative(allele):
    """
    Since allele is list: returns if 'alt' in allele list
    """
    current_allele = allele

    # check case if allele is NA (float) or None
    if allele is None or allele == 'NA':
        return False
    if isinstance(allele, float):
            return False
    if type(allele) is list:
        current_allele = allele[0]
    if 'alt' in current_allele:
        return True
    return False


### Functions for interacting with enformer information
def get_enformer_class_list(enformer_variant_info):
    """
    return the enformer class list
    enformer_variant_info: [('gene_set', 'variant_type', 'enformer_class'), ]
    """
    enformer_class_list = []
    for elem in enformer_variant_info:
        enformer_class_list.append(elem[2])
    return list(set(enformer_class_list))

def get_variant_type_list(enformer_variant_info):
    """
    return the variant type list
    enformer_variant_info: [('gene_set', 'variant_type', 'enformer_class'), ]
    """
    variant_type_list = []
    for elem in enformer_variant_info:
        variant_type_list.append(elem[1])
    return list(set(variant_type_list))

def get_gene_set_list(enformer_variant_info):
    """
    return the geme set list
    enformer_variant_info: [('gene_set', 'variant_type', 'enformer_class'), ]
    """
    gene_set_list = []
    for elem in enformer_variant_info:
        gene_set_list.append(elem[0])
    return list(set(gene_set_list))


def is_element_of(column_list, pattern):
    """Check for the pattern in the column_list"""
    return pattern in column_list


def split_ids(row, id_col, separator=';'):
    """Function to split the name and create new rows while conserving all columns"""
    ids = row[id_col].split(separator)
    new_rows = []
    for id in ids:
        new_row = row.copy()
        new_row[id_col] = id
        new_rows.append(new_row)
    return pd.DataFrame(new_rows)