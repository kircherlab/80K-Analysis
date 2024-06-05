import gzip # zipped files
import re # regex for header matching
import pandas as pd
from Bio import SeqIO



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


def fasta_to_dataframe(fasta_file):
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
    return df


def get_label(header):
    """
    Returns the label of a header. Expected format: <label>:<rest_of_header>
    Can be used by applying on header column of a dataframe.
    """
    return header.split(":")[0]


def is_tested(header):
    """
    Returns for a header if it is a tested sequence.
    """
    return get_label(header) == 'cardiac_neuro_cava_random'


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
    if has_variant_info_in_header(header):
        pattern = r'([A-Z]|[0-9]+)-[0-9]+-[A-Z]-[A-Z]'
        matches = re.search(pattern, header)
        if matches:
            return matches.group()
    return "NA"


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