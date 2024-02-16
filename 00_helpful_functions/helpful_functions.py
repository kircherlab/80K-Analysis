import gzip # zipped files
import pandas as pd
from Bio import SeqIO

def read_zipped_fasta(fasta_file):
    """Read zipped fasta file with Biopython."""
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