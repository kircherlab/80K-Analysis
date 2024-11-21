# import
import pandas as pd

# load helpful functions
import sys
sys.path.append('../../00_helpful_functions')
import helpful_functions as hf


closed_chromatin_wtc11 = hf.fasta_to_dataframe('/home/kisa/coding/80K_MPRA/modeling/closed_chromatin_wtc11_ngn2_270bp.fa')
# make all chars upper case
closed_chromatin_wtc11['sequence'] = closed_chromatin_wtc11['sequence'].str.upper()
closed_chromatin_wtc11

# add column: is_open = False
closed_chromatin_wtc11['is_open'] = False

closed_chromatin_wtc11.to_csv('/home/kisa/coding/80K_MPRA/modeling/closed_chromatin_wtc11_ngn2_270bp_outcome.tsv', sep="\t", index=False)