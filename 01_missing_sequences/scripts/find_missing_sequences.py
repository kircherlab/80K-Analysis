## script finds all missing sequences of the current assignment: from notebook: /data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/notebooks/missing_sequences_per_label.ipynb

# import:
import pandas as pd

## input:
aissigned_barcodes = '/fast/work/groups/ag_kircher/MPRA/IGVF_Y1_design/experiment/results/assignment/assignIGVFDesign/assignment_barcodes.standardConfig.sorted.tsv.gz'
all_headers = '/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/experiment/standard_results/results/assignment/standardAssignIGVFDesignNoTemp/reference/all_headers.tsv'
used_for_80K = True

# all assigned sequences with barcode
assigned_seq_df = pd.read_csv(aissigned_barcodes, sep='\t', header=None)
assigned_seq_df.columns = ['barcode', 'oligo_name', 'quality', 'number of matches']

# all sequences / headers in the design
all_seq_df = pd.read_csv(all_headers, sep='\t', header=None)
all_seq_df.columns = ['oligo_name']

# for 80K analysis: check if number of unique labels is 29 and throw error if not
if used_for_80K:
    # add second column which is split at first : and contains the label (otherwise label column is assumed)
    all_seq_df['label'] = all_seq_df['oligo_name'].str.split(':').str[0]
    if len(all_seq_df['label'].unique()) != 29:
        raise ValueError('Number of unique labels is not 29')

# print the length of the data frames:
print('\n----- Summary of the input files: -------')
print('Number of assigned sequences: ', len(assigned_seq_df))
print('Number of all sequences in the design: ', len(all_seq_df))


# left join on reference by oligo_name (to check which oligos are only in reference)
all_seq_with_assignment = all_seq_df.merge(assigned_seq_df, on='oligo_name', how='left')

# get only the missing sequence names and store them in a tsv file (name and label) (barcode is na because left join did not find a match)
missing_sequences = all_seq_with_assignment[all_seq_with_assignment['barcode'].isna()][['oligo_name', 'label']]

# print the number of missing sequences
print('\n----- Summary of the missing sequences: -------')
print('Number of missing sequences: ', len(missing_sequences))


