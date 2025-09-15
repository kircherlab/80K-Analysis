import pandas as pd
## write function which prepares the enformer files (getting max and dropping columns => storing again)

def get_max_and_column_enformer_file(enformer_file_path):
    """Returns the modified enformer file with max and max column"""
    enformer_file = pd.read_csv(enformer_file_path, sep='\t')
    # compute max and id of max as mohan did in combine_tsv_vcf.py
    enformer_df = enformer_file.iloc[:,:679]
    print(enformer_df.iloc[:,5:].head())
    enformer_df['DNase_max'] = enformer_df.iloc[:,5:].values.max(axis = 1)
    enformer_df['max_col'] = enformer_df.iloc[:,5:].idxmax(axis = 1)
    # drop all columns with 'DNASE' in name
    enformer_df = enformer_df.drop(enformer_df.filter(regex='DNASE').columns, axis=1)
    return enformer_df

def write_modified_enformer_file(gene_type, variant_type, enformer_out_name):
    enformer_file_name = f'/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/{gene_type}/enformer/{gene_type}.{variant_type}_680_columns.tsv.gz'
    enformer_df = get_max_and_column_enformer_file(enformer_file_name)
    print('Write to file: ', enformer_out_name)
    enformer_df.to_csv(enformer_out_name, sep='\t', index=False)

gene_list = ['cardiac', 'cava', 'neuro', 'random']
# gene_list = ['neuro']
variant_type_list = ['ultra-rare', 'singleton']

for gene_type in gene_list:
    for variant_type in variant_type_list:
        enformer_out_name = f'/data/gpfs-1/users/kisa11_c/work/coding/MPRA/IGVF_Y1_design/design/MPRA_design/results_5k/{gene_type}/enformer/{gene_type}.{variant_type}_max_values.tsv'
        write_modified_enformer_file(gene_type, variant_type, enformer_out_name)