"""
Author: Kilian Salomon <kilian.salomon@bih-charite.de
Description:
This file validates all spdis from the metadata files it finds within the given directory: metadata files are to have the following ending (.metadata.tsv.gz) but this can be changed by the user.
"""

import ast
from importlib import reload
import os
import pandas as pd
import sys
import subprocess

sys.path.append('../notebooks/helpful_functions/')
import helpful_functions as hf
reload(hf)

# column names
col_name = 'name'
col_header = 'name'
col_sequence = 'sequence'
col_category = 'category'
col_class = 'class'
col_source = 'source'
col_ref = 'ref'
col_chr = 'chr'
col_start = 'start'
col_end = 'end'
col_strand = 'strand'
col_variant_class = 'variant_class'
col_variant_pos = 'variant_pos'
col_SPDI = 'SPDI'
col_allele = 'allele'
col_info = 'info'
my_col_ref_base = 'tmp_ref_base'
my_col_alt_base = 'tmp_alt_base'


interesting_columns = [col_name, col_sequence, col_category, col_class, col_source, col_ref,
                       col_chr, col_start, col_end, col_strand, col_variant_class, col_variant_pos, col_SPDI, col_allele, col_info]

# Usage
directory_path = '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format'
metadata_file_ending = '.metadata.tsv.gz'

def list_metadata_files_in_subdirectories(directory, metadata_file_ending):
    """Returns all the files in the subdirectories of a given directory"""
    file_list = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if metadata_file_ending in file:
                file_list.append(os.path.join(root, file))
                # print(os.path.join(root, file))
    return file_list

# Function to safely evaluate string representations of lists
def safe_eval(x):
    if pd.isna(x):
        return None
    try:
        return ast.literal_eval(x)
    except (ValueError, SyntaxError):
        return x

def verify_spdi_result(verified_spdi_file, all_spdi_df, error_file):
    """"""
    verified_spdis = pd.read_csv(verified_spdi_file, header=None,sep="\t")
    verified_spdis.columns = ['input_spdi', 'verified_spdi']
    # print the rows where a second column is not in the file
    verified_spdis['is_verified'] = verified_spdis.apply(lambda row: row['input_spdi'] == row['verified_spdi'], axis=1)
    if verified_spdis.loc[~verified_spdis['is_verified']].shape[0] > 0:
        print(f'Found problem in SPDI. The output can be found in {error_file}.')
    verified_spdis.loc[~verified_spdis['is_verified']].to_csv(error_file, header=False, index=None, sep="\t")
    # check which SPDIs are missing in the second column of the resulting file => rerun them a second time
    all_spdis = set(all_spdi_df[col_SPDI].to_list())
    matchable_spdis = set(verified_spdis['verified_spdi'].to_list())
    # plot the list or let the script run again
    not_matchable_spdis = all_spdis - matchable_spdis
    return not_matchable_spdis



# file_list = [
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/C_positive_heart_AB/C_positive_heart_AB.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Glut_Chengyu/GC_Glut_Chengyu.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/C_negative_heart_MK/C_negative_heart_MK.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Liang/GC_Liang.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_GABA_Chengyu/GC_GABA_Chengyu.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/SLEA/SLEA.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/DNase_NegControls/GC_DNase_negative_brain.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/DNase_NegControls/GC_DNase_negative_brain_shuffeled.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/DNase_NegControls/GC_DNase_negative_blood.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/DNase_NegControls/GC_DNase_negative_blood_shuffeled.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Mendelian_variants/GC_Mendelian_variants.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Vista/GC_Vista.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Cort_Chengyu/GC_Cort_Chengyu.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Hon/GC_Hon.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/C_positive_heart_MK/C_positive_heart_MK.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Selvarajan/GC_Selvarajan.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Kircher/GC_Kircher.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/DNase_PosControls/GC_DNase_positive_shuffeled.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/DNase_PosControls/GC_DNase_positive.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Atrial_fib/GC_Atrial_fib.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/neuro_controls/neuro_controls.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/MK/MK.metadata.tsv.gz']
# # identify where the wrong spdis within this group come from
# file_list = ['/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/neuro_controls/C_positive_neuron_NP.metadata.tsv.gz',
# '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/neuro_controls/C_positive_neuron_MK.metadata.tsv.gz',
# '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/neuro_controls/C_negative_heart_MK.metadata.tsv.gz',
# '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/neuro_controls/C_negative_neuron_NP.metadata.tsv.gz',
# '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/neuro_controls/C_positive_heart_MK.metadata.tsv.gz',
# '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/neuro_controls/C_positive_neuron_CD.metadata.tsv.gz',]
# file_list = ['/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/neuro_controls/C_positive_neuron_CD.metadata.tsv.gz',]

# # only the indel groups:
# file_list = [
#     '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/neuro_controls/C_positive_neuron_NP.metadata.tsv.gz',
#     '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/neuro_controls/C_negative_heart_MK.metadata.tsv.gz',
# ]

# ## second run:
# file_list= ['/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/GC_Liang/GC_Liang.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/GC_Mendelian_variants/GC_Mendelian_variants.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/GC_Hon/GC_Hon.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/GC_Mohlke/GC_Mohlke.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/C_positive_heart_CAD/C_positive_heart_CAD.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/GC_Cort_Chengyu/GC_Cort_Chengyu.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/GC_Kircher/GC_Kircher.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/cardiac_neuro_cava_random/cardiac_neuro_cava_random.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/GC_Atrial_fib/GC_Atrial_fib.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/C_negative_heart_MK/C_negative_heart_MK.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/C_positive_heart_AB/C_positive_heart_AB.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/C_positive_heart_MK/C_positive_heart_MK.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/GC_Vista/GC_Vista.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/GC_GABA_Chengyu/GC_GABA_Chengyu.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/GC_Selvarajan/GC_Selvarajan.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/neuro_controls/C_negative_neuron_NP.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/neuro_controls/C_positive_neuron_NP.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/neuro_controls/C_positive_neuron_CD.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/neuro_controls/C_positive_neuron_MK.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/neuro_controls/C_negative_neuron_MK.metadata.tsv.gz',
#  '/home/kisa/coding/80K_MPRA/MPRAOligoDesign/resources/controls/GC_Glut_Chengyu/GC_Glut_Chengyu.metadata.tsv.gz']


# # on server:
# file_list = ['/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Mohlke/GC_Mohlke.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/C_positive_heart_CAD/C_positive_heart_CAD.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/C_positive_heart_AB/C_positive_heart_AB.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Glut_Chengyu/GC_Glut_Chengyu.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/C_negative_heart_MK/C_negative_heart_MK.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Liang/GC_Liang.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_GABA_Chengyu/GC_GABA_Chengyu.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/SLEA/SLEA.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/DNase_NegControls/GC_DNase_negative_brain.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/DNase_NegControls/GC_DNase_negative_brain_shuffeled.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/DNase_NegControls/GC_DNase_negative_blood.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/DNase_NegControls/GC_DNase_negative_blood_shuffeled.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Mendelian_variants/GC_Mendelian_variants.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Vista/GC_Vista.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Cort_Chengyu/GC_Cort_Chengyu.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Hon/GC_Hon.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/C_positive_heart_MK/C_positive_heart_MK.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Selvarajan/GC_Selvarajan.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Kircher/GC_Kircher.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/DNase_PosControls/GC_DNase_positive_shuffeled.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/DNase_PosControls/GC_DNase_positive.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/GC_Atrial_fib/GC_Atrial_fib.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/neuro_controls/C_negative_neuron_MK.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/neuro_controls/C_positive_neuron_NP.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/neuro_controls/C_positive_neuron_MK.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/neuro_controls/C_negative_neuron_NP.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/neuro_controls/C_positive_neuron_CD.metadata.tsv.gz',
#  '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/MK/MK.metadata.tsv.gz']

# # on server #2:
# file_list = ['/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/MK/MK.modified.metadata.tsv.gz']

def main():
    file_list = list_metadata_files_in_subdirectories(directory_path, metadata_file_ending)
    list_columns = [col_variant_class, col_variant_pos, col_SPDI, col_allele]

    for file in file_list:
        group_name = file.split('/')[-1].split(metadata_file_ending)[0]
        print(group_name)
        # load metadata file
        metadata_file = pd.read_csv(file, sep="\t")
        # Apply the safe_eval function to the specified columns
        for col in list_columns:
            metadata_file[col] = metadata_file[col].apply(safe_eval)

        # extract the spdi column of the alternatives
        alternatives_df = metadata_file.loc[metadata_file[col_allele].apply(hf.is_alternative)].copy()
        if alternatives_df.shape[0] == 0:
            continue
        output_directory = '/home/kisa/coding/80K_MPRA/80K-Analysis/07_quality_control/notebooks/control_metadata/spdi_verification' #os.path.abspath(os.getcwd())
        output_directory = '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/design/mpra_metadata_format/SPDI_verification' #os.path.abspath(os.getcwd())
        output_file = os.path.join(output_directory, f'SPDI_{group_name}.csv')
        verified_spdi_file = os.path.join(output_directory, f'verified_SPDI_{group_name}.csv')

        if "MK/MK." in file:
            # Extract all elements from the column of lists
            SPDI_list = [item for spdi_list in alternatives_df['SPDI'] for item in spdi_list if len(spdi_list) > 0]

            # Create a new DataFrame with these elements as a single column
            SPDI_df = pd.DataFrame(SPDI_list, columns=['SPDI'])
            # Write the DataFrame to a CSV file
            SPDI_df[['SPDI']].to_csv(output_file, header=False, index=None)
        else:
            alternatives_df[['SPDI']].to_csv(output_file, header=False, index=None)
        # run the spdi batch checking script
        # Run a command and capture its output
        script_name = "/home/kisa/coding/80K_MPRA/igvf_spdi_demo_mike/spdi_batch.py"
        input_arg = f"-i {output_file}"
        command = [
            "python",
            script_name,
            "-i", output_file,
            "-t", "SPDI"
        ]
        with open(verified_spdi_file, "w") as verified_file_handle:
            result = subprocess.run(command, stdout=verified_file_handle, text=True)

        spdi_error_file = os.path.join(output_directory, f'error_SPDI_{group_name}.tsv')
        not_matchable_spdis = verify_spdi_result(verified_spdi_file, alternatives_df, spdi_error_file)

        if len(not_matchable_spdis) > 0:
            output_file_second = os.path.join(output_directory, f'second_SPDI_{group_name}.csv')
            alternatives_df.loc[alternatives_df[col_SPDI].isin(not_matchable_spdis)][['SPDI']].to_csv(output_file_second, header=False, index=None)
            command = [
            "python",
            script_name,
            "-i", output_file_second,
            "-t", "SPDI"
            ]
            verified_spdi_file_second = os.path.join(output_directory, f'second_verified_SPDI_{group_name}.csv')
            with open(verified_spdi_file_second, "w") as verified_file_handle:
                result = subprocess.run(command, stdout=verified_file_handle, text=True)
            spdi_error_file_second = os.path.join(output_directory, f'second_error_SPDI_{group_name}.tsv')
            verify_spdi_result(verified_spdi_file_second, alternatives_df, spdi_error_file_second)

if __name__ == '__main__':
    main()