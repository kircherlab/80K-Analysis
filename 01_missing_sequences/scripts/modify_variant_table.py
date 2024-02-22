# Description: This script modifies the variant region map file to match the new design fasta file. It changes the IDs of the variant region map to match the new design fasta file by matching the sequences. The output is a modified variant region map file.
# Use: python modify_variant_table.py --input-design-variant-map /data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/design/final_design/results/final_design/cardiac_neuro_cava_random/variant_region_map.tsv.gz --input-new-design /fast/groups/ag_kircher/MPRA/IGVF_Y1_design/resources/association_data/design_no_duplicates_sequence_and_header.fa --input-dup-design /data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/design/final_design/results/final_design/design.fa.gz --output-variant-map /data/gpfs-1/groups/ag_kircher/work/MPRA/IGVF_Y1_design/projects/80K_MPRA/design/variant_region_map_deduplicated.tsv.gz --verbose True

import gzip # zipped files
import pandas as pd
from Bio import SeqIO
import click


@click.command()
@click.option('--input-design-variant-map',
              'variant_region_map',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Input variant region map file')
@click.option('--input-new-design',
              'new_design_path',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Input design fasta file without duplicates (new design)')
@click.option('--input-dup-design',
              'dup_design_path',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Input design fasta file with duplicates (old design)')
@click.option('--output-variant-map',
              'output_path',
              required=True,
              type=click.Path(writable=True),
              help='Output for modified variant region map')
@click.option('--verbose',
              'verbose',
              required=False,
              default=False,
              type=bool,
              help='Verbose output (number of changed rows), default False')

def cli(variant_region_map, new_design_path, dup_design_path, output_path, verbose):

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

    
    def change_IDs(missing_rows, variant_region_map, design_df, design_df_dup, column="REF_ID"):
        """
            Iterates the missing IDs and changes the entries of the variant_region_map to the new IDs
            @param missing_rows: pandas.DataFrame (missing_rows) rows with missing IDs (subset of variant_region_map)
            @param variant_region_map: pandas.DataFrame (variant_region_map) Variant, Region, REF_ID, ALT_ID
            @param design_df: pandas.DataFrame (design_df) header and sequence
            @param design_df_dup: pandas.DataFrame (design_df_dup) header and sequence
            @param column: str (column) column to read either REF_ID or ALT_ID
        """
        change_count = 0
        
        for i, row in missing_rows.iterrows():
            change_count += 1
            ref_id = row[column]
            sequence = dup_cardiac_neuro_cava_random_df[dup_cardiac_neuro_cava_random_df["header"] == ref_id]["sequence"].values[0]
            new_ref_id = cardiac_neuro_cava_random_df[cardiac_neuro_cava_random_df["sequence"] == sequence]["header"].values[0]
            varaint_region_map.at[i, column] = new_ref_id
        return change_count, variant_region_map
        
        
    def modify_variant_region_map(variant_region_map, design_df, design_df_dup, output_path, verbose=False):
        """
        Read the variant_region_map and find the IDs not matching in the new design, change them to the IDs of the new design by matching the sequences
        @param variant_region_map: pandas.DataFrame (variant_region_map) Variant, Region, REF_ID, ALT_ID
        @param design_df: pandas.DataFrame (design_df) header and sequence
        @param design_df_dup: pandas.DataFrame (design_df_dup) header and sequence
        @param output_path: str (output_path) path to write the modified variant_region_map
        @param verbose: bool (verbose) print the number of changed rows
        """
        for current_id in ["REF_ID", "ALT_ID"]:
            # find the not matching IDs
            merged_variant_region_map = varaint_region_map.merge(design_df, left_on=current_id, right_on="header", how="left")
            not_in_new_design = merged_variant_region_map[merged_variant_region_map["sequence"].isna()]

            change_count, variant_region_map = change_IDs(not_in_new_design, variant_region_map, design_df, design_df_dup, column=current_id)
            if verbose:
                print("Changed %d rows using %s"% (change_count, current_id))
            
        # write variant_region_map to file
        variant_region_map.to_csv(output_path, sep="\t", index=False)
        return variant_region_map
    
    # read map and fasta files with pandas
    varaint_region_map = pd.read_csv(varaint_region_map_path, sep="\t")

    # design file
    design_df = hf.fasta_to_dataframe(new_design_path)

    # design file with duplicates
    design_df_dup = hf.fasta_to_dataframe(dup_design_path)

    # modify and write the variant region map
    modify_variant_region_map(varaint_region_map, design_df, design_df_dup, output_path, verbose=verbose)