# example call: python rename_region_bed_2_design_fa.py --input-bed /path/to/region.bed --output-file /path/to/region.bed
import click
import pandas as pd
import gzip # zipped files
from Bio import SeqIO


@click.command()
@click.option('--input-file', type=click.Path(exists=True), required=True)
@click.option('--output-file', type=click.Path(readable=True), required=True)
@click.option('--file-type', default='bed', help='Type of file to rename (bed, variant_map, design_fasta).')

def cli(input_file, output_file, file_type):
    def rename_column(name):
        """
        Rename regions to be consistent with the naming in the final design fasta
        """
        if name.startswith('cardiac_neuro_cava_random'):
            return name.replace(',', '~')
        if name.startswith('GC_Selvarajan'):
            return name.replace(',', '~')
        if name.startswith('GC_Atrial_fib'):
            return name.replace(',', '~')
        if name.startswith('GC_Kircher'):
            return name.replace(',', '~')
        if name.startswith('GC_Mohlke'):
            return name.replace(',', '~')
        if name.startswith('GC_Mendelian_variants'):
            return name.replace('>', '*')
        else: return name


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


    def write_fasta(sequence_df, output_path, header=['header', 'sequence']):
        """
        Write the fasta file with the header and sequence
        """
        with open(output_path, 'w') as f:
            for index, row in sequence_df.iterrows():
                f.write('>' + row[header[0]] + '\n' + row[header[1]] + '\n')
        return True


    input_df = pd.read_csv(input_file, sep='\t', header=None)

    if file_type == 'design_fasta':
        # read design fasta file and rename the name column ("," => "~"; ">" => "*")
        input_df = fasta_to_dataframe(input_file)
        input_df['header'] = input_df['header'].apply(rename_column)
        # write file
        return write_fasta(input_df, output_file)

    if file_type == 'bed':
        # read bed file like csv and rename the name column ("," => "~") (only for tested regions)
        input_df.columns = ['region_chr', 'region_start', 'region_end', 'region_name', 'region_score', 'region_strand']

        input_df['region_name'] = input_df['region_name'].apply(rename_column)

    if file_type == 'variant_map':
        # read variant_map and rename all columns with "," to "~", ">" to "*"
        # columns: for 80K ID Region  REF_ID  ALT_ID
        input_df.columns = ['ID', 'Region', 'REF_ID', 'ALT_ID']
        # if input_df['ID'].nunique() != input_df.shape[0]:
        #     raise ValueError("ID column has duplicates")
        input_df['ID'] = input_df['ID'].apply(rename_column)
        # if input_df['ID'].nunique() != input_df.shape[0]:
        #     raise ValueError("ID column has duplicates")
        input_df['Region'] = input_df['Region'].apply(rename_column)
        input_df['REF_ID'] = input_df['REF_ID'].apply(rename_column)
        input_df['ALT_ID'] = input_df['ALT_ID'].apply(rename_column)
        # MPRAsnakeflow expects ID, REF, ALT columns
        mprasnakeflow_df = input_df[['ID', 'REF_ID', 'ALT_ID']]
        mprasnakeflow_df.columns = ['ID', 'REF', 'ALT']
        # write to file (output_file.split(".")[0] + "_mprasnakeflow_input.tsv.gz")
        mprasnakeflow_df.to_csv(output_file.split(".")[0] + "_mprasnakeflow_input.tsv.gz", sep='\t', header=True, index=False, compression='gzip')

    # write file
    input_df.to_csv(output_file, sep='\t', header=False, index=False, compression='gzip')


if __name__ == '__main__':
    cli()