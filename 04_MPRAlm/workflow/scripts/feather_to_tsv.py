# script for load feather file and save as tsv file
# Usage: python feather_to_tsv.py input.feather output.tsv

import pandas as pd
import click


@click.command()
@click.option('--input-feather',
              'input_feather',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Input feather table')
@click.option('--output-tsv',
              'output_tsv',
              required=True,
              type=click.Path(writable=True),
              help='Output for tsv table')

def cli(input_feather, output_tsv):

    def load_feather_file(feather_file):
        df = pd.read_feather(feather_file)
        return df

    def save_tsv_file(df, tsv_file):
        df.to_csv(tsv_file, sep='\t', index=False)
    
    df = load_feather_file(input_feather)

    save_tsv_file(df, output_tsv)
       
if __name__ == '__main__':
    cli()