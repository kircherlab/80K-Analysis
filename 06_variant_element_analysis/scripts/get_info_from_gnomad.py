import os
import pandas as pd
from gnomad_db.database import gnomAD_DB
import click
import numpy
print(numpy.version.version)

# usage: python get_info_from_gnomad.py --gnomadLocation /data/cephfs-1/home/users/kisa11_c/unmirrored/projects/MPRA/IGVF_Y1_design/projects/gnomadDB/v3.1.2/ --inputPath /data/cephfs-1/home/users/kisa11_c/unmirrored/projects/MPRA/IGVF_Y1_design/projects/gnomadDB/all_80K_variants_for_gnomadDB.tsv --outputPath /data/cephfs-1/home/users/kisa11_c/unmirrored/projects/MPRA/IGVF_Y1_design/projects/gnomadDB/all_80K_variants_for_gnomadDB_with_gnomad_info.tsv --infoColumns "AC, AF, AF_popmax, AF_eas, AF_nfe, AF_fin, AF_afr, AF_asj"
# usage: python get_info_from_gnomad.py --gnomadLocation /data/cephfs-1/scratch/groups/kircher/MPRA/IGVF_Y1_design/projects/gnomad_db/v3.1.2/ --inputPath /data/cephfs-1/scratch/groups/kircher/MPRA/IGVF_Y1_design/projects/gnomad_db/all_80K_variants_for_gnomadDB.tsv --outputPath /data/cephfs-1/scratch/groups/kircher/MPRA/IGVF_Y1_design/projects/gnomad_db/all_80K_variants_for_gnomadDB_with_gnomad_info.tsv --infoColumns "AC, AF, AF_popmax, AF_eas, AF_nfe, AF_fin, AF_afr, AF_asj"

@click.command()
@click.option(
    "--gnomadLocation",
    "database_location",
    required=True,
    multiple=False,
    type=click.Path(exists=True, readable=True),
    default="/data/cephfs-1/home/users/kisa11_c/unmirrored/projects/MPRA/IGVF_Y1_design/projects/gnomadDB/v3.1.2/",
    help="Path to the gnomad database location (directory)",
)
@click.option(
    "--inputPath",
    "input_path",
    required=True,
    multiple=False,
    type=click.Path(exists=True, readable=True),
    default="/data/cephfs-1/home/users/kisa11_c/unmirrored/projects/MPRA/IGVF_Y1_design/projects/gnomadDB/all_80K_variants_for_gnomadDB.tsv",
    help="Input path to the tsv file with variants to query(cols: header, chrom, pos, ref, alt)",
)
@click.option(
    "--outputPath",
    "output_path",
    required=True,
    multiple=False,
    type=click.Path(writable=True),
    default="/data/cephfs-1/home/users/kisa11_c/unmirrored/projects/MPRA/IGVF_Y1_design/projects/gnomadDB/all_80K_variants_for_gnomadDB_with_gnomad_info.tsv",
    help="Input path to the tsv file with variants to query(cols: header, chrom, pos, ref, alt)",
)
@click.option(
    "--infoColumns",
    "info_columns",
    required=False,
    multiple=False,
    type=str,
    default="AC, AF, AF_popmax, AF_eas, AF_nfe, AF_fin, AF_afr, AF_asj",
    help="Information you want to query from gnomad database (see here: https://github.com/KalinNonchev/gnomAD_DB/blob/master/gnomad_db/pkgdata/gnomad_columns.yaml)",
)

def cli(database_location, input_path, output_path, info_columns):
    # database_location = "/data/cephfs-1/home/users/kisa11_c/unmirrored/projects/MPRA/IGVF_Y1_design/projects/gnomadDB/v3.1.2/"
    db = gnomAD_DB(database_location, gnomad_version="v3")
    # input_path = os.path.join('/data/cephfs-1/home/users/kisa11_c/unmirrored/projects/MPRA/IGVF_Y1_design/projects/gnomadDB', 'cava_variants_for_gnomadDB.tsv')
    var_df = pd.read_csv(input_path, sep='\t')
    var_df = db.get_info_from_df(var_df, info_columns)

    print(var_df.head())

    var_df.to_csv(output_path, sep='\t', index=False)

    # # pass dir
    # db = gnomAD_DB(database_location, gnomad_version="v3")
    # # read in tsv of variants expected to have id, chrom, pos, ref, alt
    # var_df = pd.read_csv(input_path, sep='\t')

    # # query from dataframe AF column (see here: https://github.com/KalinNonchev/gnomAD_DB/blob/master/gnomad_db/pkgdata/gnomad_columns.yaml)
    # info_columns = "AC, AF, AF_popmax, AF_eas, AF_nfe, AF_fin, AF_afr, AF_asj" # east asian, non-finnish european, finnish, african/african-american, ashkenazi jewish
    # var_df['AF'] = db.get_info_from_df(var_df, info_columns)

    # output_path = os.path.join('/data/cephfs-1/home/users/kisa11_c/unmirrored/projects/MPRA/IGVF_Y1_design/projects/gnomadDB', 'cava_variants_for_gnomadDB_with_AF.tsv')
    # var_df.to_csv(output_path, sep='\t', index=False)

    # # query from dataframe AF and AF_popmax columns
    # # print(db.get_info_from_df(dummy_var_df, "AF, AF_popmax"))

if __name__ == "__main__":
    cli()