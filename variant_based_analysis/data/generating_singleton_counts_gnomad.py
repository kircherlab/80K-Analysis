# extract all singleton variants:
import sqlite3
import pandas as pd

# NOTE: you have to download the gnomAD SQLite database like this: BEFORE you can use this script

# from gnomad_db.database import gnomAD_DB
# download_link = " https://zenodo.org/records/11077663/files/gnomad_db.sqlite3.gz?download=1"
# output_dir = "/data/cephfs-2/unmirrored/groups/kircher/users/kisa11_c/projects/population_constraint_gnomad_singlton" # database_location
# gnomAD_DB.download_and_unzip(download_link, output_dir)

# And after the script was run you should sort the resulting file to be more efficient for downstream analysis:
# zcat singleton_variants.csv.gz | awk -F',' '{OFS="\t"; print $1, $2-1, $2}' | sort -t$'\t' -k1,1 -k2,2n > sorted_singleton_variants.bed

# bedtools intersect -a all_SCREEN_GRCh38-cCREs_without_chr.bed.gz -b sorted_singleton_variants.bed -c > all_CRE_singleton_overlaps_with_counts.bed
# NOTE: code to generate all_SCREEN_GRCh38-cCREs_without_chr.bed.gz are in the generate_element_annotation_table.ipynb notebook

# Define the path to your gnomAD SQLite database file
db_path = 'gnomad_db.sqlite3'

# The table and column names based on your schema
table_name = 'gnomad_db'
allele_count_column = 'AC'

try:
    # Connect to the database
    with sqlite3.connect(db_path) as conn:
        # Define the SQL query to select all variants where the Allele Count (AC) is 1.0.
        # Use 1.0 because the AC column type is REAL (a floating-point number).
        query = f"SELECT * FROM {table_name} WHERE {allele_count_column} = 1.0;"

        # Use pandas to read the query results directly into a DataFrame
        print("Querying the database for singleton variants...")
        singleton_variants_df = pd.read_sql_query(query, conn)
        print("Query complete.")

    # Display the results
    print(f"\nFound {len(singleton_variants_df)} singleton variants.")
    print("\nFirst 5 singleton variants:")
    print(singleton_variants_df.head())

    # Optionally, save the results to a CSV file for easy access
    output_filename = 'singleton_variants.csv.gz'
    singleton_variants_df.to_csv(output_filename, index=False)
    print(f"\nSingleton variants saved to {output_filename}")

except sqlite3.OperationalError as e:
    print(f"Error: {e}")
    print("Please ensure the database file path is correct.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")