import pandas as pd
import os
from Bio import Entrez
import pickle
import sys

# useful function
def get_rsid_list(variant_list_df):
    """Concatenate all rsids within the lists of the rsid column"""
    return sum(variant_list_df["rsid"], [])

# Set your email (required by NCBI policy)
Entrez.email = "your_email@example.com"

# Function to fetch PubMed articles for an rsID
def fetch_pubmed_articles(rsid):
    query = f"{rsid}[All Fields]"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=50)  # Adjust retmax as needed
    record = Entrez.read(handle)
    handle.close()
    # Get list of PubMed IDs (PMIDs)
    pmids = record["IdList"]
    return pmids

def main():
    # Read arguments
    variant_table_path = sys.argv[1]
    chunk_index = int(sys.argv[2])
    num_chunks = int(sys.argv[3])

    # Read variant file:
    sig_variant_df = pd.read_csv(variant_table_path, sep="\t")
    sig_variant_df.columns = ['chr_a', 'start_a', 'end_a', 'id_a', 'chr_b', 'start_b', 'end_b', 'id_b']
    sig_variant_df['rsid'] = sig_variant_df['id_b'].apply(lambda elem: [elem])

    # Split rsid list into chunks
    rsid_list = get_rsid_list(sig_variant_df)
    chunk_size = len(rsid_list) // num_chunks
    rsid_chunk = rsid_list[chunk_index * chunk_size : (chunk_index + 1) * chunk_size]

    D = {}
    for rsid in rsid_chunk:
        if rsid not in D.keys():
            D[rsid] = {}
        pmid_list = fetch_pubmed_articles(rsid)
        D[rsid] = pmid_list

    # Save results for this chunk
    output_file = f"rsid_to_pmid_chunk_{chunk_index}.pkl"
    with open(output_file, "wb") as f:
        pickle.dump(D, f)

if __name__ == '__main__':
    main()

## merge:

# import pickle
# import glob

# def merge_chunks(output_file):
#     D = {}
#     for chunk_file in glob.glob("rsid_to_pmid_chunk_*.pkl"):
#         with open(chunk_file, "rb") as f:
#             chunk_data = pickle.load(f)
#             D.update(chunk_data)
#     with open(output_file, "wb") as f:
#         pickle.dump(D, f)

# merge_chunks("rsid_to_pmid_combined.pkl")
