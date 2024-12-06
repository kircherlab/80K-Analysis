import pandas as pd
import os
from Bio import Entrez
import pickle

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
    # print(record)
    # Get list of PubMed IDs (PMIDs)
    pmids = record["IdList"]
    return pmids

# Function to fetch article abstract
def fetch_abstract(pmid):
    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text")
    abstract = handle.read()
    handle.close()
    return abstract


def main():

    variant_table_path = '/home/kisa/coding/80K_MPRA/literature_for_rsids/overlapping_variants_w_1000.bed'
    variant_table_path = '/home/kisa11/work/projects/80K_MPRA/literature_rsid/overlapping_variants_w_1000.bed'
    # read variant file:
    sig_variant_df = pd.read_csv(variant_table_path, sep="\t")
    sig_variant_df.columns = ['chr_a', 'start_a', 'end_a', 'id_a', 'chr_b', 'start_b', 'end_b', 'id_b']
    sig_variant_df['rsid'] = sig_variant_df['id_b'].apply(lambda elem: [elem])

    rsid_list = get_rsid_list(sig_variant_df)

    D = {}

    for rsid in rsid_list:
        if rsid not in D.keys():
            D[rsid] = {}
        pmid_list = fetch_pubmed_articles(rsid)
        D[rsid] = pmid_list
        # for pmid in pmid_list:
        #     print(f"RSID: {rsid}, PubMed ID: {pmid}")
        #     abstract = fetch_abstract(pmid)
            # output = check_abstract(rsid, abstract)
            # print([output.condition, output.rs_mentioned, output.effect, output.celltype])
            # D[rsid][pmid] = [output.condition, output.rs_mentioned, output.effect, output.celltype]

    print("OUTPUT:")
    print(D)
    with open("rsid_to_pmid.pkl", "wb") as f:
        pickle.dump(D,f)


    #
    # variant_rsid_with_literature = variant_rsid_df.loc[variant_rsid_df["rsid"].apply(lambda rsids: rsids_with_literature(rsids, target_rsids))]

if __name__ == '__main__':
    main()