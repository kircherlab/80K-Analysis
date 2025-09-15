import pandas as pd
import os
import yaml


config_path = "/home/kisa/coding/80K_MPRA/80K-Analysis/global80K_config.yaml"
with open(config_path) as conf:
    config = yaml.load(conf, Loader=yaml.FullLoader)
    conf.close()

# load helpful functions
import sys
sys.path.append('../../00_helpful_functions')
import helpful_functions as hf

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


# read the file
df_elements_all_regions = pd.read_csv(config['files']['creating']['element_info_ngn2_undiffwtc11_202506'], sep='\t')
df_elements_all_regions = pd.read_csv(config['files']['creating']['element_info_ngn2_undiffwtc11_202507'], sep='\t')

# make the start and end to integers
df_elements_all_regions[col_start] = df_elements_all_regions[col_start].astype(int)
df_elements_all_regions[col_end] = df_elements_all_regions[col_end].astype(int)

df_elements_all_regions

overlap_1_3_df = pd.read_csv(
    '/home/kisa/coding/80K_MPRA/WTC11_ATAC_Ahituv/overlap_1_3_tested_sequences.bed',
    sep='\t',
    header=None,
    usecols=[3], # Get the 4th column (index 3), which is 'name'
    names=['name']
)

overlap_1_3_df
# Get unique names of sequences that overlap
overlapping_names_1_3 = overlap_1_3_df['name'].unique()

# Add the new column to your original DataFrame
df_elements_all_regions['is_open_NGN2_WTC11_ATAC_1_3_overlap'] = df_elements_all_regions['name'].isin(overlapping_names_1_3)


df_elements_no_alternative = df_elements_all_regions.loc[~df_elements_all_regions['name'].str.startswith("cardiac_neuro_cava_random:ALT_")].copy()
df_elements_no_alternative

df_elements_no_alternative_count = df_elements_no_alternative.groupby(["is_dCRE_NGN2", "is_npc_or_excitatory_SCREEN_cre", "is_open_NGN2_WTC11_ATAC_1_3_overlap"]).size().reset_index(name='count')
df_elements_no_alternative_count

columns_of_interest = ["is_dCRE_NGN2", "is_npc_or_excitatory_SCREEN_cre", "is_open_NGN2_WTC11_ATAC_1_3_overlap"]
df_elements_no_alternatives_active_no_screen_not_open = df_elements_no_alternative.loc[df_elements_no_alternative["is_dCRE_NGN2"] & (~df_elements_no_alternative["is_npc_or_excitatory_SCREEN_cre"]) & (~df_elements_no_alternative["is_open_NGN2_WTC11_ATAC_1_3_overlap"])].copy()
df_elements_no_alternatives_active_no_screen_not_open["SCREEN_ID"].nunique()

import pandas as pd
import requests
import json
from tqdm import tqdm


eh_screen_ids = list(df_elements_no_alternatives_active_no_screen_not_open["SCREEN_ID"].unique())
print(f"Number of EH38E IDs to query: {len(eh_screen_ids)}")

# --- 2. Define the GraphQL API Endpoint ---
FACTORBOOK_API_URL = "https://factorbook.api.wenglab.org/graphql"

# --- 3. Helper function to execute GraphQL queries ---
def run_graphql_query(url, query, variables=None):
    """
    Executes a GraphQL query against a specified URL.
    Handles JSON payload and error checking.
    """
    # Using 'json' parameter directly and an empty 'headers' dict as per user's working example
    payload = {"query": query}
    if variables:
        payload["variables"] = variables

    response = requests.post(url, json=payload, headers={}) # Matching user's working example
    response.raise_for_status() # Raise an HTTPError for bad responses (4xx or 5xx)
    return response.json()

# --- 4. Define the combined GraphQL query ---
# This query fetches both biosample Z-scores (via ccREBiosampleQuery)
# and cCRE group/maxZ scores (via cCREQuery) from the Factorbook API.
# It uses variables for accession (list of IDs) and assembly.

combined_query = """
  query topTissues($accession: [String!], $assembly: String!) {
    ccREBiosampleQuery(assembly: $assembly) {
      biosamples {
        sampleType
        displayname
        cCREZScores(accession: $accession) {
          score
          assay
          experiment_accession
        }
        name
        ontology
      }
    }
    cCREQuery(assembly: $assembly, accession: $accession) {
      accession
      group
      dnase: maxZ(assay: "DNase")
      h3k4me3: maxZ(assay: "H3K4me3")
      h3k27ac: maxZ(assay: "H3K27ac")
      ctcf: maxZ(assay: "CTCF")
      atac: maxZ(assay: "ATAC")
    }
  }
"""

# --- 5. Prepare base variables for the query ---
base_query_variables = {
    "assembly": "grch38" # Assuming GRCh38 as per your examples
}

print("\n--- Executing combined query against Factorbook API ---")

# --- 6. Execute the combined query with BATCH_SIZE = 1 for isolation ---
# This ensures we query one SCREEN_ID at a time as suggested by the user.
BATCH_SIZE = 1

all_biosample_data = [] # To store detailed biosample Z-scores
all_ccre_group_data = [] # To store cCRE group and maxZ scores

for i in tqdm(range(0, len(eh_screen_ids), BATCH_SIZE), desc="Querying Factorbook"):
    batch_ids = eh_screen_ids[i:i + BATCH_SIZE]
    s_id = batch_ids[0] # Get the single SCREEN_ID for this iteration

    query_variables = base_query_variables.copy()
    query_variables["accession"] = batch_ids # Pass the single ID as a list

    try:
        response_data = run_graphql_query(FACTORBOOK_API_URL, combined_query, query_variables)

        # Process cCREQuery results (group and maxZ scores)
        # We explicitly associate the result with the current s_id
        if "data" in response_data and response_data["data"].get("cCREQuery"):
            found_ccre_info = False
            for ccre_info in response_data["data"]["cCREQuery"]:
                if isinstance(ccre_info, dict) and ccre_info.get('accession') == s_id: # Ensure it's the data for our current s_id
                    all_ccre_group_data.append({
                        'SCREEN_ID': s_id, # Use the input s_id for mapping
                        'group': ccre_info.get('group'),
                        'dnase_maxZ': ccre_info.get('dnase'),
                        'h3k4me3_maxZ': ccre_info.get('h3k4me3'),
                        'h3k27ac_maxZ': ccre_info.get('h3k27ac'),
                        'ctcf_maxZ': ccre_info.get('ctcf'),
                        'atac_maxZ': ccre_info.get('atac')
                    })
                    found_ccre_info = True
                    break # Found the relevant cCRE info for this s_id

            if not found_ccre_info: # If cCREQuery data section was there, but no info for this s_id
                 all_ccre_group_data.append({'SCREEN_ID': s_id, 'group': None, 'error': f'No matching cCREQuery data found for {s_id}'})
        else: # If "data" or "cCREQuery" is missing/empty in response_data
            all_ccre_group_data.append({'SCREEN_ID': s_id, 'group': None, 'error': f'No cCREQuery section in response for {s_id}'})

        # Process ccREBiosampleQuery results (biosamples and Z-scores)
        # Explicitly associate with the current s_id
        biosample_data_for_s_id_found = False
        if "data" in response_data and response_data["data"].get("ccREBiosampleQuery"):
            for biosample_entry in response_data["data"]["ccREBiosampleQuery"].get("biosamples", []):
                if biosample_entry.get('cCREZScores'):
                    for zscore_entry in biosample_entry.get('cCREZScores'):
                        if isinstance(zscore_entry, dict): # Ensure score is for current s_id
                            all_biosample_data.append({
                                'SCREEN_ID': s_id, # Use the input s_id for mapping
                                'biosample_name': biosample_entry.get('name'),
                                'biosample_displayname': biosample_entry.get('displayname'),
                                'biosample_ontology': biosample_entry.get('ontology'),
                                'sample_type': biosample_entry.get('sampleType'),
                                'z_score': zscore_entry.get('score'),
                                'assay': zscore_entry.get('assay'),
                                'experiment_accession': zscore_entry.get('experiment_accession')
                            })
                            biosample_data_for_s_id_found = True

        if not biosample_data_for_s_id_found: # If biosample data section was there, but no scores for this s_id
            all_biosample_data.append({'SCREEN_ID': s_id, 'biosample_name': None, 'error': f'No valid Biosample Z-score data for {s_id}'})


        if "errors" in response_data: # Capture top-level GraphQL errors
            error_messages = [err.get('message', 'Unknown error') for err in response_data['errors']]
            print(f"\nGraphQL Errors for SCREEN_ID {s_id}: {error_messages}")
            # Add error entry if not already added by specific data processing
            if not any(d.get('SCREEN_ID') == s_id and 'error' in d for d in all_ccre_group_data):
                all_ccre_group_data.append({'SCREEN_ID': s_id, 'group': None, 'error': f'GraphQL Error: {error_messages}'})
            if not any(d.get('SCREEN_ID') == s_id and 'error' in d for d in all_biosample_data):
                all_biosample_data.append({'SCREEN_ID': s_id, 'biosample_name': None, 'error': f'GraphQL Error: {error_messages}'})


    except requests.exceptions.RequestException as e:
        # This catches HTTP errors like 400, 404, etc.
        print(f"\nNetwork error for SCREEN_ID {s_id}: {e}")
        all_ccre_group_data.append({'SCREEN_ID': s_id, 'group': None, 'error': f'Network Error: {e}'})
        all_biosample_data.append({'SCREEN_ID': s_id, 'biosample_name': None, 'error': f'Network Error: {e}'})
    except Exception as e:
        # Catch any other unexpected Python errors during processing
        print(f"\nAn unexpected error occurred for SCREEN_ID {s_id}: {e}")
        all_ccre_group_data.append({'SCREEN_ID': s_id, 'group': None, 'error': f'Unexpected Error: {e}'})
        all_biosample_data.append({'SCREEN_ID': s_id, 'biosample_name': None, 'error': f'Unexpected Error: {e}'})


df_ccre_groups = pd.DataFrame(all_ccre_group_data)
df_biosample_zscores = pd.DataFrame(all_biosample_data)

print("\n--- Raw cCRE Group Data (first 5 rows) ---")
# print(df_ccre_groups.head())
print(f"Shape: {df_ccre_groups.shape}")

print("\n--- Raw Biosample Z-Scores Data (first 5 rows) ---")
# print(df_biosample_zscores.head())
print(f"Shape: {df_biosample_zscores.shape}")


# --- 7. Process and Merge Data ---

# Merge cCRE group information with biosample Z-scores
# Note: A single SCREEN_ID might have multiple biosample Z-scores, so this merge will expand rows.
merged_df = pd.merge(df_biosample_zscores, df_ccre_groups, on='SCREEN_ID', how='left')

# Define target groups for filtering
target_groups = ["dELS", "pELS", "PLS"]

# Filter the merged DataFrame to include only cCREs in the target groups
# Ensure 'group' column exists and is not None before filtering
# Use .dropna() on 'group' to ensure no None values if they slipped through
filtered_by_group_df = merged_df[merged_df['group'].isin(target_groups)].copy()

# Identify active biosamples based on a Z-score threshold (e.g., >= 1.64 for ~95th percentile)
ZSCORE_THRESHOLD = 1.64
active_biosamples_df = filtered_by_group_df[filtered_by_group_df['z_score'].notna()]
active_biosamples_df = active_biosamples_df[active_biosamples_df['z_score'] >= ZSCORE_THRESHOLD]

# Group by original SCREEN_ID and collect unique active biosample display names
# This aggregates the biosample information back to one row per cCRE,
# listing all active biosamples for that cCRE within the target groups.
active_biosamples_per_ccre = active_biosamples_df.groupby('SCREEN_ID')['biosample_displayname'].apply(
    lambda x: list(x.dropna().unique())
).reset_index(name='Active_Biosamples_in_Target_Group')

# Create a final summary DataFrame for each original SCREEN_ID
final_summary_df = df_elements_no_alternatives_active_no_screen_not_open.copy() # Start with all original SCREEN_IDs

# Add the cCRE group information to the summary (deduplicated as it's one group per ID)
final_summary_df = pd.merge(final_summary_df, df_ccre_groups[['SCREEN_ID', 'group']].drop_duplicates(), on='SCREEN_ID', how='left')

# Add the active biosamples list
final_summary_df = pd.merge(final_summary_df, active_biosamples_per_ccre, on='SCREEN_ID', how='left')

# Fill NaN values for 'group' and 'Active_Biosamples_in_Target_Group' if no data was found or
# if the cCRE was not in the target groups / had no active biosamples.
# Check if 'group' column exists before trying to fillna, if it might not exist at all
if 'group' not in final_summary_df.columns:
    final_summary_df['group'] = None # Add it if missing
final_summary_df['group'].fillna('Not Found or Error', inplace=True)

if 'Active_Biosamples_in_Target_Group' not in final_summary_df.columns:
    final_summary_df['Active_Biosamples_in_Target_Group'] = None # Add it if missing
final_summary_df['Active_Biosamples_in_Target_Group'].fillna(value='No Active Biosamples in Target Groups', inplace=True)

print("\n--- Final Summary: EH38E IDs with their Group and Associated Active Biosamples ---")
print(final_summary_df.shape[0])

print("Number of missing active Biosamples: ", (final_summary_df['Active_Biosamples_in_Target_Group'] == "No Active Biosamples in Target Groups").sum())

final_summary_df.to_csv("/home/kisa/coding/80K_MPRA/cCREs_overlap/SCREEN_tissue/screen_summary_of_SCREEN_ids_without_marks_in_Neurons.csv.gz", sep="\t", index=None)