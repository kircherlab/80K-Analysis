## 17.02 - Optional script, took Mohans scripts and added variant type and enformer prediction result
# example row: 
# chr6	43797164	VEGFA|ENSG00000112715.26|EH38E3709352|6-43797164-G-A	G	A	4254.341	392_DNASE:CD14-positive 	1	PASS	AF=6.57168e-06;AC=1	cardiac	singleton	enformer_high
## This notebook is meant to be a module to read and filter and prepare the enformer and AF tables of mohan
# - last column: AF=...;AC=... (remove AF= and split to two column)
# - add a header
# - make it a good pandas df to merge with the header

import yaml
import os
import pandas as pd


# config 
config_path = "/home/kisa/coding/80K_MPRA/80K-Analysis/04_MPRAlm/config/config.yaml"
with open(config_path) as conf:
    config = yaml.load(conf, Loader=yaml.FullLoader)
    conf.close()

_verbose = config["verbose"]
enformer_prediction_value = config["enformer_prediction_column"] # "enformer_value"
variant_sum=5000 
higher_fraction = 0.7
lower_fraction = 0.15

## functions for this type of files
def load_df_add_columns(file_path):
    """Add column names and return dataframe"""
    enformer_af_df = pd.read_csv(file_path, sep="\t")
    enformer_af_df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "enformer_value", "max_enformer_column",  "QUAL", "FILTER", "INFO"]
    return enformer_af_df

def seperate_columns_by_pattern(df, column_with_pattern, pattern, new_column_names, drop_original_column=True):
    """Seperates columns by pattern and adds new column names"""
    extracted_df = df[column_with_pattern].str.extract(pattern)
    extracted_df.columns = new_column_names
    if drop_original_column:
        df = df.drop(column_with_pattern, axis=1)
    df = pd.concat([df, extracted_df], axis=1)
    return df

def get_AF_AC_from_INFO_column(enformer_df):
    """Adds AF and AC columns instead of INFO column"""
    df = seperate_columns_by_pattern(enformer_df, "INFO", r'AF=(?P<AF>[\d.e-]+);AC=(?P<AC>\d+)', ["AF", "AC"], drop_original_column=True)
    return df

def get_rows(df):
    """Get a specific number of rows"""
    pass

def get_higher_lower_subset(enformer_df, variant_sum=5000, higher_fraction = 0.7, lower_fraction = 0.15, column="enformer_value"):
    """Get higher 70% and lower 15% of the enformer df sorted by column such that it sums up to 5000"""
    # sort the df 
    enformer_sorted = enformer_df.sort_values(by=column, ascending=False)
    
    # get top 0.7 * 5000 rows
    top_rows = int(variant_sum * higher_fraction)
    bottom_rows = int(variant_sum * lower_fraction)
    enformer_high = enformer_sorted.iloc[:top_rows]
    
    enformer_low = enformer_sorted.iloc[len(enformer_sorted)-bottom_rows:]
    if _verbose:
        print("expected number of enformer high: ", top_rows)
        print("current number of enformer high: ", len(enformer_high))
        print("expected number of enformer low: ", bottom_rows)
        print("current number of enformer low: ", len(enformer_low))
        # check if len(enformer_low) equals bottom_rows
        if len(enformer_high) != top_rows:
            print("Top rows does not work")
        else:
            print("Top rows work")
        if len(enformer_low) == bottom_rows:
            print("bottom rows work")
        else: 
            print("bottom rows do not work")
    return enformer_high, enformer_low
    


# read as tab separated file 
for name in config["condition_names"]:
    print(f"----{name}------")
    file_name = name + "." + ".".join(config["default_name"].split(".")[1:]) # pathing the file name from the config
    enformer_af_path = os.path.join(config["enformer_results_dir"], file_name)
    print(enformer_af_path)
    # load table and add columns
    enformer_af_df = load_df_add_columns(enformer_af_path)  
    # print(enformer_af_df.head())
    ## change info column AF=...;AC=... into two columns (AF, AC)
    enformer_af_ac_df = get_AF_AC_from_INFO_column(enformer_af_df)
    # sort by ac (to test)
    enformer_sorted = enformer_af_ac_df.sort_values(by=enformer_prediction_value, ascending=False)
    
    # take num * fraction first / last elements
    enformer_high, enformer_low = get_higher_lower_subset(enformer_df=enformer_af_ac_df, column=enformer_prediction_value)
    #! TODO: enformer_high set is broken (intitial table has a length of 2500 (expected 5000)) check if you downloaded the correct one
    # write enformer_low
    low_number = ""
    high_number = ""
    if _verbose:
        # expected number _ real number
        low_number = "_" + str(int(variant_sum * lower_fraction)) + "_" + str(len(enformer_low))
        high_number = "_" + str(int(variant_sum * higher_fraction)) + "_" + str(len(enformer_high))
    output_low = ".".join(file_name.split(".")[:-1]) + f"{low_number}_low.tsv"
    output_high = ".".join(file_name.split(".")[:-1]) + f"{high_number}_high.tsv"
    enformer_low.to_csv(os.path.join(config["output_dir"], output_low), sep="\t", header=True, index=False)
    enformer_high.to_csv(os.path.join(config["output_dir"], output_high), sep="\t", header=True, index=False)

    # break # debug
    
# todo: add main function