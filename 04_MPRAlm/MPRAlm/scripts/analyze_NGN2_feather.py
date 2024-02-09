import pandas as pd 
import requests
import yaml 
import matplotlib.pyplot as plt 
import pysam
import os
from pysam import VariantFile

# load config file
config_path = "/home/kisa/coding/80K_MPRA/80K-Analysis/04_MPRAlm/MPRAlm/config/config.yaml"
with open(config_path, "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)


toptable_mpralm = pd.read_feather(config["files"]["toptable"])
weights_mpralm = pd.read_feather(config["files"]["weights"])

# print(weights_mpralm)
# print(toptable_mpralm)

## 
print(f"number of tested variants: {len(toptable_mpralm)}")
print(toptable_mpralm.head())
# [21102 rows x 7 columns]
# Index(['logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B', 'variant_id'], dtype='object')
# summary of dataframe
print(toptable_mpralm.describe())
p_sig_toptable = toptable_mpralm[toptable_mpralm["adj.P.Val"] < 0.05]
print(len(p_sig_toptable)) # 150


log2_toptable = toptable_mpralm[toptable_mpralm["logFC"] > 1]
print("Log2 higher than 1: ", len(log2_toptable)) # 7 

log2_low_toptable = toptable_mpralm[toptable_mpralm["logFC"] < -1]
print("Log2 lower than -1: ", len(log2_low_toptable)) # 4



# log2 foldchange distribution histogram from significant variants
plt.hist(toptable_mpralm["logFC"], bins=100)
plt.savefig(config["files"]["log2_all_hist"])
plt.clf()


plt.hist(p_sig_toptable["logFC"], bins=100)
plt.savefig(config["files"]["log2_p_sig_hist"])
plt.clf()

# # check overlapping variant_ids between p_sig and log2

# # are variants from low log2 in p_sig?
low_log2_sig = log2_low_toptable[log2_low_toptable["adj.P.Val"] < 0.05] # 4
print("Significantly lower expression: ", len(low_log2_sig))
print(low_log2_sig.variant_id.tolist())
# 'cardiac_neuro_cava_random:DISC1|ENSG00000162946.23|EH38E2873486|1-231734633-A-G', 
# genecard: https://www.genecards.org/cgi-bin/carddisp.pl?gene=DISC1&keywords=DISC1 
# associated with Schizophrenia
# 'cardiac_neuro_cava_random:R3HDML|ENSG00000101074.5|EH38E2114680|20-44320090-C-T', 
# gene card: https://www.genecards.org/cgi-bin/carddisp.pl?gene=R3HDML&keywords=R3HDML 
# predicted to enable peptidase inhibitor activity (Peptidase activity will increase if this gene is downregulated)
# 'cardiac_neuro_cava_random:IRX3|ENSG00000177508.12|EH38E3181492|16-54316677-A-C', 
# gene card: https://www.genecards.org/cgi-bin/carddisp.pl?gene=IRX3&keywords=IRX3
# early neural deveolpment; TF; restrict generation of motor neurons to the apporprioate region of the neural tube (UniProtKB)
# 'cardiac_neuro_cava_random:NRXN1|ENSG00000179915.25|EH38E1997520|2-51080053-G-C'
# gene cards: https://www.genecards.org/cgi-bin/carddisp.pl?gene=NRXN1&keywords=NRXN1
# important for synaptic contact: membrane protein

# # are variant ids from log2 in p_sig?
high_log2_sig = log2_toptable[log2_toptable["adj.P.Val"] < 0.05] # 7
# print(high_log2_sig.variant_id.tolist())
# # ['cardiac_neuro_cava_random:TRIO|ENSG00000038382.23|EH38E2358394|5-14408059-A-G', 'cardiac_neuro_cava_random:TRIO|ENSG00000038382.23|EH38E3627037|5-14259916-C-T', 'cardiac_neuro_cava_random:ANKZF1|ENSG00000163516.14|EH38E2076021|2-219267001-T-C', 'cardiac_neuro_cava_random:SLC1A2|ENSG00000110436.13|EH38E1532674|11-35289966-C-T', 'cardiac_neuro_cava_random:CHI3L1|ENSG00000133048.13|EH38E2859237|1-203197736-G-C', 'cardiac_neuro_cava_random:LMNA|ENSG00000160789.24|EH38E1387716|1-156152803-G-A', 'cardiac_neuro_cava_random:GATA4|ENSG00000136574.19|EH38E2610941|8-11642458-G-A']
# # make vcf file from variant ids (chr - pos - ref - alt)
# vcf_prepare = [id.split("|")[-1].replace("-", "\t") for id in high_log2_sig.variant_id.tolist()]
# print(vcf_prepare)

# TRIO (neuro) (ASD SFARI) (NPCs control option) https://docs.google.com/spreadsheets/d/1UGT3mE4l_I7s351m7SvfaVhezCbAqHsT/edit#gid=1945784398
# status	gene-symbol	gene-name	ensembl-id	chromosome	genetic-category	gene-score	syndromic	eagle	number-of-reports																
# 9	TRIO	Trio Rho guanine nucleotide exchange factor	ENSG00000038382	5	Rare Single Gene Mutation, Syndromic, Functional	1	0		38																
# Gene cards: https://www.genecards.org/cgi-bin/carddisp.pl?gene=TRIO&keywords=TRIO
# GDP to GTP exchange factor => reorganization of actin (important for cell growth)
# associated with Intelllectual Developmental Disorder

# ANKZF1 (neuro) (GnomAD: https://gnomad.broadinstitute.org/gene/ENSG00000163516?dataset=gnomad_r4)
# looked in all sheets cardiac (https://docs.google.com/spreadsheets/d/1tyZmnJwkbtq_h1xn2ZGLiEMUMpVJ4-_V/edit#gid=795471055)
# looked in all sheets neurological (https://docs.google.com/spreadsheets/d/1UGT3mE4l_I7s351m7SvfaVhezCbAqHsT/edit#gid=1825322252)
# Gene cards: https://www.genecards.org/cgi-bin/carddisp.pl?gene=ANKZF1&keywords=ANKZF1
# ANKZF1 (Ankyrin Repeat And Zinc Finger Peptidyl TRNA Hydrolase 1) -> role in ribosome associated quality control (RQC) pathway
# Associated with Van Maldergem Syndrome 1 (https://www.malacards.org/card/van_maldergem_syndrome_1) fetal disease (onset at birth), neurological disease

# SLC1A2 (neuro) (GnomAD: https://gnomad.broadinstitute.org/gene/ENSG00000110436?dataset=gnomad_r4)
# Gene card: (https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLC1A2&keywords=SLC1A2) Gene encodes for solute transporter proteins (neurotransmitter transport (Glutamate) from synapsis) associated with Deveolpmental and Epileptic Encephalopathy
# Why does this variant has a high log2 fold change? This neurotransmitter transport is important in neural cells and therefore might be activated


# CHI3L1 chitinase 3 like 1 (GnomAD: https://gnomad.broadinstitute.org/gene/ENSG00000133048?dataset=gnomad_r4)
# Gene cards: https://www.genecards.org/cgi-bin/carddisp.pl?gene=CHI3L1&keywords=CHI3L1 
# Glycoprotein member and is thought to play a role in the process of inflammation and tissue remodeling. Disease associated (Schizophrenia)
# Why does this variant has a high log2 fold change?

# LMNA lamin A/C (GnomAD: https://gnomad.broadinstitute.org/gene/ENSG00000160789?dataset=gnomad_r4)
# Gene cards: https://www.genecards.org/cgi-bin/carddisp.pl?gene=LMNA&keywords=LMNA
# The lamin family of proteins make up the matrix next to the nucleus and are **highly conserved** in evolution

# GATA4 gene card: https://www.genecards.org/cgi-bin/carddisp.pl?gene=GATA4
# zinc-finger transcription factor (gene is involved in embryogenesis and myocardial differentiation)

# Go term enrichment analysis with these 6 different genes (GO: https://geneontology.org/) biological process, molecular function, cellular component
# ENSG00000038382.23
# ENSG00000163516.14
# ENSG00000110436.13
# ENSG00000133048.13
# ENSG00000160789.24
# ENSG00000136574.19)



# Expectation: more significant results with negative log2FC
sign_negative_log2FC = p_sig_toptable[p_sig_toptable["logFC"] < 0]
print("Significant but negative log2FC: ", len(sign_negative_log2FC)) # 69
sign_positive_log2FC = p_sig_toptable[p_sig_toptable["logFC"] > 0]
print("Significant but positive log2FC: ", len(sign_positive_log2FC)) # 81

sign_log2FC_non_zero = p_sig_toptable[p_sig_toptable["logFC"] != 0]
print("Significant variants with non-zero log2FC", len(sign_log2FC_non_zero))
### add the AF information to the variants
# vcf file: "/home/kisa/coding/80K_MPRA/80K-Analysis/05_variant_region_list/resources/variants_5K.vcf"

variants_5K_vcf = config["files"]["variant_vcf"]

vcf_in = VariantFile(variants_5K_vcf)
found_afs_negative = 0
found_afs_positive = 0
found_af = 0

# add af to the mpralm results with non-zero log2FC
for record in vcf_in.fetch():
    for variant in sign_log2FC_non_zero.variant_id.tolist():
        if record.id in variant:
            found_af += 1
            # add the AF information to the dataframe
            sign_log2FC_non_zero.loc[sign_log2FC_non_zero.variant_id == variant, "AF"] = record.info["AF"]
            sign_log2FC_non_zero.loc[sign_log2FC_non_zero.variant_id == variant, "AC"] = record.info["AC"]
            break
print("Found: ", found_af)
print("Expected significant results: ", len(sign_log2FC_non_zero))

# remove singleton:
# sign_log2FC_non_zero = sign_log2FC_non_zero[sign_log2FC_non_zero["AC"] != 1]

def split_into_common_rare(df_with_AF_column,  common_threshold=0.05, column="AF"):
    """Splits df into rare and common at threshold"""
    common_vars = df_with_AF_column[df_with_AF_column[column] >= common_threshold]
    rare_vars = df_with_AF_column[df_with_AF_column[column] < common_threshold]
    print(f"Found {len(common_vars)} ({len(common_vars)/len(df_with_AF_column)}) common variants and {len(rare_vars)} ({len(rare_vars)/len(df_with_AF_column)}) rare variants with the common threshold {common_threshold}")
    ultra_rare = df_with_AF_column[df_with_AF_column[column] < 0.01]
    ultra_rare = ultra_rare[ultra_rare["AC"] != 1]
    print(f"Found {len(ultra_rare)} ultra-rare variants")
    print(f"Found {len(df_with_AF_column[df_with_AF_column["AC"] == 1])} singleton variants")
    return common_vars, rare_vars

def split_df_in_pos_neg(df_with_logFC_column, column="logFC"):
    """Splits a given data frame based on the logFC into positive and negative values"""
    positive_set = df_with_logFC_column[df_with_logFC_column[column] > 0]
    negative_set = df_with_logFC_column[df_with_logFC_column[column] < 0]
    print(f"Found {len(positive_set)} ({len(positive_set)/len(df_with_logFC_column)}) positive results\nFound {len(negative_set)} ({len(negative_set)/len(df_with_logFC_column)}) negative results")
    return positive_set, negative_set

## divide the variants in common and rare variants (rare = AF<5%) # 150 variants with significant adjusted p value
common_threshold = 0.05
common_vars, rare_vars = split_into_common_rare(sign_log2FC_non_zero, common_threshold)
## check how many have positive and negative results?

### common
print("Common variants")
common_pos_vars, common_neg_vars = split_df_in_pos_neg(common_vars)
### rare
print("Rare variants")
rare_pos_vars, rare_neg_vars = split_df_in_pos_neg(rare_vars)

## do the same with rare = AF<1%
common_threshold = 0.01
print(f"Change rare definition to {common_threshold}")
common_vars, rare_vars = split_into_common_rare(sign_log2FC_non_zero, common_threshold)

### common
print("Common variants")
common_pos_vars, common_neg_vars = split_df_in_pos_neg(common_vars)
### rare
print("Rare variants")
rare_pos_vars, rare_neg_vars = split_df_in_pos_neg(rare_vars)

## How many are singleton?
print("Singleton")
print(len(sign_log2FC_non_zero[sign_log2FC_non_zero["AC"] == 1]))
# sign_singleton = sign_log2FC_non_zero[sign_log2FC_non_zero["AC"] == 1]
# singleton_pos_vars, singleton_neg_vars = split_df_in_pos_neg(sign_singleton)

# significant results: 150
# Found 35 (0.23333333333333334) common variants and 115 (0.7666666666666667) rare variants with the common threshold 0.05
# Found 65 (0.43333333333333335) common variants and 85 (0.5666666666666667) rare variants with the common threshold 0.01
# (mobil)

# # add af to the mpralm results
# toptable_mpralm
# for record in vcf_in.fetch():
#     for variant in toptable_mpralm.variant_id.tolist():
#         if record.id in variant:
#             found_af += 1
#             # add the AF information to the dataframe
#             toptable_mpralm.loc[toptable_mpralm.variant_id == variant, "AF"] = record.info["AF"]
#             break
# print("Found: ", found_af) # 
# print(len(toptable_mpralm))

# # store the dataframe with AF information
# toptable_mpralm.to_csv(config["files"]["toptable_annotated"])

# # add AF to significant results
# for record in vcf_in.fetch():
#     # for each record id try to find it in the sign_negative_log2FC or sign_positive_log2FC list and add the AF information
#     # if it is not in the list, then continue but count number of found variants
#     # if it is in the list, then add the AF information to the dataframe
#     for variant in sign_negative_log2FC.variant_id.tolist():
#         if record.id in variant:
#             found_afs_negative += 1
#             # add the AF information to the dataframe
#             sign_negative_log2FC.loc[sign_negative_log2FC.variant_id == variant, "AF"] = record.info["AF"]
#             break
    
#     for variant in sign_positive_log2FC.variant_id.tolist():
#         if record.id in variant:
#             found_afs_positive += 1
#             # add the AF information to the dataframe
#             sign_positive_log2FC.loc[sign_positive_log2FC.variant_id == variant, "AF"] = record.info["AF"]
#             break

# # print number of positive and negative variants fround in the vcf file

# print("Found negative: ", found_afs_negative) # 69
# print("Found positive: ", found_afs_positive) # 81

# # store the dataframes with AF information
# sign_negative_log2FC.to_csv(config["files"]["sign_negative"])
# sign_positive_log2FC.to_csv(config["files"]["sign_positive"])


# Found:  150
# Expected significant results:  150
# Found 35 (0.23333333333333334) common variants and 115 (0.7666666666666667) rare variants with the common threshold 0.05
# Found 36 ultra-rare variants
# Found 49 singleton variants
# Common variants
# Found 16 (0.45714285714285713) positive results
# Found 19 (0.5428571428571428) negative results
# Rare variants
# Found 65 (0.5652173913043478) positive results
# Found 50 (0.43478260869565216) negative results
# Change rare definition to 0.01
# Found 65 (0.43333333333333335) common variants and 85 (0.5666666666666667) rare variants with the common threshold 0.01
# Found 36 ultra-rare variants
# Found 49 singleton variants
# Common variants
# Found 27 (0.4153846153846154) positive results
# Found 38 (0.5846153846153846) negative results
# Rare variants
# Found 54 (0.6352941176470588) positive results
# Found 31 (0.36470588235294116) negative results
# Singleton
# 49


### Analysis of enformer high and low group

# read enformer low and try to match to significant results

# config["output_dir"]
# enformer_low_path = "/home/kisa/coding/80K_MPRA/80K-Analysis/04_MPRAlm/enformer_results/high_low_tables/cardiac.ultra-rare_variants_750_750_low.tsv"
# enformer_low = pd.read_csv(enformer_low_path, sep="\t", header=True)
# print(enformer_low.head())

# significant results
# p_sig_toptable   




### Analysis of neuro vs random and cardiac vs random
