# Goal: input: path to a file with a list of gene names (given column name of a pandas dataframe (tsv)) and output omim information for each of the genes
import ast
import omim
from omim import util
from omim.db import Manager, OMIM_DATA
from collections import defaultdict
import matplotlib.pyplot as plt
import pandas as pd
import os
from collections import defaultdict

manager = Manager(dbfile=omim.DEFAULT_DB)

# # show columns
# print(util.get_columns_table())

# # show stats
# generated, table = util.get_stats_table(manager)
# print(generated)
# print(table)

# count the database
manager.query(OMIM_DATA).count()

cava_genes_from_designtable = ['PALB2', 'MORF4L1', 'BRCA2', 'BARD1', 'RAD51D', 'BRIP1', 'ATM', 'SDHD', 'CTCF', 'FARS2', 'SFPQ', 'AHDC1', 'CDH1', 'CYP2C19', 'CYP2D6', 'G6PD', 'F9', 'NF2', 'NBN', 'FH', 'CTNNB1', 'BAP1', 'TSC2', 'ZC4H2', 'SRC', 'LMNA', 'MYBPC3', 'MYH7', 'BCL10', 'MALT1', 'CARD11', 'XRCC2', 'IRF4', 'SDHA', 'SDHB', 'SDHC', 'SDHAF2', 'MAX', 'RAD51', 'ATR', 'RINT1', 'MRE11', 'RAD51B', 'FANCI', 'FANCM', 'SLX4', 'FLCN', 'NTHL1']

# sort list alphabetically
cava_genes_from_designtable.sort()


# query with key-value
# res = manager.query(OMIM_DATA, 'prefix', '*')
# res = manager.query(OMIM_DATA, 'mim_number', '600799')
res = manager.query(OMIM_DATA, 'hgnc_gene_symbol', 'BMPR2')
# res = manager.query(OMIM_DATA, 'geneMap', '%Pulmonary hypertension%', fuzzy=True)  # fuzzy query
def query_omim_for_phenotypes(gene):
    """
    Query the omim gene and return the dict of phenotypes
    """
    res = manager.query(OMIM_DATA, 'hgnc_gene_symbol', gene)
    print(gene)
    get_phenotypes = []
    for result in res.all():
        get_phenotypes.append(result.as_dict)


phenotype_dict = defaultdict(list)
for gene in cava_genes_from_designtable:
    res = manager.query(OMIM_DATA, 'hgnc_gene_symbol', gene)
    print(gene)
    for item in res.all():
        res_dict = item.as_dict
        if 'geneMap' in res_dict.keys():
            gene_map = res_dict['geneMap']
            if not gene_map:
                print("gene_map is empty")
                continue
            # Convert the string to a list
            gene_map_list = ast.literal_eval(gene_map)
            # iterate list and get the phenotype
            for gene_map_item in gene_map_list:
                print(gene_map_item['Phenotype'])
                # use default dict to store the phenotypes (key) and the gene (value)
                phenotype_dict[gene_map_item['Phenotype']].append(gene)
        else:
            print("gene_map is empty or does not exist")

# # Assume you have a function to query OMIM and return a dictionary of phenotypes
# def query_omim_for_genes(gene_list):
#     phenotype_dict = defaultdict(list)
#     for gene in gene_list:
#         phenotypes = query_omim(gene)  # This should return a list of phenotypes for each gene
#         for phenotype in phenotypes:
#             phenotype_dict[phenotype].append(gene)
#     return phenotype_dict

# Count genes per phenotype
phenotype_counts = {k: len(v) for k, v in phenotype_dict.items()}
print(phenotype_counts)

# # histogram of the phenotype_counts
# plt.hist(phenotype_counts.values(), bins=range(1, max(phenotype_counts.values()) + 1))
# plt.show()



# # # Plotting the pie chart only with phenotypes of more than 1 gene
# filtered_phenotype_counts = {k: v for k, v in phenotype_counts.items() if v > 1}
# plt.pie(filtered_phenotype_counts.values(), labels=filtered_phenotype_counts.keys(), autopct='%1.1f%%')
# plt.show()

# # determine the number of different genes included in the filtered phenotype counts
# gene_plotted = [gene for phenotype, gene_list in phenotype_dict.items() for gene in gene_list if phenotype in filtered_phenotype_counts.keys()]
# gene_plotted_set = set(gene_plotted)
# print(gene_plotted_set)
# print(len(gene_plotted_set))
# print('length of gene_list: ', len(cava_genes_from_designtable))
# Plotting the pie chart
plt.figure(figsize=(10, 6))
plt.pie(phenotype_counts.values(), labels=phenotype_counts.keys(), autopct='%1.1f%%')
plt.title("Gene Distribution by Phenotype")
plt.show()

# {'geneMap': '[
#     {"Location": "11q22.3", "Phenotype": "{Breast cancer, susceptibility to}", "Phenotype MIM number": "114480", "Inheritance": "AD, SMu", "Phenotype mapping key": "3"
#     },
#     {"Location": "11q22.3", "Phenotype": "Ataxia-telangiectasia", "Phenotype MIM number": "208900", "Inheritance": "AR", "Phenotype mapping key": "3"},
#     {"Location": "11q22.3", "Phenotype": "Lymphoma, B-cell non-Hodgkin, somatic", "Phenotype MIM number": "", "Inheritance": "", "Phenotype mapping key": "3"},
#     {"Location": "11q22.3", "Phenotype": "Lymphoma, mantle cell, somatic", "Phenotype MIM number": "", "Inheritance": "", "Phenotype mapping key": "3"},
#     {"Location": "11q22.3", "Phenotype": "T-cell prolymphocytic leukemia, somatic", "Phenotype MIM number": "", "Inheritance": "", "Phenotype mapping key": "3"}
#     ]',
#     'mim_type': 'gene',
#     'entrez_gene_id': '472',
#     'hgnc_gene_symbol': 'ATM', 'prefix': '*',
#     'references': '11238376, 825857, 12556884, 3093854, 9733514, 11418864, 10449794, 10192382, 8689683, 15829956, 17136093, 10571946, 12034743, 16958054, 23103869, 15459181, 10716718, 187527, 16799570, 9792409, 12673794, 10677309, 9050866, 15042666, 8789452, 9733515, 6504056, 10639175, 364941, 9259193, 3338800, 10550055, 2491181, 15054841, 16931761, 9600235, 8917548, 14695534, 27421701, 15758953, 11298456, 11850621, 12697903, 19153073, 642007, 10802669, 3200306, 2005780, 8101622, 8968760, 9497252, 8845835, 12915485, 7671310, 20966255, 11181576, 7671309, 8843191, 15496926, 22071889, 21960636, 23708966, 7637733, 9843217, 17554310, 10980530, 9781027, 22002603, 15834407, 15790808, 18650924, 10817650, 10910365, 14553952, 10766245, 9707615, 16141284, 12362033, 8786135, 17525332, 8755918, 7545545, 9771717, 118375, 11889466, 7671311, 8661102, 16906133, 16832357, 15665079, 9887333, 22345219, 11826028, 7792600, 8589678, 10706620, 10397742, 30549301, 15213104, 11805335, 18674748, 18483401, 12195425, 9463314, 11679583, 9334731, 10556216, 16652348, 16141325, 15174027, 3574400, 12086603, 12554677, 9443866, 10330348, 15928302, 19377469, 9450874, 8660985, 9521587, 9288106, 8797579, 9405657, 10783165, 9241281, 12540856, 8808599, 10839545, 16497931, 8672141, 8843194, 8843193, 7671296, 21160472, 9223307, 10839544, 16150740, 1551665', 'phenotypeMap': None, 'ensembl_gene_id': 'ENSG00000149311',
#     'generated': datetime.datetime(2021, 4, 14, 0, 0),
#     'mim_number': '607585',
#     'title': 'ATM SERINE/THREONINE KINASE; ATM'}