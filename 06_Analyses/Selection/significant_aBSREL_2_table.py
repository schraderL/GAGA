import os
import sys
import json
from Bio import Phylo
from io import StringIO
import argparse

parser = argparse.ArgumentParser(description='Script that retrieves significant braches from aBSREL JSON file and builds a table.')
parser.add_argument('json', type = str,
                     help = 'aBSREL JSON file')
parser.add_argument('output', type = str,
                     help = 'Name of the output table.')

parser.add_argument("Corrected", type=float,
                    help="P-value corrected for multiple testing.")
parser.add_argument("Uncorrected", type=float,
                    help="Raw p-value without correction for multiple testing.")

args = parser.parse_args()

file = args.json
output_name = args.output
p_value = args.Corrected
uncorrected = args.Uncorrected

file = open(file,'r')
file = json.load(file)
barch_attributes=file['branch attributes']['0']

# Get significant branches
significant={}
for i in barch_attributes:
    if barch_attributes[i]["Corrected P-value"]<= p_value or barch_attributes[i]["Uncorrected P-value"]<= uncorrected:
        significant[i]={"Corrected P-value":barch_attributes[i]["Corrected P-value"], 
                        "Uncorrected P-value":barch_attributes[i]["Uncorrected P-value"]}

# Get leafs (protein names) that are significant in the branches
# read tree
tree = Phylo.read(StringIO(file['input']['trees']['0']), "newick")
for branch in significant.keys():
    clade = tree.find_elements(branch)
    clade = next(clade)
    leafs = [i.name for i in clade.get_terminals()]
    significant[branch]['Leafs']=leafs
    
# Sort by corrected p-value
def take_pvalue(elem):
    return elem[1]['Corrected P-value']
significant_sorted =sorted(list(significant.items()), key= take_pvalue)


output = open(output_name, 'w')
output.write('Name\tP-value\tUncorrected P-value\tLeafs\n')
for i in range(len(significant_sorted)):
    output.write(str(significant_sorted[i][0])+'\t'+str(significant_sorted[i][1]['Corrected P-value'])+'\t'
                 +str(significant_sorted[i][1]['Uncorrected P-value'])+'\t'
                 +', '.join(significant_sorted[i][1]['Leafs'])+'\n')
output.close()
