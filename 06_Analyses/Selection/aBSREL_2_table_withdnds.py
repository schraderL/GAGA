import os
import sys
import json
from Bio import Phylo
from io import StringIO
import argparse

parser = argparse.ArgumentParser(description='Script that retrieves info and builds a table from aBSREL JSON file.')
parser.add_argument('json', type = str,
                     help = 'aBSREL JSON file')
parser.add_argument('output', type = str,
                     help = 'Name of the output table.')


args = parser.parse_args()

file = args.json
output_name = args.output


file = open(file,'r')
file = json.load(file)
barch_attributes=file['branch attributes']['0']

# Get branches
significant={}
for i in barch_attributes:
    significant[i]={"Corrected P-value":barch_attributes[i]["Corrected P-value"], 
                    "Uncorrected P-value":barch_attributes[i]["Uncorrected P-value"],
                    "Baseline MG94xREV omega ratio": barch_attributes[i]['Baseline MG94xREV omega ratio'],
                    "Full adaptive model (non-synonymous subs/site)": barch_attributes[i]['Full adaptive model (non-synonymous subs/site)'],
                    "Full adaptive model (synonymous subs/site)": barch_attributes[i]['Full adaptive model (synonymous subs/site)'],
                    "Rate classes": barch_attributes[i]['Rate classes'],
                    "omega distribution over sites": barch_attributes[i]['Rate Distributions'],
                   "omega":0,
                   "dnds":0}
    for j in barch_attributes[i]['Rate Distributions']:
        significant[i]["omega"]+=float(j[0]*j[1])
    significant[i]["dnds"]=str((float(significant[i]['Full adaptive model (non-synonymous subs/site)'])/float(significant[i]['Full adaptive model (synonymous subs/site)']))/3) # To calculate dn/ds
# Get leafs (protein names) that are significant in the branches
# read tree
tree = Phylo.read(StringIO(file['input']['trees']['0']), "newick")
for branch in significant.keys():
    clade = tree.find_elements(branch)
    clade = next(clade)
    leafs = [i.name for i in clade.get_terminals()]
    significant[branch]['Leafs']=leafs

output = open(output_name, 'w')
output.write('Name\tP-value\tUncorrected P-value\tLeafs\tBaseline omega\tOmega distribution over sites\tOmega\tRate classes\tNon-synonymous subs/site\tSynonymous subs/site\tdN/dS\n')
for i in significant:
    output.write(i+'\t'+str(significant[i]['Corrected P-value'])+'\t'
                 +str(significant[i]['Uncorrected P-value'])+'\t'+', '.join(significant[i]['Leafs'])+'\t'
                 +str(significant[i]['Baseline MG94xREV omega ratio'])+'\t'
                 +str(significant[i]['omega distribution over sites'])+'\t'+str(significant[i]['omega'])+'\t'
                 +str(significant[i]['Rate classes'])+'\t'+str(significant[i]['Full adaptive model (non-synonymous subs/site)'])+'\t'+str(significant[i]['Full adaptive model (synonymous subs/site)'])+'\t'
                 +str(significant[i]['dnds'])+'\n')
output.close()
