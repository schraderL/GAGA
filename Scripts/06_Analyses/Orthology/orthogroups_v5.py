import os
import sys
import re
import pandas as pd
import numpy as np

#module load ngs tools
#module load anaconda3/4.4.0

#usage: python3 orthogroups.py N0.tsv species.txt 0.8(Percentage of species included in the orthogroup as single copy)

# v4 multicopy now prints all 1:N (excluding single-copy)
# v5 print an additional file with all genes in the single copy orthogroups (also those with 2 or more copies in few species)

ortho_table= sys.argv[1] # orthogroup table
species= sys.argv[2] # file with species ids 
flexibility = float(sys.argv[3]) # Percentage of species included in the orthogroup as single copy/multicopy

species_list = [i.replace('\n','') for i in open(species, 'r').readlines()]
flex_len= flexibility * len(species_list)

data = pd.read_csv(ortho_table, sep='\t', dtype='string')
data_numbers= data.copy()
data_numbers = data_numbers[species_list]
first = data.columns[0]

for i in species_list:
    data_numbers[i]=data_numbers[i].str.count('\S+')
data_numbers = data_numbers.fillna(0)
c = (data_numbers.eq(1).sum(axis=1)>=flex_len)
indices1 = c[c].index.tolist()
hogs1 = data[first][indices1].tolist()

table1 = data_numbers.iloc[indices1,:].copy()
table1.insert(loc=0, column=first, value=hogs1)
table1.to_csv(species.replace(".txt","")+'_'+str(int(flexibility*100))+"percsp_orthogroups_singlecopy_counts.tsv", sep='\t', index=False)
table2= table1.copy()
for i in species_list:
    table2[i] = table2[i].apply(str)
    for j in indices1:
        if table1[i][j]==1:
            table2.at[j, i]= data[i][j]
        else:
            table2.at[j, i]= np.nan
table2.to_csv(species.replace(".txt","")+'_'+str(int(flexibility*100))+"percsp_orthogroups_singlecopy_genes.tsv", sep='\t', index=False)

table8= table1.copy()
for i in species_list:
    table8[i] = table8[i].apply(str)
    for j in indices1:
           table8.at[j, i]= data[i][j]
table8.to_csv(species.replace(".txt","")+'_'+str(int(flexibility*100))+"percsp_orthogroups_singlecopy_all_genes.tsv", sep='\t', index=False)


d = (data_numbers.gt(0).sum(axis=1)>=flex_len)
indices2 = d[d].index.tolist()
indices2 = [x for x in indices2 if x not in indices1]
hogs2 = data[first][indices2].tolist()

table3 = data_numbers.iloc[indices2,:].copy()
table3.insert(loc=0, column=first, value=hogs2)
table3.to_csv(species.replace(".txt","")+'_'+str(int(flexibility*100))+"percsp_orthogroups_multicopy_counts.tsv", sep='\t', index=False)
table4= table3.copy()
for i in species_list:
    table4[i] = table4[i].apply(str)
    for j in indices2:
        table4.at[j, i]= data[i][j]

table4.to_csv(species.replace(".txt","")+'_'+str(int(flexibility*100))+"percsp_orthogroups_multicopy_genes.tsv", sep='\t', index=False)

hogs3 = [hog for hog in data[first].tolist() if hog not in hogs1+hogs2]
table5 = data_numbers[~data_numbers.index.isin(indices1+indices2)].copy()
table5.insert(loc=0, column=first, value=hogs3)
table5.to_csv(species.replace(".txt","")+'_'+str(int(flexibility*100))+"percsp_orthogroups_unassigned_counts.tsv", sep='\t', index=False)
table6 = data[~data.index.isin(indices1+indices2)].copy()
table7 = table6[[first]+species_list].copy()

table7.to_csv(species.replace(".txt","")+'_'+str(int(flexibility*100))+"percsp_orthogroups_unassigned_genes.tsv", sep='\t', index=False)

