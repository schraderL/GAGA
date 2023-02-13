import os
import sys
from Bio import SeqIO
import pandas as pd

# Run analysis of the gene family database and retrieve the rare sequences in another file for future analysis. If the sequences in rare are not good, then they have to be removed from the gene family db.

# Input are the excel with gene families "Data/gene_families.xlsx": Edit line 11
# Second input is a directory with the sequences in fasta "Data/Gene_families/" and with format GENENAME_db.fasta: Edit line 14

df = pd.read_excel('Data/gene_families.xlsx') 
gene_families= [df['Gene family'][i].replace(" ", "") for i in range(len(df))]
# Get the gene families fasta database file:
gene_families_db = ['Data/Gene_families/'+i+'_db.fasta' for i in gene_families]

# Get mean of the proteins family
for db in gene_families_db:
    records = list(SeqIO.parse(db, "fasta"))
    mean=0
    for i in range(len(records)):
        mean+=len(records[i].seq)
    mean=mean/(len(records))

    # Get a file with rare sequences (being rare if they are smaller or higher then 60% of the average) and a file with good sequences
    rare_sequences=[]
    for i in range(len(records)):
        if len(records[i].seq)<mean*0.6 or len(records[i].seq)>mean/0.6:
            rare_sequences.append(records[i])

    rare_file = db.replace("_db.fasta","_rare_db.fasta")   
    SeqIO.write(rare_sequences, rare_file, "fasta")
