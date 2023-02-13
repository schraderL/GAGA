#a How to run: $ python3 run_blast.py query database evalue
import os
import os.path
from os import path
import sys
import pandas as pd


############ INPUTS ############
db = sys.argv[1]
query = sys.argv[2]
evalue = sys.argv[3]
threads = sys.argv[4]

# NOTE: be sure that the database is decompressed: gzip -dc file.fasta.gz 

# Generate blast database from fasta file if is not already created
if str(path.isfile(db+'.pdb'))=='False':
    os.system("makeblastdb -dbtype prot -in %s" %db)

#output name = gene family name + _blast_output (for example OBP_blast_otuput)
output_name=  query.replace("_db.fasta","")+"_blast_output.txt"

# Produce tabular blastp output
os.system("blastp -query %s -db %s -out %s -evalue %s -outfmt '6 std qlen slen' -num_threads %s"% (query, db, output_name, evalue, threads))
