import os
import sys
import numpy as np

genome = sys.argv[1]
blast = sys.argv[2]
inter = sys.argv[3]
prot_avg = sys.argv[4]
pept_avg =  int(float(sys.argv[5]))
output = sys.argv[6]
parsed_blast = np.loadtxt(blast, dtype=object)
parsed_domain = np.loadtxt(inter, dtype=object)
gene_family = blast.replace("_parsed_blast.txt","")

# function that checks if two ranges overlap
def overlap(a,b):
    return a[0] <= b[0] <= a[1] or b[0] <= a[0] <= b[1]

parsed_blast= np.ndarray.tolist(parsed_blast)
parsed_domain= np.ndarray.tolist(parsed_domain)

# Check if splits from blastp are errors or not by checking if they overlap
for i in range(len(parsed_blast)):
    for j in range(len(parsed_blast)):
        if j!=i:
            if parsed_blast[i][0][:-1]==parsed_blast[j][0][:-1]:
                if overlap((int(float(parsed_blast[i][1])),int(float(parsed_blast[i][2]))),(int(float(parsed_blast[j][1])),int(float(parsed_blast[j][2])))):
                    parsed_blast[i][0]="delete"
for i in parsed_blast:
    if i[0]=="delete":
        parsed_blast.remove(i)
for i in range(len(parsed_blast)):
    if parsed_blast[i][0][-7:-1] == '_split':
        only_split=0
        for j in range(len(parsed_blast)):
            if parsed_blast[i][0][:-1]==parsed_blast[j][0][:-1]:
                only_split+=1
        if only_split==1:
            parsed_blast[i][0]=parsed_blast[i][0][:-7]
parsed_blast= np.asarray(parsed_blast)

splits= [(item[0],pos)  for pos, item in enumerate(parsed_blast) if item[0][-7:-1] == '_split']

# loop through domain proteins
for i in range(len(parsed_domain)):
    # loop through the splited proteins in blast
    for j in splits:
        # if the proteins names are the same and the proteins overlap
        if parsed_domain[i][0]==j[0][:-7] and overlap((int(parsed_domain[i][1]),int(parsed_domain[i][2])),(int(parsed_blast[j[1]][1]),int(parsed_blast[j[1]][2]))):
            parsed_domain[i][0]=j[0]

parsed_domain= np.asarray(parsed_domain)


# Merge the parsed blast and Interproscan outputs  
merged={}
for i in parsed_domain:
    if i[0] not in merged.keys():
        merged[i[0]]={"start":i[1], "stop":i[2], "lenght":i[3], "source":i[4]}
    else: # Obtain the longest protein from the two outputs
        if int(merged[i[0]]["start"])>int(i[1]):
            merged[i[0]]["start"]=i[1]
        if int(merged[i[0]]["stop"])<int(i[2]):
            merged[i[0]]["stop"]=i[2]
        merged[i[0]]["source"]=i[4]

for i in parsed_blast:
    if i[0] not in merged.keys():
        merged[i[0]]={"start":i[1], "stop":i[2], "lenght":i[3], "source":i[4]}
    else: # Obtain the longest protein from the two outputs
        if int(merged[i[0]]["start"])>int(i[1]):
            merged[i[0]]["start"]=i[1]
        if int(merged[i[0]]["stop"])<int(i[2]):
            merged[i[0]]["stop"]=i[2]
        merged[i[0]]["source"]="Pfam_blastp"
        
# Add label to the final table: Chimera, Complete or Parcial
for i in merged:
    # if the start position is lower than the average peptide length then start position is 1
    if int(merged[i]["start"])<pept_avg:
        merged[i]["start"]='1'
    if abs(int(merged[i]["lenght"])-int(merged[i]["stop"]))<pept_avg:
        merged[i]["stop"]=merged[i]["lenght"]
    # If two times the alignment length is lower than the subject length then it is a possible chimera
    if (int(merged[i]["stop"])-int(merged[i]["start"])+2)*2 < int(merged[i]["lenght"]):
        merged[i]['label']="chimera"
    # If the protein length is less than 60% of the gene family proteins average length 
    elif (int(merged[i]["stop"])-int(merged[i]["start"])+2)< (float(prot_avg) * 0.6):
        merged[i]['label']="parcial"
    else:
        merged[i]['label']="complete"

delete=[]
for i in merged:
    if "delete" in i:
        delete.append(i)
for i in delete:
    del merged[i]


from Bio import SeqIO
records = list(SeqIO.parse(genome, "fasta"))
for j in merged:
    if merged[j]["start"]!="1":
        if j[-7:-1]=='_split':
            for i in range(len(records)):
                if records[i].id == j[:-7]:
                    for k in range(pept_avg):
                        if records[i].seq[int(merged[j]["start"])-pept_avg+k]=="M":
                            merged[j]["start"]=str(int(merged[j]["start"])-pept_avg+k+1)
                            break
        else:
            for i in range(len(records)):
                if records[i].id == j:
                    for k in range(pept_avg):
                        if records[i].seq[int(merged[j]["start"])-pept_avg+k]=="M":
                            merged[j]["start"]=str(int(merged[j]["start"])-pept_avg+k+1)
                            break
        
# Generate table with relevant proteins and their information
f = open(gene_family+"_merged.txt", "w")
for i in merged:
    f.write(i)
    f.write("\t")
    f.write(str(merged[i]["start"]))
    f.write("\t")
    f.write(str(merged[i]["stop"]))
    f.write("\t")
    f.write(str(merged[i]["lenght"]))
    f.write("\t")
    f.write(str(merged[i]["source"]))
    f.write("\t")
    f.write(str(merged[i]["label"]))
    f.write("\n")
f.close()

# Generate a fasta file with the sequences of the proteins of interest

f = open(gene_family+"_parsed_sequences.fasta", "w")

for j in merged:
    if j[-7:-1]=='_split':
        for i in range(len(records)):
            if records[i].id == j[:-7]:
                f.write(">")
                f.write(j)
                f.write("\n")
                f.write(str(records[i].seq[int(merged[j]["start"])-1:int(merged[j]["stop"])]))
                f.write("\n")
    else:
        for i in range(len(records)):
            if records[i].id == j:
                f.write(">")
                f.write(j)
                f.write("\n")
                f.write(str(records[i].seq[int(merged[j]["start"])-1:int(merged[j]["stop"])]))
                f.write("\n")
f.close()

# Merge identified sequences with the gene family protein db good
records1 = list(SeqIO.parse(gene_family+"_parsed_sequences.fasta", "fasta"))
records2 = list(SeqIO.parse(gene_family+"_db.fasta", "fasta"))
sequences = records1+records2
SeqIO.write(sequences, output, "fasta")
