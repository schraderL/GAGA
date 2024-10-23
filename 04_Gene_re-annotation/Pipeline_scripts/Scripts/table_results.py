import os
import os.path
from os import path
import sys
import pandas as pd

Genome = sys.argv[1]
xlsx = sys.argv[2]
genome_dir = sys.argv[3]
# Get family names from 'gene_families.xlsx'
df = pd.read_excel(xlsx) 
# Store Gene families:
gene_families= [df['Gene family'][i].replace(" ", "") for i in range(len(df))]
average= [str(df['Average number of genes'][i]) for i in range(len(df))]
trimmed= [str(df['Trimmed output'][i]) for i in range(len(df))]
Bitacora= [df['Bitacora'][i] for i in range(len(df))]
with open(genome_dir+"/table_results.txt", "w") as fp:
    fp.write("Gene/Gene family\tStep_1_complete\tStep_1_partial\tStep_1_chimera\tStep_1_total\tNumber of annotated genes identified\tNumber of putative not annotated genes\tTotal number of identified genes (Annotated + Genomic)\tTotal number of identified genes clustering identical sequences\tTotal number of identified genes clustering highly identical sequences and proteins shorter than 30 aa\tAverage number in ant genomes\tTRIMMED\n")
    fp.close()
for i in range(len(gene_families)):
    merged = genome_dir+"/gene_families_pipeline/"+gene_families[i]+"/"+gene_families[i]+"_merged.txt"
    step1 = len(open(merged).readlines())
    complete = open(merged).read().count("complete")
    chimeras = open(merged).read().count("chimera")
    parcial = open(merged).read().count("parcial")
    # check if the output file from that Gene family exists
    result= genome_dir+"/gene_families_pipeline/"+gene_families[i]+"/Result/"+Genome+"_genecounts_summary.txt"
    if str(path.isfile(result))=='True':
        f=open(result)
        lines=f.readlines()
        gf, n1, n2, n3, n4, n5 = lines[1].split()
        # Look if the expected amount of proteins per gene/gene family is a single number or a range e.g. 1 or 0-20
        if len(average[i].split("-"))==2: # if the expected is a range, we'll choose the step that's inside the range having preference with Step2
            av=average[i].split("-")
            if int(av[0]) <= int(n5) <= int(av[1]) and int(av[0]) <= int(step1) <= int(av[1]):
                if int(n5)>=int(step1) and Bitacora[i]!="No":
                    get_from="Step 2"
                else:
                    get_from="Step 1"
            elif int(av[0]) <= int(n5) <= int(av[1]) and Bitacora[i]!="No":
                get_from="Step 2"
            elif int(av[0]) <= int(step1) <= int(av[1]):
                get_from="Step 1"
            else:
                get_from="Review"
        else: # if the expected is a single number, then we will choose between step1 or step2 depending on which is closer to the expected
            av=int(average[i])
            if av==1 and int(step1)!=1 and int(n5)!=1:
                get_from="Review"
            elif min([int(step1),int(n5)], key=lambda x:abs(x-av))==int(step1):
                if abs(int(step1)-av)>10:
                    get_from="Review"
                else:
                    get_from="Step 1"
            else:
                if abs(int(n5)-av)>10:
                    get_from="Review"
                elif Bitacora[i]=="No":
                    get_from="Review"
                else:
                    get_from="Step 2"
        with open(genome_dir+"/table_results.txt", "a") as fp:
            fp.write(gf+"\t"+str(complete)+"\t"+str(parcial)+"\t"+str(chimeras)+"\t"+str(step1)+"\t"+n1+"\t"+n2+"\t"+n3+"\t"+n4+"\t"+n5+"\t"+average[i]+"\t"+get_from+"\t"+trimmed[i]+"\n")
            fp.close()

