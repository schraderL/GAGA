###########################################################
#                                                         #
#     Run gene annotation pipeline Step 1 + Bitacora      #
#                                                         #
###########################################################
import os
import sys
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Pipeline for gene family re-annotation accross hundreds of genomes.')
parser.add_argument('pipeline_dir', type = str,
                     help = 'Main directory of the pipeline that contains the README.md.')
parser.add_argument('gene_families_info', type = str,
                     help = 'Excel file containing all the gene families information.')
parser.add_argument('gene_families_db', type = str,
                     help = 'Directory containing the query protein databases (GENEFAMILY-NAME_db.fasta), where the “GENEFAMILY-NAME” label is a gene family name from the Excel file. The addition of ”_db” to the database name with its proper extension is mandatory. ej: PATH/TO/DB')
parser.add_argument('proteome', type = str,
                     help = 'File with predicted proteins in FASTA format.')
parser.add_argument('interpro', type = str,
                     help = 'File with predicted domains from InterPro in TSV format.')
parser.add_argument('gff', type = str,
                     help = 'File with structural annotations in GFF3 format.')
parser.add_argument('genome', type = str,
                     help = 'File with genomic sequences in FASTA format.')
parser.add_argument('bitacora', type = str,
                     help = 'Path to script runBITACORA_command_line.sh.')
parser.add_argument('name', type = str,
                     help = 'Name/ID of the genome ej. GAGA-0001.')
parser.add_argument('out_dir', type = str,
                     help = 'Output directory.')
parser.add_argument("threads", type=int,
                    help="Number of threads to be used.")

args = parser.parse_args()

path = args.pipeline_dir
proteome = args.proteome
interpro = args.interpro
genome_directory = args.out_dir
gff = args.gff
genome = args.genome
run_bitacora = args.bitacora
genome_name = args.name
threads = args.threads
excel_table = args.gene_families_info
gf_db = args.gene_families_db

# Get family names from 'gene_families.xlsx'
df = pd.read_excel(excel_table) 

# Store Gene families:
gene_families= [df['Gene family'][i].replace(" ", "") for i in range(len(df))]

# Get if it is possible to run blast for each gene family
blast= [df['Blast'][i] for i in range(len(df))]

# Get the gene families fasta database file:
gene_families_db = [gf_db+'/'+i+'_db.fasta' for i in gene_families]
# Copy the gene family db in a folder inside the proteome
for i in range(len(gene_families_db)):
    if blast[i]=='Yes':
        new = genome_directory+'/gene_families_pipeline/'+gene_families[i]
        os.system("mkdir -p %s && cp %s %s" % (new, gene_families_db[i], new))
        gene_families_db[i] = new + '/' + gene_families[i]+'_db.fasta'
    else:
        new = genome_directory+'/gene_families_pipeline/'+gene_families[i]
        os.system("mkdir -p %s && echo > %s" % (new, new+'/' + gene_families[i] + '_db.fasta'))
        gene_families_db[i] = new + '/' + gene_families[i]+'_db.fasta'

# Get InterPro and Pfam domains
InterPro = [df['InterPro Domain'][i] for i in range(len(df))]
Pfam = [df['Pfam domain'][i] for i in range(len(df))]

# Get Average length of the proteins of each gene family and peptide average length
prot_avg_len = [float(df['Average protein length'][i]) for i in range(len(df))]
pept_avg_len = [float(df['Peptide start length'][i]) for i in range(len(df))]

# Get e-values for each gene family to run blast
evalues = [float(df['E-value'][i]) for i in range(len(df))]
# Get if it is necessary to run bitacora for each gene family
Bitacora= [df['Bitacora'][i] for i in range(len(df))]

# run blast for every gene family if able
for i in range(len(gene_families)):
    if blast[i]=='Yes':
        os.system("python3 %s %s %s %s %s" % (path+'/Scripts/run_blast.py' ,proteome, gene_families_db[i], evalues[i], threads))
        name = gene_families_db[i].replace("_db.fasta","_blast_output.txt")
        # check if blast output is empty then create blank file
        if os.stat(name).st_size == 0:
            with open(gene_families_db[i].replace("_db.fasta","_parsed_blast.txt"), 'w') as fp:
                fp.write('delete\t1\t10\t10\tblastp\n')
                fp.write('delete\t1\t10\t10\tblastp\n')
                fp.close()
        else:
            outname = gene_families_db[i].replace("_db.fasta","")
            os.system("perl %s %s %s %s" % (path+'/Scripts/get_blastp_parsed_newv2.pl' ,name, outname, evalues[i]))
        with open(gene_families_db[i].replace("_db.fasta","_parsed_blast.txt"), 'a') as fp:
            fp.write('delete\t1\t10\t10\tblastp\n')
            fp.write('delete\t1\t10\t10\tblastp\n')
            fp.close()
    else:# create a blank file called GeneFamilyName_parsed_blast.txt
        with open(gene_families_db[i].replace("_db.fasta","_parsed_blast.txt"), 'w') as fp:
            fp.write('delete\t1\t10\t10\tblastp\n')
            fp.write('delete\t1\t10\t10\tblastp\n')
            fp.close()

# read InterproScan output
tsv = pd.read_csv(interpro, sep='\t', header=None, names = list(range(0,14)))
# Check if there is Pfam domain, if not check if there is InterPro domain. Also, if there are more than one ID then retrieve all of them.
for i in range(len(InterPro)):
    name= gene_families_db[i].replace("_db.fasta","_parsed_domain.txt")
    # Check if there is Pfam domain and Interpro
    if str(Pfam[i])[0]=='P' and str(InterPro[i])[0]=='I':
        # Obtain a table with all the domains
        l=[]
        if "&&" not in Pfam[i].split():
            for j in Pfam[i].split():
                if j[0]=='P':
                    l.append(tsv.loc[tsv[4] == j])
            new_table = pd.concat(l)[[0, 6, 7, 2, 3]]
        if "&&" in Pfam[i].split():
            for j in Pfam[i].split():
                if j[0]=='P':
                    l.append(tsv.loc[tsv[4] == j].drop_duplicates(subset=[0], keep='first'))
            new_table = pd.concat(l)
            new_table = new_table[new_table.duplicated(subset=[0], keep=False)][[0, 6, 7, 2, 3]]
        
        ll=[]
        if "&&" not in InterPro[i].split():
            for j in InterPro[i].split():
                if j[0]=='I':
                    ll.append(tsv.loc[tsv[11] == j])
            new_table2 = pd.concat(ll)[[0, 6, 7, 2, 3]]
        if "&&" in InterPro[i].split():
            for j in InterPro[i].split():
                if j[0]=='I':
                    ll.append(tsv.loc[tsv[11] == j].drop_duplicates(subset=[0], keep='first'))
            new_table2 = pd.concat(ll)
            new_table2 = new_table2[new_table2.duplicated(subset=[0], keep=False)][[0, 6, 7, 2, 3]]
        if len(new_table2)!=0 or len(new_table)!=0:
            pd.concat([new_table,new_table2]).to_csv(path_or_buf= name ,header=None, index=None, sep='\t', mode='a')
        else:  # if no matches then write a blank file
            with open(name, 'w') as fp: 
                fp.write('delete\t1\t10\t10\tPfam\n')
                fp.write('delete\t1\t10\t10\tPfam\n')
                fp.close() 
    # if there is Pfam and no Interpro:
    elif str(Pfam[i])[0]=='P' and str(InterPro[i])[0]!='I':
        # Obtain a table with all the pfam domains
        l=[]
        if "&&" not in Pfam[i].split():
            for j in Pfam[i].split():
                if j[0]=='P':
                    l.append(tsv.loc[tsv[4] == j])
            new_table = pd.concat(l)[[0, 6, 7, 2, 3]]
        else:
            for j in Pfam[i].split():
                if j[0]=='P':
                    l.append(tsv.loc[tsv[4] == j].drop_duplicates(subset=[0], keep='first'))
            new_table = pd.concat(l)
            new_table = new_table[new_table.duplicated(subset=[0], keep=False)][[0, 6, 7, 2, 3]]
        if len(new_table)!=0:
            new_table.to_csv(path_or_buf= name ,header=None, index=None, sep='\t', mode='a')
        else:  # if no matches then write a blank file
            with open(name, 'w') as fp: 
                fp.write('delete\t1\t10\t10\tPfam\n')
                fp.write('delete\t1\t10\t10\tPfam\n')
                fp.close()       
    # if there is no Pfam but Interpro:
    elif str(InterPro[i])[0]=='I' and str(Pfam[i])[0]!='P':
        # Obtain a table with all the domains
        l=[]
        if "&&" not in InterPro[i].split():
            for j in InterPro[i].split():
                if j[0]=='I':
                    l.append(tsv.loc[tsv[11] == j])
            new_table = pd.concat(l)[[0, 6, 7, 2, 3]]
        else:
            for j in InterPro[i].split():
                if j[0]=='I':
                    l.append(tsv.loc[tsv[11] == j].drop_duplicates(subset=[0], keep='first'))
            new_table = pd.concat(l)
            new_table = new_table[new_table.duplicated(subset=[0], keep=False)][[0, 6, 7, 2, 3]]
        if len(new_table)!=0:
            new_table.to_csv(path_or_buf= name ,header=None, index=None, sep='\t', mode='a')
        else:  # if no matches then write a blank file
            with open(name, 'w') as fp: 
                fp.write('delete\t1\t10\t10\tPfam\n')
                fp.write('delete\t1\t10\t10\tPfam\n')
                fp.close() 
    else: # otherwise create a blank file
        with open(name, 'w') as fp:
            fp.write('delete\t1\t10\t10\tPfam\n')
            fp.write('delete\t1\t10\t10\tPfam\n')
            fp.close()
    with open(name, 'a') as fp:
        fp.write('delete\t1\t10\t10\tPfam\n')
        fp.write('delete\t1\t10\t10\tPfam\n')
        fp.close()  


# Merge blast and interproscan results into a table and retrieve the fasta sequences
for i in range(len(gene_families)):
    blast= gene_families_db[i].replace("_db.fasta","_parsed_blast.txt")
    inter=  gene_families_db[i].replace("_db.fasta","_parsed_domain.txt")
    prot_avg= prot_avg_len[i]
    pept_avg= pept_avg_len[i]
    directory = gene_families_db[i].replace(gene_families[i]+"_db.fasta", "Result")
    os.system("mkdir -p %s" % (directory))
    output = directory+"/"+gene_families[i]+'_db.fasta'
    os.system("python3 %s %s %s %s %s %s %s" % (path+'/Scripts/merge_blast_domain.py', proteome, blast, inter, prot_avg, pept_avg, output))

# Produce GFF of the step1 results
for i in range(len(gene_families)):
    merged = gene_families_db[i].replace("_db.fasta","_merged.txt")
    output = gene_families_db[i].replace("_db.fasta","")
    os.system("perl %s %s %s %s" % (path+'/Scripts/get_annot_genes_gff_v2.pl', gff, merged, output))
    # Encode proteins and CDS from the generated GFFs
    os.system("perl %s %s %s %s" % (path+'/bitacora-master/Scripts/gff2fasta_v3.pl', genome, output+'_annot_genes.gff3', output+'_gff'))
    os.system("perl %s %s %s %s" % (path+'/bitacora-master/Scripts/gff2fasta_v3.pl', genome, output+'_annot_genes_trimmed.gff3', output+'_gfftrimmed'))
    os.system("perl %s %s %s %s" % (path+'/bitacora-master/Scripts/gff2fasta_v3.pl', genome, output+'_annot_genes.gff3', output+'_gff'))
    os.system("perl %s %s %s %s" % (path+'/bitacora-master/Scripts/gff2fasta_v3.pl', genome, output+'_annot_genes_trimmed.gff3', output+'_gfftrimmed'))

# Produce HMM profile of the gene family db (good) + the parsed sequences
for i in range(len(gene_families)):
    fasta = gene_families_db[i].replace(gene_families[i]+"_db.fasta","Result/"+gene_families[i]+"_db.fasta")
    alignment = fasta.replace("_db.fasta","_alignment.txt")
    profile = fasta.replace(".fasta",".hmm")
    os.system("mafft --auto %s > %s" % (fasta, alignment))
    os.system("hmmbuild %s  %s" % (profile, alignment))
    os.system("rm %s" % (alignment))
    if os.stat(fasta).st_size == 0:
        Bitacora[i]="No"

# Run bitacora if necessary or not if the database is empty
for i in range(len(gene_families)):
    if Bitacora[i]=="full":
        directory = gene_families_db[i].replace(gene_families[i]+"_db.fasta","Result")
        os.chdir("%s" % (directory))
        os.system("bash %s -m %s -q %s -g %s -f %s -p %s -n %s -t %s -e %s" % (run_bitacora, 'full',directory, genome, gff, proteome, genome_name, threads, evalues[i]))
    elif Bitacora[i]=="genome":
        directory = gene_families_db[i].replace(gene_families[i]+"_db.fasta","Result")
        os.chdir("%s" % (directory))
        os.system("bash %s -m %s -q %s -g %s -n %s -t %s -e %s" % (run_bitacora,'genome' ,directory, genome, genome_name, threads, evalues[i]))
        with open(genome_name+"_genecounts_genomic_proteins.txt", 'r') as fp:
            lines = fp.read().split('\n')
            vari = lines[1].split()
            fp.close()
        with open(genome_name+"_genecounts_summary.txt", 'w') as fp:
            fp.write('Bitacora runned in "genome" mode\n')
            fp.write(gene_families[i]+"\t"+vari[1]+"\t0\t"+vari[1]+"\t"+vari[2]+"\t"+vari[3]+"\n")
            fp.close()
    elif Bitacora[i]=="protein":
        directory = gene_families_db[i].replace(gene_families[i]+"_db.fasta","Result")
        os.chdir("%s" % (directory))
        os.system("bash %s -m %s -q %s -p %s -n %s -t %s -e %s" % (run_bitacora, 'protein', directory, proteome, genome_name, threads, evalues[i]))
        with open(genome_name+"_genecounts_annotated_proteins.txt", 'r') as fp:
            lines = fp.read().split('\n')
            vari = lines[1].split()
            fp.close()
        with open(genome_name+"_genecounts_summary.txt", 'w') as fp:
            fp.write('Bitacora runned in "genome" mode\n')
            fp.write(gene_families[i]+"\t"+vari[1]+"\t0\t"+vari[1]+"\t"+vari[1]+"\t"+vari[4]+"\n")
            fp.close()
    else:
        directory = gene_families_db[i].replace(gene_families[i]+"_db.fasta","Result")
        os.chdir("%s" % (directory))
        with open(genome_name+"_genecounts_summary.txt", 'w') as fp:
            fp.write('Bitacora not runned\n')
            fp.write(gene_families[i]+"\t0\t0\t0\t0\t0\n")
            fp.close()
        
os.chdir("%s" % (path))

os.system("python3 %s %s %s %s" % (path+'/Scripts/table_results.py', genome_name, excel_table, genome_directory))

# Reorganize output files
for i in range(len(gene_families)):
    directory = gene_families_db[i].replace(gene_families[i]+"_db.fasta", "Result")
    step2 = directory.replace("Result","Step2_bitacora")
    step1 = directory.replace("Result","Step1")
    os.system("mv %s %s" % (directory, step2))
    os.system("mkdir -p %s" % (step1))
    os.system("mkdir -p %s" % (step1+'/temporary_files'))
    pblast = gene_families_db[i].replace("db.fasta", "parsed_blast.txt")
    pdomain = gene_families_db[i].replace("db.fasta", "parsed_domain.txt")
    oblast = gene_families_db[i].replace("db.fasta", "blast_output.txt")
    os.system("mv %s %s" % (pblast, step1+'/temporary_files'))
    os.system("mv %s %s" % (pdomain, step1+'/temporary_files'))
    if str(os.path.isfile(oblast))=='True':
        os.system("mv %s %s" % (oblast, step1+'/temporary_files'))
    os.system("mv %s %s" % (gene_families_db[i].replace("db.fasta", "*"), step1))
