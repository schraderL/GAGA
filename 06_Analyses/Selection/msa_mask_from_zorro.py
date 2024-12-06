import os
import sys
from Bio import SeqIO

#Usage: python msa_mask_from_zorro.py MSA.fa OUTOUT_ZORRO.txt threshold output_name

multiple_seq_alg = sys.argv[1]
zorro_file = sys.argv[2]
threshold = sys.argv[3]
output = sys.argv[4]

x1 = [float(i.replace('\n','')) for i in open(zorro_file,'r').readlines()]
alg = list(SeqIO.parse(multiple_seq_alg, "fasta"))
masked_file = open(output, 'w')

for gene in range(len(alg)):
    masked_seq = ''
    for pos in range(len(str(alg[gene].seq))):
        if x1[pos]>float(threshold):
            masked_seq+=str(alg[gene].seq)[pos]
        else:
            masked_seq+='N'
    masked_file.write('>'+str(alg[gene].id)+'\n')
    masked_file.write(masked_seq+'\n')
masked_file.close()
