#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn=4
#PBS -l mem=40gb
#PBS -l walltime=6:00:00

# Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

### Here follows the user commands:
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes

# Load all required modules for the job
#First load modules
module load ngs
module load tools

# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

GENOME="/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/GAGA_all_final_assemblies_softmasked/GAGA-0001_SLR-superscaffolder_final_dupsrm_filt.softMasked.fasta"
OUTDIR="GAGA-0001" # Directory must not exist before, or it will fail (unless it is used --overwrite)

mkdir $OUTDIR
cd $OUTDIR

## Getting annotations from BITACORA already run in the re-annotation pipeline

## Getting annotations from BITACORA already run in the re-annotation pipeline
#cp /home/projects/ku_00039/people/igngod/Gene-annotation-pipeline-main/Data/Genomes/$OUTDIR/gene_families_pipeline/CD36/Step2_bitacora/CD36/CD36_genomic_and_annotated_proteins_trimmed_idseqsclustered.gff3 Bitacora.gff3
cat /home/projects/ku_00039/people/igngod/Gene-annotation-pipeline-main/Data/Genomes/$OUTDIR/gene_families_pipeline/CD36/Step2_bitacora/CD36/Intermediate_files/CD36_annot_genes.gff3 /home/projects/ku_00039/people/igngod/Gene-annotation-pipeline-main/Data/Genomes/$OUTDIR/gene_families_pipeline/CD36/Step2_bitacora/CD36/Intermediate_files/CD36_genomic_genes.gff3 > Bitacora.gff3


######## Generate a GFF3 and protein file

perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $GENOME Bitacora.gff3 Bitacora
sed s/X*$// Bitacora.pep.fasta > Bitacora.pep.fasta.tmp
mv Bitacora.pep.fasta.tmp Bitacora.pep.fasta

######## Run Interpro in the protein set
module purge
module load ngs tools
module load anaconda3/4.4.0
module load jdk/19
module load perl
module load signalp/4.1g
#module load interproscan/5.51-85.0
module load interproscan/5.52-86.0

f="Bitacora.pep.fasta"

interproscan.sh -i $f -t p -goterms -iprlookup -cpu 4     ## COMMENTED BECAUSE IT IS ALREADY RUN!!!!!


######## Run blast with CD36 and SNMPs to obtain names
module load ncbi-blast/2.11.0+

blastp -query Bitacora.pep.fasta -db /home/projects/ku_00039/people/joeviz/OR_annotation/Chemo_db/CD36_db.fasta -outfmt "6 std qlen slen" -out Bitacora.pep.fasta.CD36blast.txt -num_threads 4


######## Run script to rename the gff3 and generate the protein file and summary table

perl /home/projects/ku_00039/people/joeviz/OR_annotation/run_classification_bitacora_CD36SNMP.pl Bitacora.gff3 $OUTDIR $GENOME


    
