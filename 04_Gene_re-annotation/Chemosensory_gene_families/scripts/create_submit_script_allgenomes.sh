#!/bin/bash

# Script to create a job array to run the pipeline in each genome. We use the GAGA-ID as filename, and the specific paths need to be modified in the script.
# Note that this example is set with the CSP

FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/GAGA_all_final_assemblies_softmasked/*fasta # Folder containing the genome assemblies in fasta

i=0

for f in $FILES
do
#    FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+).fasta/){print "$1"}' $f)"
	FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*(GAGA-\d\d\d\d)\S+/ || $ARGV[0] =~ /.*(NCBI-\d\d\d\d)\S+/ || $ARGV[0] =~ /.*(OUT-\d\d\d\d)\S+/){print "$1"}' $f)"
#    mkdir -p $FILENAME

    echo '#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn=4
#PBS -l mem=40gb
#PBS -l walltime=6:00:00' > submit_script_$FILENAME.sh

echo '
# Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

### Here follows the user commands:
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes

# Load all required modules for the job
#First load modules
module load ngs tools
module load anaconda3/4.4.0
module load jdk/19
module load perl
module load signalp/4.1g
#module load interproscan/5.51-85.0
module load interproscan/5.52-86.0
module load ncbi-blast/2.11.0+


# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

GENOME="'"$f"'"
OUTDIR="'"$FILENAME"'" # Directory must not exist before, or it will fail (unless it is used --overwrite)

mkdir $OUTDIR
cd $OUTDIR

# Getting the annotation from BITACORA that are run the re-annotation pipeline
cp /home/projects/ku_00039/people/igngod/Gene-annotation-pipeline-main/Data/Genomes/$OUTDIR/gene_families_pipeline/CSP/Step2_bitacora/CSP/CSP_genomic_and_annotated_proteins_trimmed_idseqsclustered.gff3 Bitacora_raw.gff3
# Using untrimmed proteins
#cat /home/projects/ku_00039/people/igngod/Gene-annotation-pipeline-main/Data/Genomes/$OUTDIR/gene_families_pipeline/CSP/Step2_bitacora/CSP/Intermediate_files/CSP_annot_genes.gff3 /home/projects/ku_00039/people/igngod/Gene-annotation-pipeline-main/Data/Genomes/$OUTDIR/gene_families_pipeline/CSP/Step2_bitacora/CSP/Intermediate_files/CSP_genomic_genes.gff3 > Bitacora_raw.gff3


######## Generate a GFF3 and protein file

sed 's/split/separated/g' Bitacora_raw.gff3 > Bitacora_raws.gff3

perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $GENOME Bitacora_raws.gff3 Bitacora_raws
sed 's/X*$//' Bitacora_raws.pep.fasta > Bitacora_raws.pep.fasta.tmp
mv Bitacora_raws.pep.fasta.tmp Bitacora_raws.pep.fasta


# Identify and curate chimeric sequences

blastp -subject Bitacora_raws.pep.fasta -query /home/projects/ku_00039/people/joeviz/OR_annotation/OR_db/CSP_db.fasta -out Bitacora_raws.pep.fasta.iblast.txt -outfmt "6 std qlen slen" -evalue 1e-5

#perl /home/projects/ku_00039/people/joeviz/programs/bitacora_modftrimlength/Scripts/get_blastp_parsed_newv2.pl Bitacora_raws.pep.fasta.iblast.txt Bitacora_raws.pep.fasta.iblast 1e-5
#perl /home/projects/ku_00039/people/joeviz/programs/bitacora_modftrimlength/Scripts/get_blast_hmmer_combined.pl Bitacora_raws.pep.fasta.iblastblastp_parsed_list.txt Bitacora_raws.pep.fasta.iblast

perl /home/projects/ku_00039/people/joeviz/OR_annotation/get_blastp_parsed_newv2_nmlencut.pl Bitacora_raws.pep.fasta.iblast.txt Bitacora_raws.pep.fasta.iblast 50 1e-5
perl /home/projects/ku_00039/people/joeviz/OR_annotation/get_blast_hmmer_combined.pl Bitacora_raws.pep.fasta.iblastblastp_parsed_list.txt Bitacora_raws.pep.fasta.iblast

perl /home/projects/ku_00039/people/joeviz/OR_annotation/get_fullproteinlist_curatingchimeras_trimsingle.pl Bitacora_raws.pep.fasta Bitacora_raws.pep.fasta.iblast_combinedsearches_list.txt
perl /home/projects/ku_00039/people/joeviz/programs/bitacora_modftrimlength/Scripts/get_annot_genes_gff_v2.pl Bitacora_raws.gff3 $GENOME Bitacora_raws.pep.fasta.iblast_combinedsearches_full_fixed_list.txt Bitacora

sed s/separated/split/g Bitacora_annot_genes_trimmed.gff3 > Bitacora.gff3


perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $GENOME Bitacora.gff3 Bitacora
sed 's/X*$//' Bitacora.pep.fasta > Bitacora.pep.fasta.tmp
mv Bitacora.pep.fasta.tmp Bitacora.pep.fasta

######## Run Interpro in the protein set

f="Bitacora.pep.fasta"

interproscan.sh -i $f -t p -goterms -iprlookup -cpu 4     ## COMMENTED BECAUSE IT IS ALREADY RUN!!!!!


######## Run blast with ORs to obtain names

blastp -query Bitacora.pep.fasta -db /home/projects/ku_00039/people/joeviz/OR_annotation/OR_db/CSP_db.fasta -outfmt "6 std qlen slen" -out Bitacora.pep.fasta.CSPblast.txt -num_threads 4


######## Run script to rename the gff3 and generate the protein file and summary table

perl /home/projects/ku_00039/people/joeviz/OR_annotation/run_OR_classification_bitacora_CSP.pl Bitacora.gff3 $OUTDIR $GENOME


    ' >> submit_script_$FILENAME.sh


echo "qsub submit_script_$FILENAME.sh"


done




