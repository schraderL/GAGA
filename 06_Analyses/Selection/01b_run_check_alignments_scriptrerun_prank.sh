#!/bin/bash

## Check if all alignments exist in ALNFOLDER with the extension, and are non zero size (-s instead of -f)

NTHREADS=1

#FILES=/home/projects/ku_00039/people/joeviz/orthology_alignments/test_orthogroups_singlecopyv2/*pep.fasta
#FILES=/home/projects/ku_00039/people/joeviz/orthology_alignments/single_copy_orthogroups_speciestree/orthogroups_seqs/*pep.fasta
FILES=/home/projects/ku_00039/people/joeviz/orthology_alignments/all_orthogroups/orthogroups_seqs/*pep.fasta

#ALNFOLDER="protein_alignments/protein_codon_alignments"
#ALNFOLDER="protein_alignments/protein_alignments"
ALNFOLDER="codon_alignments"
FEXTENSION=".cds.aln"
#FEXTENSION=".pep.aln"
#FEXTENSION=".cds.pal2nal.aln"

#OUTDIR=single_copy_speciestree

#mkdir -p $OUTDIR


i=1
e=0

for f in $FILES
do
    FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+)\.pep.fasta/){print "$1"}' $f)"
    FILECDS="$(perl -e 'if ($ARGV[0] =~ /(.*\/\S+)\.pep.fasta/){print "$1".".cds.fasta"}' $f)"
#    mkdir -p $FILENAME


#	if [  -f "$ALNFOLDER/$FILENAME$FEXTENSION" ]; 
#	then
#	    echo "$FILENAME exists."
#	else
#		echo "$FILENAME\n"
#	fi

	if [ ! -s "$ALNFOLDER/$FILENAME$FEXTENSION" ]; then

		echo '#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn='"$NTHREADS"'
#PBS -l mem=120gb
#PBS -l walltime=192:00:00

# Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

### Here follows the user commands:
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes

# Load all required modules for the job
module load tools
module load ngs
#module load anaconda3/4.4.0
#module load mafft/7.453
#module load muscle/3.8.425
module load gcc/11.1.0
module load prank/170427

# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

prank -d='"$FILECDS"' -o=codon_alignments/'"$FILENAME"'.cds.aln -codon -F
mv codon_alignments/'"$FILENAME"'.cds.aln.best.fas codon_alignments/'"$FILENAME"'.cds.aln

	' >> submit_run_aln_$FILENAME.sh
	echo "qsub submit_run_aln_$FILENAME.sh"

	fi

done




