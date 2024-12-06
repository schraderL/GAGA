#!/bin/bash

## Maximum array length is 2000, so grouping multiple commands per run

# Number of threads per alignment
NTHREADS=1

# Create array job to submit

#FILES=/home/projects/ku_00039/people/joeviz/orthology_alignments/single_copy_orthogroups_speciestree/orthogroups_seqs/*cds.fasta
FILES=/home/projects/ku_00039/people/joeviz/orthology_alignments/all_orthogroups/orthogroups_seqs/*cds.fasta


#OUTDIR=single_copy_speciestree

#mkdir -p $OUTDIR
mkdir -p codon_alignments
#mkdir -p protein_alignments
#mkdir -p protein_codon_alignments


i=1
e=0

for f in $FILES
do
    FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+)\.cds.fasta/){print "$1"}' $f)"
#    FILECDS="$(perl -e 'if ($ARGV[0] =~ /(.*\/\S+)\.pep.fasta/){print "$1".".cds.fasta"}' $f)"
#    mkdir -p $FILENAME


	echo "#!/bin/sh
	" >> run_aln_$i.sh


	if [ $e -gt 14 ]
	then   

		# PRANK
	    echo "prank -d=$f -o=codon_alignments/$FILENAME.cds.aln -codon -F
		mv codon_alignments/$FILENAME.cds.aln.best.fas codon_alignments/$FILENAME.cds.aln
		" >> run_aln_$i.sh

		# Prot aln
		# MAFFT
#	   	echo "mafft --thread $NTHREADS --auto $f > protein_alignments/$FILENAME.pep.aln
#	   	" >> run_aln_$i.sh

		# Tcoffee
#		echo "t_coffee $f -mode mcoffee -multi_core jobs
#		perl /programs/clustal2fasta.pl < $f\_proteins.aln > protein_alignments/$FILENAME\.pep.aln
#		" >> run_aln_$i.sh
		#Muscle
#		echo "muscle -in $f -out protein_alignments/$FILENAME\.pep.aln
#		" >> run_aln_$i.sh

		# CDS aln from prot aln
		#pal2nal
#		echo "perl /home/projects/ku_00039/people/joeviz/programs/pal2nal.pl protein_alignments/$FILENAME.pep.aln $FILECDS -output fasta > protein_codon_alignments/$FILENAME.cds.pal2nal.aln
#		" >> run_aln_$i.sh	

#		echo "
#	    " >> run_aln_$i.sh

	    i=$((i+1))
	    e=0

	else

	    # PRANK
	    echo "prank -d=$f -o=codon_alignments/$FILENAME.cds.aln -codon -F
		mv codon_alignments/$FILENAME.cds.aln.best.fas codon_alignments/$FILENAME.cds.aln
		" >> run_aln_$i.sh


	    e=$((e+1))

	fi


done


    echo '#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn='"$NTHREADS"'
#PBS -l mem=60gb
#PBS -l walltime=96:00:00' > submit_run_aln_all.sh

    echo "#PBS -t 1-$i" >> submit_run_aln_all.sh

echo '
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

bash run_aln_${PBS_ARRAYID}.sh

    ' >> submit_run_aln_all.sh





