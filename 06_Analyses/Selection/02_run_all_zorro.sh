#!/bin/bash


## Create a job array to run ZORRO

## Maximum array length is 2000, so grouping multiple command per run


# Number of threads per alignment
NTHREADS=1

# Create array job to submit

#FILES=/home/projects/ku_00039/people/joeviz/orthology_alignments/test_alns/codon_alignments/*fas
#FILES=/home/projects/ku_00039/people/joeviz/orthology_alignments/single_copy_orthogroups_speciestree/codon_alignments/codon_alignments/*cds.aln
#FILES=/home/projects/ku_00039/people/joeviz/orthology_alignments/single_copy_orthogroups_speciestree/protein_alignments/protein_codon_alignments_pal2nal/*aln
FILES=/home/projects/ku_00039/people/joeviz/orthology_alignments/all_orthogroups/codon_alignments/codon_alignments/*cds.aln


OUTDIR=zorro_codon_alignments
#OUTDIR=zorro_protein_alignments
#OUTDIR=zorro_protein_codon_alignments

mkdir -p $OUTDIR

i=1
e=0

for f in $FILES
do
    FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+)\.aln/){print "$1"}' $f)"
 #   FILECDS="$(perl -e 'if ($ARGV[0] =~ /(.*\/\S+)\.pep.fasta/){print "$1".".cds.fasta"}' $f)"


	echo "#!/bin/sh
	" >> $OUTDIR/run_zorro_$i.sh


#	if [ $e -gt 3 ]
	if [ $e -gt 14 ]
	then   

#/home/projects/ku_00039/people/joeviz/programs/zorro_linux_x86_64 -sample in > out

		# ZORRO
	   	echo "/home/projects/ku_00039/people/joeviz/programs/zorro_linux_x86_64 -sample $f > $OUTDIR""/$FILENAME""_zorro.txt
	   	" >> run_zorro_$i.sh

	    i=$((i+1))
	    e=0

	else

		# ZORRO
	   	echo "/home/projects/ku_00039/people/joeviz/programs/zorro_linux_x86_64 -sample $f > $OUTDIR""/$FILENAME""_zorro.txt
	   	" >> run_zorro_$i.sh

	    e=$((e+1))

	fi


done


    echo '#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn='"$NTHREADS"'
#PBS -l mem=20gb
#PBS -l walltime=96:00:00' > submit_run_zorro_all.sh

    echo "#PBS -t 1-$i" >> submit_run_zorro_all.sh

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
module load anaconda3/4.4.0
module load fasttree/2.1.11
module load perl/5.24.0


# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

bash run_zorro_${PBS_ARRAYID}.sh

    ' >> submit_run_zorro_all.sh





