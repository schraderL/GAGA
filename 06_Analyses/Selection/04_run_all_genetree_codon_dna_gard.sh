#!/bin/bash


# Number of threads per alignment
NTHREADS=2

# Create jobs to submit


#### !!! Use codon alignments without stop codons or iqtree will stop

FILES=/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Hyphy_gard/partition_files_v2/*cds.nonstop.aln

OUTDIR=partition_files_v2

#mkdir -p $OUTDIR


for f in $FILES
do
    FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+)\.cds/){print "$1"}' $f)"
#    FILECDS="$(perl -e 'if ($ARGV[0] =~ /(.*\/\S+)\.pep.fasta/){print "$1".".cds.fasta"}' $f)"
#    mkdir -p $FILENAME

    echo '#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn='"$NTHREADS"':thinnode
#PBS -l mem=40gb
#PBS -l walltime=12:00:00' > submit_run_iqtree_"$FILENAME".sh

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
module load iqtree/2.1.3

# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here


    ' >> submit_run_iqtree_"$FILENAME".sh

	# iqtree
	 echo "iqtree2 -s $f -B 1000 -T $NTHREADS -m MFP --prefix $OUTDIR/$FILENAME
	" >> submit_run_iqtree_"$FILENAME".sh

	echo "qsub submit_run_iqtree_"$FILENAME".sh" >> run_submit_all.sh




done



