#!/bin/bash


# Number of threads per run
NTHREADS=4

# Create job array to submit

#FILES=/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Hyphy_gard/partition_files_v2/*.cds.nonstop.aln
FILES=/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Cleaned_alignments/run_all_withgard/hmmcleaner_50_aln/*cds.nonstop_hmmclean_gapfilt.aln


#TREEFOLDER=/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Hyphy_gard/partition_files_v2/
TREEFOLDER=/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Cleaned_alignments/run_all_withgard/hmmcleaner_50_aln/
# !!! Change line 32 with the treefile extension

#OUTDIR=absrel
OUTDIR=absrel_allog

mkdir -p $OUTDIR

COUNTER=1

for f in $FILES
do
    FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+)\.cds/){print "$1"}' $f)"
#    FILECDS="$(perl -e 'if ($ARGV[0] =~ /(.*\/\S+)\.pep.fasta/){print "$1".".cds.fasta"}' $f)"
#    mkdir -p $FILENAME

    #TREEFILE="$TREEFOLDER""/""$FILENAME"".treefile"
    #TREEFILE="$TREEFOLDER""/""$FILENAME"".hmmclean_gapfilt.treefile"
    TREEFILE="$TREEFOLDER""/""$FILENAME"".hmmclean_gapfilt.treefile"
    if [ ! -s "$TREEFILE" ]; then
        echo "$TREEFILE not found or empty in tree folder for $f, skipping..."
        continue
    fi



   echo '#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn='"$NTHREADS"':thinnode
#PBS -l mem=16gb
#PBS -l walltime=99:00:00' > submit_hyphy_"$FILENAME".sh


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
module load fasttree/2.1.11
module load perl/5.24.0
module load openmpi/gcc
module load hyphy/2.5.38

# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

    ' >> submit_hyphy_"$FILENAME".sh


    echo "hyphy CPU=$NTHREADS absrel --alignment $f --tree $TREEFILE --output $OUTDIR/$FILENAME""_ABSREL.json > $OUTDIR/$FILENAME""_absrel.output " >> submit_hyphy_"$FILENAME".sh


    i=$((i+1))

    echo "qsub submit_hyphy_$FILENAME.sh" >> submit_all_hyphy.sh

    COUNTER=$((COUNTER+1))

    if [[ "$COUNTER" -gt 100 ]]; then
        COUNTER=1
        echo "sleep 60s" >> submit_all_hyphy.sh
    fi   

done


