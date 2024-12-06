#!/bin/bash

### Modified script to sleep for a minute each 100 jobs, so the queue has time to process the jobs. 
## Maximum array length is 2000, so grouping multiple commands per run ## NO
## Single job per HOG
# Time set to 150 hours per job

module purge
module load ngs tools
module load anaconda3/4.4.0

# Number of threads per alignment
NTHREADS=1

# Create job array to submit

#TRAITFILE=/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/Worker_sterility/Species_selection_Worker_ovary_clades_all.txt
#TRAITFILE=/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/Worker_polymorphism/Species_selection_Polymorphism_clades_all.txt
TRAITFILE=$1

# Input folder with the codon alignments
FILES=/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Cleaned_alignments/run_all_withgard/hmmcleaner_50_aln/*cds.nonstop_hmmclean_gapfilt.aln

# Input folder with the trees for each alignment, edit line 41 with the extension for each tree
TREEFOLDER=/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Cleaned_alignments/run_all_withgard/hmmcleaner_50_aln/

# Output directory
OUTDIR=relax

mkdir -p $OUTDIR

COUNTER=1

for f in $FILES
do
    FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+)\.cds/){print "$1"}' $f)"
#    FILECDS="$(perl -e 'if ($ARGV[0] =~ /(.*\/\S+)\.pep.fasta/){print "$1".".cds.fasta"}' $f)"
#    mkdir -p $FILENAME

    #TREEFILE="$TREEFOLDER""/""$FILENAME"".treefile"
    TREEFILE="$TREEFOLDER""/""$FILENAME"".hmmclean_gapfilt.treefile"    
    if [ ! -s "$TREEFILE" ]; then
        echo "$TREEFILE not found or empty in tree folder for $f, skipping..."
        continue
    fi

    # Create tree for relax

    RELAXTREEFILE="$OUTDIR""/""$FILENAME""_relax.treefile"
    python3 /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/label_relax_fromcladetable.py $TREEFILE $TRAITFILE $RELAXTREEFILE
#    python3 /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/label_relax_fromcladetable_winternalnodes_usingclades.py $TREEFILE $TRAITFILE $RELAXTREEFILE


   echo '#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn='"$NTHREADS"':thinnode
#PBS -l mem=8gb
#PBS -l walltime=96:00:00' > submit_hyphy_"$FILENAME".sh


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
#module load hyphy/2.5.29
#module load hyphy/2.5.38
module load hyphy/2.5.50

# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

    ' >> submit_hyphy_"$FILENAME".sh


    echo "hyphy CPU=$NTHREADS relax --starting-points 5 --alignment $f --tree $RELAXTREEFILE --output $OUTDIR/$FILENAME""_RELAX.json --test test --reference reference --models Minimal > $OUTDIR/$FILENAME""_relax.output " >> submit_hyphy_"$FILENAME".sh
 

    i=$((i+1))

    echo "qsub submit_hyphy_$FILENAME.sh" >> submit_all_hyphy.sh

    COUNTER=$((COUNTER+1))

    if [[ "$COUNTER" -gt 100 ]]; then
        COUNTER=1
        echo "sleep 60s" >> submit_all_hyphy.sh
    fi    

done




