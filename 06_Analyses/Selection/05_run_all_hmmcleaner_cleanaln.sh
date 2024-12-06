#!/bin/bash

#### Script to submit hmmcleaner and clean alignment gaps before running absrel

### Add the sleep step between submitted jobs not to saturate the queue
# First run hmmcleaner
# Second use the perl script to keep only codons shared with more than 50% of the species


## Maximum array length is 2000, so grouping multiple commands per run

# Number of threads per alignment
NTHREADS=1

# Create job array to submit

FILES=/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Hyphy_gard/partition_files_v2/*.cds.nonstop.aln

#TREEFOLDER=/home/projects/ku_00039/people/joeviz/Suz/ortholog_alignments/all/codon_aln_trees/
#TREEFOLDER=/home/projects/ku_00039/people/joeviz/Suz/ortholog_alignments/all_2batch/codon_aln_trees/
#TREEFOLDER=/home/projects/ku_00039/people/joeviz/Suz/ortholog_alignments/all_2batch/Hyphy/gard_partitions/codon_aln_trees
#TREEFOLDER=/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/codon_dna_trees/codon_dna_aln_trees/
TREEFOLDER=/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Hyphy_gard/partition_files_v2/
#TREEFOLDER=/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Cleaned_alignments/testfiles/

OUTDIR=hmmcleaner_aln
OUTDIRFILT=hmmcleaner_50_aln 

mkdir -p $OUTDIR
mkdir -p $OUTDIRFILT

COUNTER=1


for f in $FILES
do
    FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+)\.cds/){print "$1"}' $f)"
#    FILECDS="$(perl -e 'if ($ARGV[0] =~ /(.*\/\S+)\.pep.fasta/){print "$1".".cds.fasta"}' $f)"
#    mkdir -p $FILENAME
    NAMENOEXT="$(perl -e 'if ($ARGV[0] =~ /(.*\/\S+)\.aln/){print "$1"}' $f)"

    TREEFILE="$TREEFOLDER""/""$FILENAME"".treefile"
    if [ ! -s "$TREEFILE" ]; then
        echo "$TREEFILE not found or empty in tree folder for $f, skipping..."
        continue
    fi

#    PREVABSRELOUTPUT="$OUTDIR/$FILENAME""_absrel.output"
#    PREVABSRELJSON="$OUTDIR/$FILENAME""_ABSREL.json"
 
#    if [ -s "$PREVABSRELJSON" ]; then
#        continue         
#    fi    


   echo '#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn='"$NTHREADS"':thinnode
#PBS -l mem=16gb
#PBS -l walltime=24:00:00' > submit_cleaning_"$FILENAME".sh


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
module load perl/5.24.0
module load hmmer/3.3.2
module load newick-utils/1.6

cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)
#/home/projects/ku_00039/people/joeviz/perl5/bin/HmmCleaner.pl

# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

    ' >> submit_cleaning_"$FILENAME".sh


#   echo "hyphy CPU=$NTHREADS absrel --alignment $f --tree $TREEFILE --output $OUTDIR/$FILENAME""_ABSREL.json > $OUTDIR/$FILENAME""_absrel.output " >> submit_cleaning_"$FILENAME".sh
    echo "perl /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Cleaned_alignments/translating_seqs_v3_onlyframe1.pl $f $OUTDIR/$FILENAME.pep.aln" >> submit_cleaning_"$FILENAME".sh
    echo "/home/projects/ku_00039/people/joeviz/perl5/bin/HmmCleaner.pl $OUTDIR/$FILENAME.pep.aln --large > $OUTDIR/$FILENAME.stout" >> submit_cleaning_"$FILENAME".sh
    echo "/home/projects/ku_00039/people/joeviz/perl5/bin/transferCleaner.pl $f -log=$OUTDIR/$FILENAME.pep_hmm.log --large -delchar - " >> submit_cleaning_"$FILENAME".sh
#    echo "mv $f"".stdout $OUTDIR" >> submit_cleaning_"$FILENAME".sh
    echo "mv $NAMENOEXT""_cleaned.ali $OUTDIR/$FILENAME.cds.nonstop_hmmcleanraw.ali" >> submit_cleaning_"$FILENAME".sh
    echo "sed 's/\*/-/g' $OUTDIR/$FILENAME.cds.nonstop_hmmcleanraw.ali > $OUTDIR/$FILENAME.cds.nonstop_hmmclean.ali" >> submit_cleaning_"$FILENAME".sh
    echo "perl /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Cleaned_alignments/clean_aln_fromhmmcleaner.pl $OUTDIR/$FILENAME"".cds.nonstop_hmmclean.ali $TREEFOLDER $OUTDIRFILT 0.50" >> submit_cleaning_"$FILENAME".sh
#    echo "mv $OUTDIRFILT""/*_rooted.treefile $OUTDIRROOT" >> submit_cleaning_"$FILENAME".sh

    echo "qsub submit_cleaning_$FILENAME.sh" >> run_submit_all_cleaning.sh

    COUNTER=$((COUNTER+1))

    if [[ "$COUNTER" -gt 100 ]]; then
        COUNTER=1
        echo "sleep 60s" >> run_submit_all_cleaning.sh
    fi   

done




