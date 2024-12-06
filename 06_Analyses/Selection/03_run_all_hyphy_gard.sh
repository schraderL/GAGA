#!/bin/bash

# Script to generate a job array to run HyPhy GARD

## Maximum array length is 2000, so grouping multiple commands per run

# Number of threads per alignment
NTHREADS=2

# Create job array to submit

FILES=/home/projects/ku_00039/people/joeviz/orthology_alignments/all_orthogroups/codon_alignments/codon_alignments_qualfiltered/*cds.nonstop.aln


#OUTDIR=protein_aln_trees
#OUTDIR=protein_aln_trees_trimalaut
#OUTDIR=protein_codon_aln_trees
OUTDIR=gard

mkdir -p $OUTDIR


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
#PBS -l walltime=120:00:00' > submit_hyphy_"$FILENAME".sh


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


	# hyphy GARD -- mode Faster could be run if it takes to much with Normal
#    echo "hyphy gard --alignment $f --type codon --mode Faster --output $OUTDIR/$FILENAME""_gard.json --output-lf $OUTDIR/$FILENAME"".best-gard > $OUTDIR/$FILENAME""_gard.output " >> submit_hyphy_"$FILENAME".sh
    echo "hyphy gard --alignment $f --mode Faster --output $OUTDIR/$FILENAME""_gard.json --output-lf $OUTDIR/$FILENAME"".best-gard > $OUTDIR/$FILENAME""_gard.output " >> submit_hyphy_"$FILENAME".sh
#    echo "hyphy gard --alignment $f --type codon --mode Normal --output $OUTDIR/$FILENAME""_gard.json --output-lf $OUTDIR/$FILENAME"".best-gard > $OUTDIR/$FILENAME""_gard.output " >> submit_hyphy_"$FILENAME".sh
#    echo "hyphy gard --alignment $f --mode Normal --output $OUTDIR/$FILENAME""_gard.json --output-lf $OUTDIR/$FILENAME"".best-gard > $OUTDIR/$FILENAME""_gard.output " >> submit_hyphy_"$FILENAME".sh

    echo "hyphy /home/projects/ku_00039/people/joeviz/programs/hyphy-analyses/extract-partitions/extract-partitions.bf --msa $OUTDIR/$FILENAME"".best-gard --output $OUTDIR/$FILENAME""_parts --extension cds.nonstop.fullpart.aln ENV=\"DATA_FILE_PRINT_FORMAT=9\" " >> submit_hyphy_"$FILENAME".sh
#hyphy /home/projects/ku_00039/people/joeviz/programs/hyphy-analyses/extract-partitions/extract-partitions.bf --msa $gardnex --output $outdir/$id\_parts  --extension cds.nonstop.fullpart.aln ENV=\"DATA_FILE_PRINT_FORMAT=9

    i=$((i+1))

    echo "qsub submit_hyphy_$FILENAME.sh" >> run_submit_all_hyphy.sh


done




