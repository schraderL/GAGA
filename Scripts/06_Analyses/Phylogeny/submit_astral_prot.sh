#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00039 -A ku_00039
### Job name (comment out the next line to get the name of the script used as the job name)
##PBS -N test
### Output files (comment out the next 2 lines to get the job name used instead)
##PBS -e test.err
##PBS -o test.log
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=40
### Memory
#PBS -l mem=170gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here 12:00:00, 12 hours)
#PBS -l walltime=24:00:00
### Forward X11 connection (comment out if not needed)
##PBS -X
  
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
module load java/1.8.0
module load astral-mp/5.15.4
module load newick-utils/1.6


# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

cat ../protein_aln_trees_trimal/*treefile > singlecopy_protein_alltrees.nwk

nw_ed singlecopy_protein_alltrees.nwk 'i & b <= 10' o > singlecopy_protein_alltrees_bbcollapsed.nwk

astral-mp -i singlecopy_protein_alltrees.nwk -o singlecopy_protein_alltrees.astral.tree -T 40 > singlecopy_protein_alltrees.astral.out 2> singlecopy_protein_alltrees.astral.err
astral-mp -i singlecopy_protein_alltrees_bbcollapsed.nwk -o singlecopy_protein_alltrees_bbcollapsed.astral.tree -T 40 > singlecopy_protein_alltrees_bbcollapsed.astral.out 2> singlecopy_protein_alltrees_bbcollapsed.astral.err
#astral-mp -i All_trees_trimmed_protein_alignments_bbcollapsed.nwk -o All_trees_trimmed_protein_alignments_bbcollapsed.astral.tree -T 40 > All_trees_trimmed_protein_alignments_bbcollapsed.astral.out 2> All_trees_trimmed_protein_alignments_bbcollapsed.astral.err

#astral-mp -i singlecopy_proteincodon_alltrees_bbcollapsed.nwk -o singlecopy_proteincodon_alltrees_bbcollapsed.tree -T 40 > singlecopy_proteincodon_alltrees_bbcollapsed.astral.out 2> singlecopy_proteincodon_alltrees_bbcollapsed.astral.err
#astral-mp -i singlecopy_90perc_strict_protein_alltrees_bbcollapsed.nwk -o singlecopy_90perc_strict_protein_alltrees_bbcollapsed.astral.tree -T 40 > singlecopy_90perc_strict_protein_alltrees_bbcollapsed.astral.out 2> singlecopy_90perc_strict_protein_alltrees_bbcollapsed.astral.err
#astral-mp -i singlecopy_90perc_strictexact_protein_alltrees_bbcollapsed.nwk -o singlecopy_90perc_strictexact_protein_alltrees_bbcollapsed.astral.tree -T 40 > singlecopy_90perc_strictexact_protein_alltrees_bbcollapsed.astral.out 2> singlecopy_90perc_strictexact_protein_alltrees_bbcollapsed.astral.err






