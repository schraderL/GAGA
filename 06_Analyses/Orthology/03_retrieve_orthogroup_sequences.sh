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
#PBS -l mem=180gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here 12:00:00, 12 hours)
#PBS -l walltime=120:00:00
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
module load ngs tools
module load anaconda3/4.4.0

# Code
# Input N0_onehogcol.tsv is the N0.tsv file removing the second, third and fourth columns, therefore containing a column with the HOG id, and the genes for all species included in hte orthology assessment

perl get_orthogroup_sequences.pl N0_onehogcol.tsv all_orthogroups

# Other examples
# perl get_orthogroup_sequences.pl species_all_80percsp_orthogroups_singlecopy_genes.tsv single_copy

