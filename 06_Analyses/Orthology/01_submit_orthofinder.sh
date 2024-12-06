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
module load perl/5.24.0
module load ncbi-blast/2.11.0+
module load mcl/14-137
module load fastme/2.1.5
module load dlcpar/1.0
module load orthofinder/2.5.4
module load mafft/7.453
module load fasttree/2.1.11
module load iqtree/2.1.3

alias iqtree='/services/tools/iqtree/2.1.3/bin/iqtree2'

# Code

# Folder DB contains the final annotated proteins (one representative isoform per gene) with name GAGA-ID.fasta

orthofinder -f DB -S diamond_ultra_sens -t 40 -a 1 -M msa -y -o output_msa -s GAGA_dated_phylogeny_newick.tre


