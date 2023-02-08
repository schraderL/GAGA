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
#PBS -l nodes=1:ppn=1
### Memory
#PBS -l mem=80gb
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
module load perl
module load miniconda3/4.8.5

conda activate ntJoin

module load samtools/1.10
module load bedtools/2.28.0
 
# Code

### First, the input genomes must have just one sequence per line; use script "get_fasta_1seqperline.pl genome" to create the required files
### Copy or create softlink to the current directory!!!
genome2="GAGA-0392_supernova.fasta" 
genome1="GAGA-0392_masurca.fasta"

out=GAGA-392_ntjoin_scaffolded



ntJoin assemble target=$genome1 target_weight=1 references=$genome2 reference_weights=2 no_cut=True k=32 w=500 G=5000 prefix=$out



