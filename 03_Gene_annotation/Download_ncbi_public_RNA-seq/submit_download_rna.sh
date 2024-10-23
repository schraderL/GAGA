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
#PBS -l nodes=1:ppn=1:thinnode
### Memory
#PBS -l mem=10gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 12 hours)
#PBS -l walltime=500:00:00
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
module load ascp/3.9.6
module load sratoolkit/2.10.7 
module load pigz/2.3.4
 
# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in separate script and execute from here
#<your script>

#fastq-dump --defline-qual "+" --split-e --outdir . --defline-seq '@$ac-$si/$ri' SRR7894261
#fastq-dump --defline-qual "+" --split-e --outdir . --defline-seq '@$ac-$si/$ri' SRR7894262

#perl run_download_sra.pl SRA_data_RNA-seq_tab_test.tsv
perl run_download_sra.pl SRA_data_RNA-seq_tab.tsv
