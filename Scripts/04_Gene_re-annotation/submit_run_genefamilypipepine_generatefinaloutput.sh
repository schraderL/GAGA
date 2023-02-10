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
#PBS -l nodes=1:ppn=5
### Memory
#PBS -l mem=40gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here 12:00:00, 12 hours)
#PBS -l walltime=24:00:00
### Forward X11 connection (comment out if not needed)
##PBS -X

# Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

### Here follows the user commands:
# Define number of processors
NPROCS='wc -l < $PBS_NODEFILE'
echo This job has allocated $NPROCS nodes

# Load all required modules for the job
module load tools
module load ngs
module load anaconda3/4.4.0
module load openjdk/15.0.1
module load perl

# Code

OUTDIR="/home/projects/ku_00039/people/joeviz/gene_family_pipeline/run_pipeline_Ectatomma/run_pipeline/Ectatomma_final_annotations" # Output directory
mkdir -p $OUTDIR
cd $OUTDIR

###### EDIT INPUT ######
path="/home/projects/ku_00039/people/joeviz/gene_family_pipeline/run_pipeline_Ectatomma/Gene-annotation-pipeline"
table="/home/projects/ku_00039/people/joeviz/gene_family_pipeline/run_pipeline_Ectatomma/Gene_reannotation_full_table.xlsx"
pipeline_output_directory="/home/projects/ku_00039/people/joeviz/gene_family_pipeline/run_pipeline_Ectatomma/run_pipeline/Ectatomma_test"
genome="/home/projects/ku_00039/people/joeviz/GAGA_genomes/Colab_Ectatomma/repeat_annot/Ectatomma_ruidum_Hifigenome.nextpolish_final.softMasked.fasta"
output_directory="$OUTDIR"
############################

python3 "$path"/Scripts/gen_final_output.py $path $genome $table $pipeline_output_directory $output_directory
