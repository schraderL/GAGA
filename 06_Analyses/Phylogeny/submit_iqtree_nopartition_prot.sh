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
#PBS -l mem=1500gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here 12:00:00, 12 hours)
#PBS -l walltime=366:00:00
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
module load anaconda3/4.4.0
module load perl
module load trimal/1.4.1
module load muscle/3.8.425
module load iqtree/2.1.3
#module load iqtree/1.6.12

 
# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

#python /home/projects/ku_00039/people/joeviz/programs/BUSCO_Phylogenomics/BUSCO_phylogenomics.py -d /home/projects/ku_00039/people/joeviz/busco/GAGA_genomes/run_all_final_assemblies_buscov512_euk -l eukaryota_odb10 -t 40 --supermatrix -psc 80 -o busco_phylogeny_supermatrix_psc80
#iqtree -s SUPERMATRIX.aln -bb 1000 -alrt 1000 -nt AUTO -ntmax " + str(threads) + "

ALN=/home/projects/ku_00039/people/joeviz/orthology_trees/single_copy_orthologs/species_trees/concatenated_matrix/concatenated_matrix_singlecopy_90perc_strictexact_protein.fasta

iqtree2 -s $ALN --prefix output_singlecopy_90perc_strictexact_protein_nopartition_b1000_modelfinder -B 1000 -alrt 1000 -m MFP -msub nuclear -T 40 --safe 

