#!/bin/sh
### Account information
#PBS -W group_list=ku_00039 -A ku_00039
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=1
### Memory
#PBS -l mem=120gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here 12:00:00, 12 hours)
#PBS -l walltime=72:00:00

  
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
module load bwa/0.7.17
module load samtools/1.9
module load bedtools/2.26.0
module load anaconda3/2020.07


## Code
# Run in the same folder after finishing submit_allhic_pipeline.sh

# Load already mapped bam file
OUTCLEANBAM="sample.clean.bam"

# AllHiC dir from: https://github.com/tangerzhang/ALLHiC
ALLHICDIR="/home/projects/ku_00039/people/joeviz/hic_scaffolding/allhic/ALLHiC/"
export PATH=$PATH:"$ALLHICDIR"/bin:"$ALLHICDIR"/scripts



#Heatmap Plot for assembly assessment
#(a) get group length; Note: script can be found here (https://github.com/tangerzhang/my_script/blob/master/getFaLen.pl)

perl "$ALLHICDIR"/scripts/getFaLen.pl -i groups.asm.fasta -o len.txt

grep 'sample.clean.counts_GATC' len.txt > chrn.list
sort -nr -k 2 len.txt > len_sort.txt
awk '{if ($2>800000) print}' len_sort.txt > len_sort_filt800k.txt

#(b) plotting; Note: only keep chromosomal level assembly for plotting.

#"$ALLHICDIR"/bin/ALLHiC_plot $OUTCLEANBAM groups.agp chrn.list 500k pdf
"$ALLHICDIR"/bin/ALLHiC_plot $OUTCLEANBAM groups.agp len_sort_filt800k.txt 500k pdf



