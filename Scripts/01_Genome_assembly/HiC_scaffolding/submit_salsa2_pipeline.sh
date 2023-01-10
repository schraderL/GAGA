#!/bin/sh
### Account information
#PBS -W group_list=ku_00039 -A ku_00039
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=40
### Memory
#PBS -l mem=180gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here 12:00:00, 12 hours)
#PBS -l walltime=200:00:00
  
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
module load boost/1.74.0


# Code


## Set input data and program paths

# Genome assembly to scaffold using the Hi-C data
GENOME="/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Polished_assemblies_10batch/GAGA-0014_nextpolish_correct.fasta"
# Hi-C reads
READ1="/home/projects/qku_00039/people/joeviz/data/GAGA-0014/Hi-C/CL100158752_L01_read_1.fq.gz"
READ2="/home/projects/ku_00039/people/joeviz/data/GAGA-0014/Hi-C/CL100158752_L01_read_2.fq.gz"

# Directory to salsa and scripts
SALSADIR="/home/projects/ku_00039/people/joeviz/hic_scaffolding/salsa2/SALSA-2.2/"
# Scripts from https://github.com/ArimaGenomics/mapping_pipeline
HICDIR="/home/projects/ku_00039/people/joeviz/hic_scaffolding/salsa2/hic-pipeline/"
export PATH=$PATH:$HICDIR

#First, index the contigs

bwa index $GENOME
samtools faidx $GENOME


#Now, align the R1 and R2 reads to your contigs separately, filtering the output through the filter_chimeras.py script to remove experimental artifacts from the alignments:

bwa mem -t 40 $GENOME $READ1 | samtools view -@ 40 -bh - | "$HICDIR"/filter_chimeras.py - > r1.bam
bwa mem -t 40 $GENOME $READ2 | samtools view -@ 40 -bh - | "$HICDIR"/filter_chimeras.py - > r2.bam


#Combine the r1 and r2 files using combine_ends.py, and then fix mates, sort, and remove PCR duplicates:

"$HICDIR"/combine_ends.py r1.bam r2.bam | samtools fixmate -@ 40 -m - - | samtools sort -@ 40 - | samtools markdup -@ 40 -r - combined.bam


##Run SALSA2
#SALSA2 takes alignments in bed rather than bam format, sorted by read name:

bamToBed -i combined.bam | sort -k 4 > combined.bed

#Now, run SALSA. 

module unload anaconda3/2020.07
module load miniconda3/4.8.5
#conda create -n python2 python=2.7
#conda install -c prometeia networkx  #v1.11
conda activate python2

python2.7 "$SALSADIR"/run_pipeline.py -a $GENOME -l "$GENOME".fai -b combined.bed -e GATC -o salsa -m yes


# Generate hic visualization

module load java/1.8.0

bash /home/projects/ku_00039/people/joeviz/hic_scaffolding/salsa2/SALSA/convert.sh salsa


