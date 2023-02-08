#!/bin/sh
### Account information
#PBS -W group_list=ku_00039 -A ku_00039
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=40
### Memory
#PBS -l mem=40gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here 12:00:00, 12 hours)
#PBS -l walltime=96:00:00

  
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

# Input data

FIRSTGENOME="/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Polished_assemblies_6batch/GAGA-0026_nextpolish.fasta" # Genome to scaffold
READ1="/home/projects/ku_00039/people/joeviz/data/Hi-C_data/GAGA-0026/CL100158547_L02_read_1.fq.gz" # Hi-C reads
READ2="/home/projects/ku_00039/people/joeviz/data/Hi-C_data/GAGA-0026/CL100158547_L02_read_2.fq.gz" # Hi-C reads
CHRNUMBER="26" # Number of expected chromosomes


# tmp files names
INBAM="combined.bam" # output bam from mapping Hi-C data into the initial genome assembly
GENOME="GAGA-0026_assembly_hic_corrected.fasta" # Genome assembly after being corrected using the Hi-C mapping

# output bam
OUTCLEANBAM="sample.clean.bam" # output bam from mapping Hi-C data into the Hi-C corrected genome assembly
OUTCLEAN="sample.clean"

#Use precomputed bam file, skip then the mapping step in lines: 62-71
#OUTCLEANBAM="/home/projects/ku_00039/people/joeviz/hic_scaffolding/allhic/GAGA-0346/allhic_salsabam_salsacleanedassembly/combined.bam"
#OUTCLEAN="/home/projects/ku_00039/people/joeviz/hic_scaffolding/allhic/GAGA-0346/allhic_salsabam_salsacleanedassembly/combined"

# Script directory from AllHiC: 
ALLHICDIR="/home/projects/ku_00039/people/joeviz/hic_scaffolding/allhic/ALLHiC/"
export PATH=$PATH:"$ALLHICDIR"/bin:"$ALLHICDIR"/scripts
# Scripts from https://github.com/esrice/hic-pipeline
HICDIR="/home/projects/ku_00039/people/joeviz/hic_scaffolding/salsa2/hic-pipeline/"
export PATH=$PATH:$HICDIR


#First, index the contigs
bwa index $FIRSTGENOME
samtools faidx $FIRSTGENOME

#Now, align the R1 and R2 reads to your contigs separately, filtering the output through the filter_chimeras.py script to remove experimental artifacts from the alignments:
bwa mem -t 40 $FIRSTGENOME $READ1 | samtools view -@ 40 -bh - | "$HICDIR"/filter_chimeras.py - > r1.bam
bwa mem -t 40 $FIRSTGENOME $READ2 | samtools view -@ 40 -bh - | "$HICDIR"/filter_chimeras.py - > r2.bam

#Combine the r1 and r2 files using combine_ends.py, and then fix mates, sort, and remove PCR duplicates:
"$HICDIR"/combine_ends.py r1.bam r2.bam | samtools fixmate -@ 40 -m - - | samtools sort -@ 40 - | samtools markdup -@ 40 -r - $INBAM



# Index bam file

samtools index $INBAM

# Run corrector, to split putative mis-joined scaffolds using Hi-C mapping data

"$ALLHICDIR"/bin/ALLHiC_corrector -m $INBAM -r $FIRSTGENOME -o $GENOME -t 40


# Index the Hi-C corrected assembly

bwa index $GENOME
samtools faidx $GENOME


#Now, align the R1 and R2 reads to your contigs separately, filtering the output through the filter_chimeras.py script to remove experimental artifacts from the alignments:

bwa mem -t 40 $GENOME $READ1 | samtools view -@ 40 -bh - | "$HICDIR"/filter_chimeras.py - > r1.bam
bwa mem -t 40 $GENOME $READ2 | samtools view -@ 40 -bh - | "$HICDIR"/filter_chimeras.py - > r2.bam


#Combine the r1 and r2 files using combine_ends.py, and then fix mates, sort, and remove PCR duplicates:

"$HICDIR"/combine_ends.py r1.bam r2.bam | samtools fixmate -@ 40 -m - - | samtools sort -@ 40 - | samtools markdup -@ 40 -r - $OUTCLEANBAM

#############################


#Partition contigs into user pre-defined groups; restriction sites (-e) and number of clusters (-k) should be modified accordingly

"$ALLHICDIR"/bin/ALLHiC_partition -b $OUTCLEANBAM -r $GENOME -e GATC -k $CHRNUMBER


#Extract CLM file and counts of restriction sites

"$ALLHICDIR"/bin/allhic extract $OUTCLEANBAM $GENOME --RE GATC


#Optimize the ordering and orientation for each clustered group

for i in $(seq 1 $CHRNUMBER); 
	do 
		
		"$ALLHICDIR"/bin/allhic optimize "$OUTCLEAN".counts_GATC."$CHRNUMBER"g"$i".txt "$OUTCLEAN".clm

	done


#Get chromosomal level assembly

"$ALLHICDIR"/bin/ALLHiC_build $GENOME


#Heatmap Plot for assembly assessment
#(a) get group length; Note: script can be found here (https://github.com/tangerzhang/my_script/blob/master/getFaLen.pl)

perl "$ALLHICDIR"/scripts/getFaLen.pl -i groups.asm.fasta -o len.txt

grep 'sample.clean.counts_GATC' len.txt > chrn.list
sort -nr -k 2 len.txt > len_sort.txt
awk '{if ($2>800000) print}' len_sort.txt > len_sort_filt800k.txt

#(b) plotting; Note: only keep chromosomal level assembly for plotting.

#"$ALLHICDIR"/bin/ALLHiC_plot $OUTCLEANBAM groups.agp chrn.list 500k pdf
"$ALLHICDIR"/bin/ALLHiC_plot $OUTCLEANBAM groups.agp len_sort_filt800k.txt 500k pdf



