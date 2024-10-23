#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn=40
#PBS -l mem=160gb
#PBS -l walltime=24:00:00

# Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

### Here follows the user commands:
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes

# Load all required modules for the job
#First load modules
module load ngs
module load tools
module load miniconda3/4.10.3


# Then activate conda happy
conda activate happy


# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

GENOME="/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/GAGA_all_final_assemblies/GAGA-0001_SLR-superscaffolder_final_dupsrm_filt.fasta"
OUTDIR="GAGA-0001" # Directory must not exist before, or it will fail (unless it is used --overwrite)

HMMDIR="path/to/RefHMMs_ant_GRs_0.40dist" # HMM files
ORANNOTATION="/home/projects/ku_00039/people/joeviz/OR_annotation/run_all_GAGA_v4/$OUTDIR/$OUTDIR\_ABCENTH_clean_OR_renamed_all.pep.fasta"  # OR annotations


HAPpy --threads 40 --annotator ABCENTH --hmm_dir $HMMDIR --genome $GENOME --output_dir $OUTDIR

cd $OUTDIR

######## Generate a GFF3 and protein file

module load genometools/1.5.10

perl /home/projects/ku_00039/people/joeviz/OR_annotation/scripts/get_abcenth_gtf_corrected.pl ABCENTH.gtf $GENOME
gt gtf_to_gff3 -o ABCENTH.gff3 ABCENTH_corrected.gtf
gt gff3 -sort -tidy -retainids -o ABCENTH_clean.gff3 ABCENTH.gff3

perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $GENOME ABCENTH_clean.gff3 ABCENTH_clean
sed s/X*$// ABCENTH_clean.pep.fasta > ABCENTH_clean.pep.fasta.tmp
mv ABCENTH_clean.pep.fasta.tmp ABCENTH_clean.pep.fasta


######## Run Interpro in the protein set
module purge
module load ngs tools
module load anaconda3/4.4.0
#module load openjdk/15.0.1
module load jdk/19
module load perl
module load signalp/4.1g
#module load interproscan/5.51-85.0
module load interproscan/5.52-86.0

f="ABCENTH_clean.pep.fasta"

interproscan.sh -i $f -t p -goterms -iprlookup -cpu 40



######## Run blast with ORs to obtain names
module load ncbi-blast/2.11.0+


blastp -query ABCENTH_clean.pep.fasta -db /home/projects/ku_00039/people/joeviz/OR_annotation/Chemo_db/ORco_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.ORcoblast.txt -num_threads 40
blastp -query ABCENTH_clean.pep.fasta -db /home/projects/ku_00039/people/joeviz/OR_annotation/Chemo_db/GR2_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.GRblast.txt -num_threads 40 -max_target_seqs 5

cat /home/projects/ku_00039/people/joeviz/OR_annotation/Chemo_db/OR_db.fasta $ORANNOTATION > OR_masAbcenth_db.fasta
makeblastdb -in OR_masAbcenth_db.fasta -dbtype prot
blastp -query ABCENTH_clean.pep.fasta -db OR_masAbcenth_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.ORblast.txt -num_threads 4 -max_target_seqs 5

blastp -query ABCENTH_clean.pep.fasta -db /home/projects/ku_00039/people/joeviz/OR_annotation/Chemo_db/GR_dmel_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.GRdmelblast.txt -num_threads 4 -max_target_seqs 5



######## Run script to rename the gff3 and generate the protein file and summary table

perl /home/projects/ku_00039/people/joeviz/OR_annotation/scripts/run_GR_classification_abcenth_GR_nonamefilter.pl ABCENTH_clean.gff3 $OUTDIR $GENOME


    
