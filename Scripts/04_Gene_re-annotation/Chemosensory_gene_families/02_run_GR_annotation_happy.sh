#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn=40
#PBS -l mem=140gb
#PBS -l walltime=12:00:00

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

GENOME="/home/projects/ku_00039/people/joeviz/GAGA_annotations/Repeat_masked_assemblies/GAGA_annotations/GAGA-0001/GAGA-0001_SLR-superscaffolder_final_dupsrm_filt.softMasked.fasta"
OUTDIR="GAGA-0001" # Directory must not exist before, or it will fail (unless it is used --overwrite)

PROTSEQS="GR_ant_genewise_db.fasta"
ORANNOTATION="/home/projects/ku_00039/people/joeviz/OR_annotation/run_all_GAGA_v4/GAGA-0001/GAGA-0001_ABCENTH_clean_OR_renamed_all.pep.fasta"

HAPpy --threads 40 --genome $GENOME --output_dir $OUTDIR --protein_seqs $PROTSEQS

cd $OUTDIR

######## Generate a GFF3 and protein file

module load genometools/1.5.10

cp genewise/all_genewise_predictions.gff .
#perl /home/projects/ku_00039/people/joeviz/OR_annotation/get_genewise_gtf_corrected.pl all_genewise_predictions.gff $GENOME
perl /home/projects/ku_00039/people/joeviz/OR_annotation/get_genewise_gtf_corrected_extrafilterfirstexon.pl all_genewise_predictions.gff $GENOME

gt gtf_to_gff3 -force -o all_genewise_predictions_corrected_unclean.gff3 all_genewise_predictions_corrected.gtf
gt gff3 -force -sort -tidy -retainids -o all_genewise_predictions_corrected.gff3 all_genewise_predictions_corrected_unclean.gff3

perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $GENOME all_genewise_predictions_corrected.gff3 all_genewise_predictions_corrected
sed s/X*$// all_genewise_predictions_corrected.pep.fasta > all_genewise_predictions_corrected.pep.fasta.tmp
mv all_genewise_predictions_corrected.pep.fasta.tmp all_genewise_predictions_corrected.pep.fasta

## convert cds to exon, and get orf to make genes start with M

#sed s/CDS/exon/g all_genewise_predictions_corrected.gff3 > all_genewise_predictions_corrected_cdstoexon.gff3
perl /home/projects/ku_00039/people/joeviz/OR_annotation/get_genewise_cds_to_exon.pl all_genewise_predictions_corrected all_genewise_predictions_corrected_cdstoexon.gff3
gt cds -force -matchdescstart -minorflen 10 -startcodon yes -seqfile $GENOME -o all_genewise_predictions_corrected_cdstoexon_tocds.gff3 all_genewise_predictions_corrected_cdstoexon.gff3
grep -v 'exon' all_genewise_predictions_corrected_cdstoexon_tocds.gff3 > Genewise.gff3


perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $GENOME Genewise.gff3 Genewise
sed s/X*$// Genewise.pep.fasta > Genewise.pep.fasta.tmp
mv Genewise.pep.fasta.tmp Genewise.pep.fasta


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

f="Genewise.pep.fasta"

interproscan.sh -i $f -t p -goterms -iprlookup -cpu 40    


######## Run blast with ORs to obtain names
module load ncbi-blast/2.11.0+


blastp -query Genewise.pep.fasta -db /home/projects/ku_00039/people/joeviz/OR_annotation/Chemo_db/ORco_db.fasta -outfmt "6 std qlen slen" -out Genewise.pep.fasta.ORcoblast.txt -num_threads 4
#blastp -query Genewise.pep.fasta -db /home/projects/ku_00039/people/joeviz/OR_annotation/Chemo_db/OR_db.fasta -outfmt "6 std qlen slen" -out Genewise.pep.fasta.ORblast.txt -num_threads 4 -max_target_seqs 5
blastp -query Genewise.pep.fasta -db /home/projects/ku_00039/people/joeviz/OR_annotation/Chemo_db/GR2_db.fasta -outfmt "6 std qlen slen" -out Genewise.pep.fasta.GRblast.txt -num_threads 4 -max_target_seqs 5

cat /home/projects/ku_00039/people/joeviz/OR_annotation/Chemo_db/OR_db.fasta $ORANNOTATION > OR_masAbcenth_db.fasta
makeblastdb -in OR_masAbcenth_db.fasta -dbtype prot
blastp -query Genewise.pep.fasta -db OR_masAbcenth_db.fasta -outfmt "6 std qlen slen" -out Genewise.pep.fasta.ORblast.txt -num_threads 4 -max_target_seqs 5

blastp -query Genewise.pep.fasta -db /home/projects/ku_00039/people/joeviz/OR_annotation/Chemo_db/GR_dmel_db.fasta -outfmt "6 std qlen slen" -out Genewise.pep.fasta.GRdmelblast.txt -num_threads 4 -max_target_seqs 5


######## Run script to rename the gff3 and generate the protein file and summary table

perl /home/projects/ku_00039/people/joeviz/OR_annotation/scripts/run_GR_classification_happy_GR.pl Genewise.gff3 $OUTDIR $GENOME


    
