#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn=4
#PBS -l mem=60gb
#PBS -l walltime=6:00:00

# Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

### Here follows the user commands:
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes

# Load all required modules for the job
module load ngs
module load tools
module load anaconda3/4.4.0
module load perl
module load hmmer/3.2.1
module load java/1.8.0


# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

GENOME="/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/GAGA_all_final_assemblies/GAGA-0001_SLR-superscaffolder_final_dupsrm_filt.fasta"
OUTDIR="GAGA-0001" # Directory must not exist before, or it will fail (unless it is used --overwrite)

ORANNOTATION="/home/projects/ku_00039/people/joeviz/OR_annotation/run_all_GAGA_v4/GAGA-0001/GAGA-0001_ABCENTH_clean_OR_renamed_all.pep.fasta"

mkdir $OUTDIR
cd $OUTDIR


# Get previous ABCENTH and genewise annotation GFF


cp /home/projects/ku_00039/people/joeviz/OR_annotation/Chemosensory/run_GR/run_abcenth_all/$OUTDIR/"$OUTDIR"_ABCENTH_clean_GR_renamed_all_nofragment.gff3 ABCENTH_clean_GR_renamed_all_nofragment.gff3
cp /home/projects/ku_00039/people/joeviz/OR_annotation/Chemosensory/run_GR/run_abcenth_all/$OUTDIR/"$OUTDIR"_GRs.txt ABCENTH_GRs.txt

cp /home/projects/ku_00039/people/joeviz/OR_annotation/Chemosensory/run_GR/run_genewise_all_ok/$OUTDIR/"$OUTDIR"_Genewise_GR_renamed_all.gff3 Genewise_GR_renamed_all.gff3


######## Get sugar and conserved GR receptors, and combine with genewise (priority for conserved GRs in ABCENTH)

perl /home/projects/ku_00039/people/joeviz/OR_annotation/get_conserved_GRs.pl ABCENTH_GRs.txt ABCENTH_clean_GR_renamed_all_nofragment.gff3 

# It will generate: ABCENTH_GR_all_renamed.gff3 and ABCENTH_GR_Dmelconserved_renamed.gff3

java -jar /home/projects/ku_00039/people/joeviz/programs/GeMoMa-1.7.1/GeMoMa-1.7.1.jar CLI CompareTranscripts p=Genewise_GR_renamed_all.gff3 a=ABCENTH_GR_Dmelconserved_renamed.gff3 outdir=gemoma_outdir_Dmel > gemoma.out 2> gemoma.err

perl /home/projects/ku_00039/people/joeviz/OR_annotation/get_combined_nr_gff.pl Genewise_GR_renamed_all.gff3 ABCENTH_GR_Dmelconserved_renamed.gff3 gemoma_outdir_Dmel/comparison.tabular Genewise_AbcenthDmel_combined.gff3


######## Combine both GFF3 genewise and ABCENTH (priority for gene wise)

java -jar /home/projects/ku_00039/people/joeviz/programs/GeMoMa-1.7.1/GeMoMa-1.7.1.jar CLI CompareTranscripts p=ABCENTH_GR_all_renamed.gff3 a=Genewise_AbcenthDmel_combined.gff3 outdir=gemoma_outdir_abcenth > gemoma.out 2> gemoma.err

perl /home/projects/ku_00039/people/joeviz/OR_annotation/get_combined_nr_gff.pl ABCENTH_GR_all_renamed.gff3 Genewise_AbcenthDmel_combined.gff3 gemoma_outdir_abcenth/comparison.tabular GenewiseAbcenth.gff3



######## Generate the GFF3 and protein file

perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $GENOME GenewiseAbcenth.gff3 GenewiseAbcenth
sed s/X*$// GenewiseAbcenth.pep.fasta > GenewiseAbcenth.pep.fasta.tmp
mv GenewiseAbcenth.pep.fasta.tmp GenewiseAbcenth.pep.fasta

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

f="GenewiseAbcenth.pep.fasta"

interproscan.sh -i $f -t p -goterms -iprlookup -cpu 4     ## COMMENTED BECAUSE IT IS ALREADY RUN!!!!!


######## Run blast with ORs to obtain names
module load ncbi-blast/2.11.0+


blastp -query GenewiseAbcenth.pep.fasta -db /home/projects/ku_00039/people/joeviz/OR_annotation/Chemo_db/ORco_db.fasta -outfmt "6 std qlen slen" -out GenewiseAbcenth.pep.fasta.ORcoblast.txt -num_threads 4
blastp -query GenewiseAbcenth.pep.fasta -db /home/projects/ku_00039/people/joeviz/OR_annotation/Chemo_db/GR2_db.fasta -outfmt "6 std qlen slen" -out GenewiseAbcenth.pep.fasta.GRblast.txt -num_threads 4 -max_target_seqs 5

cat /home/projects/ku_00039/people/joeviz/OR_annotation/OR_db/OR_db.fasta $ORANNOTATION > OR_masAbcenth_db.fasta
makeblastdb -in OR_masAbcenth_db.fasta -dbtype prot
blastp -query GenewiseAbcenth.pep.fasta -db OR_masAbcenth_db.fasta -outfmt "6 std qlen slen" -out GenewiseAbcenth.pep.fasta.ORblast.txt -num_threads 4 -max_target_seqs 5

blastp -query GenewiseAbcenth.pep.fasta -db /home/projects/ku_00039/people/joeviz/OR_annotation/Chemo_db/GR_dmel_db.fasta -outfmt "6 std qlen slen" -out GenewiseAbcenth.pep.fasta.GRdmelblast.txt -num_threads 4 -max_target_seqs 5


######## Run script to rename the gff3 and generate the protein file and summary table

perl /home/projects/ku_00039/people/joeviz/OR_annotation/scripts/run_GR_classification_happy_GR.pl GenewiseAbcenth.gff3 $OUTDIR $GENOME


    
