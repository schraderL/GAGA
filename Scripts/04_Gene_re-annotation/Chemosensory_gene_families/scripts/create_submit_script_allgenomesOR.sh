#!/bin/bash

# Script to create a job array to run the pipeline in each genome. We use the GAGA-ID as filename, and the specific paths need to be modified in the script.
# Note that this example is set with the ORs

FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/GAGA_all_final_assemblies_oldv/*fasta

i=0

for f in $FILES
do
#    FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+).fasta/){print "$1"}' $f)"
	FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*(GAGA-\d\d\d\d)\S+/ || $ARGV[0] =~ /.*(NCBI-\d\d\d\d)\S+/ || $ARGV[0] =~ /.*(OUT-\d\d\d\d)\S+/){print "$1"}' $f)"
#    mkdir -p $FILENAME

    echo '#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn=40
#PBS -l mem=160gb
#PBS -l walltime=24:00:00' > submit_abcenth_$FILENAME.sh

echo '
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
conda activate happyv2


# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

GENOME="'"$f"'"
OUTDIR="'"$FILENAME"'" # Directory must not exist before, or it will fail (unless it is used --overwrite)

HAPpy --threads 40 --annotator ABCENTH --hmm_dir /home/projects/ku_00039/people/joeviz/OR_annotation/RefHMMs_HsSrObCfSi_0.45dist --genome $GENOME --output_dir $OUTDIR

cd $OUTDIR

######## Generate a GFF3 and protein file

module load genometools/1.5.10

perl /home/projects/ku_00039/people/joeviz/OR_annotation/get_abcenth_gtf_corrected.pl ABCENTH.gtf $GENOME
gt gtf_to_gff3 -o ABCENTH.gff3 ABCENTH_corrected.gtf
gt gff3 -sort -tidy -retainids -o ABCENTH_clean.gff3 ABCENTH.gff3

perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $GENOME ABCENTH_clean.gff3 ABCENTH_clean
sed 's/X*$//' ABCENTH_clean.pep.fasta > ABCENTH_clean.pep.fasta.tmp
mv ABCENTH_clean.pep.fasta.tmp ABCENTH_clean.pep.fasta


######## Run Interpro in the protein set
module purge
module load ngs tools
module load anaconda3/4.4.0
module load openjdk/15.0.1
module load perl
module load signalp/4.1g
#module load interproscan/5.51-85.0
module load interproscan/5.52-86.0

f="ABCENTH_clean.pep.fasta"

/services/tools/interproscan/5.47-82.0/interproscan.sh -i $f -t p -goterms -iprlookup -cpu 40



######## Run blast with ORs to obtain names
module load ncbi-blast/2.11.0+


blastp -query ABCENTH_clean.pep.fasta -db /home/projects/ku_00039/people/joeviz/OR_annotation/OR_db/ORco_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.ORcoblast.txt -num_threads 40
blastp -query ABCENTH_clean.pep.fasta -db /home/projects/ku_00039/people/joeviz/OR_annotation/OR_db/OR_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.ORblast.txt -num_threads 40 -max_target_seqs 5
blastp -query ABCENTH_clean.pep.fasta -db /home/projects/ku_00039/people/joeviz/OR_annotation/OR_db/GR2_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.GRblast.txt -num_threads 40 -max_target_seqs 5


######## Run script to rename the gff3 and generate the protein file and summary table

perl /home/projects/ku_00039/people/joeviz/OR_annotation/run_OR_classification.pl ABCENTH_clean.gff3 $OUTDIR $GENOME


    ' >> submit_abcenth_$FILENAME.sh


echo "qsub submit_abcenth_$FILENAME.sh"


done




