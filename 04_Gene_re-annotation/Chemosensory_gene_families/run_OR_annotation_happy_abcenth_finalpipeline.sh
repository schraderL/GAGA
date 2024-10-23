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
#PBS -l mem=160gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here 12:00:00, 12 hours)
#PBS -l walltime=24:00:00
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
#First load modules
module load ngs
module load tools
module load miniconda3/4.10.3


# Then activate conda happy
conda activate happyv2


# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

GENOME="/home/projects/ku_00039/people/joeviz/GAGA_genomes/GAGA-0001_genome.fasta"
OUTDIR="GAGA-0001" # Directory must not exist before, or it will fail (unless it is used --overwrite)

HAPpy --threads 40 --annotator ABCENTH --hmm_dir /home/projects/ku_00039/people/joeviz/OR_annotation/RefHMMs_HsSrObCfSi_0.45dist --genome $GENOME --output_dir $OUTDIR

cd GAGA-0001 # Change same as OUTDIR

######## Generate a GFF3 and protein file

## Deprecated steps
#cat ABCENTH.gtf | awk -F'[\t; ]' '{a[$10] = a[$10] "\n" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\tID=" $10 ":CDS-" $4 ";Parent=" $10 "-RA"}!starts[$10]||$4 < starts[$10]{starts[$10] = $4}!stops[$10]||$5>stopts[$10]{stops[$10]=$5}{loc[$10]=$1;source[$10]=$2;strand[$10]=$7}{for (i in a) print loc[i] "\t" source[i] "\tgene\t" starts[i] "\t" stops[i] "\t.\t" strand[i] "\t.\tID=" i "\n" loc[i] "\t" source[i] "\tmRNA\t" starts[i] "\t" stops[i] "\t.\t" strand[i] "\t.\tID=" i "-RA;Parent=" i gensub("CDS","exon","G",a[i]) a[i]}' > ABCENTH.gff3
#perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $GENOME ABCENTH.gff3 ABCENTH
## Sean command fails

#module load gffread/0.12.4
#gffread -E ABCENTH.gtf -o- > ABCENTH_gffread.gff3


module load genometools/1.5.10
gt gtf_to_gff3 -o ABCENTH.gff3 ABCENTH.gtf
gt gff3 -sort -tidy -retainids -o ABCENTH_clean.gff3 ABCENTH.gff3

#perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $GENOME test.gff3 test
#perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $GENOME testclean.gff3 testclean
perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $GENOME ABCENTH_clean.gff3 ABCENTH_clean
sed 's/X*$//' ABCENTH_clean.pep.fasta > ABCENTH_clean.pep.fasta.tmp
mv ABCENTH_clean.pep.fasta.tmp ABCENTH_clean.pep.fasta



### Validation of the annotated ORs, comment from here if you would like to skip the next commands (Interpro and blast) 

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


blastp -query ABCENTH_clean.pep.fasta -db Chemo_db/ORco_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.ORcoblast.txt -num_threads 40
blastp -query ABCENTH_clean.pep.fasta -db Chemo_db/OR_db/OR_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.ORblast.txt -num_threads 40 -max_target_seqs 5
blastp -query ABCENTH_clean.pep.fasta -db Chemo_db/OR_db/GR2_db.fasta -outfmt "6 std qlen slen" -out ABCENTH_clean.pep.fasta.GRblast.txt -num_threads 40 -max_target_seqs 5


######## Run script to rename the gff3 and generate the protein file and summary table


perl /home/projects/ku_00039/people/joeviz/OR_annotation/run_OR_classification.pl ABCENTH_clean.gff3 $OUTDIR $GENOME















