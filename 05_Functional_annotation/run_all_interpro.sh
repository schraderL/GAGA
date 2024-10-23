#!/bin/bash

# Create array job to submit
## Edit also outdir in interpro command: currently to interpro_batch3

#FILES=/home/projects/ku_00039/people/joeviz/GAGA_annotations/Monomorium_pharaonis/*pep
#FILES=$(cat ../protein_list_GAGA_all_norna_nompha_20210624.txt)
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_files/*/*pep
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_files_batch4/GAGA_annotations/*/*annotation_representative_repfilt.pep.fasta
#FILES=/home/projects/ku_00039/people/joeviz/interpr/testset.fasta
FILES=/home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_GAGA_annotations/orthofinder_input/*fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/new_outgroups/final_annotations/*representative.gff.pep

i=0

for f in $FILES
do
#	FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+).pep/){print "$1"}' $f)"
#	FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+)\/(\S+).fasta/){print "$1\_$2"}' $f)"
#	FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+)\/(\S+)\/(\S+)\/(\S+).fasta/){print "$1\_$2\_$3\_$4"}' $f)"
#	FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+)\/(\S+)\/(\S+)\/(\S+).pep/){print "$1\_$2\_$3\_$4"}' $f)"  ### If file ends in .pep
#	mkdir -p $FILENAME

#	GAGAID="$(perl -e 'if ($ARGV[0] =~ /.*(GAGA-\d\d\d\d)\S+.fasta/ || $ARGV[0] =~ /.*(NCBI-\d\d\d\d)\S+.fasta/ || $ARGV[0] =~ /.*(OUT-\d\d\d\d)\S+.fasta/){print "$1"}' $f)"
#	GAGAID="$(perl -e 'if ($ARGV[0] =~ /.*(GAGA-\d\d\d\d)\S+.pep/ || $ARGV[0] =~ /.*(NCBI-\d\d\d\d)\S+.pep/ || $ARGV[0] =~ /.*(OUT-\d\d\d\d)\S+.pep/){print "$1"}' $f)" ### .pep version instead of fasta
	GAGAID="$(perl -e 'if ($ARGV[0] =~ /.*(GAGA-\d\d\d\d)\S*.fasta/ || $ARGV[0] =~ /.*(NCBI-\d\d\d\d)\S*.fasta/ || $ARGV[0] =~ /.*(OUT-\d\d\d\d)\S*.fasta/){print "$1"}' $f)"

    i=$((i+1))

    echo "#!/bin/sh
    mkdir -p $GAGAID
    cd $GAGAID
interproscan.sh -i $f -t p -goterms -iprlookup -d /home/projects/ku_00039/people/joeviz/functional_annotation/interpro_finalannot/$GAGAID -cpu 40
    " > run_interpro_$i.sh

done


    echo '#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn=40
#PBS -l mem=160gb
#PBS -l walltime=24:00:00' > submit_run_interpro_all.sh

    echo "#PBS -t 1-$i" >> submit_run_interpro_all.sh

echo '
# Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

### Here follows the user commands:
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes

# Load all required modules for the job
module load ngs tools
module load anaconda3/4.4.0
module load jdk/19
module load perl
module load signalp/4.1g
module load interproscan/5.52-86.0


# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

bash run_interpro_${PBS_ARRAYID}.sh

    ' >> submit_run_interpro_all.sh






