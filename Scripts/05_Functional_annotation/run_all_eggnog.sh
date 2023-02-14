#!/bin/bash

# Create array job to submit EGGNOG in all genome annotations

FILES=/home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_GAGA_annotations/orthofinder_input/*fasta # Protein dataset

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

	tmpdir=\"\$(pwd)\"\"/temp\"
	mkdir -p \$tmpdir
	curdir=\"\$(pwd)\"
    scrdir=\"\$(pwd)\"\"/scratchdir\"
    mkdir -p \$scrdir
    
    emapper.py --cpu 40 -i $f --itype proteins --dbmem --data_dir /services/tools/eggnog-mapper/2.1.6/data/ --output out_$GAGAID --output_dir \$curdir --temp_dir \$tmpdir --scratch_dir \$scrdir -m diamond --evalue 0.001 --score 40 --pident 20 --query_cover 20 --subject_cover 20 --dmnd_ignore_warnings --tax_scope auto --target_orthologs all --go_evidence non-electronic --report_orthologs
    " > run_eggnog_$i.sh

done


    echo '#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn=40
#PBS -l mem=160gb
#PBS -l walltime=12:00:00' > submit_run_eggnog_all.sh

    echo "#PBS -t 1-$i" >> submit_run_eggnog_all.sh

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
#module load anaconda3/4.4.0
module load anaconda3/2021.05
module load hmmer/3.2.1
#module load diamond/2.0.6
module load diamond/2.0.13
#module load eggnog-mapper/2.0.8
module load eggnog-mapper/2.1.6


# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

bash run_eggnog_${PBS_ARRAYID}.sh

    ' >> submit_run_eggnog_all.sh






