#!/bin/bash

# Create array job to submit

#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/stLFR_assemblies_withmasurca/newmasurca/*fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/stLFR_assemblies_withmasurca/*.fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/PacBio_nonpolished_assemblies_allfinalbatchv2/*fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/stLFR_assemblies_2batch/*fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Polished_assemblies_9batch/*fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/tmpcheck/*fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Polished_assemblies_SLRscaffolder_4batch/*fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Polished_assemblies_SLRscaffolder/*fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Polished_assemblies_thirdbatch/*fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Polished_assemblies/*fasta
FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/GAGA_all_final_assemblies/*fasta
QUERYF=/home/projects/ku_00039/people/joeviz/bitacora_barcoding/barcoding_nt_DB

i=0

for f in $FILES
do
	FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+).fa/){print "$1"}' $f)"
	mkdir -p $FILENAME
	#cd $FILENAME

	#qsub -pe ompi24 1 -cwd -V -j y -N runbt_$FILENAME -b y /users-d1/jvizueta/bitacora/runBITACORA_command_line.sh -m genome -q $QUERYF -g $f -b T 

	i=$((i+1))

    echo "#!/bin/sh
cd $FILENAME
/home/projects/ku_00039/people/joeviz/bitacora_barcoding/bitacora_nt/runBITACORA_command_line.sh -m nt -q $QUERYF -p $f -b T -l 360 -e 1e-20
cd ..
    " > run_$i.sh

	#cd ..

done


    echo '#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn=1:thinnode
#PBS -l mem=16gb
#PBS -l walltime=6:00:00' > submit_run_all_nt.sh

	echo "#PBS -t 1-$i" >> submit_run_all_nt.sh

echo '
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
module load anaconda3/4.4.0
module load perl
#module load ncbi-blast/2.10.0+ #it yields inconsisten results when multithreading and therefore it is best to use 2.2 version
module load ncbi-blast/2.2.31+
module load hmmer/3.2.1
module load java/1.8.0


# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

bash run_${PBS_ARRAYID}.sh

    ' >> submit_run_all_nt.sh

