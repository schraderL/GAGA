#!/bin/bash

# Number of threads per job
NTHREADS=1

# Create array job to submit

#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/*/*fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Polished_assemblies_slavemaker_subproject/*fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Polished_assemblies_thirdbatch/*fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Polished_assemblies_8batch/*fasta
FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/tmpcheck/*_final.fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/stLFR_assemblies_2batch/*fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Fix_duplicated_busco_assemblies/assemblies_tocheck_busco/*fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Polished_assemblies_stLFRscaffolding/*fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Polished_assemblies_SLRscaffolder_3batch/*fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Final_PacBio_assemblies/*fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Final_PacBio_assemblies/temp/newbatch8/*final.fasta
#FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/Chr_level_assemblies/new_2batch_final/*final.fasta

i=0

for f in $FILES
do
    FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+).fasta/){print "$1"}' $f)""_dupsrm.fasta"
    FILENAMEL="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+).fasta/){print "$1"}' $f)""_dupsrm.log"
    FOLDERNAME="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+).fasta/){print "$1"}' $f)"
    mkdir -p $FOLDERNAME

    i=$((i+1))

    echo "#!/bin/sh
cd $FOLDERNAME
funannotate clean -i $f -o $FILENAME > $FILENAMEL
    " > run_job_$i.sh

done


    echo '#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn='"$NTHREADS"'
#PBS -l mem=60gb
#PBS -l walltime=24:00:00' > submit_run_all.sh

    echo "#PBS -t 1-$i" >> submit_run_all.sh

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
module load genemark-es/4.62
module load signalp/4.1c
module load funannotate/1.8.3


# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

bash run_job_${PBS_ARRAYID}.sh

    ' >> submit_run_all.sh





