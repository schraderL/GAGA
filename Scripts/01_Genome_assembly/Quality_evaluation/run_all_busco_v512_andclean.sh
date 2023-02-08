#!/bin/bash

# Number of threads per busco
NTHREADS=20

# Create array job to submit

FILES=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/GAGA_all_final_assemblies_softmasked/*fasta

i=0

for f in $FILES
do
    FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+).fasta/){print "$1"}' $f)"
#    mkdir -p $FILENAME

    i=$((i+1))

    echo "#!/bin/sh
perl /home/projects/ku_00039/people/joeviz/scripts/N50Stat_v3.pl -i $f
perl /home/projects/ku_00039/people/joeviz/scripts/rename_assembly_scaffold.pl $f
/services/tools/busco/5.1.2/bin/busco --config /home/projects/ku_00039/people/joeviz/busco/config.ini -i $f -c $NTHREADS -m geno --auto-lineage-euk --augustus --augustus_species fly --offline -o busco_$FILENAME
rm -rf busco_$FILENAME/run_*
rm -rf busco_$FILENAME/logs
rm -rf busco_$FILENAME/auto_lineage
rm -rf busco_$FILENAME/blast_db
    " > run_busco_$i.sh

done


    echo '#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn='"$NTHREADS"'
#PBS -l mem=80gb
#PBS -l walltime=48:00:00' > submit_run_busco_all.sh

    echo "#PBS -t 1-$i" >> submit_run_busco_all.sh

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

module load lp_solve/5.5.2.5
module load suitesparse/5.7.2
module load augustus/3.3.3

module load prodigal/2.6.3

module load java/1.8.0
module load pasta/20200406
module load sepp/4.3.10
module load metaeuk/4-a0f584d

#module load busco/4.0.6
module load busco/5.1.2


export AUGUSTUS_CONFIG_PATH="/home/projects/ku_00039/people/joeviz/busco/augustus_config"


# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

bash run_busco_${PBS_ARRAYID}.sh

    ' >> submit_run_busco_all.sh





