#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn=20
#PBS -l mem=180gb
#PBS -l walltime=96:00:00


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
module load bwa/0.7.16a


# Code

# First fake input assembly
# intalled in: /home/projects/ku_00039/people/joeviz/programs/SLR-superscaffolder/SLR-superscaffolder_installed/ 
#YOUR-INSTALL-DIR/bin/FakeSOAPContig < your-contig-sequence-file 1>your-prefix.contig 2>your-prefix.ContigIndex

#Then: Prepare the conf.ini file
#3rd: Generate the pipeline and working folder
#YOUR-INSTALL-DIR/scaffold/prepare.sh ./your-conf.ini

cd slr_output
./run.sh >log_pipeline_run2try 2>&1 

    
