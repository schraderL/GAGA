#!/bin/bash

# Script to create the scripts to run CAFE and submit it in Computerome cluster
# First edit the input, parameters and command to run each step
# Usage: bash run_cafe_analyses.sh

#### INPUT ####

# Table with gene family counts
INTABLE="/home/projects/ku_00039/people/joeviz/CAFE/runCAFE_withoutgroups/N0_GeneCounts_nostLFR_sfilt20_nolow.tsv" # Table with gene family counts filtering orthogroups with high variance
#INTABLE="/home/projects/ku_00039/people/joeviz/CAFE/runCAFE_withoutgroups/N0_GeneCounts_nostLFR_noTEs.tsv" # Table with gene family counts without filtering 

# Dated tree (species names need to be the same as in the table)
INTREE="/home/projects/ku_00039/people/joeviz/CAFE/runCAFE_withoutgroups/GAGA_dated_phylogeny_newick_woutgroup_DatesPeters2017_CurrBio_nostLFR.tre" # Ultrametric tree

# Tree with the number of lambdas to estimate across the tree. Comment for global lambda model
#LAMBDATREE="/home/projects/ku_00039/people/joeviz/CAFE/runCAFE_withoutgroups/GAGA_dated_phylogeny_newick_woutgroup_DatesPeters2017_CurrBio_nostLFR.tre" # Tree with lambda for multiple lambda model

# Output directory to run CAFE
DIR="CAFE_lambda_witherrormodel" # Directory to create and run CAFE

# Output name prefix for CAFE output
OUTPUT="singlelambdapoisson" # Output prefix for each CAFE run


#### PARAMETERS ####

NTHREADS=8 # Number of threads for each CAFE run
NUMRUNS=20 # Number of independent runs of CAFE to ensure convergence
NTIME=48 # Number of hours in qsub script, use 150 for long runs estimating lambda


#### CAFE commands - Comment/Uncomment the COMMAND line to run each step ####

# Step 1 - Global lambda model using a poisson distribution 
#COMMAND="/home/projects/ku_00039/people/joeviz/CAFE/CAFE/CAFE5/bin/cafe5 -i $INTABLE -t $INTREE -c $NTHREADS -p " # CAFE command, single lambda poisson

# Step 2 - Global lambda model using a poisson distribution, with the input table filtering genes with high variance. Edit the INTABLE (line 10). 
COMMAND="/home/projects/ku_00039/people/joeviz/CAFE/CAFE/CAFE5/bin/cafe5 -i $INTABLE -t $INTREE -c $NTHREADS -p " # CAFE command, single lambda poisson

# Step 3 - Global lambda model using a poisson distribution and estimating an error model. Using the input table filtering genes with high variance. Edit the INTABLE (line 10).  
#COMMAND="/home/projects/ku_00039/people/joeviz/CAFE/CAFE/CAFE5/bin/cafe5 -i $INTABLE -t $INTREE -c $NTHREADS -p -e" # CAFE command, with error model

# Step 4 - Multiple lambda model using a poisson distribution and the error model estimated in the best run of Step 3. Using the input table filtering genes with high variance, Edit the INTABLE (line 10) and set the file LAMBDATREE (line 17) containing the tree with multiple lambdas. 
# Edit also the file in -e (in the line below) which should contain the error model estimated in Step 3
#COMMAND="/home/projects/ku_00039/people/joeviz/CAFE/CAFE/CAFE5/bin/cafe5 -i $INTABLE -t $INTREE -c $NTHREADS -y $LAMBDATREE -p -e/home/projects/ku_00039/people/joeviz/CAFE/runCAFE/2set_nostlfr_withoutlep/CAFE_singlelambda_poisson_sfilt20_nofiltlow_errormodel_addruns/singlelambdapoisson_run3/Base_error_model.txt " # CAFE command, multiple lambda with error model already calculated

# Step 5 - Fixed lambda(s) based on the best fitting model, a poisson distribution and error model estimated in Step 3. Using the input table containing ALL orthogroups. Edit the INTABLE (line 10). 
# Edit also the file in -e which should contain the error model estimated in Step 3
# Edit the lambda value(s) in -l  containing the lamba(s) estimated in the best fitting model
# Command if the best model is a global lambda (it used -l to indicate the lambda value)
#COMMAND="/home/projects/ku_00039/people/joeviz/CAFE/CAFE/CAFE5/bin/cafe5 -i $INTABLE -t $INTREE -c $NTHREADS -p -l 0.0048104885261371 -e/home/projects/ku_00039/people/joeviz/CAFE/runCAFE/2set_nostlfr_withoutlep/CAFE_singlelambda_poisson_sfilt20_nofiltlow_errormodel_addruns/singlelambdapoisson_run3/Base_error_model.txt" # CAFE command, global fixed lambda with error model already calculated
# Command if the best model is a multiple lambda (it uses -m to indicate the multiple lambda values corresponding the order to the lambda number if the phylogeny in -y). Set the file LAMBDATREE (line 17) 
#COMMAND="/home/projects/ku_00039/people/joeviz/CAFE/CAFE/CAFE5/bin/cafe5 -i $INTABLE -t $INTREE -c $NTHREADS -m 0.003845644880218,0.0050413215441222 -y $LAMBDATREE -p -e/home/projects/ku_00039/people/joeviz/CAFE/runCAFE/2set_nostlfr_withoutlep/CAFE_singlelambda_poisson_sfilt20_nofiltlow_errormodel_addruns/singlelambdapoisson_run3/Base_error_model.txt " # CAFE command, multiple  fixed lambdas with error model already calculated


# Other models and parameters for CAFE

# Fixing lambda without using error model. Note in our GAGA data: the error model had better fit than without using the error model.
#COMMAND="/home/projects/ku_00039/people/joeviz/CAFE/CAFE/CAFE5/bin/cafe5 -i $INTABLE -t $INTREE -c $NTHREADS -p -l 0.0005" # CAFE command fix lambda

# Global lambda model using a poisson distribution and two discrete gamma rate (-k 2) categories. Note in our GAGA data: The gamma model did not converge even in the global lambda model (simplest).
#COMMAND="/home/projects/ku_00039/people/joeviz/CAFE/CAFE/CAFE5/bin/cafe5 -i $INTABLE -t $INTREE -c $NTHREADS -p -k 2 " # CAFE command, gamma model

# Multiple lambda model using a poisson distribution without error model. Note in our GAGA data: the error model had better fit than without using the error model.
#COMMAND="/home/projects/ku_00039/people/joeviz/CAFE/CAFE/CAFE5/bin/cafe5 -i $INTABLE -t $INTREE -c $NTHREADS -y $LAMBDATREE -p " # CAFE command, multiple lambda


#### Code ####

mkdir -p $DIR
cd $DIR

i=1

until [ $i -gt $NUMRUNS ] # Iterate to create files for each independent run
do

    OUTNAME="$OUTPUT""_run$i"


   echo '#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn='"$NTHREADS"':thinnode
#PBS -l mem=64gb
#PBS -l walltime='"$NTIME"':00:00' > submit_cafe_"$OUTNAME".sh


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
#module load anaconda3/4.4.0
module load anaconda3/2022.10
#module load gcc/12.1.0
module load gcc/7.4.0

# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

    ' >> submit_cafe_"$OUTNAME".sh


    echo "$COMMAND -o $OUTNAME > $OUTNAME""_stdout.txt 2> $OUTNAME""_stderr.txt" >> submit_cafe_"$OUTNAME".sh
 

    echo "qsub submit_cafe_$OUTNAME.sh" >> submit_all_cafe.sh


    ((i=i+1))

done





