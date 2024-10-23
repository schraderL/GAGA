#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn=4
#PBS -l mem=16gb
#PBS -l walltime=48:00:00


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
#module load openjdk/15.0.1
module load perl
module load ncbi-blast/2.2.31+
module load hmmer/3.2.1
module load mafft/7.453
module load java/1.8.0


# This is where the work is done

OUTDIR="/home/projects/ku_00039/people/joeviz/gene_family_pipeline/run_pipeline_Ectatomma/run_pipeline" # Output directory
mkdir -p $OUTDIR
cd $OUTDIR

### Input files
path="/home/projects/ku_00039/people/joeviz/gene_family_pipeline/run_pipeline_Ectatomma/Gene-annotation-pipeline"
table="/home/projects/ku_00039/people/joeviz/gene_family_pipeline/run_pipeline_Ectatomma/Gene_reannotation_full_table.xlsx"
gene_fam_db="/home/projects/ku_00039/people/joeviz/gene_family_pipeline/run_pipeline_Ectatomma/Gene_family_database"
genome_name="Ectatomma"
output_directory="/home/projects/ku_00039/people/joeviz/gene_family_pipeline/run_pipeline_Ectatomma/run_pipeline/Ectatomma"
proteome="/home/projects/ku_00039/people/joeviz/GAGA_genomes/Colab_Ectatomma/Gene_annotation_final_repfilt/Ectatomma_ruidum_final_annotation_representative_repfilt.pep.fasta"
gff="/home/projects/ku_00039/people/joeviz/GAGA_genomes/Colab_Ectatomma/Gene_annotation_final_repfilt/Ectatomma_ruidum_final_annotation_representative_repfilt.gff3"
genome="/home/projects/ku_00039/people/joeviz/GAGA_genomes/Colab_Ectatomma/repeat_annot/Ectatomma_ruidum_Hifigenome.nextpolish_final.softMasked.fasta"
interpro="/home/projects/ku_00039/people/joeviz/GAGA_genomes/Colab_Ectatomma/Gene_annotation_final_repfilt/Interpro_out/Ectatomma_ruidum_final_annotation_representative_repfilt.pep.fasta.tsv"
run_bitacora="/home/projects/ku_00039/people/joeviz/programs/bitacora/runBITACORA_command_line.sh" # Make sure to edit this script to include the path to BITACORA scripts and GeMoMa jar file
num_threads=4
############################

python3 "$path"/Scripts/run_analysis.py $path $table $gene_fam_db $proteome $interpro $gff $genome $run_bitacora $genome_name $output_directory $num_threads


