#!/bin/bash

## Maximum array length is 2000, so grouping multiple commands per run

# Error from translatorX: Length (681) of amino acid seq from GAGA-0365_TrufIR101 does not correspond to 1/3 of nucelotide seq (682)

# Number of threads per alignment
NTHREADS=4

# Create array job to submit

#FILES=/home/projects/ku_00039/people/joeviz/orthology_alignments/test_orthogroups_singlecopyv2/*pep.fasta
#FILES=/home/projects/ku_00039/people/joeviz/orthology_alignments/single_copy_orthogroups_speciestree/orthogroups_seqs/*pep.fasta
#FILES=/home/projects/ku_00039/people/joeviz/orthology_alignments/all_orthogroups/orthogroups_seqs/*pep.fasta
FILES=/home/projects/ku_00039/people/joeviz/orthology_alignments/orthogroups_with_outgroups/single_copy_100pc/*pep.fasta

#OUTDIR=single_copy_speciestree

#mkdir -p $OUTDIR
#mkdir -p codon_alignments
mkdir -p protein_alignments
#mkdir -p protein_codon_alignments
#mkdir -p protein_aln_trees

OUTDIR=protein_aln_trees_trimal
mkdir -p protein_aln_trees_trimal

i=1
e=0

for f in $FILES
do
    FILENAME="$(perl -e 'if ($ARGV[0] =~ /.*\/(\S+)\.pep.fasta/){print "$1"}' $f)"
    FILECDS="$(perl -e 'if ($ARGV[0] =~ /(.*\/\S+)\.pep.fasta/){print "$1".".cds.fasta"}' $f)"
#    mkdir -p $FILENAME


	echo "#!/bin/sh
	" >> run_aln_$i.sh


#	if [ $e -gt 1 ]
#	then   

#	    echo "iqtree2 -s $f -bb 1000 -T $NTHREADS -msub nuclear

		# PRANK
#	    echo "/users-d1/jvizueta/programs/prank_debian7/prank -d=$outdir/$file\_cds_codon3multiple.fasta -o=$outdir/$file\_cds_aligned.fasta -codon -noxml -notree -F
#		echo "mv $outdir/$file\_cds_aligned.fasta.2.fas $outdir/$file\_cds_aligned.fasta
#		echo "mv $outdir/$file\_*dnd $outdir/$file\_*xml $outdir/$file\_*.1* $outdir/temp/

		# Prot aln
		# MAFFT
	   	echo "mafft-linsi --thread $NTHREADS --auto $f > protein_alignments/$FILENAME.pep.aln
	   	" >> run_aln_$i.sh

		# Tcoffee
#		echo "t_coffee $f -mode mcoffee -multi_core jobs
#		perl /programs/clustal2fasta.pl < $f\_proteins.aln > protein_alignments/$FILENAME\.pep.aln
#		" >> run_aln_$i.sh
		#Muscle
#		echo "muscle -in $f -out protein_alignments/$FILENAME\.pep.aln
#		" >> run_aln_$i.sh

		# CDS aln from prot aln
#		echo "perl /home/projects/ku_00039/people/joeviz/programs/translatorx_vLocal.pl -i $FILECDS -o protein_codon_alignments/$FILENAME.cds.aln -a protein_alignments/$FILENAME.pep.aln
#		mv protein_codon_alignments/$FILENAME.cds.aln.nt_ali.fasta protein_codon_alignments/$FILENAME.cds.aln
#		rm -rf protein_codon_alignments/$FILENAME*html protein_codon_alignments/$FILENAME*ali.*
#		" >> run_aln_$i.sh
#		#pal2nal
#		echo "perl /home/projects/ku_00039/people/joeviz/programs/pal2nal.pl protein_alignments/$FILENAME.pep.aln $FILECDS -output fasta > protein_codon_alignments/$FILENAME.cds.pal2nal.aln
#		" >> run_aln_$i.sh	

		#trimal protein alignment
		echo "trimal -in protein_alignments/$FILENAME.pep.aln -out protein_alignments/$FILENAME.pep.trimalgap.aln -gappyout
		trimal -in protein_alignments/$FILENAME.pep.aln -out protein_alignments/$FILENAME.pep.trimalaut.aln -automated1
		" >> run_aln_$i.sh	

		#tree
		echo "perl /home/projects/ku_00039/people/joeviz/orthology_trees/rename_fasta_to_gagaid.pl protein_alignments/$FILENAME.pep.trimalgap.aln "$OUTDIR/$FILENAME"".renamed.aln"
		iqtree2 -s $OUTDIR/$FILENAME"".renamed.aln -B 1000 -T $NTHREADS -m MFP -msub nuclear --prefix $OUTDIR/$FILENAME
		" >> run_aln_$i.sh

#		echo "
#	    " >> run_aln_$i.sh

	    i=$((i+1))
	    e=0

#	else

#	    echo "iqtree2 -s $f -bb 1000 -T $NTHREADS -msub nuclear
	    # MAFFT
#	   	echo "mafft --thread $NTHREADS --auto $f > protein_alignments/$FILENAME.pep.aln
#	   	" >> run_aln_$i.sh
#		echo "perl /home/projects/ku_00039/people/joeviz/programs/translatorx_vLocal.pl -i $FILECDS -o protein_codon_alignments/$FILENAME.cds.aln -a protein_alignments/$FILENAME.pep.aln
#		mv protein_codon_alignments/$FILENAME.cds.aln.nt_ali.fasta protein_codon_alignments/$FILENAME.cds.aln
#		rm -rf protein_codon_alignments/$FILENAME*html protein_codon_alignments/$FILENAME*ali.*
#		" >> run_aln_$i.sh
#		echo "perl /home/projects/ku_00039/people/joeviz/programs/pal2nal.pl protein_alignments/$FILENAME.pep.aln $FILECDS -output fasta > protein_codon_alignments/$FILENAME.cds.pal2nal.aln
#		" >> run_aln_$i.sh	
#		#trimal protein alignment
#		echo "trimal -in protein_alignments/$FILENAME.pep.aln -out protein_alignments/$FILENAME.pep.trimalgap.aln -gappyout
#		trimal -in protein_alignments/$FILENAME.pep.aln -out protein_alignments/$FILENAME.pep.trimalaut.aln -automated1
#		" >> run_aln_$i.sh	
#		#tree
#		echo "perl /home/projects/ku_00039/people/joeviz/orthology_trees/rename_fasta_to_gagaid.pl protein_alignments/$FILENAME.pep.trimalgap.aln "$OUTDIR/$FILENAME"".renamed.aln"
#		iqtree2 -s $OUTDIR/$FILENAME"".renamed.aln -B 1000 -T $NTHREADS -m MFP -msub nuclear --prefix $OUTDIR/$FILENAME
#		" >> run_aln_$i.sh


#	    e=$((e+1))

#	fi


done


    echo '#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn='"$NTHREADS"':thinnode
#PBS -l mem=64gb
#PBS -l walltime=48:00:00' > submit_run_aln_all.sh

    echo "#PBS -t 1-$i" >> submit_run_aln_all.sh

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
module load mafft/7.453
module load muscle/3.8.425
module load trimal/1.4.1
#module load prank/170427
module load iqtree/2.1.3

# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

bash run_aln_${PBS_ARRAYID}.sh

    ' >> submit_run_aln_all.sh





