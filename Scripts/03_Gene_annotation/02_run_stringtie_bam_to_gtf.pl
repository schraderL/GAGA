#!/usr/bin/perl
use strict;
use warnings;

#########################################################
#####
#
# Create script to submit stringtie
#
#####

# Usage
# perl run_stringtie_bam_to_gtf folder_with_bams(folder obtained from running run_gemoma_annotation.pl ) 


## Input variables

my $idlist = $ARGV[0];

my $gemomadir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/all/run_hic_batch/";
#my $gemomadir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/all/run_firstbatch_goodmapping_gemomanompha/"; # Dir with gemoma output
#my $gemomadir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/all/run_all_norna_nompha/"; # Dir with gemoma output
#my $gemomadir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/all/run_secondbatch_subproject_goodmapping_gemomanompha_estesi"; # Dir with gemoma output

#my $assemblydir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Repeat_masked_assemblies/GAGA_annotations/"; # Folder with unzipped soft-masked genome assemblies classified in folders (as in the directory GAGA_annotations in ERDA)
#my $rnadir = "/home/projects/ku_00039/people/joeviz/RNA-seq/GAGA-IDs";

#my $gemomainputdir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Dataset/Whole_genomes/"; # Folder containing the folders with the annotations used in gemoma

my $gagaidlist = "/home/projects/ku_00039/people/joeviz/GAGA_species_list.txt";


## Start

my $numid = "0";
open(File, "<", $idlist);
while(<File>){
	chomp;
	my $gagaid = $_;
	next if ($gagaid !~ /\S+/);
	$gagaid =~ s/\s+//;


	# Retrieve GFF from gemoma
	my $gfffile = "$gemomadir\/$gagaid\/gemoma_output/final_annotation.gff";

	## Check if GFF is well read
	system ("wc -l $gfffile > testdb.out 2> /dev/null");
	my $dbout = 0;
	open (File2, "<", "testdb.out");
	while(<File2>){
		chomp;
		my $line2 = $_;
		next if ($line2 !~ /\S+/);	
		if ($line2 =~ /(\d+) /){
			$dbout = $1;
		} else {
			die "ERROR in reading gemoma gff from $gagaid in $gfffile, testdb.out in $line2\n";
		}
	}
	close File2;

	if ($dbout == 0){
		print "ERROR in $gagaid gemoma gff,  $gfffile is empty!\n";	
		next;	
	}

	# Create script to run Stringtie

	$numid++;
	open ("Outscript", ">", "run_stringtie_$numid\.sh");
	print Outscript "\#\!\/bin\/bash\n\n";
	#print Outscript "cd $gagaid\n\n";
	

	## Find bam files
	system ("ls $gagaid\/\*\.bam > inbam.tmp 2> /dev/null");
	#print "ls $gagaid\/\*\.bam > inbam.tmp 2> /dev/null\n";
	my $bamfilenum = "0";
	my $allstringtout = "";
	open (Tmpfile, "<", "inbam.tmp");
	while (<Tmpfile>){
		chomp;
		my $nline = $_;
		next if ($nline !~ /\S+/);		
		if ($nline =~ /(\S+)Aligned\.sortedByCoord\.out\.bam/){
			my $rnafile = $1;
			my $bamfile = $nline;
			my $logfile = "$rnafile"."Log.final.out";
			$bamfilenum++;

			my $sampleid = "$gagaid\_"."$bamfilenum";
			my $stringout = "$rnafile\_stringtie.gtf";

			$allstringtout .= "$stringout ";

			print Outscript "stringtie -p 8 -G $gfffile -o $stringout -l $sampleid $bamfile\n";

		}
	}
	close Tmpfile;
	system ("rm inbam.tmp");

	if ($bamfilenum < 1){
		print "ERROR: No BAM file found for $gagaid\n";
		next;
	} elsif ($bamfilenum == 1) {
		#print Outscript "\ncp $allstringtout $gagaid\/$gagaid\_stringtie\_allmerged.gtf\n"; Run also merge
		print Outscript "\nstringtie --merge -G $gfffile -o $gagaid\/$gagaid\_stringtie\_allmerged.gtf $allstringtout\n";		


	} else {
		print Outscript "\nstringtie --merge -G $gfffile -o $gagaid\/$gagaid\_stringtie\_allmerged.gtf $allstringtout\n";

	}

	close Outscript

}
close File;

## Create submission script

open ("Submitscript", ">", "submit_run_stringtie_all\.sh");
print Submitscript "#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn=8
#PBS -l mem=100gb
#PBS -l walltime=24:00:00
#PBS -t 1-$numid

# Go to the directory from where the job was submitted (initial directory is \$HOME)
echo Working directory is \$PBS_O_WORKDIR
cd \$PBS_O_WORKDIR

### Here follows the user commands:
# Define number of processors
NPROCS=`wc -l < \$PBS_NODEFILE`
echo This job has allocated \$NPROCS nodes

# Load all required modules for the job

module load tools
module load ngs
module load anaconda3/4.4.0
module load perl
module load stringtie/2.1.5

# This is where the work is done

bash run_stringtie_\$\{PBS_ARRAYID\}.sh


";


