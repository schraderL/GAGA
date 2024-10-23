#!/usr/bin/perl
use strict;
use warnings;

#########################################################
#####
#
# Filter repeats from annotation file
#
#####

# Usage
# perl run_gene_annotation_filtering_all.pl 

#perl run_gene_annotation_filtering.pl input.gff input.proteins.fasta input.genome input.genome.hardmasked input.repeat.gff GAGA-ID
#perl run_gene_annotation_filtering.pl /home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_files/GAGA-0014/GAGA-0014.annotation.merge.renamed.fixed.representative.gff3 /home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_files/GAGA-0014/GAGA-0014.annotation.merge.renamed.fixed.representative.pep /home/projects/ku_00039/people/joeviz/GAGA_annotations/Repeat_masked_assemblies/GAGA_annotations/GAGA-0014/GAGA-0014_nextpolish_correct_final_dupsrm_filt.softMasked.fasta /home/projects/ku_00039/people/joeviz/GAGA_annotations/Repeat_masked_assemblies/GAGA_annotations/GAGA-0014/GAGA-0014_nextpolish_correct_final_dupsrm_filt.repeatMasked.fasta /home/projects/ku_00039/people/joeviz/GAGA_annotations/Repeat_masked_assemblies/GAGA_annotations/GAGA-0014/GAGA-0014_nextpolish_correct_final_dupsrm_filt.repeats.gff GAGA-0014


## Input variables

#my $gagaid = "$ARGV[0]";
#my $annotdir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_files"; # Folder with unzipped soft-masked genome assemblies classified in folders (as in the directory GAGA_annotations in ERDA)
my $annotdir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_files_batch4/GAGA_annotations/";
my $repeatdir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Repeat_masked_assemblies/GAGA_annotations/"; # Folder with unzipped soft-masked genome assemblies classified in folders (as in the directory GAGA_annotations in ERDA)


## Start

system ("ls $annotdir/*/*representative.gff3 > infiles.txt");

my $numid = "0";
open(File, "<", "infiles.txt");
while(<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	$line =~ s/\s+//;

	my $ingff = "$line";
	my $inprots = "";
	my $ingenome = "";
	my $ingenomemasked = "";
	my $inrepeatgff = ""; 
	my $gagaid = "";
	my $infullgff = "";

	if ($line =~ /.*(GAGA-\d\d\d\d)\S+/ || $line =~ /.*(NCBI-\d\d\d\d)\S+/ || $line =~ /.*(OUT-\d\d\d\d)\S+/){
		$gagaid = $1;
	} else {die "Can't find GAGA-ID in $line\n";}   


	## Find prot fasta
	system ("ls $annotdir\/$gagaid\/\*representative*pep > inassembly.tmp 2> /dev/null");
	open (Tmpfile, "<", "inassembly.tmp");
	while (<Tmpfile>){
		chomp;
		my $nline = $_;
		next if ($nline !~ /\S+/);		
		if ($nline =~ /pep/){
			$inprots = $nline;
		}
	}
	close Tmpfile;
	system ("rm inassembly.tmp");

	if ($inprots !~ /pep/){
		print "ERROR: Skipping $gagaid Protein fasta could not be found, Maybe it is gzipped?, in $annotdir\/$gagaid\/\n";
		die;
	}


	## Find masked assembly and repeat gff
	system ("ls $repeatdir\/$gagaid\/\* > inassembly.tmp 2> /dev/null");
	open (Tmpfile, "<", "inassembly.tmp");
	while (<Tmpfile>){
		chomp;
		my $nline = $_;
		next if ($nline !~ /\S+/);		
		if ($nline =~ /repeatMasked\.fasta$/){
			$ingenomemasked = $nline;
		} elsif ($nline =~ /softMasked\.fasta$/){
			$ingenome = $nline;
		} elsif ($nline =~ /repeats\.gff$/){
			$inrepeatgff = $nline;
		}
	}
	close Tmpfile;
	system ("rm inassembly.tmp");

	if ($inrepeatgff !~ /gff/){
		print "ERROR: Skipping $gagaid; $inrepeatgff GFF could not be found, Maybe it is gzipped?, in $repeatdir\/$gagaid\/\n";
		die;
	}

	## Find full gff
	system ("ls $annotdir\/$gagaid\/\*fixed.gff3 > inassembly.tmp 2> /dev/null");
	open (Tmpfile, "<", "inassembly.tmp");
	while (<Tmpfile>){
		chomp;
		my $nline = $_;
		next if ($nline !~ /\S+/);		
		if ($nline =~ /gff/){
			$infullgff = $nline;
		}
	}
	close Tmpfile;
	system ("rm inassembly.tmp");

	if ($infullgff !~ /gff3/){
		print "ERROR: Skipping $gagaid full GFF3 could not be found, Maybe it is gzipped?, in $annotdir\/$gagaid\/\n";
		die;
	}

	# Create script to run filtering script

	$numid++;
	open ("Outscript", ">", "run_$numid\.sh");
	print Outscript "\#\!\/bin\/bash\n\n";
	print Outscript "mkdir -p $gagaid\ncd $gagaid\n\n";

	# Run script
#	print Outscript "perl /home/projects/ku_00039/people/joeviz/GAGA_annotations/Repeat_filtering/final_filtering/run_gene_annotation_filtering_withgff.pl $ingff $inprots $ingenome $ingenomemasked $inrepeatgff $gagaid\n";
#	print Outscript "perl /home/projects/ku_00039/people/joeviz/GAGA_annotations/Repeat_filtering/final_filtering/run_gene_annotation_filtering_onlygff.pl $ingff $inprots $ingenome $ingenomemasked $inrepeatgff $gagaid\n";
#	print Outscript "perl /home/projects/ku_00039/people/joeviz/GAGA_annotations/Repeat_filtering/final_filtering/run_gene_annotation_filtering_onlyfullgff.pl $ingff $inprots $ingenome $ingenomemasked $inrepeatgff $gagaid\n";
	print Outscript "perl /home/projects/ku_00039/people/joeviz/GAGA_annotations/Repeat_filtering/final_filtering/run_gene_annotation_filtering.pl $ingff $inprots $ingenome $ingenomemasked $inrepeatgff $gagaid $infullgff\n";

}
close File;

print "\n";


## Create submission script

open ("Submitscript", ">", "submit_run_filtrep_all\.sh");
print Submitscript "#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn=4
#PBS -l mem=20gb
#PBS -l walltime=1:00:00
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
module load genemark-es/4.62
module load signalp/4.1c
module load funannotate/1.8.3

# This is where the work is done

bash run_\$\{PBS_ARRAYID\}.sh


";



