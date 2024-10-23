#!/usr/bin/perl
use strict;
use warnings;

#########################################################
#####
#
# Create scripts to submit GeMoMa annotation for GAGA assemblies without RNA-seq, using the annotation from the closes relative with RNA-seq
#
#####

# Usage
# perl run_gemoma_annotation.pl GAGA-ID_list_nornabatch_withcloserrelative.txt


## Input variables

my $idlist = $ARGV[0]; # Two columns, GAGA to annotate, and close relative
my $assemblydir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Repeat_masked_assemblies/GAGA_annotations/"; # Folder with unzipped soft-masked genome assemblies classified in folders (as in the directory GAGA_annotations in ERDA)
#my $rnadir = "/home/projects/ku_00039/people/joeviz/RNA-seq/GAGA-IDs";

my $gemomainputdir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Dataset/Whole_genomes/"; # Folder containing the folders with the annotations used in gemoma
my $annotdir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_files_batch3/GAGA_annotations/"; # Directory containing the annotations from gaga genomes annotated with RNA-seq to use in gemoma

my $gagaidlist = "/home/projects/ku_00039/people/joeviz/GAGA_species_list.txt";


## Start

# Read species list and prefix
my %shortname;
open(Gagafile, "<", $gagaidlist);
while(<Gagafile>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	$line =~ s/\n//g; $line =~ s/\r//g;
	my @subl = split (/\t/, $line);
	$shortname{$subl[0]} = $subl[3];
}
close Gagafile;
#

print "Processing GAGA IDs:\n";

my $numid = "0";
open(File, "<", $idlist);
while(<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
#	$line =~ s/\s+//;
	my $gagaid = ""; my $closeid = "";
	if ($line =~ /(\S+)\t(\S+)/){
		$gagaid = $1;
		$closeid = $2;
	}


	## Find fasta
	system ("ls $assemblydir\/$gagaid\/$gagaid\*softMasked.fasta > inassembly.tmp 2> /dev/null");
	my $assembly = "";
	open (Tmpfile, "<", "inassembly.tmp");
	while (<Tmpfile>){
		chomp;
		my $nline = $_;
		next if ($nline !~ /\S+/);		
		if ($nline =~ /fasta/){
			$assembly = $nline;
		}
	}
	close Tmpfile;
	system ("rm inassembly.tmp");

	if ($assembly !~ /fasta/){
		print "ERROR: Skipping $gagaid Assembly could not be found, Maybe it is gzipped?, in $assemblydir\/$gagaid\/\n";
		next;
	}

	## Find close relative GFF3 and fasta
	system ("ls $annotdir\/$closeid\/$closeid\*annotation_representative_repfilt.gff3 > ingff.tmp 2> /dev/null");
	my $gff = "";
	open (Tmpfile, "<", "ingff.tmp");
	while (<Tmpfile>){
		chomp;
		my $nline = $_;
		next if ($nline !~ /\S+/);		
		if ($nline =~ /gff3/){
			$gff = $nline;
		}
	}
	close Tmpfile;
	system ("rm ingff.tmp");

	if ($gff !~ /gff3/){
		print "ERROR: Skipping $gagaid GFF3 could not be found, Maybe it is gzipped?, in $annotdir\/$closeid\/\n";
		next;
	}

	system ("ls $assemblydir\/$closeid\/$closeid\*softMasked.fasta > inassembly.tmp 2> /dev/null");
	my $closeassembly = "";
	open (Tmpfile, "<", "inassembly.tmp");
	while (<Tmpfile>){
		chomp;
		my $nline = $_;
		next if ($nline !~ /\S+/);		
		if ($nline =~ /fasta/){
			$closeassembly = $nline;
		}
	}
	close Tmpfile;
	system ("rm inassembly.tmp");

	if ($closeassembly !~ /fasta/){
		print "ERROR: Skipping $gagaid Assembly could not be found, Maybe it is gzipped?, in $assemblydir\/$closeid\/\n";
		next;
	}




	# Create script to run GeMoMa

	$numid++;
	open ("Outscript", ">", "run_$numid\.sh");
	print Outscript "\#\!\/bin\/bash\n\n";
	print Outscript "mkdir -p $gagaid\ncd $gagaid\n\n";



	# Run GeMoMa
	my $prefixname = "$shortname{$gagaid}\_";

	print Outscript "\nmkdir -p gemoma_output

# Input files
cds_parts_dm=\"$gemomainputdir\/D_mel/ncbigffreformatted/cds-parts.fasta\"
assig_file_dm=\"$gemomainputdir\/D_mel/ncbigffreformatted/assignment.tabular\"

cds_parts_nv=\"$gemomainputdir\/N_vit/ncbigffreformatted/cds-parts.fasta\"
assig_file_nv=\"$gemomainputdir\/N_vit/ncbigffreformatted/assignment.tabular\"

cds_parts_am=\"$gemomainputdir\/A_mel/ncbigffreformatted/cds-parts.fasta\"
assig_file_am=\"$gemomainputdir\/A_mel/ncbigffreformatted/assignment.tabular\"

cds_parts_ob=\"$gemomainputdir\/O_biroi/ncbigffreformatted/cds-parts.fasta\"
assig_file_ob=\"$gemomainputdir\/O_biroi/ncbigffreformatted/assignment.tabular\"

cds_parts_cf=\"$gemomainputdir\/C_flo/ncbigffreformatted/cds-parts.fasta\"
assig_file_cf=\"$gemomainputdir\/C_flo/ncbigffreformatted/assignment.tabular\"

cds_parts_tc=\"$gemomainputdir\/T_cas/ncbigffreformatted/cds-parts.fasta\"
assig_file_tc=\"$gemomainputdir\/T_cas/ncbigffreformatted/assignment.tabular\"

cds_parts_mp=\"$gemomainputdir\/Mpha/cds-parts.fasta\"
assig_file_mp=\"$gemomainputdir\/Mpha/assignment.tabular\"


#Gemoma
java -Xms20G -Xmx170G -jar /home/projects/ku_00039/people/joeviz/programs/GeMoMa-1.7.1/GeMoMa-1.7.1.jar CLI GeMoMaPipeline threads=6 t=$assembly \\
	s=pre-extracted i=obir c=\$cds_parts_ob a=\$assig_file_ob \\
 	s=pre-extracted i=dmel c=\$cds_parts_dm a=\$assig_file_dm \\
	s=pre-extracted i=nvit c=\$cds_parts_nv a=\$assig_file_nv \\
	s=pre-extracted i=amel c=\$cds_parts_am a=\$assig_file_am \\
	s=pre-extracted i=cflo c=\$cds_parts_cf a=\$assig_file_cf \\
	s=pre-extracted i=tcas c=\$cds_parts_tc a=\$assig_file_tc \\
	s=own i=$closeid a=$gff g=$closeassembly \\
	tblastn=false outdir=gemoma_output p=true pc=true AnnotationFinalizer.r=SIMPLE AnnotationFinalizer.p=$prefixname AnnotationFinalizer.n=false \\
	m=/services/tools/mmseqs2/release_12-113e3/bin/
	\n"; 


	print "$gagaid OK\n";
}
close File;

print "\n";


## Create submission script

open ("Submitscript", ">", "submit_run_gemoma_all\.sh");
print Submitscript "#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn=6
#PBS -l mem=180gb
#PBS -l walltime=96:00:00
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
module load gcc
module load star/2.7.2b
module load samtools/1.12
module load ncbi-blast/2.2.31+
module load java/1.8.0
module load mmseqs2/release_12-113e3

# This is where the work is done

bash run_\$\{PBS_ARRAYID\}.sh


";



