#!/usr/bin/perl
use strict;
use warnings;

#########################################################
#####
#
# Create scripts to submit GeMoMa annotation for GAGA assemblies with RNA-seq
#
#####

# Usage
# perl run_gemoma_annotation.pl GAGA-ID_all.txt


## Input variables

my $idlist = $ARGV[0]; # List with GAGA-IDs to generate the scripts
my $assemblydir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Repeat_masked_assemblies/GAGA_annotations/"; # Folder with unzipped soft-masked genome assemblies classified in folders (as in the directory GAGA_annotations in ERDA). Expected file name to finish with softMasked.fasta
#my $rnadir = "/home/projects/ku_00039/people/joeviz/RNA-seq/GAGA-IDs";
my $rnadir = "/home/projects/ku_00039/people/joeviz/RNA-seq_1June/GAGA-IDs";

my $gemomainputdir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Dataset/Whole_genomes/"; # Folder containing the folders with the annotations used in gemoma

my $gagaidlist = "/home/projects/ku_00039/people/joeviz/GAGA_species_list.txt"; # Table with information about GAGA species, id and names


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
	my $gagaid = $_;
	next if ($gagaid !~ /\S+/);
	$gagaid =~ s/\s+//;

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

	## Find RNA-seq reads
	system ("ls $rnadir\/$gagaid\/\*gz $rnadir\/$gagaid\/\*\/\*gz > inrna.tmp 2> /dev/null");
	my $paired = "";
	my %secondpair;
	my $unpaired = "";
	open (Tmpfile, "<", "inrna.tmp");
	while (<Tmpfile>){
		chomp;
		my $nline = $_;
		next if ($nline !~ /\S+/);		
		if ($nline =~ /(\S+)_1\.fq\.gz/ || $nline =~ /(\S+)_1\.fastq\.gz/ || $nline =~ /(\S+)\_R1\S*\.fastq\.gz/ || $nline =~ /(\S+)\_R1\S*\.fq\.gz/){
			$paired .= "$nline ";
		}
		elsif ($nline =~ /(\S+_)2(\.fq\.gz)/ || $nline =~ /(\S+_)2(\.fastq\.gz)/ || $nline =~ /(\S+\_R)2(\S*\.fastq\.gz)/ || $nline =~ /(\S+\_R)2(\S*\.fq\.gz)/){
			my $firstpair = "$1"."1$2";
			if ($paired =~ /$firstpair /){
				$secondpair{$firstpair} = $nline;
			}
			else {
				print "ERROR: No read pair could be found for $nline. Expected file $firstpair. Skipping this RNA-seq, check the files.\n";
				next;
			}
		}
		elsif ($nline =~ /(\S+)\.fq\.gz/ || $nline =~ /(\S+)\.fastq\.gz/){
			$unpaired .= "$nline ";
		}
	}
	close Tmpfile;
	system ("rm inrna.tmp");

	if ($paired !~ /gz/ && $unpaired !~ /gz/){
		print "ERROR, skipping $gagaid No RNA-seq data could be found in $rnadir\/$gagaid\/\n";
		next;
	}

	# Create script to process and map RNA-seq, and run GeMoMa

	$numid++;
	open ("Outscript", ">", "run_$numid\.sh");
	print Outscript "\#\!\/bin\/bash\n\n";
	print Outscript "mkdir -p $gagaid\ncd $gagaid\n\n";

	# Index genome
	print Outscript "mkdir -p $gagaid\_starindex\n";
	print Outscript "STAR --runThreadN 40 --runMode genomeGenerate --genomeSAindexNbases 13 --genomeDir $gagaid\_starindex --genomeFastaFiles $assembly\n\n";

	# Process RNA-seq

	my $bamfiles = "";


### Paired reads

	my @subpaired = split (/\s/, $paired);
	foreach my $fpair (@subpaired){
		my $spair = $secondpair{$fpair};
		my $rname = "";
		if ($fpair =~ /.*\/(\S+)\.f/){
			$rname = $1;
		} else {
			print "ERROR: Can't find read name in $fpair Skipping...\n";
			next;
		}

		# Verify RNA-seq gzip integrity 

		print Outscript "if gzip -t $fpair; then
    echo '$fpair is OK'
else 
    echo '$fpair is corrupt'
    exit 1
fi\n";
		print Outscript "if gzip -t $spair; then
    echo '$spair is OK'
else 
    echo '$spair is corrupt'
    exit 1
fi\n\n";


		# No trimming

		# Map RNA-seq
		print Outscript "STAR --genomeDir $gagaid\_starindex --runThreadN 40 --readFilesIn $fpair $spair --readFilesCommand zcat --outFileNamePrefix $rname --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 183273429364 --outSAMstrandField intronMotif --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonical --alignIntronMin 20\n";

		print Outscript "samtools index $rname"."Aligned.sortedByCoord.out.bam\n\n";

		$bamfiles .= "ERE.m=$rname"."Aligned.sortedByCoord.out.bam ";

	}


### Unpaired reads

	my @subunpaired = split (/\s/, $unpaired);
	foreach my $unread (@subunpaired){
		my $rname = "";
		if ($unread =~ /.*\/(\S+)\.f/){
			$rname = "$1";
		} else {
			print "ERROR: Can't find read name in $unread Skipping...\n";
			next;
		}

		# Verify RNA-seq gzip integrity 

		print Outscript "if gzip -t $unread; then
    echo '$unread is OK'
else 
    echo '$unread is corrupt'
    exit 1
fi\n\n";


		# No trimming

		# Map RNA-seq
		print Outscript "STAR --genomeDir $gagaid\_starindex --runThreadN 40 --readFilesIn $unread --readFilesCommand zcat --outFileNamePrefix $rname --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 183273429364 --outSAMstrandField intronMotif --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonical --alignIntronMin 20\n";

		print Outscript "samtools index $rname"."Aligned.sortedByCoord.out.bam\n\n";

		$bamfiles .= "ERE.m=$rname"."Aligned.sortedByCoord.out.bam ";

	}


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
	r=MAPPED ERE.s=FR_UNSTRANDED $bamfiles \\
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
#PBS -l nodes=1:ppn=40
#PBS -l mem=180gb
#PBS -l walltime=72:00:00
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



