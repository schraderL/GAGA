#!/usr/bin/perl
use strict;
use warnings;

#########################################################
#####
#
# Download raw read data from NCBI SRA using a table as input with the sra id and gaga-id
#
#####

# Usage
# perl run_download_sra.pl SRA_data_RNA-seq_tab.tsv 


## Load all required modules for the job when submitting to the queue
#module load tools
#module load ngs
#module load ascp/3.9.6
#module load sratoolkit/2.10.7 
#module load pigz/2.3.4

## Start

open(File, "<", $ARGV[0]);
while(<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	my @subl = split (/\t/, $line);

	next if ($subl[3] !~ /X/); #Skip if not labelled to download 
	my $gagaid = $subl[0];
	my $sra = $subl[8];
	my $outname = $subl[4];

	system ("mkdir -p $gagaid");
	system ("fastq-dump --defline-qual \"+\" --split-e --outdir \. --defline-seq \'\@\$ac-\$si\/\$ri\' $sra");

	my $downloadedfiles = "0";
	system ("ls $sra\* > tmpfiles.txt");
	open (Tmpfile, "<", "tmpfiles.txt");
	while (<Tmpfile>){
		chomp;
		my $line2 = $_;
		next if ($line2 !~ /\S+/);
		if ($line2 =~ /(\S+)\_1.fastq/){
			system ("gzip $1\_1.fastq");
			system ("mv $1\_1.fastq.gz $gagaid\/$outname\_1.fastq.gz");
			$downloadedfiles++;
		} elsif ($line2 =~ /(\S+)\_2.fastq/){
			system ("gzip $1\_2.fastq");
			system ("mv $1\_2.fastq.gz $gagaid\/$outname\_2.fastq.gz");
			$downloadedfiles++;
		} elsif ($line2 =~ /(\S+)\.fastq/){
			system ("gzip $1\.fastq");
			system ("mv $1\.fastq.gz $gagaid\/$outname\.fastq.gz");
			$downloadedfiles++;
		}

	}
	close Tmpfile;
	if ($downloadedfiles == 0){
		print "Warning: No file downloaded for $gagaid $sra\n"
	}

}
close File;


