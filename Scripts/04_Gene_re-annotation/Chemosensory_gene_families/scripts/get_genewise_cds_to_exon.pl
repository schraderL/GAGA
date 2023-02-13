#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# Parse genewise GFF3 and convert CDS to exon to find best orfs; Note: it will skip pseudogenes

# usage: perl get_genewise_gtf_corrected.pl all_genewise_predictions.gff genome.fasta

my $prefix = $ARGV[0];
my $output = $ARGV[1];
my $gfffile = "$prefix\.gff3";
my $protfile = "$prefix\.pep.fasta";

my ($name, $line);

my %fasta;
open(File, "<", $protfile);
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	

	if ($line =~ />(\S+)/){
		$name = $1;
	} else {
		$fasta{$name} .= "$line";
	}
}
close File;

my $pseudogenes = "";
my $pseudocount = "0";
my $framecount = "0";

foreach my $prot (sort keys %fasta){
	my $length = length ($fasta{$prot});
	my $seq = uc ($fasta{$prot});

	# Avoid error of substr outside of string; short genes will be filtered out at the end anyway
	next if ($length < 40);

#	my $middleseq = substr ($seq, 5, -10);
	my $middleseq = substr ($seq, 12, -15);	
	if ($middleseq =~ /X/){
		my $lenx = () = $middleseq =~ /X/g;
		if ($lenx > 5){
			# Frameshit error 
			$framecount++;
#			print "Erroneus frameshit in $prot, number of X is $lenx\n$seq\n";
		} else {
			$pseudocount++;
			$pseudogenes .= "$prot ";
#			print "Pseudogene in $prot, number of X is $lenx\n$seq\n";			
		}
	}
}

print "Pseudogenes skipped: $pseudocount\n";
print "Erroneus frameshift found: $framecount\n";

# Write and edit the gff

open (Results, ">", "$output");

open (GFFfile , "<", $gfffile); 
while (<GFFfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);

	if ($line =~ /^#/){
		print Results "$line\n";
		next;
	}
	
	my @subline = split (/\t/, $line);

	if ($subline[2] =~ /CDS/){
		my $genename = "";
		my $rest = "";
		if ($subline[8] =~ /Parent=([^;]+)(\S*)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
			$rest = $2;
		}
		else {die "ERROR in run_OR_classification.pl: It fails detecting Parent ID in $line\n";}

		# Printing pseudogenes as they are
		if ($pseudogenes =~ /$genename /){
			print Results "$line\n";
		} else { # Converting CDS to exon
			print Results "$subline[0]\t$subline[1]\texon\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\t$subline[8]\n";
		}
		
	}
	elsif ($subline[2] =~ /mRNA/){
		my $genename = "";
		my $rest = "";
		if ($subline[8] =~ /ID=([^;]+)\;Parent=[^;]+(\;\S*)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
			$rest = $2;
		}
		else {print "ERROR in run_OR_classification.pl: It fails detecting ID in $line\n";}

		# Printing pseudogenes as they are
		if ($pseudogenes =~ /$genename /){
			print Results "$line\n";
		} else {
			print Results "$line\n";
		}
	}
	elsif ($subline[2] =~ /gene/){
		my $genename = "";
		my $rest = "";
		if ($subline[8] =~ /ID=([^;]+)\;(\S*)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
			$rest = $2;
		}
		else {print "ERROR in run_OR_classification.pl: It fails detecting ID in $line\n";}

		# Printing pseudogenes as they are
		if ($pseudogenes =~ /$genename /){
			print Results "$line\n";
		} else {
			print Results "$line\n";
		}
	}
}
close GFFfile;
close Results;


### END

