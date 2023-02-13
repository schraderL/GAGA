#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# usage: perl get_abcenth_gtf_corrected.pl ABCENTH.gtf genome.fasta

#  It corrects errors in the Abcenth gtf, such as cds longer than the actual scaffold length, and negative starting values

my ($line, $name, $nameout);
my $gff = $ARGV[0];
my $genome = $ARGV[1];


# Reading Protein fasta

my %fasta;
my %flength;
open(File, "<", $genome);
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

foreach my $seq (sort keys %fasta){
	my $length = length ($fasta{$seq});
	$flength{$seq} = "$length";
}


# Reading GTF
open (Results, ">", "ABCENTH_corrected.gtf");

open (GFFfile , "<", $gff); 
while (<GFFfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);

	if ($line =~ /^#/){
		print Results "$line\n";
	}

	my @subline = split (/\t/, $line);

	if ($subline[2] =~ /CDS/){
		my $ini = $subline[3];
		my $end = $subline[4];
		my $scaf = $subline[0];

		my $ok = "1";
		if ($ini < 1){
			$ini = 1;
			$ok++;
		} 
		if ($end > $flength{$scaf}){
			$end = $flength{$scaf};
			$ok++;
		}

		if ($ok == 1){
			print Results "$line\n";
		} else {
			print Results "$subline[0]\t$subline[1]\t$subline[2]\t$ini\t$end\t$subline[5]\t$subline[6]\t$subline[7]\t$subline[8]\n";
			print "Warning: Found erroneus positions in $line\n";
		}


	} else {
		print "Warning: Line in ABCENTH.gtf does not contain CDS: $line\n";
		print Results "$line\n";
	}
}
close GFFfile;
close Results;


### END

