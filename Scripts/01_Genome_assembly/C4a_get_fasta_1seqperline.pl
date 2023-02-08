#!/usr/bin/perl
use strict;
use warnings;

# Create a fasta with 1 sequence per line required for ntJoin

#usage: perl get_fasta_1seqperline.pl infasta

my ($name, $line, $nameout, $contig, $frame, $name2);
my (%nrfa, %inicio, %fin, %fastacutted, %fastanames);
my $skip = 0;
my $ids = "";

my $file= "$ARGV[0]";
my %fasta;
open (Fasta , "<", $file);
while (<Fasta>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ /^>(\S+)/) {
		$name = $1;
	} else {
		$fasta{$name} .= $line;
	}
}
close Fasta;

open (Rest, ">", "$ARGV[0]\_1seqperline.fasta");

foreach my $key (sort keys %fasta){
	print Rest ">$key\n$fasta{$key}\n";
}
#close Res;
close Rest;


