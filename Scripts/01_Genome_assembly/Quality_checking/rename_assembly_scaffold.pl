#!/usr/bin/perl
use strict;
use warnings;

# Remove scaffolds sorted by length

#usage: perl get_fasta_lengthstats.pl infasta

my ($name, $line, $nameout, $contig, $frame, $name2);
my (%nrfa, %inicio, %fin, %fastacutted, %fastanames);
my $skip = 0;
my $ids = "";

my $out = "";

#if ($ARGV[0] =~ /(\S+)(GAGA-\d\d\d\d)(\S+)(\S+).fasta/){
#	$out = "$1"."$2"."final_$3";
#}
if ($ARGV[0] =~ /(\S+).fasta/){
	$out = $1;
} else {
	$out = $ARGV[0];
}


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
		$fasta{$name} .= uc($line); # Change all sequences to upper case nucleotides
	}
}
close Fasta;

#open (Res, ">", "$ARGV[0]\_scaflonger100Kb.fasta");
open (Rest, ">", "$out\_length.txt");
open (Resren, ">", "$out\_final.fasta");

foreach my $key (sort keys %fasta){
	my $length = length ($fasta{$key});
	print Rest "$key\t$length\n";
	if ($length >= 100000){
#		print Res ">$key\n$fasta{$key}\n";
	}
}
#close Res;
close Rest;

system ("sort -n -r -k 2 $out\_length.txt > $out\_length_sorted.txt");

my $n = 1;

open (File , "<", "$out\_length_sorted.txt");
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my @subl = split (/\t/, $line);
	print Resren ">Scaffold$n\n$fasta{$subl[0]}\n";
	$n++;
}
close Fasta;




