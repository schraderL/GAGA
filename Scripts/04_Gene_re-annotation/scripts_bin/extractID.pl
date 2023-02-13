#!/usr/bin/perl -w
use strict;

if(@ARGV != 2) {
	warn "#Usage: perl $0 <annotpipeline.gff3> <prefix> \n";
	exit;
}

my $file = shift;
my $prefix = shift;
open IN, "<$file" or die "failed to open $file:$!\n";
open OUT1, ">$prefix.mRNA.lst" or die "failed:$!\n";
open OUT2, ">$prefix.gene.lst" or die "failed:$!\n";
while (<IN>) {
	if(/\smRNA\s/)	{
		if(/($prefix\_g\d+_i\d+)/) {
			print OUT1 "$1\t";
		}else {
			print OUT1 "NA\t";
		}
		if(/ID=([^;]+);/) {
			print OUT1 "$1\n";
		}

		if(/($prefix\_g\d+)/) {
			print OUT2 "$1\t";
		}else {
			print OUT2 "NA\t";
		}
		if(/ID=(\S+)_$prefix/) {
			print OUT2 "$1\n";
		}else {
			print OUT2 "NA\n";
		}
	}
}
close (IN);
