#!/usr/bin/perl -w
use strict;
die "Usage: <genewise reuslt>\n" unless @ARGV == 1;

open IN, $ARGV[0];
while (<IN>) {
	if (/^Bits/) {
		my $line = <IN>;
		my ($score, $gene) = (split /\s+/, $line)[0,1];
		print "$gene\t$score\n";
	}
}
close IN;
