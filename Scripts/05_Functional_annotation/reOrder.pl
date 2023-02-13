#!/usr/bin/perl -w
use strict;

my $file = shift;

my %hash;
open IN, "<$file" or die "failed to open $file:$!\n";
while (<IN>) {
	chomp;
	my @t=split /\t/;
	if(/\sgene\s/) {
		my $id = $1 if(/ID=([^;]+);/);
	}elsif(/\smRNA\s/) {

	}elsif(/\sCDS\s/) {

	}
}
close (IN);
