#!/usr/bin/perl -w
use strict;

if(@ARGV !=2) {
	warn "#Usage: perl $0 <gene_list> <cluster> \n";
	exit;
}

my $file1 = shift;
my $file2 = shift;

my %hash;
open IN, "<$file1" or die "failed to open $file1:$!\n";
while (<IN>) {
	chomp;
	$hash{$_} = 1;
}
close (IN);

open IN, "<$file2" or die "failed to open $file2:$!\n";
while (<IN>) {
	chomp;
	my @t = split /\s+/;
	shift @t;
	shift @t;
	my @out;
	for my $i(0..$#t) {
		if(exists $hash{$t[$i]}) {
			my $out = $t[$i];
			push @out,$out;
		}
	}
	if (scalar @out >=1) {
		print join("\t",@out),"\n";
	}else {
		print $t[0],"\n";
	}
}
close (IN);

