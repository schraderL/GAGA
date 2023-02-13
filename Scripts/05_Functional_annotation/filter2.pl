#!/usr/bin/perl -w
use strict;

if(@ARGV != 2) {
	warn "#Usage: perl $0 <representative> <gff> \n";
	exit;
}

my $file1 = shift;
my $file2 = shift;

my %hash;
open IN, "<$file1" or die "failed: $!\n";
while (<IN>) {
	chomp;
	my @t=split /\s+/;
	$hash{$t[0]} = 1;
}
close (IN);

open IN, "<$file2" or die "failed:$!\n";
while (<IN>) {
	if(/\sgene\s/) {
		my $id = $1 if(/ID=([^;]+);/);
		if(exists $hash{$id}) {
			print;
		}else {
			print STDERR $_;
		}
	}else {
		print;
	}
}
close (IN);
