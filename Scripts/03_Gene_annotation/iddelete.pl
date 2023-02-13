#!/usr/bin/perl -w
use strict;

my %hash;

my $file1 = shift;
my $file2 = shift;

open IN,"<$file1" or die "failed to open $file1:$!\n";
while (<IN>) {
	chomp;
	my @t=split /\s+/;
	$hash{$t[0]} = 1;
}
close (IN);

open IN, "<$file2" or die "failed to open $file2:$!\n";
while (<IN>) {
	if(/#/) {
		print;
		next;
	}
	my $id = $1 if(/ID=([^;]+);/ || /Parent=([^;]+);/);
	next if(exists $hash{$id});
	print;
}
close (IN);
