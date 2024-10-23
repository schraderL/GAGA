#!/usr/bin/perl -w
use strict;

my $file = shift;

open IN, "<$file" or die "failed to open $file:$!\n";
while (<IN>) {
	chomp;
	my $gaga_id = $_;
	`ls $gaga_id/*/*.gff3 | grep -v CombinedEvidences > tmp`;
	open F, "<tmp" or die "failed: tmp\n";
	while (<F>) {
		my $count = 0;

		chomp;
		my $gff_file = $_;
		my @t=split /\//;
		open F2, "<$gff_file" or die "$!\n";
		while (<F2>) {
			$count++ if(/\sgene\s/);
		}
		close (F2);
		print "$t[0]\t$t[1]\t$count\n";
	}
	close (F);
}
close (IN);
