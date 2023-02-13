#!/usr/bin/perl -w
use strict;

if(@ARGV != 2) {
	warn "#Usage: perl $0 <Overlap.lst> <prefix>\n";
	exit;
}

my $file = shift;
my $prefix = shift;

my $overlap_size = 0;
my $OR_id   = '';
my $gene_id = '';
my $mRNA_id = '';
my $OR_id_for = '';
my $out_id = '';
open IN, "<$file" or die "failed to open $file:$!\n";
<IN>;
while (<IN>) {
	my @t=split /\t/;
	if(/(\S+)\s+($prefix\_g\d+)/) {

		$OR_id   = $1;
		$gene_id = $2;
		$mRNA_id = $t[1];
	}
	if($OR_id eq $OR_id_for) {
		if($t[10] >= $overlap_size) {
			$overlap_size = $t[10];
			$out_id = $mRNA_id;
		}
	
	}elsif($OR_id ne $OR_id_for) {
		$overlap_size = 0;
		print "$OR_id_for\t$out_id\n";
		$out_id = $mRNA_id;

	}

	$OR_id_for = $OR_id;
}
close (IN);

print "$OR_id_for\t$out_id\n";
