#!/usr/bin/perl -w
use strict;

if(@ARGV != 2) {
	warn "#Usage: perl $0 <gff> <prefix>\n";
	exit;
}

my $file = shift;
my $prefix = shift;

my %length;
readFile2($file,\%length);

open IN, "<$file" or die "failed to open $file:$!\n";
while (<IN>) {
	chomp;
	if(/\smRNA\s/) {
		my $type = "GeMoMa";
		if(/GeMoMa/) {
			$type = "GeMoMa";
		}elsif(/StringTie/) {
			$type = "StringTie";
		}elsif(/AUGUSTUS/){
			$type = "AUGUSTUS";
		}else{
			$type = "JoelAdd";
		}
		my $tag="normal";
		my @t=split /\t/;
		my $id = $1 if(/ID=([^;]+);/);
		my $gene_id = $1 if(/Parent=([^;]+);/);
		if($id =~ /split/i) {
			$tag = "split";
		}elsif($id =~ /^$prefix\_g/) {
			$tag = "normal";
		}else {
			$tag = "curated";
		}
		if(exists $length{$id}) {
			print "$gene_id\t$id\t$tag\t$length{$id}\t$type\n";
		}else {
			print STDERR  "$id\n";
		}
	}
}

#############
#######
sub readFile2 {
	my ($file,$p) = @_;

	open IN, "<$file" or die "failed to open $file:$!\n";
	while (<IN>) {
		chomp;
		if(/\sCDS\s/) {
			my @t=split /\t/;
			my $id = $1 if(/Parent=([^;]+);/);
			$p->{$id} += ($t[4]-$t[3]+1);
		}
	}
	close (IN);
}
