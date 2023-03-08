#!/usr/bin/perl
use strict;

die "\nFunction: make simple name for swissprot/trembl protein ID\n\nUsage: perl simplify_swissprot.pl uniprot_sprot.fasta > uniprot_sprot.fasta.simple\n\n" if(@ARGV==0);

my $uniprot_file = shift;

open(IN, $uniprot_file) || die ("can not open $uniprot_file\n");
$/=">"; <IN>; $/="\n";
while (<IN>) {
	chomp;
	#sp|Q01525|14332_ARATH 14-3-3-like protein GF14 omega OS=Arabidopsis thaliana GN=GRF2 PE=2 SV=2
	my ($id,$desc) = ($2,$3) if(/^(sp|tr)\|([^\|]+)\|(.+)/);	
	$/=">";
	my $seq = <IN>;
	chomp $seq;
	$/="\n";

	print ">$id  $desc\n$seq";
}
close(IN);
