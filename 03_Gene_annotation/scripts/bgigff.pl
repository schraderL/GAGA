#!/usr/bin/perl -w
use strict;

while (<>) {
	next if(/^#/);
	chomp;
	my @t=split /\t/;
	if($t[2] eq "mRNA" || $t[2] eq "CDS") {
		if(/;$/) {
			print $_,"\n";
		}else {
			print "$_;\n";
		}
	}
}
