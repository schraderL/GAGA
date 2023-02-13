#!/usr/bin/perl -w
use strict;

my $file = shift;

open IN, "<$file" or die "failed to open $file:$!\n";
my $gene_id = '';
my $mRNA_id = '';
while (<IN>) {
    chomp;
    if(/\sgene\s/) {
        $gene_id = $1 if(/ID=([^;]+);/);
    }elsif(/\smRNA\s/) {
	my $id = $1 if(/ID=([^;]+);/);
        if(/Parent/) {
            if(/Parent=(\w+)_i\d+;/) {
		my @t=split /\t/;
		$t[8] = "ID=$id;Parent=$gene_id;";
		$_ = join("\t",@t);
		print STDERR "$gene_id\n";
	    }
        }else {
            $_ .= "Parent=$gene_id;";
        }
#        $gene_id = '';
    }
    print $_,"\n";
}
close (IN);
