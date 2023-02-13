#!/usr/bin/perl -w
use strict;

my $file = shift;
open IN, "<$file" or die "failed:$!\n";
my $parent_id = "";
while (<IN>) {
    chomp;
    my @t=split /\t/;
    my $id;
#    my $parent_id;
    if ($t[2] eq "mRNA") {
        $id = $1 if(/ID=([^;]+);/);
#        $parent_id = $1 if(/Parent=([^;]+);/);
        $t[8] = "ID=$id;gene_name=$id;";
        print join("\t",@t),"\n";
    }elsif($t[2] eq "CDS") {
        $id = $1 if(/Parent=([^;]+);/);
        $t[8] = "Parent=$id;gene_name=$id;";
        print join("\t",@t),"\n";
    }
}
close (IN);
