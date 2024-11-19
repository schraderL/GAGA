#!/usr/bin/perl -w
use strict;

if(@ARGV != 1) {
    warn "#Usage: perl $0 <net> \n";
    exit;
}

my $file = shift;

open IN, "<$file" or die "failed to open $file:$!\n";
my $query_id="";
while (<IN>) {
    chomp;
    if(/net\s/) {
        my @t=split /\s+/;
        $query_id = $t[1];
    }elsif(/fill\s/) {
        s#^\s+##;
        my @t=split /\s+/;
        print "$query_id\t$t[1]\t",$t[1]+$t[2]-1,"\t$t[3]\t$t[5]\t",$t[5]+$t[6]-1,"\t$t[4]\t$t[2]\t$t[6]\n";
    }
}
close (IN);
