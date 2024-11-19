#!/usr/bin/perl -w
use strict;

my $file1 = shift;
my $file2 = shift;

my %sizes;
open IN, "<$file1" or die "failed to open $file1:$!\n";
while (<IN>) {
    chomp;
    my @t=split /\t/,$_;
    $sizes{$t[0]} = $t[1];
}
close (IN);

open IN, "<$file2" or die "failed to open $file2:$!\n";
while (<IN>) {
    my @t=split /\s+/,$_;
    if($t[6] eq "+") {
        print "$t[0]\t$t[1]\t$t[2]\t$t[3]\t$t[4]\t$t[5]\t",$t[2]-$t[1]+1,"\t",$t[5]-$t[4]+1,"\n";
    }elsif($t[6] eq "-") {
        ($t[4],$t[5]) = (($sizes{$t[3]}-$t[5]+1),($sizes{$t[3]}-$t[4]+1));

        $t[6] = "+";
        print "$t[0]\t$t[1]\t$t[2]\t$t[3]\t$t[4]\t$t[5]\t",$t[2]-$t[1]+1,"\t",$t[5]-$t[4]+1,"\n";        
    }
}
close (IN);
