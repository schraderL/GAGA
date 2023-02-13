#!/usr/bin/perl -w
use strict;

if(@ARGV != 2) {
	warn "#Usage: perl $0 <representative> <gff> \n";
	exit;
}

my $file1 = shift;
my $file2 = shift;

my %hash;
open IN, "<$file1" or die "failed to open $file1:$!\n";
while (<IN>) {
    chomp;
    my @t=split /\s+/;
    $hash{$t[1]} = 1;
}
close (IN);

open IN, "<$file2" or die "failed to open $file2:$!\n";
my $tag = 0;
while (<IN>) {
    if(/\sgene\s/) {
        print;
        $tag = 0;
    }elsif(/\smRNA\s/) {
        my $id = $1 if(/ID=([^;]+);/);
        if(exists $hash{$id}) {
            $tag = 1;
            print;
        }else {
            $tag = 0;
        }
    }elsif(/\sCDS\s/) {
        print if($tag);
    }else {
        print if($tag);
    }
}
close (IN);
