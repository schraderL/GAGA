#!/usr/bin/perl -w
use strict;

my $file1 = shift;
my $file2 = shift;

my %hash;
open IN, "<$file1" or die "failed to open $file1:$!\n";
while (<IN>) {
    chomp;
    my @t=split /\s+/;
    $hash{$t[0]} = $t[1];
}
close (IN);

open IN, "<$file2" or die "failed to open $file2:$!\n";
my $num=0;
my $id = '';
my $new_id = '';
print "##gff-version 3\n";
while (<IN>) {
    next if(/#/);
    chomp;
    my @t=split /\t/;

    if($t[2] eq "gene") {
        $id = $1 if($t[8] =~ /ID=(\S+)/);
        $new_id = $hash{$id};
        $t[8] = "ID=$new_id;";
        $num=0;
        print "###\n";
        print "#Gene: $new_id\n";
    }elsif($t[2] eq "mRNA") {
        $num++;

	my $tmp_id = $1 if($t[8] =~ /ID=([^;]+);/);

        my $new_id2 = "$new_id\_i$num";
        $t[8] = "ID=$new_id2;";
	print STDERR "$tmp_id\t$new_id2\n";
    }elsif($t[2] eq "CDS") {
        my $new_id2 = "$new_id\_i$num";
        $t[8] = "Parent=$new_id2;";
    }else {
        my $new_id2 = "$new_id\_i$num";
        $t[8] = "Parent=$new_id2;";
    }
    if(/StringTie/) {
        $t[5] = ".";
    }
    print join("\t",@t),"\n";
}
close (IN);
