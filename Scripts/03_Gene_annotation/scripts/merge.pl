#!/usr/bin/perl -w
use strict;

open IN, "<$ARGV[0]" or die "failed:$!\n";
my %pen;
while (<IN>) {
    chomp;
	next if(/#/);
    my @t=split /\s+/;
    if(!exists $pen{$t[1]}{score}) {
        $pen{$t[1]}{id} = $t[0]; 
        $pen{$t[1]}{score} = $t[3];
    }elsif($t[3] > $pen{$t[1]}{score} ) {
        $pen{$t[1]}{id} = $t[0]; 
        $pen{$t[1]}{score} = $t[3];
    }
}
close (IN);

open IN, "<$ARGV[1]" or die "failed:$!\n";
while (<IN>) {
    chomp;
	next if(/#/);
    my @t=split /\s+/;
    my $id = $1 if(/ID=(\S+);/ || /Parent=(\S+);/);
    if(exists $pen{$id}) {
        $t[8] .="gene_name=$pen{$id}{id};";
        print join("\t",@t),"\n";
    }else {
#        print STDERR $id,"\n";
    }
}
