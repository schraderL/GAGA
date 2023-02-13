#!/usr/bin/perl -w
use strict;

my $file = shift;

open IN, "<$file" or die "failed to open $file:$!\n";
my$CDS_i = 0;
while (<IN>) {
    if(/\s+CDS\s+/) {
        chomp;
        $CDS_i++;
        my $CDS="CDS$CDS_i";

        my @t=split /\t/;
        my $id = $1 if($t[8] =~ /Parent=([^;]+);/);

	if($t[8]=~/CDS\d+/) {
	    print $_,"\n";
	}elsif($t[8]=~/^ID/) {
            s#$id#$id.$CDS#;
            print $_,"\n";
        }else {
            $t[8] = "ID=$id.$CDS;Parent=$id;";
            print join("\t",@t),"\n";
        }
    }else {
        print;
        $CDS_i=0;
    }
}
close (IN);
