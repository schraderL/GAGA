#!/usr/bin/perl -w
use strict;

if (@ARGV != 1) {
	warn "#Usage: perl $0 <Len> \n";
	exit;
}

my $file2 = shift;

my %id;
open IN, "<$file2" or die "failed to open $file2:$!\n";
while (<IN>) {
    chomp;
    my @t=split/\s+/;
    $id{$t[0]}{$t[1]}{len} = $t[3];
    if($t[2] eq "split") {
	$id{$t[0]}{$t[1]}{tag} = "split";
    }else{
        $id{$t[0]}{$t[1]}{tag} = $t[4];
    }
}
close (IN);

foreach my $gene_id (keys %id) {
    my $long_len = 0;
    my $out_id = '';
    my $out_tag = 0;
    foreach my $mrna_id (keys %{$id{$gene_id}}) {
        my $len = $id{$gene_id}{$mrna_id}{len};
        my $tag = $id{$gene_id}{$mrna_id}{tag};
	if($tag eq "split") {
		print "$gene_id\t$mrna_id\t$len\tsplit\n";
		next;
	}

        #    $tag = $tag eq "normal" ? 0 : 1;
        if($tag eq "JoelAdd") {
            $tag = 5;
        }elsif($tag eq "GeMoMa") {
            $tag = 4;
        }else{
	    $tag = 3;
	}

        if($tag > $out_tag) {
            $long_len = $len;
            $out_id = $mrna_id;
            $out_tag = $tag;
        }elsif($tag == $out_tag && $len > $long_len) {
            $long_len = $len;
            $out_id = $mrna_id;
            $out_tag = $tag;
        }
    }
    print "$gene_id\t$out_id\t$long_len\t";
    my $output = $out_tag == 5 ? "curated" : "normal";
    print "$output\n";
}
