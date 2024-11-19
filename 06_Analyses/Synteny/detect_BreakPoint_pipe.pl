#!/usr/bin/perl -w
use strict;

if(@ARGV != 3) {
	warn "#Usage: perl $0 <spec_id1> <spec_id2> <spec1_spec2_anc> \n";
	exit;
}

my $species1 = shift;
my $species2 = shift;
my $spec1_spec2_anc = shift;

`sort  -k2,2 -k3n ./$species1-$species2/synteny_blocks.lst > ./$species1-$species2/synteny_blocks.lst.sort`;
`perl /run/media/dell/data/User/xiongzj/GAGA_project/23.OrthoFinder/run_135/bin/sbs_filterExtremities.pl ./$species1-$species2/synteny_blocks.lst.sort ~/GAGA_project/23.OrthoFinder/run_135/05.synteny/Length/genesSTE.$species2.list.bz2.len ~/GAGA_project/23.OrthoFinder/run_135/05.synteny/Length/genesSTE.$species1.list.bz2.len > ./$species1-$species2/synteny_blocks.lst.bra`;
