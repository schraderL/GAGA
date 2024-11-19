#!/usr/bin/perl -w
use strict;

if(@ARGV != 3) {
	warn "#Usage: perl $0 <spec_id1> <spec_id2> <spec1_spec2_anc> \n";
	exit;
}

my $species1 = shift;
my $species2 = shift;
my $spec1_spec2_anc = shift;

`mkdir $species1-$species2` if(!-d "$species1-$species2");

`/run/media/dell/data/public/software/conda/miniconda3/envs/py2/bin/python2.7 /run/media/dell/data/public/software/PhyDiag/PhylDiag-master/src/phylDiag.py ~/GAGA_project/23.OrthoFinder/run_135/05.synteny/genes/genesSTE.$species2.list.bz2 ~/GAGA_project/23.OrthoFinder/run_135/05.synteny/genes/genesSTE.$species1.list.bz2 ~/GAGA_project/23.OrthoFinder/run_135/04.run/output3/ancGenomes/generic-workflow/ancGenome.$spec1_spec2_anc.list.CAR --no-imr > ./$species1-$species2/synteny_blocks.lst`;
