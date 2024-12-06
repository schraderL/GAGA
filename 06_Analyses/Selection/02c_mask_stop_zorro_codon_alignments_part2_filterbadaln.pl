#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# Part 2, use the obtained table with Zorro scores and filter out all bad quality alignments. It copies the good-quality MSAs to a folder to continue the analyses. 

# usage: perl mask_stop_zorro_codon_alignments.pl

## Load computerome2 modules
#module load ngs tools
#module load anaconda3/4.4.0


my $cdsdir = "/home/projects/ku_00039/people/joeviz/orthology_alignments/all_orthogroups/codon_alignments/codon_alignments/";
my $zorrofile = "Zorro_average_scores.txt"; # Generated in the first script

my $outdir = "codon_alignments_qualfiltered";
system ("mkdir -p $outdir");


my ($line, $name);

open (Scores, ">", "Zorro_average_scores_qualfiltered_good.txt");
print Scores "OG\tAverage score\tNumber of sequences\tAverage number of unaligned positions\tNumber of positions check2\tNumber of positions\tNumber of positions > 9\tNumber of positions > 6\tNumber of positions > 4\tNumber of positions < 1\n";

open (Scoresbad, ">", "Zorro_average_scores_qualfiltered_bad.txt");
print Scoresbad "OG\tAverage score\tNumber of sequences\tAverage number of unaligned positions\tNumber of positions check2\tNumber of positions\tNumber of positions > 9\tNumber of positions > 6\tNumber of positions > 4\tNumber of positions < 1\n";


# Reading Zorro table

open(File, "<", $zorrofile);
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	next if ($line =~ /Average score/); # Skip first line

	my @subl = split (/\t/, $line);

	my $id = $subl[0];

	# Quality filter

	if ($subl[2] < 17){ # Skip alignments with less than 17 sequences
		print Scoresbad "$line\n";
		next;
	}
=head
	if ($subl[1] < 4){ # Skip alignments with averaga zorro score lower than 4
		print Scoresbad "$line\n";
		next;
	}
=cut
=head
	my $filtalnlen = $subl[4]*0.1;
	if ($subl[8] < $filtalnlen){ # Skip alignments with less than 10% of the total alignment with score below 4
		print Scoresbad "$line\n";
		next;
	}
=cut	
	my $filtlen = $subl[3]*0.4;
	if ($subl[8] < $filtlen){ # Skip alignments with less than 40% of original sequence avergae length with score below 4
		$filtlen = $subl[3]*0.3; # Additional to keep some more orthogroups with more than 80 sequences, and alignment above 4 in more than 30% of the original sequence length
		if ($subl[2] < 81 || $subl[8] < $filtlen){
			$filtlen = $subl[3]*0.20; # Additional to keep some more orthogroups with more than 130 sequences, and alignment above 4 in more than 25% of the original sequence length
			if ($subl[2] < 130 || $subl[8] < $filtlen){
				print Scoresbad "$line\n";
				next;
			}
		}
	}
	# Save good quality alignments

	print Scores "$line\n";
	system ("cp $cdsdir\/$id\.cds\.nonstop.aln $outdir\/");


}
close File;
close Scores;
close Scoresbad;



