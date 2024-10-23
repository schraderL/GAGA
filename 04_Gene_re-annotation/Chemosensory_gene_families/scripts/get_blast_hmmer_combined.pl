#!/usr/bin/perl
use strict;
use warnings;

# Usage: 	perl Scripts/get_blast_hmmer_combined.pl outname\_allsearches_list.txt outname

my %mdom;

open (File, "<", $ARGV[0]);
my %genes;
while (<File>){
	chomp;
	my $line = $_;
	my @subline = split(/\s/, $line);
	my $gene = "";
	if ($subline[0] =~ /(\S+)\_split\d/){
		$gene = $1;
	} else {
		$gene = $subline[0];
	}
	push (@{$genes{$gene}}, $line);

	if ($subline[4] =~ /hmmer/){
		$mdom{$gene} = $subline[5];
	}

}
close File;
open (Results, ">", "$ARGV[1]\_combinedsearches_list.txt");

foreach my $gene (sort keys %genes){
	my $fragcomp = "-";	
	my (@ini, @fin);
	$ini[0] = "999999999999999999";
	$fin[0] = "1"; 
	my $method = "";
	foreach my $list (@{$genes{$gene}}){
		my @subline = split(/\s/, $list);

		if ($subline[4] =~ /blastp/){
			unless ($fragcomp =~ /complete/){
				$fragcomp = $subline[1];
			}
			unless ($method =~ /blastp/){
				$method .= "blastp";
			}
		} elsif ($subline[4] =~ /hmmer/){
			unless ($method =~ /hmmer/){
				$method .= "hmmer";
			}
		}

		my $lini = int($subline[2]); # New lines to edit the start and end position to allow overlapping positions
		my $lend = int($subline[3]);

		my $n = 0;
		my $extrahit = 0;
		foreach my $i (@ini){
			my $f = $fin[$n];
#			my $f2 = $f + 10;
#			my $i2 = $i - 10;
			my $f2 = $f;
			my $i2 = $i;
				if ($i >= $lini && $f <= $lend){
				$ini[$n] = $lini;
				$fin[$n] = $lend;
			}
			elsif ($i <= $lini && $f < $lend && $f2 >= $lini){
				$fin[$n] = $lend;
			}
			elsif ($i > $lini && $f >= $lend && $i2 <= $lend){
				$ini[$n] = $lini;
			}
			elsif ($f2 < $lini) {
				$extrahit++;
				my $dif = int(($lini - $f2 - 10)/2);
				if ($dif > 0 && $dif < 50){
					$lini = $lini - $dif;
					$fin[$n] = $f2 + $dif;
				}
			}
			elsif ($i2 > $lend) {
				$extrahit++;
				my $dif = int(($i2 - $lend - 10)/2);
				if ($dif > 0 && $dif < 50){
					$lend = $lend + $dif;
					$ini[$n] = $i2 - $dif;
				}
			}
			$n++;
		}
		if ($extrahit >= $n) {
				$ini[$n] = $lini;
				$fin[$n] = $lend;				
		}

	}

	my $hits = scalar(@ini);
	if ($hits == 1){
		if ($method =~ /hmmer/){
			print Results "$gene $fragcomp $ini[0] $fin[0] $method $mdom{$gene}\n";
		} else {
			print Results "$gene $fragcomp $ini[0] $fin[0] $method 0\n";
		}	
	}
	else {
		my $n = 0;
		foreach my $i (@ini){
			my $f = $fin[$n];
			my $nn = $n+1;
			if ($method =~ /hmmer/){
				my $numsplit = scalar (@ini);
				if ($numsplit >= $mdom{$gene}){
					print Results "$gene\_split$nn $fragcomp $ini[$n] $fin[$n] $method 1\n";
				} else {
					print Results "$gene\_split$nn $fragcomp $ini[$n] $fin[$n] $method $mdom{$gene}\n";
				}
			} else {
				print Results "$gene\_split$nn $fragcomp $ini[$n] $fin[$n] $method 0\n";
			}				
			$n++;
		}
	}

}
close Results;







