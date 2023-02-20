#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# Filtering gene families with high variance for CAFE analysis, and also those present in very few species. In addition, filter transposable elements that may remain in the annotations.

# usage: perl filter_genefams_forcafe.pl
# Edit input and output files within the script

## Input ##
my $inputtable = "N0_GeneCounts_nostLFR.tsv"; # Input table containing gene counts that is used in CAFE
#my $inputtable = "N0_GeneCounts_nostLFR_wLept.tsv";

## Output ##
my $output = "N0_GeneCounts_nostLFR_sfilt20_nolow.tsv";
#my $output = "N0_GeneCounts_nostLFR_wLept_sfilt20_nolow.tsv";

## Parameters ##
my $filtdif = "20"; # Maximum difference between minimum and maximum gene count to retain the orthogroup. I used 20 to achieve convergence in CAFE
my $filtminsp = "0"; # Minimum number of species that contain genes in an orthogroup in order to retain it. No need to apply this filter. 


# Get transposable element proteins from the annotations:
# grep -i 'transposable' /home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_functional_annotation_wKEGG/Orthogroups/Orthogroups_functional_annotation_summary.tsv > list_TE_tofilt.txt
# grep -i 'transposase' /home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_functional_annotation_wKEGG/Orthogroups/Orthogroups_functional_annotation_summary.tsv >> list_TE_tofilt.txt
# awk '{print $1}' list_TE_tofilt.txt | sort | uniq > list_TE_tofilt_hogid.txt


## Script

my $TEtofilt = "";
open(File, "<", "list_TE_tofilt_hogid.txt");
while(<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);	
	$TEtofilt .= " $line ";
}
close File;

open (Results, ">", $output);
open (Resultssummary, ">", "$output\_summary.txt");

my @header; my $headercount = 0;
my %sptotalhogcount; my %spfilthogcount;
my $totalhogs = 0; my $filthogs = 0; my $filthogste = 0; my $filthogshighn = 0; my $filthogslown = 0; 

open(File, "<", "$inputtable");
while(<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);	
	$line =~ s/\r//g; # remove this new line from excel
	my @subl = split (/\t/, $line);
	my $desc = shift(@subl);
	my $hog = shift(@subl);

	if ($headercount == 0){ # header line
		@header = @subl;
		print Results "$line\t\n";
		print Resultssummary "$line\t\n";
		$headercount++;
	} else {
		$totalhogs++;
		my $max = 0; my $min = 99999999999; my $spcount = 0;
		my $spheadercount = 0;
		foreach my $value (@subl) {
			next if ($value !~ /\S+/); # Avoid last column without data

			if ($value > $max){
				$max = $value;
			}
			if ($value < $min){
				$min = $value;
			}
			if ($value > 0){
				$spcount++;
			}

			$sptotalhogcount{$header[$spheadercount]} += $value; # Save and sum total number of genes per species
			$spheadercount++;
		}

		my $dif = $max - $min;

		if ($dif > $filtdif){
			$filthogshighn++;
		} elsif ($spcount < $filtminsp) {
			$filthogslown++;
		} elsif ($TEtofilt =~ / $hog /) {
			$filthogste++;
		} else { # Retain the HOG
			$filthogs++;
			print Results "$line\t\n";

			$spheadercount = 0;
			foreach my $value (@subl) {
				next if ($value !~ /\S+/); # Avoid last column without data
				$spfilthogcount{$header[$spheadercount]} += $value; # Save and sum total number of genes per species
				$spheadercount++;
			}
		}

	}


}
close File;

my $totalsphigh = 0; my $totalsplow = 9999999999; my $totalspsum = 0; my $totalspcount = 0;
my $filtsphigh = 0; my $filtsplow = 9999999999; my $filtspsum = 0; my $filtspcount = 0;

print Resultssummary "Total\tTotal\t";
foreach my $sp (@header){
	print Resultssummary "$sptotalhogcount{$sp}\t";
	if ($sptotalhogcount{$sp} > $totalsphigh){
		$totalsphigh = $sptotalhogcount{$sp};
	}
	if ($sptotalhogcount{$sp} < $totalsplow){
		$totalsplow = $sptotalhogcount{$sp};
	}
	$totalspsum += $sptotalhogcount{$sp};
	$totalspcount++;
}
print Resultssummary "\nFilt\tFilt\t";
foreach my $sp (@header){
	print Resultssummary "$spfilthogcount{$sp}\t";
	if ($spfilthogcount{$sp} > $filtsphigh){
		$filtsphigh = $spfilthogcount{$sp};
	}
	if ($spfilthogcount{$sp} < $filtsplow){
		$filtsplow = $spfilthogcount{$sp};
	}
	$filtspsum += $spfilthogcount{$sp};	
	$filtspcount++;
}
print Resultssummary "\n\n\n";

my $totalspave = $totalspsum/$totalspcount;
my $filtspave = $filtspsum/$filtspcount;


print Resultssummary "Species with highest number of genes in all orthogroups (Original without filtering): $totalsphigh\n";
print Resultssummary "Species with lowest number of genes in all orthogroups: $totalsplow\n";
print Resultssummary "Average genes per species in all orthogroups: $totalspave\n\n";
print Resultssummary "Species with highest number of genes in retained orthogroups (after filtering): $filtsphigh\n";
print Resultssummary "Species with lowest number of genes in all orthogroups: $filtsplow\n";
print Resultssummary "Average genes per species in all orthogroups: $filtspave\n\n\n";


print Resultssummary "Total number of orthogroups: $totalhogs\nNumber of orthogroups retained after filtering: $filthogs\n";
print Resultssummary "Orthogroups filtered with high variance ($filtdif as the maximum difference between the maximum and minimum count): $filthogshighn\n";
print Resultssummary "Orthogroups filtered that are present in less than $filtminsp species: $filthogslown\n";
print Resultssummary "Orthogroups filtered annotated as TE proteins: $filthogste\n\n";

# print table with total genes per species, and second line with total genes after filtering
# Then print some summaries such as species with lowest genes, highest n genes, average before/after filtering
# print also number of orthogroups filtered because of TE, high variance, or low species presence

close Results;
close Resultssummary;

