#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# Script to get the sequences from the alignments

# usage: perl get_orthogroup_sequences_fromaln.pl ortologtable(Format 1 column OG, rest are genes) dir_with_aln extension outputdir
# perl get_orthogroup_sequences_fromaln.pl /home/projects/ku_00039/people/joeviz/orthofinder/run_allGAGA_final_annotations/script_orthogroups_subset/species_subset_80percsp_orthogroups_multicopy_genes.tsv all_orthogroups/orthogroups_seqs/ pep.fasta test_subset_fromaln/


my $intable = "$ARGV[0]";
my $alndir = "$ARGV[1]";
my $extension = "$ARGV[2]"; # Ej: .pep.fasta
my $outdir = "$ARGV[3]";

system ("mkdir -p $outdir");


my ($line, $name);


# Reading fasta

my %protfasta;
system ("ls $alndir\/\*$extension > tmp_protlist.txt");
open(File, "<", "tmp_protlist.txt");
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	

	my $filenodir = "";
	if ($line =~ /.*\/(\S+)/){
		$filenodir = $1;
	} else {
		$filenodir = "$line";
	}

	my $hogid = "";
	if ($filenodir =~ /([^\.]+)/){
		$hogid = $1;
	} else {
		die "Can't find hog id (expected HOGXXXX.extension) in $line\n";
	}  

	open(Filef, "<", $line);
	while(<Filef>){
		chomp;
		my $line2 = $_;
		if ($line2 =~ />(\S+)/){
			$name = $1;
		} else {
			$protfasta{$hogid}{$name} .= "$line2";
		}
	}
	close Filef;
}
close File;
system ("rm tmp_protlist.txt");



# Read orthology table

my $headerl = "0";
my @colheader;

open(File, "<", $intable);
while(<File>){
	chomp;
	$line = $_;
	$line =~ s/\r//g;
	next if ($line !~ /\S+/);
	my @subl = split (/\t/, $line);

	if ($headerl == 0){
		foreach my $col (@subl){
			push (@colheader, $col);
		}

	} else {
		my $og = "";
		if ($subl[0] =~ /N\d+\.(\S+)/){
			$og = $1;
		} else {
			$og = $subl[0];
		}
			

		open (Resultsp, ">", "$outdir\/$og\.pep.fasta");
		open (Resultsc, ">", "$outdir\/$og\.cds.fasta");

		my $i = "0";
		foreach my $col (@subl){
			if ($i == 0){
				$i++;
				next;
			} else {
				my $gagaid = $colheader[$i];
				my @subgene = (split /\, /, $col);
				foreach my $gene (@subgene){
					my $genename = "$gagaid\_$gene";
					if (exists $protfasta{$og}{$genename}){
						my $seq = $protfasta{$og}{$genename};
						print Resultsp ">$genename\n$seq\n";
					} else {
						die "Can't find $gagaid $gene from $og in the protein file\n";
					}	
				}

				$i++;
			}

		}

		close Resultsp;
		close Resultsc;
	}

	$headerl++;
}
close File;


