#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);


# usage: perl get_orthogroup_sequences.pl ortologtable(Format 1 column OG, rest are genes) outputdir

# Inputs
my $intable = "$ARGV[0]";
my $outdir = "$ARGV[1]";
my $proteomedir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_GAGA_annotations/"; # Folder containing protein and cds sequences, it will look a folder that is the GAGA-ID (i.e.: GAGA-0001), and use the representative.pep.fasta and representative.cds.fasta files (edit line 24 and 57 with other file extensions)


# Code

system ("mkdir -p $outdir");


my ($line, $name);


# Reading Protein fasta

my %protfasta;
system ("ls $proteomedir\/\*/\*\_representative.pep.fasta > tmp_protlist.txt");
open(File, "<", "tmp_protlist.txt");
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	

	my $gagaid = "";
	if ($line =~ /(GAGA\-\d\d\d\d)/ || $line =~ /(NCBI\-\d\d\d\d)/ || $line =~ /(OUT\-\d\d\d\d)/ || $line =~ /(OUTG\-\d\d\d\d)/){
		$gagaid = $1;
	} else {
		die "Can't find GAGA id in $line\n";
	}  

	open(Filef, "<", $line);
	while(<Filef>){
		chomp;
		my $line2 = $_;
		if ($line2 =~ />(\S+)/){
			$name = $1;
		} else {
			$protfasta{$gagaid}{$name} .= "$line2";
		}
	}
	close Filef;
}
close File;
system ("rm tmp_protlist.txt");


# Reading CDS fasta

my %cdsfasta;
system ("ls $proteomedir\/\*/\*\_representative.cds.fasta > tmp_protlist.txt");
open(File, "<", "tmp_protlist.txt");
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	

	my $gagaid = "";
	if ($line =~ /(GAGA\-\d\d\d\d)/ || $line =~ /(NCBI\-\d\d\d\d)/ || $line =~ /(OUT\-\d\d\d\d)/ || $line =~ /(OUTG\-\d\d\d\d)/){
		$gagaid	= $1;
	} else {
		die "Can't find GAGA id in $line\n";
	}  

	open(Filef, "<", $line);
	while(<Filef>){
		chomp;
		my $line2 = $_;
		if ($line2 =~ />(\S+)/){
			$name = $1;
		} else {
			$cdsfasta{$gagaid}{$name} .= "$line2";
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
					if (exists $protfasta{$gagaid}{$gene}){
						my $seq = $protfasta{$gagaid}{$gene};
						print Resultsp ">$gagaid\_$gene\n$seq\n";
					} else {
						die "Can't find $gagaid $gene from $og in the protein file\n";
					}	

					if (exists $cdsfasta{$gagaid}{$gene}){
						my $seq = uc($cdsfasta{$gagaid}{$gene});
						my $nstopseq = "";
						my $multseq = "";

						# Check that the cds is multiple of 3
						my $lengthcheck = length($seq)/3;
						if ($lengthcheck =~ /\.33/){
							#print Results "$line2"."NN\n";	
							$multseq = "$seq"."NN";
							print "$gagaid $gene in $og is not multiple of 3! Fixed in the output sequences for $og\n";	# Print seq id to check
						} elsif ($lengthcheck =~ /\.66/) {
							#print Results "$line2"."N\n";
							$multseq = "$seq"."N";
							print "$gagaid $gene in $og is not multiple of 3! Fixed in the output sequences for $og\n";	# Print seq id to check
						} else {
							#print Results "$line2\n";
							$multseq = "$seq";
						}


						# Delete stop codon at the end of the cds
						my $length = length($multseq);
						for (my $n = 0; $n < ($length - 2); $n += 3) {
							my $codon = substr ($multseq, $n, 3);
							if ($codon =~ /TAA/ || $codon =~ /TAG/ || $codon =~ /TGA/ || $codon =~ /\S\SN/){ # also discard codon finished in Ns
#								print Nonstop "---";
								if ($n >= ($length - 3)){ # Only discard stop if it is the last one
									$nstopseq .= "";	
								} else {
									$nstopseq .= "$codon";	
								}
							}
							else {
#								print Nonstop "$codon";
								$nstopseq .= "$codon";
							}
						}

						print Resultsc ">$gagaid\_$gene\n$nstopseq\n";
					} else {
						die "Can't find $gagaid $gene from $og in the cds file\n";
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


