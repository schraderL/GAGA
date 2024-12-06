#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# Uncomment line 151 to retrieve a masked alignment

# usage: perl mask_stop_zorro_codon_alignments.pl

## Load computerome2 modules
#module load ngs tools
#module load anaconda3/4.4.0


my $cdsdir = "/home/projects/ku_00039/people/joeviz/orthology_alignments/all_orthogroups/codon_alignments/codon_alignments/";
my $zorrodir = "/home/projects/ku_00039/people/joeviz/orthology_alignments/all_orthogroups/codon_alignments/zorro_codon_alignments/";

my $zorroscore = "4"; # Zorro score used to mask codons

#system ("mkdir -p $outdir");


my ($line, $name);

open (Scores, ">", "Zorro_average_scores.txt");
print Scores "OG\tAverage score\tNumber of sequences\tAverage number of unaligned positions\tNumber of positions check2\tNumber of positions\tNumber of positions > 9\tNumber of positions > 6\tNumber of positions > 4\tNumber of positions < 1\n";

# Reading Protein fasta

#my %protfasta;
system ("ls $cdsdir\/\*\.cds.aln > tmp_protlist.txt");
open(File, "<", "tmp_protlist.txt");
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	

	my $id = "";
	if ($line =~ /.*\/(\S+)\.cds/){
		$id = $1;
	} else {
		die "Can't find species id in $line\n";
	}  

	my %fasta;
	my $ognumseqs = "0"; my $ogalnlength = "0";
	my $ogunalignedlength = "0";
	open(Filef, "<", $line);
	while(<Filef>){
		chomp;
		my $line2 = $_;
		if ($line2 =~ />(\S+)/){
			$name = $1;
			$ognumseqs++;
		} else {
			$fasta{$name} .= "$line2";
			if ($ognumseqs == 1){ # only count aln length for the first sequence 
				$ogalnlength += length($line2);
			}
			my $nogapseq = $line2;
			$nogapseq =~ s/\-//g;
			$ogunalignedlength += length($nogapseq);
		}
	}
	close Filef;
#=h
	open (Results, ">", "$cdsdir\/$id\.cds.nonstop.aln");

	foreach my $gene (keys %fasta){
		my $seq = uc($fasta{$gene});
		my $nstopseq = "";
		my $multseq = "";

		my $lengthcheck = length($seq)/3;
		if ($lengthcheck =~ /\.33/){
			#print Results "$line2"."NN\n";	
			$multseq = "$seq"."NN";
			print "$id $gene is not multiple of 3!\n";	# Print seq id to check
		} elsif ($lengthcheck =~ /\.66/) {
			#print Results "$line2"."N\n";
			$multseq = "$seq"."N";
			print "$id $gene is not multiple of 3!\n";	# Print seq id to check
		} else {
			#print Results "$line2\n";
			$multseq = "$seq";
		}

		# Delete stop codonz
		my $length = length($multseq);
		for (my $n = 0; $n < ($length - 2); $n += 3) {
			my $codon = substr ($multseq, $n, 3);
			if ($codon =~ /TAA/ || $codon =~ /TAG/ || $codon =~ /TGA/ || $codon =~ /\S\SN/){ # also discard codon finished in Ns
#				print Nonstop "NNN";
				$nstopseq .= "NNN";	
#				if ($n >= ($length - 3)){ # Only discard stop if it is the last one
#					$nstopseq .= "";	
#				} else {
#					$nstopseq .= "$codon";	
#				}
			}
			else {
#				print Nonstop "$codon";
				$nstopseq .= "$codon";
			}
		}

		print Results ">$gene\n$nstopseq\n";

	}

	close Results;
#=cut
	# Now mask the alignment based on zorro scores

	my $zorrofile = "$zorrodir\/$id\.cds_zorro.txt";

	my $total = 0; my $pos = 0;	my $posh9 = 0; my $posh6 = 0; my $posh4 = 0; my $posl1 = 0;
	open (Mask, "<", $zorrofile);
	while (<Mask>){
		chomp;
		my $line2 = $_;
		next if ($line2 !~ /\S+/);
		unless ($line2 =~ /Killed/){
			$pos++;
			$total += $line2;	
			if ($line2 >= 9){
				$posh9++;
			}
			if ($line2 >= 6){
				$posh6++;
			}
			if ($line2 >= 6){
				$posh4++;
			}
			if ($line2 < 1){
				$posl1++;
			}
		}
	}	
	close Mask;

	my $media = "";

	if ($pos == 0){
		# Empty zorro file
		print "Warning: $id with $ognumseqs sequences, and alignment length of $ogalnlength; zorro file is empty\n";
		$media = "0";
	} else {
		$media = $total/$pos;
		system ("python3 /home/projects/ku_00039/people/joeviz/orthology_alignments/msa_mask_from_zorro.py $cdsdir\/$id\.cds.nonstop.aln $zorrofile $zorroscore $cdsdir\/$id\.cds.nonstop.zorromasked.aln ");
		#print ("python3 /home/projects/ku_00039/people/joeviz/Suz/ortholog_alignments/msa_mask_from_zorro.py $cdsdir\/$id\.cds.nonstop.aln $zorrofile 5 $cdsdir\/$id\.cds.nonstop.zorromasked.aln \n");

	}

	my $avgnumpos =int ($ogunalignedlength/$ognumseqs + 0.5);
	print Scores "$id\t$media\t$ognumseqs\t$avgnumpos\t$ogalnlength\t$pos\t$posh9\t$posh6\t$posh4\t$posl1\n";


}
close File;
close Scores;
system ("rm tmp_protlist.txt");



