#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# usage: perl get_genewise_gtf_corrected.pl all_genewise_predictions.gff genome.fasta

#  It filters bad predictions in gene wise
# It also fiters duplicated cds, it keeps the ones with best score

my ($line, $name, $nameout);
my $gff = $ARGV[0];
my $genome = $ARGV[1];

my $mincds = "130"; # Length to filter short cds (in nt) # Only the inital exon
my $cdsdist = "260"; # Min distance between cds to filter short cds

# Reading Protein fasta

my %fasta;
my %flength;
open(File, "<", $genome);
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	

	if ($line =~ />(\S+)/){
		$name = $1;
	} else {
		$fasta{$name} .= "$line";
	}
}
close File;

foreach my $seq (sort keys %fasta){
	my $length = length ($fasta{$seq});
	$flength{$seq} = "$length";
}


my %gffcds;
my %gffscafcds;

# Reading GTF
open (Results, ">", "all_genewise_predictions_corrected.gtf");

open (GFFfile , "<", $gff); 
while (<GFFfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);

	if ($line =~ /^#/){
		print Results "$line\n";
	}

	my @subline = split (/\t/, $line);

	if ($subline[2] =~ /CDS/){

		my $genename = "";
		if ($subline[8] =~ /(transcript_id|Transcript).([^;]+)/){  # Reading GTF formats.
			$genename = $2;

		}
		else {die "ERROR in get_genewise_gtf_corrected.pl: It does not recognize Parent ID in the gff in: $line\n";}


		if ($subline[6] =~ /\+/){ # Forward strand
			#OK
		} elsif ($subline[6] =~ /\-/){
			#OK
		} else {die "ERROR in get_genewise_gtf_corrected.pl: No forward/reverse in: $line\n";}


		if ($subline[7] =~ /\d+/){ # Forward strand
			#OK
		} else {die "ERROR in get_genewise_gtf_corrected.pl: No frame in CDS in: $line\n";}



		# Correcting positions lower or higher than scaffold length
		my $ini = $subline[3];
		my $end = $subline[4];

		if ($subline[4] > $subline[3]){
			#OK
		} else { # reversed positions
			$ini = $subline[4];
			$end = $subline[3];			
		}

		my $scaf = $subline[0];
		my $dif = $end - $ini;

		my $ok = "1";
		if ($ini < 1){
			$ini = 1;
			$ok++;
		} 
		if ($end > $flength{$scaf}){
			$end = $flength{$scaf};
			$ok++;
		}

		if ($ok == 1){
			#print Results "$line\n";
		} else {
			#print Results "$subline[0]\t$subline[1]\t$subline[2]\t$ini\t$end\t$subline[5]\t$subline[6]\t$subline[7]\t$subline[8]\n";
			print "Warning: Found erroneus positions in $line\n";
		}

		next if ($subline[5] < 20); ## Checking if filtering by this score helps 

		my $complexname = "$genename\_coordinate\_$subline[0]\_$ini\_$end";

		if ($dif > 21){ # Avoid saving cds really short
#		 	push ( @{$gffcds{$genename}{$ini}}, join ("\t", $subline[0],$subline[1],$subline[2],$ini,$end,$subline[5],$subline[6],$subline[7],$subline[8]));
#			push ( @{$gffscafcds{$subline[0]}{$genename}}, join ("\t" , $ini, $end, $subline[5])) ;

#			print Results "$subline[0]\t$subline[1]\t$subline[2]\t$ini\t$end\t$subline[5]\t$subline[6]\t$subline[7]\t$subline[8]\n"; ### Testing first prints

			###### Filter duplicated CDS

			my $repeated = "0";

			foreach my $prevcds (sort keys %gffcds){
				my @prevsubl = split (/\t/, $gffcds{$prevcds});
				my $pini = $prevsubl[3];
				my $pend = $prevsubl[4];
				my $psca = $prevsubl[0];
				my $psco = $prevsubl[5];

				if ($psca =~ /^$subline[0]$/){
					if ($ini >= $pini && $ini < $pend){ # Overlapping CDS
						if ($psco >= $subline[5]){
							#OK - keep previous
							$repeated++;
						} else {
							delete $gffcds{$prevcds};
	#						$gffcds{$complexname} = "$subline[0]\t$subline[1]\t$subline[2]\t$ini\t$end\t$subline[5]\t$subline[6]\t$subline[7]\t$subline[8]";
						}
					} elsif ($pini >= $ini && $pini < $end){ # Overlapping CDS
						if ($psco >= $subline[5]){
							#OK - keep previous
							$repeated++;
						} else {
							delete $gffcds{$prevcds};
	#						$gffcds{$complexname} = "$subline[0]\t$subline[1]\t$subline[2]\t$ini\t$end\t$subline[5]\t$subline[6]\t$subline[7]\t$subline[8]";
						}
					}
				}
			}
			
			if ($repeated == 0){
				$gffcds{$complexname} = "$subline[0]\t$subline[1]\t$subline[2]\t$ini\t$end\t$subline[5]\t$subline[6]\t$subline[7]\t$subline[8]";
			}

		}


	} else {
		#print "Warning: Line in input gtf does not contain CDS: $line\n";
		#print Results "$line\n";
		next; # Skipping intron lines
	}
}
close GFFfile;


foreach my $complexname (sort keys %gffcds){

	my $genename = "";
	if ($complexname =~ /(\S+)\_coordinate/){
		$genename = $1;
	} else {
		die "Can't find gene name after coordinate in $complexname\nGFF line $gffcds{$complexname}\n\n";
	}

	my @subl = split (/\t/, $gffcds{$complexname});

#	print Results "$gffcds{$complexname}\n";

#	push ( @{$gffscafcds{$genename}}, $gffcds{$complexname}) ;

	$gffscafcds{$genename}{$subl[6]}{$subl[3]} = $gffcds{$complexname};

}


# Now process the gff cds hits, and discard lost cds in different chain, short cds more far than 500bp AND shorter than 40aa = 120 nt; 


foreach my $gene (sort keys %gffscafcds){
	foreach my $strand (sort keys %{$gffscafcds{$gene}}){
		if ($strand =~ /\+/){

			# Check and discard if exists cds from the same gene in the other strand
			my $opstrand = "-";
			if (exists $gffscafcds{$gene}{$opstrand}){
				#die "Found opposite strands for gene $gene\n";
				my $cdcount = keys %{$gffscafcds{$gene}{$strand}};
				my $opcdcount = keys %{$gffscafcds{$gene}{$opstrand}};
				if ($cdcount >= $opcdcount){
					delete $gffscafcds{$gene}{$opstrand};
				} else{
					delete $gffscafcds{$gene}{$strand};
					next;
				}
			}
			#	

			my $prevrow = "";
			my $firstrow = "0";

			foreach my $ini (sort { $a <=> $b } keys %{$gffscafcds{$gene}{$strand}}){
				my $cdsline = $gffscafcds{$gene}{$strand}{$ini};
				if ($firstrow == 0){
					$prevrow = $cdsline;
					$firstrow++;
				} else {
					my @subl = split (/\t/, $cdsline);
					my @prevsubl = split (/\t/, $prevrow);

					my $prevlen = $prevsubl[4] - $prevsubl[3];
					my $distance = $subl[3] - $prevsubl[4];

#					if ($distance > $cdsdist && $prevlen <= $mincds){
					if ($distance > $cdsdist && $prevlen <= $mincds && $firstrow == 1){ # Only filter first exons, and extend 30 nt to the start of the next

						# Skip line
#						$prevrow = $cdsline;
						my $newini = $subl[3] - 30 + $prevsubl[7] + $subl[7];
						$prevrow = "$subl[0]\t$subl[1]\t$subl[2]\t$newini\t$subl[4]\t$subl[5]\t$subl[6]\t0\t$subl[8]";

					} else {
						print Results "$prevrow\n";
						$prevrow = $cdsline;
					}

					$firstrow++;

				}

			}

			# Print last cds 
			print Results "$prevrow\n";


		} elsif ($strand =~ /\-/){

			# Check and discard if exists cds from the same gene in the other strand
			my $opstrand = "+";
			if (exists $gffscafcds{$gene}{$opstrand}){
				#die "Found opposite strands for gene $gene\n";
				my $cdcount = keys %{$gffscafcds{$gene}{$strand}};
				my $opcdcount = keys %{$gffscafcds{$gene}{$opstrand}};
				if ($cdcount >= $opcdcount){
					delete $gffscafcds{$gene}{$opstrand};
				} else{
					delete $gffscafcds{$gene}{$strand};
					next;
				}
			}
			#	

			my $prevrow = "";
			my $firstrow = "0";

			foreach my $ini (sort { $b <=> $a } keys %{$gffscafcds{$gene}{$strand}}){
				my $cdsline = $gffscafcds{$gene}{$strand}{$ini};
				if ($firstrow == 0){
					$prevrow = $cdsline;
					$firstrow++;
				} else {
					my @subl = split (/\t/, $cdsline);
					my @prevsubl = split (/\t/, $prevrow);

					my $prevlen = $prevsubl[4] - $prevsubl[3];
					my $distance = $prevsubl[3] - $subl[4];

#					if ($distance > $cdsdist && $prevlen <= $mincds){
					if ($distance > $cdsdist && $prevlen <= $mincds && $firstrow == 1){ # Only filter first exons, and extend 30 nt to the start of the next

						# Skip line
#						$prevrow = $cdsline;
						my $newini = $subl[4] + 30 - $prevsubl[7] - $subl[7];
						$prevrow = "$subl[0]\t$subl[1]\t$subl[2]\t$subl[3]\t$newini\t$subl[5]\t$subl[6]\t0\t$subl[8]";

					} else {
						print Results "$prevrow\n";
						$prevrow = $cdsline;
						
					}

				}

			}

			# Print last cds 
			print Results "$prevrow\n";




		} else {
			die "Can't find strand in $gene $strand\n";
		}

	}

}


#foreach my $scaf (sort keys %gffscafcds){
#	foreach my $gene (sort keys %{$gffscafcds{$scad}}){




#	}

#}





close Results;


### END

