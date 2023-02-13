#!/usr/bin/perl -w
use strict;
use Data::Dumper;

if (@ARGV != 2) {
	warn "#Usage: perl $0 <gff1> <gff2> \n";
	exit;
}

my $gff1file = shift;
my $gff2file = shift;

my (%gff1,%gff2);
readGff($gff1file,\%gff1);
readGff($gff2file,\%gff2);

#print Dumper \%gff1;

foreach my $chr (keys %gff1) {
	foreach my $id1 (keys %{$gff1{$chr}}) {
		my ($gene1_start,$gene1_end) = @{$gff1{$chr}{$id1}{'mRNA'}};
		my $score1 = $gff1{$chr}{$id1}{'score'};
		my $gene1_len = $gene1_end-$gene1_start+1;
		my @cds1 = @{$gff1{$chr}{$id1}{'CDS'}};
		if (exists $gff2{$chr}) {
			foreach my $id2 (keys %{$gff2{$chr}}) {
				my ($gene2_start,$gene2_end) = @{$gff2{$chr}{$id2}{'mRNA'}};
				my $score2 = $gff2{$chr}{$id2}{'score'};
				my @cds2 = @{$gff2{$chr}{$id2}{'CDS'}};
				my $gene2_len = $gene2_end-$gene2_start+1;
				next if ($gene1_end < $gene2_start or $gene2_end < $gene1_start);
				my $overlapSize = overlapSize($gene1_start,$gene1_end,$gene2_start,$gene2_end);
				my $pencen1 = int($overlapSize*100/$gene1_len)/100;
				my $pencen2 = int($overlapSize*100/$gene2_len)/100;
				
				print "$chr\tmRNA\t$id1\t$gene1_start\t$gene1_end\t$gene1_len\t$score1\t$pencen1\t$id2\t$gene2_start\t$gene2_end\t$gene2_len\t$score2\t$pencen2\t";
				my ($cds1_len,$cds2_len);
				my $cds_overlapSize;
				my $cds_tag = 1;
				for(my $j=0;$j<@cds1;) {
					my $cds1_start = $cds1[$j];
					my $cds1_end   = $cds1[$j+1];
					$cds1_len += $cds1_end-$cds1_start+1;
					$j+=2;
					for(my $i=0;$i<@cds2;)  {
						my $cds2_start = $cds2[$i];
						my $cds2_end   = $cds2[$i+1];
						$cds2_len += $cds2_end-$cds2_start+1 if ($cds_tag);
						$i+=2;

#						next if ($cds1_end < $cds2_start or $cds2_end < $cds1_start);
						if (($cds1_start >= $cds2_start and $cds1_start <= $cds2_end) or ($cds2_start >= $cds1_start and $cds2_start <= $cds1_end)) {
							$cds_overlapSize += overlapSize($cds1_start,$cds1_end,$cds2_start,$cds2_end);
						}
#						my $pencen1 = int ($cds_overlapSize*10000/$cds1_len)/100;
#						my $pencen2 = int ($cds_overlapSize*10000/$cds2_len)/100;
#						print "$chr\tCDS",(($i-2))/2 +1,"\t$id1\t$cds1_start\t$cds1_end\t\t$pencen1\t$id2\t$cds2_start\t$cds2_end\t$cds_overlapSize\t$pencen2\n";
					}
					$cds_tag = 0;
				}
				if (defined $cds_overlapSize) {
#					print "$chr\tCDS\t$id1\t$cds1_len\t$cds_overlapSize\t",int($cds_overlapSize*100/$cds1_len)/100,"\t.\t.\t$id2\t$cds2_len\t$cds_overlapSize\t",int($cds_overlapSize*100/$cds2_len)/100,"\n";
					print "CDS\t$cds_overlapSize\t$cds1_len\t",int($cds_overlapSize*100/$cds1_len)/100,"\t$cds2_len\t",int($cds_overlapSize*100/$cds2_len)/100,"\n";
				}else {
					print "\n";
				}
			}
		}
	}
}
###################################
########### sub routine ###########
sub overlapSize {
	my ($block1_start,$block1_end,$block2_start,$block2_end) = @_;
	my $combine_start = $block1_start <= $block2_start ? $block2_start : $block1_start;
	my $combine_end   = $block1_end   <= $block2_end   ? $block1_end   : $block2_end;

	my $overlapSize = $combine_end - $combine_start+1;
#	return ($combine_start,$combine_end);
	return $overlapSize;
}


sub readGff {
	my $file = shift;
	my $p = shift;

	open IN, "<$file" or die "failed to open $file: $!\n";
	while (<IN>) {
		chomp;
		my $id;
		my @t = split /\t/;

		next if(/^#/);
		next if($t[2] ne "mRNA" and $t[2] ne "CDS");
		if ($t[8] =~ /ID=([^;]+)/ or $t[8] =~ /Parent=([^;]+)/) {
			$id = $1;
		}
		if (/\s+mRNA\s+/) {
			@{$p->{$t[0]}{$id}{'mRNA'}} = ($t[3],$t[4]);
			$p->{$t[0]}{$id}{'score'} = $t[6];
		}elsif (/\s+CDS\s+/) {
			push @{$p->{$t[0]}{$id}{'CDS'}},($t[3],$t[4]);
		}
	}
	close (IN);
}
