#!/usr/bin/perl -w
use strict;
use Data::Dumper;

if (@ARGV != 1) {
	warn "#Usage: perl $0 <gff3> \n";
	exit;
}

my $file = shift;

my %gff3;
open IN, "<$file" or die "failed to open $file: $!\n";
while (<IN>) {
	chomp;
	my @t = split /\t/;
	my $gene_name = $1 if ($t[8]=~/gene_name=([^;]+);/);
	my $trans_id  = $1 if ($t[8]=~/ID=([^;]+);/ || $t[8]=~/Parent=([^;]+);/);
	if (defined $gene_name && $t[2] eq 'mRNA' ) {
		$t[7] = '.';
		my $trans_id  = $1 if ($t[8]=~/ID=([^;]+);/);
		push @{$gff3{$gene_name}{mRNA}}, [@t];
	}elsif (defined $gene_name && $t[2] eq 'CDS') {
		my $trans_id  = $1 if ($t[8]=~/Parent=([^;]+);/);
		push @{$gff3{$gene_name}{$trans_id}{CDS}},[@t];
	}elsif (defined $gene_name && $t[2] eq 'exon') {
		$t[7] = '.';
		my $trans_id  = $1 if ($t[8]=~/Parent=([^;]+);/);
		push @{$gff3{$gene_name}{$trans_id}{exon}},[@t];
	}elsif (defined $gene_name && $t[2] eq 'start_codon') {
		$t[7] = 0; 
		my $trans_id  = $1 if ($t[8]=~/Parent=([^;]+);/);
		@{$gff3{$gene_name}{$trans_id}{start_codon}} = @t;
	}elsif (defined $gene_name && $t[2] eq 'stop_codon') {
		$t[7] = 0;
		my $trans_id  = $1 if ($t[8]=~/Parent=([^;]+);/);
		@{$gff3{$gene_name}{$trans_id}{stop_codon}} = @t;
	}elsif (defined $gene_name && $t[2] eq 'five_prime_UTR') {
		$t[7] = '.';
		my $trans_id  = $1 if ($t[8]=~/Parent=([^;]+);/);
		push @{$gff3{$gene_name}{$trans_id}{five_prime_UTR}}, [@t];
	}elsif (defined $gene_name && $t[2] eq 'three_prime_UTR') {
		$t[7] = '.';
		my $trans_id  = $1 if ($t[8]=~/Parent=([^;]+);/);
		push @{$gff3{$gene_name}{$trans_id}{three_prime_UTR}}, [@t];
	}else {
		print STDERR "$_\n";
	}
}
close (IN);

#print Dumper \%gff3;
#exit;

print "##gff-version 3\n";
foreach my $gene (keys %gff3) {
	print "#Gene: $gene\n";
	my @genes = @{$gff3{$gene}{mRNA}};
	my @genes_temp1 = sort{$a->[3] <=> $b->[3]} @genes;
	my @genes_temp2 = sort{$b->[4] <=> $a->[4]} @genes;

	my $gene_start = $genes_temp1[0]->[3];
	my $gene_end   = $genes_temp2[0]->[4];

	my $scaff  = $genes_temp1[0]->[0];
	my $strand = $genes_temp1[0]->[6];
	print join("\t",$scaff,'StringTie','gene',$gene_start,$gene_end,".",$strand,".","ID=$gene"),"\n";
	foreach (@genes) {
		my @mRNA = @{$_};
		$mRNA[8] =~ s/gene_name/Parent/;
		print join("\t",@mRNA),"\n";

		my $trans_id = $1 if ($mRNA[8]=~/ID=([^;]+);/);

		my @CDS  = @{$gff3{$gene}{$trans_id}{CDS}};
		my @EXON = @{$gff3{$gene}{$trans_id}{exon}} if (exists $gff3{$gene}{$trans_id}{exon});
		my @START_CODON = @{$gff3{$gene}{$trans_id}{start_codon}} if (exists $gff3{$gene}{$trans_id}{start_codon});
		my @STOP_CODON = @{$gff3{$gene}{$trans_id}{stop_codon}}   if (exists $gff3{$gene}{$trans_id}{stop_codon});
		my @UTR5 = @{$gff3{$gene}{$trans_id}{five_prime_UTR}}     if (exists $gff3{$gene}{$trans_id}{five_prime_UTR});
		my @UTR3 = @{$gff3{$gene}{$trans_id}{three_prime_UTR}}    if (exists $gff3{$gene}{$trans_id}{three_prime_UTR});


#		if (defined $UTR5[8]) {
            foreach my $utr5 (@UTR5) {
                my @out = @{$utr5};
                $out[8] =~ s/gene_name=([^;]+);//;
                print join("\t",@out),"\n";
            }
#			$UTR5[8] =~ s/gene_name=([^;]+);//;
#			print join("\t",@UTR5),"\n";
#		}
		if (defined $START_CODON[8]) {
			$START_CODON[8] =~ s/gene_name=([^;]+);//;
			print join("\t",@START_CODON),"\n";
		}
		foreach my $cds(@CDS) {
			my @out = @{$cds};
			$out[8] =~ s/gene_name=([^;]+);//;
			print join("\t",@out),"\n";
		}
		if (defined $STOP_CODON[8]) {
			$STOP_CODON[8] =~ s/gene_name=([^;]+);//;
			print join("\t",@STOP_CODON),"\n";
		}
#		if (defined $UTR3[8]) {
            foreach my $utr3 (@UTR3) {
                my @out = @{$utr3};
                $out[8] =~ s/gene_name=([^;]+);//;
                print join("\t",@out),"\n";
            }
#			$UTR3[8] =~ s/gene_name=([^;]+);//;
#			print join("\t",@UTR3),"\n";
#		}
		if (defined $EXON[0]) {
			foreach my $exon(@EXON) {
				my @out = @{$exon};
				$out[8] =~ s/gene_name=([^;]+);//;
				print join("\t",@out),"\n";
			}
		}
	}
	print "###\n";
}
