#!/usr/bin/perl
use strict;
use warnings;
use File::Basename qw(dirname basename);
use Data::Dumper;

if (@ARGV < 2){
	print "Usage:\n\tperl $0 <input_gff> <genome_fa> [err_cutoff,default=1] >output_pseudo_gff\n";
	exit;
}

my $get_gene = "perl /nas/GAG_01A/assembly/yanglinfeng/BGI_BC_GAG/Annotation/common_bin/getGene.pl ";
my $cds2aa = "perl /nas/GAG_01A/assembly/yanglinfeng/BGI_BC_GAG/Annotation/common_bin/cds2aa.pl ";

my ($file1, $file2, $err_cutoff) = @ARGV;
$err_cutoff ||=1;

my %gff;
open (IN, $file1) || die $!;
while (<IN>){
	next if (/^#/);
	chomp;
	my @c = split;
	my $id = $1 if ($c[8] =~ /=(\S+?);/);
	if ($c[2] eq 'mRNA'){
		my $shift = 0;
		$shift = $1 if ($c[8] =~/Shift=(\d+);/ );
		$c[8] =~ s/Shift=(\d+);//;
		$gff{$id}{shift} = $shift;
	}elsif ($c[2] eq 'CDS'){
		$gff{$id}{cds_num} ++;
	}
	push @{$gff{$id}{out}}, [@c];
}
close IN;

my $dirname = dirname ($file1);
my $basename = basename ($file1, '.gff');

`$get_gene $file1 $file2 >./$basename.cds`;
`$cds2aa -check ./$basename.cds | awk \'\$4 == 0\' >./$basename.inner_stop`;

my %list;
open (IN, "./$basename.inner_stop") || die $!;
while (<IN>){
	next if (/^Id/);
	my @c=split;
	$list{$c[0]} = 1;
}
close IN;

foreach my $id (keys %gff){
	my $inner_stop = ($list{$id}) ? $list{$id} : 0;
	if ($inner_stop + $gff{$id}{shift} >= $err_cutoff){
		foreach (@{$gff{$id}{out}}){
			$_->[8] .= "Shift=$gff{$id}{shift};inner_stop_codon=$inner_stop;" if ($_->[2] eq 'mRNA');
			my $tmp = join "\t", @$_;
			print "$tmp\n";
		}
	}
}

