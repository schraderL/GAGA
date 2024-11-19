#!/usr/bin/perl -w
use strict;
use Data::Dumper;

if(@ARGV != 3) {
	warn "#Usage: perl $0 <sbs> <car_len> <gene_len> \n";
	exit;
}

my $sbs_file      = shift;
my $car_lenfile   = shift;
my $gene_lenfile  = shift;

my %block_size;
readBlock($sbs_file,\%block_size);

my %length_chr1;
my %genes1;
readGeneLength($car_lenfile,\%length_chr1,\%genes1);

my %length_chr2;
my %genes2;
readGeneLength($gene_lenfile,\%length_chr2,\%genes2);

open IN, "<$sbs_file" or die "failed to open $sbs_file:$!\n";
<IN>;
my $line = <IN>;
my @cache;
my @arr = split /\t/,$line;
my @gene1 = split /\s+/,$arr[12];
my @gene2 = split /\s+/,$arr[13];
push @cache,[$arr[1],$arr[9],$gene1[0],$gene2[0],$arr[0]];
while (<IN>) {
	chomp;
	my @arr=split /\t/;
	next if($block_size{$arr[0]} < 10);
	my @gene1 = split /\s+/,$arr[12];
	my @gene2 = split /\s+/,$arr[13];


	push @cache,[$arr[1],$arr[9],$gene1[0],$gene2[0],$arr[0]];

	if($length_chr1{$arr[1]} >=10 && $length_chr2{$arr[9]} >=10) {

		if($cache[0][4] ne $arr[0]) {

			if($cache[0][1] ne $arr[9] && $cache[0][0] ne $arr[1]) {
				;
			}elsif ($genes2{$cache[0][3]} <=3 && $genes2{$cache[0][3]} >= $length_chr2{$cache[0][1]}-3) {
				;
			}elsif ($genes1{$cache[0][2]} <=3 && $genes1{$cache[0][2]} >= $length_chr1{$cache[0][0]}-3) {
				;
			}else {
				print "$cache[0][4]-$arr[0]\t$cache[0][0]\t$arr[1]\t$cache[0][1]\t$arr[9]\n";
			}
		}
	}else{
        warn "$_\n";
    }
	shift @cache;
}
close (IN);

sub readCarLength {
	my ($file,$p) = @_;

	open IN, "<$file" or die "failed to open $file:$!\n";
	while (<IN>) {
		chomp;
		my @t=split /\s+/;
		$p->{$t[0]} = $t[1];
	}
	close (IN);
}

sub readGeneLength {
	my ($file,$p,$p2) = @_;

	open IN, "<$file" or die "failed to open $file:$!\n";
	while (<IN>) {
		chomp;
		my @t=split /\s+/;
		$p->{$t[0]} = $t[1];
		$p2->{$t[2]} = $t[3];
	}
	close (IN);
}

sub readBlock {
	my ($file,$p) = @_;

	open IN, "<$file" or die "failed to open $file:$!\n";
	while (<IN>) {
		chomp;
		my @t=split /\s+/;
		$p->{$t[0]}++;
	}
	close (IN);
}
