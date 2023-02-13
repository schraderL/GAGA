#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

if (@ARGV < 1){
	print "Usage:\n\tperl $0 <gff> [mRNA | CDS]>output\n";
	exit;
}

my $file = shift;
my $flag = shift;

$flag ||= 'CDS';

my %gff;
open (IN, $file) || die $!;
while (<IN>){
	chomp;
	next if (/^#/);
	my @c = split /\t/, $_;
	my $id = $1 if ($c[8] =~ /=(\S+?);/);
	next unless ($c[2] =~ /$flag/i );
	$gff{$id}{$flag} += ($c[4] - $c[3] + 1);
}
close IN;

foreach my $id (sort keys %gff){
	print "$id\t$gff{$id}{$flag}\n";
}

