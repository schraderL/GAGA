#!/us/bin/perl
use strict;
use warnings;
use Data::Dumper;

if (@ARGV!=2){
	print "Usage:\n\tperl $0 <gff> <threshold> >output_gff\n";
	exit;
}

my $file = shift;
my $threshould = shift;

my %gff;
open (IN, $file) || die;
while (<IN>){
	chomp;
	next if (/^#/);
	my @c = split /\t/, $_;
	my $id = $1 if ($c[8] =~ /=(\S+?);/);
	if($c[2] eq 'mRNA'){
		$gff{$id}{value} = $c[5];
		$gff{$id}{Shift}=$1 if($c[8]=~/Shift=(\d+);/);
	}
	$gff{$id}{out} .= "$_\n";
}
close IN;
#print Dumper %gff;
#exit;

foreach my $id (sort keys %gff){
	next if ($gff{$id}{value} < $threshould);
	print $gff{$id}{out};
}

