#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $Bestn;
GetOptions(
	"best:i"=>\$Bestn
);
$Bestn ||= 1;
die "perl $0 <*.m8> [-best 1]" if @ARGV<1;
my $m8=shift;
open (IN,$m8) or die "$!";
my %hash;my $index;
while (<IN>) {
	my @arr=split /\t/;
	unless (exists $hash{$arr[0]}) {
		$hash{$arr[0]}=1;
		$index=0;
	}
	if (exists $hash{$arr[0]}) {
		$index++;
	}
	if ($index <= $Bestn) {
		print "$_";
	}
}
close IN;
