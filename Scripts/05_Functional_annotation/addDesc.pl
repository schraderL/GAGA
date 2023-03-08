#!/usr/bin/perl
#Author: liaoxinhui@genomics.org.cn
#Date:  Fri Jul 10 04:30:29 CST 2015


use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use Cwd 'abs_path';
use File::Basename;

my ($input, $annot);
GetOptions(
	"input:s" => \$input,
	"annot:s" => \$annot
);

if (!$input || !$annot) {
	die "usage:perl $0 -input /indir/*.xls -annot annotation.xls\n";
}

my @files = glob("$input");
foreach my $f (@files) {die "Error: can't find $f!\n" unless (-f $f)}

my %desc = ();
open ANNOT,$annot or die $!;
chomp(my $header = <ANNOT>); $header =~ s/^(\S+)\s+//;
while (<ANNOT>) {
	next if (/^#/ || /^\s*$/);
	chomp; my @a = split /\s+/,$_,2;
	$desc{$a[0]} = $a[1];
}
close ANNOT;

foreach my $f (@files) {
	open IN,$f or die $!;
	chomp(my $title = <IN>);
	my $out = "$title\t$header\n";
	while (<IN>) {
		next if (/^#/ || /^\s*$/);
		chomp; my @a = split /\s+/,$_,2;
		die "Error: can't find $a[0] in $annot\t$f!\n" unless ($desc{$a[0]});
		$out .= "$_\t$desc{$a[0]}\n";
	}
	close IN;
	open OUT,">$f" or die $!;
	print OUT $out;
	close OUT;
}
