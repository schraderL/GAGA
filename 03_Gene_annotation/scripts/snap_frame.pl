#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper; 
use File::Basename qw(dirname basename);

if (@ARGV < 1){
	print "Usage:\n\tperl $0 <input_gff> >output_gff\n";
	exit;
}

my $file = shift;
my %gff;
open(IN, $file) || die $!;
while (<IN>){
	chomp;
	next if (/^#/);
	my @c = split;
	($c[3], $c[4]) = ($c[4], $c[3]) if ($c[3] > $c[4]);
	my $id = $1 if ($c[8] =~ /=(\S+?);/);
	push @{$gff{$c[0]}{$id}}, [@c];
}
close IN;

foreach my $chr (keys %gff){
	foreach my $id (keys %{$gff{$chr}}){
		my $temp = join "\t", @{$gff{$chr}{$id}[0]};
		$temp  .= "\n";
		shift @{$gff{$chr}{$id}};
		my @table = sort {$a->[3] <=> $b->[3]} @{$gff{$chr}{$id}};
		@table = reverse @table if ($table[0][6] eq '-');
		my $shift = 0;
		foreach my $p (@table){
			$p->[7] = $shift;
			$temp .= join "\t", @$p;
			$temp .= "\n";
			my $length = $p->[4] - $p->[3] + 1 + $shift;
			$shift = $length % 3;
		}
		print $temp;
	}
}
