#!/usr/bin/perl -w
use strict;

if(@ARGV !=2) {
    warn "#Usage: perl $0 <blast> <swissprot_pep> \n";
    exit;
}

my $blastOut   = shift;
my $trembl_pep = shift;

my %Out;
open IN, "<$blastOut" or die "failed to open $blastOut: $!\n";
while (<IN>) {
	chomp;
	my @t = split /\s+/;
	push @{$Out{$t[1]}},$t[0];
}
close (IN);

open IN, "<$trembl_pep" or die "failed to open $trembl_pep: $!\n";
while (<IN>) {
	chomp;
	if (/^>(\S+)/) {
		my $id = $1;
		if (exists $Out{$id}) {
			foreach my $name (@{$Out{$id}}) {
				print "$name\t$_\n";
			}
		}
	}
}
close (IN);
