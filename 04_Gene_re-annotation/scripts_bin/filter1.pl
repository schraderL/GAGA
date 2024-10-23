#!/usr/bin/perl -w
use strict;

if (@ARGV != 2)  {
	warn "#Usage: perl $0 <Overlap.lst> <prefix> \n";
	exit;
}
my $file = shift;
my $prefix = shift;

open IN, "<$file" or die "failed to open $file:$!\n";
<IN>;
while (<IN>) {
	my @t=split /\s+/;
	if(/$prefix[a-zA-Z1-9]/)	{
		print $_;
	}
}
close (IN);
