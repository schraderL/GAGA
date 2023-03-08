#!/usr/bin/perl
#Author: liaoxinhui@genomics.org.cn
#Date:  Wed Jul  8 16:44:53 CST 2015


use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use Cwd 'abs_path';
use File::Basename;

die "usage: perl $0 [in.iprscan.tsv] [out.xls]\n" unless (@ARGV == 2);

open IN,"$ARGV[0]" or die $!;
open OUT,">$ARGV[1]" or die $!;
print OUT "Query_id\tSubject_id\tSubject_DB\tQuery_start\tQuery_end\tE_value\tSubject_annotation\n";

while (<IN>) {
	next if (/^#/ || /^\s*$/);
	chomp; my @a = split /\t/;
#	$a[0] =~ s/_\d+$//;
	$a[5] = ($a[5] =~ /^\s*$/) ? "--" : $a[5];
	print OUT "$a[0]\t$a[4]\t$a[3]\t$a[6]\t$a[7]\t$a[8]\t$a[5]\n";
}
close IN;
close OUT;
