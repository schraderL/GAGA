#!/usr/bin/perl
use strict;

## Modified: add @a to split function, so that it works well for new perl (>5.11)
## Author: fangqi@genomics.cn
## 20180809

my %filter;
open IN,$ARGV[0] or die "$!";
while(<IN>)
{
	next if(/^Query_id/);
	chomp; #split;
    my @a = split /\s+/;
	my $align_rate1=abs(($a[3]-$a[2])/$a[1]);
	my $align_rate2=abs(($a[7]-$a[6])/$a[5]);
	if($align_rate1>=0.3 && $align_rate2>=0.3 && ($a[3]-$a[2])>90)
	{
		$filter{$a[0]}=1;	
	}
}
close IN;

open IN,$ARGV[1] or die "$!";
open OUT,">$ARGV[1].filter_gene" or die "$!";
$/=">";<IN>;$/="\n";
while(<IN>)
{
		my $name=$1 if(/(\S+)/);

		$/=">";
		my $seq=<IN>;
		chomp $seq;
		$/="\n";

		unless(exists $filter{$name})
		{
				print OUT ">$name\n$seq";
		}
}
close IN;
close OUT;
