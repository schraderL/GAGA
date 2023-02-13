#!/usr/bin/perl
use strict;
my %filter;
open IN,$ARGV[0] or die "$!";
while(<IN>)
{
	next if(/^Query_id/);
	chomp;split;
	my $align_rate1=abs(($_[3]-$_[2])/$_[1]);
	my $align_rate2=abs(($_[7]-$_[6])/$_[5]);
	if($align_rate1>=0.8 || $align_rate2>=0.8)
	{
		$filter{$_[4]}=1;	
	}
}
close IN;

open IN,$ARGV[1] or die "$!";
open OUT,">$ARGV[1].filter_contaminate" or die "$!";
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
