#!/usr/bin/perl
use strict;
die "Usage:perl $0 <.classified>" if(@ARGV==0); 
open IN,$ARGV[0] or die "$!";
open OUT,">$ARGV[0].unknown" or die "$!";
open KW,">$ARGV[0].known" or die "$!";
$/=">";<IN>;$/="\n";
while(<IN>)
{
		my $name=$_ ;

		$/=">";
		my $seq=<IN>;
		chomp $seq;
		$/="\n";

		if($name=~/Unknown/)
		{
				print OUT ">$name$seq";
		}
		else
		{
				print KW ">$name$seq";
		}
		
}
close IN;
close OUT;
close KW;
