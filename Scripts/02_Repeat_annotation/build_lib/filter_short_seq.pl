#!/usr/bin/perl
use strict;

open IN,$ARGV[0] or die "$!";
my $flag=$ARGV[1];##### Please give the name of repeat software you used! #####;
open OUT,">$ARGV[0].100bp" or die "$!";
$/=">";<IN>;$/="\n";
while(<IN>)
{
		my $name=$1 if(/(\S+)/);

		$/=">";
		my $seq=<IN>;
		chomp $seq;
		$seq=~s/\s//g;
		$/="\n";

		my $length=length $seq;
		if($length>=100)
		{
				print OUT ">$flag"."_$name\n";
				for(my $i=0;$i<=$length;$i+=70)
				{
						my $sub=substr($seq,$i,70);
						print OUT "$sub\n";
				}
		}
}
close IN;
close OUT;
