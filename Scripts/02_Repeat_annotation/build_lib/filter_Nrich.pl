#!/usr/bin/perl
use strict;
open IN,$ARGV[0] or die "$!";
open OUT,">$ARGV[0].tmp" or die "$!";
$/=">";<IN>;$/="\n";
while(<IN>)
{
		my $head=$_;
		
		$/=">";
		my $seq=<IN>;
		chomp $seq;
		$seq=~s/\s//g;
		$/="\n";
		
		my $tmp=$seq;
		my $length=length $seq;
		my $N_num=($tmp=~s/N/X/ig);
		if($N_num/$length<0.05)
		{
				print OUT ">$head";
				for (my $i=0;$i<$length;$i+=70)
				{
						my $sub=substr($seq,$i,70);
						print OUT "$sub\n";
				}
		}
		
}
close IN;
close OUT;
if($ARGV[1]=~/LTR/){
`sed 's/Unknown/LTR/' $ARGV[0].tmp > $ARGV[0].tmp2`;
`rm -f $ARGV[0] $ARGV[0].tmp`;
`mv $ARGV[0].tmp2 $ARGV[0]`;
}
else{
		`rm -f $ARGV[0]`;
		`mv $ARGV[0].tmp $ARGV[0]`;
	}
