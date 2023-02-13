#!/usr/bin/perl
use strict;
my (%filter,%len);
die "Usage : perl $0 <fa><.m8> " if(@ARGV==0);
use FindBin qw($Bin $Script);
my $get_chr_len = "$Bin/get_chr_len.pl";
`perl $get_chr_len  $ARGV[0]`;
#`perl /nas/GAG_02/pengchunfang/GACP/GACP-7.0/03.repeat_finding/denovo_predict/bin/get_chr_len.pl   $ARGV[0]`;
open IN,"$ARGV[0].chr.len" or die "$!";
while(<IN>)
{
		chomp;split;
		$len{$_[0]}=$_[1];
}
close IN;

open IN,$ARGV[1] or die "$!";
while(<IN>)
{
		chomp;split;
		next if($_[0] eq $_[1] || $_[2]<80);
		my $align_len=$_[3];
		my $rate1=$align_len/$len{$_[0]};
		my $rate2=$align_len/$len{$_[1]};
		if($rate1>=0.5 && $rate2>=0.5)
		{
				if($len{$_[0]}<=$len{$_[1]})
				{
						$filter{$_[0]}=1;
				}
				else
				{
						$filter{$_[1]}=1;
				}
		}
}
close IN;

open FA,$ARGV[0] or die "$!";
open OUT,">$ARGV[0].unredundance" or die "$!";
$/=">";<FA>;$/="\n";
while(<FA>)
{
		my $name=$1 if(/(\S+)/);

		$/=">";
		my $seq=<FA>;
		chomp $seq;
		$/="\n";
		
		unless (exists $filter{$name})
		{
				print OUT ">$name\n$seq";
		}
}
close FA;
close OUT;
