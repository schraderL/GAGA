#!usr/bin/perl
use strict;

=head1 description

this program is to caculate the length of each chromosome of given species

=head1 example

perl get_chr_len.pl <.fa>

=cut
die `pod2text $0` if(@ARGV==0);
my $fa=shift;
my (%len,$name,$length,$row,$col);
open IN,$fa or die "$!";
open OUT,">$fa.chr.len" or die "$!";
while(<IN>)
{
	chomp;
	if(/^>(\S+)/){$name=$1;$length=0;$row=0;}
	$/="(>\w+)";
	$length+=length $_ unless($row==0);
	$row++;
	$/="\n";
	$len{$name}=$length;
}
foreach my $chr (sort keys %len)
{
	print OUT "$chr\t$len{$chr}\n";
}
close IN;
close OUT;
#`/opt/blc/self-software/bin/msort  -krn2 $fa.chr.len > $fa.chr.len.sort`;
`/share/app/msort-20160805/bin/msort -krn2 $fa.chr.len > $fa.chr.len.sort`;
