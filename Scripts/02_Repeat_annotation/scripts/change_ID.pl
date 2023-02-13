#!/usr/bin/perl

$infile = shift;
open (IN,$infile) or die $!;
open (OUT , ">" . "$infile.denovo.gff") or die $!;

while (chomp ($line = <IN>))
{
	$line =~ s/ID=TE/ID=Denovo_TE/g;
	print OUT "$line\n";
}
close (IN);
close (OUT);
