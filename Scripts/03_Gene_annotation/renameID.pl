#!/usr/bin/perl 

use strict;

die "Usage:
	perl $0 <gff> <new name of ID>\n" if(@ARGV <2);

my $file=shift;
my $name=shift;

my %gff;
open IN,$file or die $!;
while(<IN>){
	next if /^\s|#/;
	chomp;
	my @a=split /\s+/;
	my $id;
	if($a[2] eq 'mRNA' && $a[8]=~/ID=([^;]+)/){
		@{$gff{$a[0]}{$1}{mRNA}}=@a;
	}elsif($a[8]=~/Parent=([^;]+)/){
		push @{$gff{$a[0]}{$1}{CDS}},[@a];
	}
}
close IN;

my $num="00000";
for my $scaff(keys %gff){
	foreach my $id (sort {@{$gff{$scaff}{$a}{mRNA}}[3] <=> @{$gff{$scaff}{$b}{mRNA}}[3]} keys %{$gff{$scaff}}) {
		$num++;
		$gff{$scaff}{$id}{mRNA}[8]="ID=$name\_$num;";
		print STDERR "$id\t","$name\_$num\n";
		print join("\t",@{$gff{$scaff}{$id}{mRNA}})."\n";
		for my $cds(@{$gff{$scaff}{$id}{CDS}}){
			$cds->[8]="Parent=$name\_$num;";
			print join("\t",@{$cds})."\n";
		}
	}
}
