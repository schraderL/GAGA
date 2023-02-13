#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

#Changed line 61 for 62, to print the trimmed positions

my ($line, $name, $nameout);

my $protfile = $ARGV[0];
my $listfile = $ARGV[1];
my $prefix = "";
if ($listfile =~ /(\S+)\_list.txt/){
	$prefix = $1;
} else {die "Can't find prefix on list file\n";}

my $outfile = "$prefix\_full_fixed_list.txt";

# Reading Protein fasta

my %fasta;
open(File, "<", $protfile);
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	

	if ($line =~ />(\S+)/){
		$name = $1;
	} else {
		$fasta{$name} .= "$line";
	}
}
close File;

my %lengthseq;
foreach my $prot (sort keys %fasta){
	my $length = length ($fasta{$prot});
	$lengthseq{$prot} = $length;
}


# Reading list and printing output

open (Results, ">", "$outfile");

open(File, "<", $listfile);
while(<File>){
	chomp;
	my $line = $_;
	$line =~ s/\r//g;
	next if ($line !~ /\S+/);
	my @subl = split (/\s/, $line);

	if ($line =~ /(\S+)\_split/){
		my $gene = $1;
		print Results "$line\n";
		delete($fasta{$gene});
	} else {
#		print Results "$subl[0] annot 1 $lengthseq{$subl[0]} blastp 0\n";
		print Results "$line\n";
		delete($fasta{$subl[0]});
	}
}
close File;

foreach my $gene (sort keys %fasta){
	print Results "$gene annot 1 $lengthseq{$gene} noblastp 0\n";
}


close Results;

