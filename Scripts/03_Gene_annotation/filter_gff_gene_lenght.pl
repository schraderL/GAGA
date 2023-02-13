#!/usr/bin/perl

use strict;
use File::Basename qw(basename dirname);
use Getopt::Long;

my $help="\n\t--threshold\tleast length of cds [default=500]\n\t--exons\tthe # of exons want to filter out [default=1]\n\t--help\tshow this help inf\n";
$help.="\n\tThe 'threshold' is the total length of exons in a mRNA,\n\tfilting which short then it. default=500\n";

my ($thre,$exons,$Help);
GetOptions(
	"threshold:i"=>\$thre,
	"exons:i"=>\$exons,
	"help"=>\$Help,
);
$thre ||= 500;
$exons ||= 1;

die"Usage:\n\tperl $0 <gff> [options] > outfile\n\t$help\n" if(@ARGV < 1 || $Help);

my $file=shift;
my $file_neme=basename($file);

my %gff;
read_gff($file,\%gff);

for my $name(sort keys %gff){
	for my $id(keys %{$gff{$name}}){
        if($gff{$name}{$id}{len} < $thre && $gff{$name}{$id}{cds_num} <= $exons) {
            ;
        }else {
            print $gff{$name}{$id}{gene};
        }
#		print $gff{$name}{$id}{gene} unless($gff{$name}{$id}{len} < $thre && $gff{$name}{$id}{cds_num} <= $exons);
	}
}

#################################################
#################################################
#################################################
sub read_gff{
	my ($file,$gff)=@_;
	open(IN,$file)||die "can not open $file\n";
	while(<IN>){
		next if(/^#/);
		chomp;
		my @a=split /\t+/;
		next if($a[2] eq 'gene');
		($a[3],$a[4])=($a[4],$a[3]) if($a[3]>$a[4]);
		my ($name,$id);
		$name=$a[0];
		if($a[2] eq 'mRNA'){
			$id=$1 if($a[8]=~/ID=([^;]+)/);
#			$a[8]="ID=$id;\n";
		}elsif($a[2] eq "CDS"){
			$id=$1 if($a[8]=~/Parent=([^;\s]+)/);
#			$a[8]="Parent=$id;\n";
			$$gff{$name}{$id}{cds_num}++;
			$$gff{$name}{$id}{len}+=(($a[4]-$a[3])+1);
		}
		$$gff{$name}{$id}{gene}.=join("\t",@a)."\n";
	}
	close IN;
}
