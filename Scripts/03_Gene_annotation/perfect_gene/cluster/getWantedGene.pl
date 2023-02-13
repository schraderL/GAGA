#!/usr/bin/perl

use strict;
use File::Basename qw(basename dirname);
use FindBin;
use Data::Dumper;

my $help=<<USAGE;

get what you wanted genes based hcluster , len and list files.

Liu Shiping liushping\@genomics.org.cn
2010-5-19

USAGE:
	perl $0 <list> <XX.hcluster> <XX.len> 

USAGE

die "$help\n" if @ARGV <3;

my $file_list=shift;
my $hcluster=shift;
my $Len=shift;

my $name=basename $hcluster;

my %len;
read_len($Len,\%len);

#print Dumper \%len;
#print "\nlen \n";

my %list;
read_list($file_list,\%list);

#print Dumper \%list;
#print "\nlist\n";

open OUT1,">$hcluster.reference" or die $!;
open OUT2,">$hcluster.all.list" or die $!;
open IN,$hcluster or die $!;
while(<IN>){
	next if /^\s|#/;
	chomp;
	my @a=split /\t+/;
	shift @a;
	shift @a;
	my @is;
	my @score;
	
	## init
	for(my $i=0;$i<@a;$i++){
		$score[$i][$i]=$len{$a[$i]};
		for(my $j=0;$j<@a;$j++){
			next if($i==$j);
			if($list{$a[$i]}{$a[$j]}==1){
				$score[$i][$j]=1;
			}else{
				$score[$i][$j]=0;
			}
		}
	}

	## find begin
	while(1){
		my $max=$score[0][0];
		my $fl=0;
		for(my $i=1;$i<@a;$i++){
			if($score[$i][$i]>$max){
				$max=$score[$i][$i];
				$fl=$i;
			}
		}
		last if($max==0);
		push @is,$fl;
		$score[$fl][$fl]=0;
		for(my $j=0;$j<@a;$j++){
			$score[$j][$j]=0 if($score[$fl][$j]==1);
		}
		
	}
	print OUT1 join("\t",@a)."\n";
	print OUT2 join("\n",@a)."\n";
	for(my $i=0;$i<@is;$i++){
		$a[$is[$i]]=undef;
	}
	for my $aa(@a){
		print "$aa\n" if($aa);
	}	
}
close IN;
close OUT1;
close OUT2;


####################3
sub read_len{
	my ($file,$hash)=@_;
	open IN,$file or die $!;
	while(<IN>){
		next if /^\s|#/;
		if(/(\S+)\s+(\S+)/){
			$$hash{$1}=$2;
		}
	}
	close IN;
}

sub read_list{
	my ($file,$hash)=@_;
	open IN,$file or die $!;
	while(<IN>){
		next if /^\s|#/;
		if(/(\S+)\s+(\S+)/){
			$$hash{$1}{$2}=1;
			$$hash{$2}{$1}=1;
		}
	}
	close IN;
}
