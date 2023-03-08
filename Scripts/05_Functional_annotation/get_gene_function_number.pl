#!/usr/bin/perl
use strict;
use warnings;

my $file1=shift;
my $file2=shift;
my $file3=shift;
my $file4=shift;
my $file5=shift;

my %ID;

open (IN1,$file1) || die "$!";
while(<IN1>){
	my @c=split;
	my $id_name=$c[0];
	$ID{$id_name}=$id_name;
}
close IN1;

open (IN2,$file2) || die "$!";
while(<IN2>){
	my @c=split;
	my $id_name=$c[0];
	$ID{$id_name}=$id_name;
}
close IN2;

open (IN3,$file3) || die "$!";
while(<IN3>){
	 my @c=split;
	 my $id_name=$c[0];
	 $ID{$id_name}=$id_name;
 }
close IN3;

open (IN4,$file4) || die "$!";
while(<IN4>){
	my @c=split;
	my $id_name=$c[0];
	$ID{$id_name}=$id_name;
}
close IN3;

open (IN5,$file5) || die "$!";
while(<IN5>){
	my @c=split;
	my $id_name=$c[0];
	$ID{$id_name}=$id_name;
}
close IN5;

 my @key=keys %ID;
 my $number=@key;
 print "$number\n";

			      
