#!/usr/bin/perl

=head1 Name

get_hcluster_redundance.pl  --  get the redundant elements from hcluster result

=head1 Description

There is two types of representative gene:
(1) the one with biggest length
(2) the one with most relations

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-6-17
  Note:

=head1 Usage
  perl get_hcluster_redundance.pl [options] <overlap_dist_file> <hcluster_file>
  --type <str>  max_length or max_relation, the type of representative gene, default length
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple



=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

##get options from command line into variables and set default values
my ($Type);
my ($Verbose,$Help);
GetOptions(
	"type:s"=>\$Type,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Type ||= "length";
die `pod2text $0` if (@ARGV < 2 || $Help);

my $overlap_file = shift; ##result file of findOverlap_CalcuDist.pl
my $hcluster_file = shift; ##result file of hcluster.pl

my $output;
my $output_reference;
my $output_all_list;
my %relation_number;
my %gene_length;

open IN,$overlap_file || die "fail $overlap_file";
while (<IN>) {
	my @t = split /\t/;
	next if($t[0] eq $t[1]);
	
	if($Type eq "max_relation"){
		$relation_number{$t[0]}++;
		$relation_number{$t[1]}++;
	}
	
	if($Type eq "max_length"){
		$gene_length{$t[0]} = $t[6];
		$gene_length{$t[1]} = $t[10];
	}
}
close IN;

open IN,$hcluster_file || die "fail $hcluster_file";
while (<IN>) {
	chomp;
	my @t = split /\t/;
	my $cluster_id = shift @t; ## remove Cluster_id
	my $element_num = shift @t; ## reomve element number
	my %temp;
	foreach my $gene_id (@t) {
		$temp{$gene_id} = $relation_number{$gene_id}  if($Type eq "max_relation");
		$temp{$gene_id} = $gene_length{$gene_id}  if($Type eq "max_length");
	}
	my @sorted = sort {$temp{$b} <=> $temp{$a}} keys %temp;
	$output_all_list .= join("\n",@sorted)."\n";
	$output_reference .= join("\t",@sorted)."\n";
	my $represent = shift @sorted; ## remove the representative gene
	foreach  (@sorted) {
		$output .= "$_\t$cluster_id\n";
	}
}
close IN;

print  $output;

open OUT,">$hcluster_file.reference" || die "fail $hcluster_file.reference";
print OUT $output_reference;
close OUT;

open OUT,">$hcluster_file.all.list" || die "fail $hcluster_file.all.list";
print OUT $output_all_list;
close OUT;