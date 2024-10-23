#!/usr/bin/perl

=head1 Name

gene_ontology_parser.pl  --  parse and convert the GO flatfiles

=head1 Description

this program read GO flat files, and get GO iterms: (GO id, GO description)
generate a class file and an alias file.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Author: Zheng Hancheng, zhenghch@genomics.org.cn
  Version: 1.0,  Date: 2008-5-22


=head1 Usage

  perl gene_ontology_parser.pl <ontology_file>
  --verbose   output verbose information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl ../bin/gene_ontology_parser.pl function.ontology 
  perl ../bin/gene_ontology_parser.pl component.ontology 
  perl ../bin/gene_ontology_parser.pl process.ontology

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

my ($Verbose,$Help);
GetOptions(
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);

my $ontology_file = shift;

my $class_file="$ontology_file.class";
my $alias_file="$ontology_file.alias";
my %GO;
my %alias;

##read the ontology file
my ($first_class,$second_class);
open IN,$ontology_file || die "fail $ontology_file";
while (<IN>) {
	chomp;
	my $rank = length($1) if(/^(\s*)/);
	next unless ($rank != 0 && /^\s*%|</);
	my ($go_str,$go_id);
	my @temp = split /;/;
	$go_str = $temp[0];
	$go_str =~ s/[%|<]//;
	$go_str =~ s/\\//g;
	$go_str =~ s/^\s+//;
	$go_str =~ s/\s+$//;

	$temp[1]=~s/(%|<) .*//g;	    ##some lines like (GO:0033592 % single-stranded RNA binding)
	my @go_id = split /, /,$temp[1];
	s/\s+//g foreach @go_id;		##the go_id may has some blank, delete it.
	my $go_id=$go_id[0];
	if (@go_id>1){
		$alias{$go_id}=[@go_id[1..@go_id-1]];
	}	

	if($rank == 1){
		$first_class = $go_str;
	}elsif($rank == 2){
		$second_class = $go_str;
		$GO{$first_class}{$second_class}{$go_id} = $go_str;
	}elsif($rank > 2){
		warn "empty GO id: $go_id " if(!$go_id);
		$GO{$first_class}{$second_class}{$go_id} = $go_str;
	}
}
close IN;


##generate the class file
my $output;
foreach my $first_class (sort keys %GO) {
	my $first_p = $GO{$first_class};
	foreach my $second_class (sort keys %$first_p) {
		next if($second_class =~ /obsolete/); ##remove the obsolete class
		my $second_p = $first_p->{$second_class};
		foreach my $go_id (sort keys %$second_p) {
			my $go_str = $second_p->{$go_id};
			$output.= "$first_class\t$second_class\t$go_id\t$go_str\n";
		}
	}
}
open OUT,">$class_file" or die "fail $class_file:$!";
print OUT $output;
close OUT;


##generate the alias file
my $output;
for  my $go_id (sort keys %alias){
	$output.="$go_id\t";
	$output.="$_\t" foreach (sort @{$alias{$go_id}});
	$output.="\n";
}
open OUT,">$alias_file" or die "fail $alias_file:$!";
print OUT $output;
close OUT; 
