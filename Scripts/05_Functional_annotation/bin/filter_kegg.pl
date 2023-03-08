#!/usr/bin/perl

=head1 Name

filter_kegg.pl  -- divide into prokaryote and eukaryote classes

=head1 Description

First divide the genes into prokaryote(Archaea/Bacteria) and eukaryote(Eukaryota) classes

There is three lineages, each with a set of model species:
     Archaea     51
     Bacteria	609
     Eukaryota	137

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-5-21
  Note:

=head1 Usage
  
  perl filter_kegg.pl <genome_file> <seq_pep_dir>
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  ##generate two files: kegg_eukaryote_clean.fa and kegg_prokaryote_clean.fa
  perl filter_kegg.pl genome seq_pep

=cut


use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;


my ($verbose,$help);
GetOptions(
	"verbose"=>\$verbose,
	"help"=>\$help
);
die `pod2text $0` if (@ARGV<2 || $help);

my ($genome,$seq_pep) = @ARGV;
my %Eukaryote;

read_genome($genome,\%Eukaryote);


divide_kegg($seq_pep,\%Eukaryote);

####################################################
################### Sub Routines ###################
####################################################

sub divide_kegg{
	my ($directory,$eukaryote_p) = @_;
	my @files = glob("$directory/*.pep");
	open ALL, ">kegg_all_clean.fa" || die "fail kegg_all_clean.fa";
	open PRO, ">kegg_prokaryote_clean.fa" || die "fail kegg_prokaryote_clean.fa";
	open EUK, ">kegg_eukaryote_clean.fa" || die "fail kegg_eukaryote_clean.fa";
	foreach my $pep_file (@files) {
		open(IN, $pep_file) || die ("can not open $pep_file\n");
		$/=">"; <IN>; $/="\n";
		while (<IN>) {
			my $head = $_;
			my $name = $1 if($head =~ /^(\S+)/);
			my $spec = $1 if($name =~ /^(\w+):/);
			$/=">";
			my $seq = <IN>;
			chomp $seq;
			$/="\n";
			
			next unless(/\s+K\d+\s+/); ##filter those don't have KO id
			print ALL ">".$head.$seq;
			if (exists $eukaryote_p->{$spec} && $eukaryote_p->{$spec} eq "Eukaryota") {
				print EUK ">".$head.$seq;
			}else{
				print PRO ">".$head.$seq;
			}
	
		}
		close(IN);
	}
	close ALL;
	close PRO;
	close EUK;
}

sub read_genome{
	my ($file,$eukaryote_p) = @_;
	$/ = "///";
	open IN,$file || die "fail open $file";
	while (<IN>) {
		my $species = $1 if(/ENTRY\s+(\S+)/);
		my $lineage = $1 if(/LINEAGE\s+(\w+);/);
		$eukaryote_p->{$species} = $lineage if($species);
	}
	close IN;
	$/ = "\n";
}

