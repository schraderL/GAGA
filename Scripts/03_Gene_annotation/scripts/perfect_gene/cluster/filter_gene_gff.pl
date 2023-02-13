#!/usr/bin/perl

=head1 Name

filter_on_genomic_level.pl  -- filter the gene gff file 

=head1 Description

this program filter the gene gff file based on cds length and exon number

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-4-19

=head1 Usage
  
  perl remove_redundance_on_genomic_level.pl <mRNA.gff>
  --min_cds_len <num>   set the minimum cds length, default 150bp
  --min_exon_num <num>  set the minimum exon number, default 1
  --outdir <str>   set the result directory, default="./" 
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

 perl ../bin/filter_gene_gff.pl ../input/rice_KOME_cDNA_5000.fa.blat.sim4.gff

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

my ($Verbose,$Help,$Outdir,$Min_cds_len,$Min_exon_num);
GetOptions(
	"min_cds_len:i"=>\$Min_cds_len,
	"min_exon_num:i"=>\$Min_exon_num,
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Min_cds_len ||= 150;
$Min_exon_num ||= 1;
$Outdir ||= "./";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $gff_file = shift;


$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

my %Info;
parse_gff($gff_file,\%Info);
##print Dumper \%Info;

open (IN,$gff_file) || die ("fail open $gff_file\n");
while (<IN>) {
	chomp;
	s/^\s+//;
	my @t = split(/\t/);
	my $qname = $2 if($t[8] =~ /^(ID|Parent)=([^;]+);*/);
	
	next if(!exists $Info{$qname} || $Info{$qname}{cds_len} < $Min_cds_len || $Info{$qname}{exon_num} < $Min_exon_num);

	print $_."\n";

}
close(IN);



##filter gff based on Min_cds_len and Min_exon_num
####################################################
sub parse_gff{
	my $file=shift;
	my $ref=shift;

	open (IN,$file) || die ("fail open $file\n");
	while (<IN>) {
		chomp;
		s/^\s+//;
		my @t = split(/\t/);
		next if(@t != 9);
		my $tname = $t[0];
		next if($t[2] ne 'CDS');
		my $qname = $1 if($t[8] =~ /^Parent=([^;]+);*/);
		my $start = $t[3];
		my $end = $t[4];

		$ref->{$qname}{cds_len} += $end - $start + 1;
		$ref->{$qname}{exon_num} ++;
	}
	close(IN);
	
}

