#!/usr/bin/perl

=head1 Name

filter_on_genomic_level.pl  -- remove redundance of a block position file

=head1 Description

This pipeline takes a gff or table format file as input, remove the redundance and 
get a purified block set. 

The table format is:
gene	chr	strand	start	end

do hcluster on each chromosome, this can save memory and increase the speed significantly.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 3.0,  Date: 2008-6-23

=head1 Usage

  perl filter_on_genomic_level.pl [options] <mRNA.gff> <cds.len>
  --format <str>   gff or table, set the file format, default gff
  --outdir <str>   set the result directory, default="./" 
  --verbose   output running progress information to screen  
  --cutoff   default 100bp
  --help      output help information to screen  

=head1 Exmple

 perl ../bin/filter_on_genomic_level.pl rice_KOME_cDNA_5000.fa.blat.sim4.gff

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

my ($Format);
my ($Verbose,$Help,$Outdir);
my $Cutoff;
GetOptions(
	"format:s"=>\$Format,
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"cutoff:i"=>\$Cutoff,
	"help"=>\$Help
);
$Format ||= "gff";
$Outdir ||= "./";
$Cutoff ||=100;
die `pod2text $0` if (@ARGV == 0 || $Help);

my $gff_file = shift;
my $cds_len_file = shift;

my $gff_file_basename = basename($gff_file);

$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

#my %config;
#parse_config("$Bin/config.txt",\%config);

`perl $Bin/findOverlap_CalcuDist.pl --format $Format --length $Cutoff $gff_file > $Outdir/$gff_file_basename.dist`;

##cut the dist file into chromosomes and do hcluster for each chromosome
my %Content;
open IN, "$Outdir/$gff_file_basename.dist" || die "fail $Outdir/$gff_file_basename.dist";
while (<IN>) {
	my @t = split /\t/;
	my $chr = $t[4];
	push @{$Content{$chr}},$_;
}
close IN;


my $cut_dir = "$Outdir/$gff_file_basename.dist.cut";
`rm -r $cut_dir` if(-d $cut_dir);
`rm $Outdir/$gff_file_basename.dist.hcluster ` if(-f "$Outdir/$gff_file_basename.dist.hcluster");
`rm $Outdir/$gff_file_basename.dist.hcluster.redundance` if(-f "$Outdir/$gff_file_basename.dist.hcluster.redundance");
`rm $Outdir/$gff_file_basename.dist.hcluster.reference` if(-f "$Outdir/$gff_file_basename.dist.hcluster.reference");
`rm $Outdir/$gff_file_basename.dist.hcluster.all.list` if(-f "$Outdir/$gff_file_basename.dist.hcluster.all.list");

mkdir($cut_dir);

foreach my $chr (sort keys %Content) {
	my $chr_p =  $Content{$chr};
	my $dist_file = "$cut_dir/$gff_file_basename.dist.$chr";
	open OUT,">$dist_file" || die "fail";
	foreach my $line (@$chr_p) {
		print OUT $line;
	}
	close OUT;
	
	`perl $Bin/hcluster.pl -type min -stop 0.5 $dist_file > $dist_file.hcluster`;
	`perl $Bin/getWantedGene.pl  $dist_file $dist_file.hcluster $cds_len_file > $dist_file.hcluster.redundance`;
	`cat $dist_file.hcluster >> $Outdir/$gff_file_basename.dist.hcluster`;
	`cat $dist_file.hcluster.redundance >> $Outdir/$gff_file_basename.dist.hcluster.redundance`;
	`cat $dist_file.hcluster.reference >> $Outdir/$gff_file_basename.dist.hcluster.reference`;
	`cat $dist_file.hcluster.all.list >> $Outdir/$gff_file_basename.dist.hcluster.all.list`;
#	`rm $dist_file $dist_file.hcluster $dist_file.hcluster.redundance $dist_file.hcluster.reference $dist_file.hcluster.all.list `;
}


`perl $Bin/fishInWinter.pl -bf table -ff $Format -except  $Outdir/$gff_file_basename.dist.hcluster.redundance $gff_file > $Outdir/$gff_file_basename.nr.gff` if(-e "$Outdir/$gff_file_basename.dist.hcluster.redundance");
`cp $gff_file $Outdir/$gff_file_basename.nr.gff` unless(-e "$Outdir/$gff_file_basename.dist.hcluster.redundance");

`perl $Bin/fishInWinter.pl -bf table -ff $Format -except  $Outdir/$gff_file_basename.dist.hcluster.all.list $gff_file > $Outdir/$gff_file_basename.uncluster.gff` if(-e "$Outdir/$gff_file_basename.dist.hcluster.all.list");
`touch $Outdir/$gff_file_basename.uncluster.gff` unless(-e "$Outdir/$gff_file_basename.dist.hcluster.all.list");


####################################################
################### Sub Routines ###################
####################################################



##parse the software.config file, and check the existence of each software
####################################################
sub parse_config{
	my $conifg_file = shift;
	my $config_p = shift;
	
	my $error_status = 0;
	open IN,$conifg_file || die "fail open: $conifg_file";
	while (<IN>) {
		if (/(\S+)\s*=\s*(\S+)/) {
			my ($software_name,$software_address) = ($1,$2);
			$config_p->{$software_name} = $software_address;
			if (! -e $software_address){
				warn "Non-exist:  $software_name  $software_address\n"; 
				$error_status = 1;
			}
		}
	}
	close IN;
	die "\nExit due to error of software configuration\n" if($error_status);
}



## the one with most relations is the representative gene, the others are take as redundance.
sub get_redundance_from_hcluster {
	my $overlap_file = shift;
	my $hcluster_file = shift;
	my $redundance_file = shift;
	
	my $output;
	my %relation_number;

	open IN,$overlap_file || die "fail $overlap_file";
	while (<IN>) {
		my @t = split /\t/;
		next if($t[0] eq $t[1]);
		$relation_number{$t[0]}++;
	}
	close IN;

	open IN,$hcluster_file || die "fail $hcluster_file";
	while (<IN>) {
		my @t = split /\s+/;
		shift @t; ## remove Cluster_id:
		my %temp;
		foreach my $gene_id (@t) {
			$temp{$gene_id} = $relation_number{$gene_id};
		}
		my @sorted = sort {$temp{$b} <=> $temp{$a}} keys %temp;
		shift @sorted; ## remove the representative gene
		foreach  (@sorted) {
			$output .= $_."\n";
		}
	}
	close IN;

	open OUT, ">$redundance_file" || die "fail $redundance_file";
	print OUT $output;
	close OUT;
}
