#!/usr/bin/perl

=head1 Name
	perfect_gene.pl
	liushiping@genomics.cn

=head1 Command-line Option
	--sco	set the score of genes.dafault frac=95
	--start	set the max numbers of the amino acid extending the start.default star=10
	--stop	set the max numbers of the amino acid extending the end.default end=10
	--outdir
	--help
	
=head1 Usage
	perl perfect_gene.pl [--sco] [--start] [--stop] [--outdir] <genome_fa> <gff> [<gff> ...] 

=head1 Example
	perl perfect_gene.pl --sco 99 --start 5 --stop 5 seqence.fa genewise1.gff genewise2.gff genewise3.gff

=cut



use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

my ($Score,$Start,$Stop,$Outdir,$Help);
GetOptions(
	"sco:n"=>\$Score,
	"start:s"=>\$Start,
	"stop:s"=>\$Stop,
	"outdir:s"=>\$Outdir,
	"help"=>\$Help
);
$Score ||= 95;
$Start ||= 10;
$Stop ||= 10;
$Outdir ||= "./";
$Outdir =~ s/\/$//;

die `pod2text $0` if ($Help || @ARGV < 2);

my ($file_genome,@file_gff)=@ARGV;
my $file_genome_name=basename($file_genome);

my $filter_gff		="$Bin/filter_gff_percent-lenght.pl";
my $mult_exon		="$Bin/mult-exon_err.pl";
my $fishInWinter	="$Bin/fishInWinter.pl";
my $add_start_stop	="$Bin/add_start_stop.pl";
my $clustergff		="$Bin/clustergff.pl.new.pl";

for(my $i=0;$i < @file_gff;$i++){
	my $file_gff_sub=$file_gff[$i];
	my $file=basename($file_gff_sub);
	`perl $filter_gff $file_gff_sub $Score > $Outdir/$file.$Score.gff`;
	`perl $mult_exon $Outdir/$file.$Score.gff $file_genome > $Outdir/$file.$Score.stop`;
	`perl $fishInWinter --bf gff --ff gff --except $Outdir/$file.$Score.stop $Outdir/$file.$Score.gff > $Outdir/$file.$Score.stop.un`;
	`perl $add_start_stop $file_genome $Outdir/$file.$Score.stop.un $Start $Stop > $Outdir/$file.$Score.stop.un.$Start\_$Stop.gff`;
}

`cat $Outdir/*$Start\_$Stop.gff > $Outdir/$file_genome_name.gff`;
`perl $clustergff $Outdir/$file_genome_name.gff`;
