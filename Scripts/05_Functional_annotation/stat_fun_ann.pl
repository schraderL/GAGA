#!/usr/bin/perl
=head1 Name

	stat_fun_ann.pl -- just use for stat function annotation result.

=head1 Version

	Author: zhouheling (zhouheling@genomics.org.cn)
	Version: 1.0    Date: 2010-08-2

=head1 Usage

	perl stat_fun_ann.pl [options] *.pep

	-Interpro		stat Interpro result
	-KEGG			stat KEGG result
	-Swissprot		stat Swissprot result
	-TrEMBL			stat TrEMBL result

	-verbose		output verbose information to screen
	-help			output help information to screen

=head1 Exmple

	perl /nas/GAG_02/zhouheling/GACP-8.0/02.function_annotation/auto_fun_ann/bin/stat_fun_ann.pl -Interpro -KEGG -Swissprot -TrEMBL *.pep

=cut


use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use lib "$Bin/../../common_bin";
use GACP qw(parse_config);

my ($Interpro,$KEGG,$Swissprot,$TrEMBL,$verbose,$help);
my ($Total,$Annotated,$Swissprot_num,$TrEMBL_num,$KEGG_num,$InterPro_num,$GO_num,$Unanotated) =(0,0,0,0,0,0,0);

GetOptions(
	"Interpro"=>\$Interpro,
	"KEGG"=>\$KEGG,
	"Swissprot"=>\$Swissprot,
	"TrEMBL"=>\$TrEMBL,
	"verbose"=>\$verbose,
	"help"=>\$help,
);

my $pep_file = shift;
my $pep_name = basename($pep_file);
my $species_name=$1 if ($pep_name=~/^([\w-]+)\./);

my $config_file = "$Bin/../../config.txt";
my $gene_function_number = parse_config($config_file,"get_gene_function_number");

open(OUT , ">" .  "$species_name\.function.statistics.xls") or die $!;

print OUT "\tNumber\tPercent(%)\n";

########## stat all gene ##########
$Total = `grep -c ">" $pep_name`;
chomp($Total);
print OUT "Total\t$Total\t\n";

########## stat Interpro result ##########
if (defined $Interpro)
{
	$InterPro_num = $1 if (`wc -l InterPro.result` =~ /(\d+)/); chomp($InterPro_num);
	$GO_num = $1 if (`wc -l GO.result` =~ /(\d+)/); chomp($GO_num);
	print OUT "InterPro\t$InterPro_num\t";
	printf OUT "%.6f\n",$InterPro_num / $Total * 100;
	print OUT "GO\t$GO_num\t";
	printf OUT "%.6f\n",$GO_num / $Total * 100;
}
########## stat KEGG result ##########
if (defined $KEGG)
{
	$KEGG_num = $1 if (`wc -l KEGG.result` =~ /(\d+)/); chomp($KEGG_num);
	$KEGG_num = $KEGG_num;
	print OUT "KEGG\t$KEGG_num\t";
	printf OUT "%.6f\n",$KEGG_num / $Total * 100;
}
########## stat Swissprot result ##########
if (defined $Swissprot)
{
	$Swissprot_num = $1 if (`wc -l Swissprot.result` =~ /(\d+)/); chomp($Swissprot_num);
	$Swissprot_num = $Swissprot_num;
	print OUT "Swissprot\t$Swissprot_num\t";
	printf OUT "%.6f\n",$Swissprot_num / $Total * 100;
}
########## stat TrEMBL result ##########
if (defined $TrEMBL)
{
	$TrEMBL_num = $1 if (`wc -l TrEMBL.result` =~ /(\d+)/); chomp($TrEMBL_num);
	$TrEMBL_num = $TrEMBL_num;
	print OUT "TrEMBL\t$TrEMBL_num\t";
	printf OUT "%.6f\n",$TrEMBL_num / $Total * 100;
}
########## stat Annotated ##########
$Annotated = `perl $gene_function_number GO.result InterPro.result KEGG.result Swissprot.result TrEMBL.result`;
chomp($Annotated);
print OUT "Annotated\t$Annotated\t";
printf OUT "%.6f\n",$Annotated / $Total * 100;
########## stat Unannotated ##########
$Unanotated = $Total - $Annotated;
print OUT "Unanotated\t$Unanotated\t";
printf OUT "%.6f\n",$Unanotated / $Total * 100;
########## finish ##########
close(OUT);
