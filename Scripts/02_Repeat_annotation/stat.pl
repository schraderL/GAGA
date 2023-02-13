#!/usr/bin/perl
=head1 Name

	stat.pl -- just use for stat repeat result.

=head1 Version

	Version: 1.0   

=head1 Usage

	perl stat.pl [options] seq_file

	-denovo			stat denovo result
	-trf			stat trf result
	-repeatmasker		stat repeatmasker result
	-proteinmask		stat proteinmask result

	-verbose		output verbose information to screen
	-help			output help information to screen

=head1 Exmple

	perl /nas/GAG_02/zhouheling/GACP-8.0/03.repeat_finding/auto_repeat/bin/stat.pl -denovo -trf -repeatmasker -proteinmask  Bm_1214.scafSeq.FG.fa
	perl /nas/GAG_02/zhouheling/GACP-8.0/03.repeat_finding/auto_repeat/bin/stat.pl -denovo -proteinmask  Bm_1214.scafSeq.FG.fa

=cut


use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use lib "$Bin/";
use GACP qw(parse_config);

my ($denovo,$trf,$repeatmasker,$proteinmask,$verbose,$help);
my ($total_len,$all_len,$all_out_trf_len,$denovo_len,$trf_len,$repeatmasker_len,$proteinmask_len,$tmp,$line) =(0,0,0,0,0,0,0,0,0);
my ($all_out_trf_DNA,$all_out_trf_LINE,$all_out_trf_SINE,$all_out_trf_LTR,$all_out_trf_Other,$all_out_trf_Satellite) = (0,0,0,0,0,0);
my ($all_out_trf_Simple_repeat,$all_out_trf_Unknown) = (0,0);
my ($denovo_DNA,$denovo_LINE,$denovo_SINE,$denovo_LTR,$denovo_Other,$denovo_Satellite,$denovo_Simple_repeat,$denovo_Unknown) = (0,0,0,0,0,0,0,0);
my ($rm_DNA,$rm_LINE,$rm_SINE,$rm_LTR,$rm_Other,$rm_Satellite,$rm_Simple_repeat,$rm_Unknown) = (0,0,0,0,0,0,0,0);
my ($pm_DNA,$pm_LINE,$pm_SINE,$pm_LTR,$pm_Other,$pm_Satellite,$pm_Simple_repeat,$pm_Unknown) = (0,0,0,0,0,0,0,0);
GetOptions(
	"denovo"=>\$denovo,
	"trf"=>\$trf,
	"repeatmasker"=>\$repeatmasker,
	"proteinmask"=>\$proteinmask,
	"verbose"=>\$verbose,
	"help"=>\$help,
);

my $seq_file = shift;
my $species_name=$1 if ($seq_file=~/^(\S+)\.fasta/);
my $config_file = "$Bin/../config.txt";
my $common_bin = "$Bin/../common_bin";
my $stat_TE = parse_config($config_file,"stat_TE2_path");
my $fastaDeal = parse_config($config_file,"fastaDeal_path");
my $svg_program = parse_config($config_file,"svg_path");
my $png_program = parse_config($config_file,"png_path");

open(OUT , ">" .  "$species_name.repeat.statistics.xls") or die $!;

########## stat all repeat ##########
`perl $fastaDeal -attribute id:len:lenwogap:gc $seq_file > $seq_file.len`;
$total_len = `awk '{(tot+=\$2)};END{print tot}' $seq_file.len`; chomp ($total_len);
print OUT "Total Gene Lenght = $total_len\n\n";

print OUT "Type\tRepeat Size\t% of genome\n";
$all_len = $1 if (`perl $stat_TE -gff all.gff -rank all` =~ /\t(\d+)/);

$all_out_trf_len = $1 if (`perl $stat_TE -gff all_without_trf.gff -rank all` =~ /\t(\d+)/);
`perl $stat_TE -gff all_without_trf.gff -rank type > all_without_trf.gff.type.tmp`;
open (IN,"all_without_trf.gff.type.tmp") or die $!;
while ($line = <IN>)
{
	chomp($line);
	$all_out_trf_DNA = $1 if ($line =~ /DNA\t(\d+)/);
	$all_out_trf_LINE = $1 if ($line =~ /LINE\t(\d+)/);
	$all_out_trf_SINE = $1 if ($line =~ /SINE\t(\d+)/);
	$all_out_trf_LTR = $1 if ($line =~ /LTR\t(\d+)/);
	$all_out_trf_Other = $1 if ($line =~ /Other\t(\d+)/);
	$all_out_trf_Unknown = $1 if ($line =~ /Unknown\t(\d+)/);
}
close(IN);

########## stat trf result ##########
if (defined $trf)
{
        $trf_len = $1 if (`perl $stat_TE -gff trf.gff -rank all` =~ /\t(\d+)/);
        print OUT "Trf\t$trf_len\t"; printf OUT "%.6f\n",$trf_len / $total_len * 100;
}

########## stat repeatmasker result ##########
if (defined $repeatmasker)
{
        $repeatmasker_len = $1 if (`perl $stat_TE -gff repeatmasker.gff -rank all` =~ /\t(\d+)/);
        print OUT "Repeatmasker\t$repeatmasker_len\t"; printf OUT "%.6f\n",$repeatmasker_len / $total_len * 100;
        `perl $stat_TE -gff repeatmasker.gff -rank type > repeatmasker.gff.type.tmp`;
	open (IN,"repeatmasker.gff.type.tmp") or die $!;
	while ($line = <IN>)
	{
		chomp($line);
		$rm_DNA = $1 if ($line =~ /DNA\t(\d+)/);
		$rm_LINE = $1 if ($line =~ /LINE\t(\d+)/);
		$rm_SINE = $1 if ($line =~ /SINE\t(\d+)/);
		$rm_LTR = $1 if ($line =~ /LTR\t(\d+)/);
		$rm_Other = $1 if ($line =~ /Other\t(\d+)/);
		$rm_Unknown = $1 if ($line =~ /Unknown\t(\d+)/);
	}
	close(IN);
}

########## stat proteinmask result ##########
if (defined $proteinmask)
{
        $proteinmask_len = $1 if (`perl $stat_TE -gff proteinmask.gff -rank all` =~ /\t(\d+)/);
        print OUT "Proteinmask\t$proteinmask_len\t"; printf OUT "%.6f\n",$proteinmask_len / $total_len * 100;
        `perl $stat_TE -gff proteinmask.gff -rank type > proteinmask.gff.type.tmp`;
	open (IN,"proteinmask.gff.type.tmp") or die $!;
	while ($line = <IN>)
	{
		chomp($line);
		$pm_DNA = $1 if ($line =~ /DNA\t(\d+)/);
		$pm_LINE = $1 if ($line =~ /LINE\t(\d+)/);
		$pm_SINE = $1 if ($line =~ /SINE\t(\d+)/);
		$pm_LTR = $1 if ($line =~ /LTR\t(\d+)/);
		$pm_Other = $1 if ($line =~ /Other\t(\d+)/);
		$pm_Unknown = $1 if ($line =~ /Unknown\t(\d+)/);
	}
}

########## stat denovo result ##########
if (defined $denovo)
{
	$denovo_len = $1 if (`perl $stat_TE -gff denovo.gff -rank all` =~ /\t(\d+)/);
	print OUT "De novo\t$denovo_len\t"; printf OUT "%.6f\n",$denovo_len / $total_len * 100;
	`perl $stat_TE -gff denovo.gff -rank type > denovo.gff.type.tmp`;
	open (IN,"denovo.gff.type.tmp") or die $!;
	while ($line = <IN>)
	{
		chomp($line);
		$denovo_DNA = $1 if ($line =~ /DNA\t(\d+)/);
		$denovo_LINE = $1 if ($line =~ /LINE\t(\d+)/);
		$denovo_SINE = $1 if ($line =~ /SINE\t(\d+)/);
		$denovo_LTR = $1 if ($line =~ /LTR\t(\d+)/);
		$denovo_Other = $1 if ($line =~ /Other\t(\d+)/);
		$denovo_Satellite = $1 if ($line =~ /Satellite\t(\d+)/);
		$denovo_Simple_repeat = $1 if ($line =~ /Simple_repeat\t(\d+)/);
		$denovo_Unknown = $1 if ($line =~ /Unknown\t(\d+)/);
	}
	close(IN);
}
#print "$total_len\t$all_len\t$all_with_trf_len\t$denovo_len\t$trf_len\t$repeatmasker_len\t$proteinmask_len\n";

print OUT "Total\t$all_len\t"; printf OUT "%.6f\n",$all_len / $total_len * 100;
print OUT "\n";

########## print denovo result ##########
print OUT "De novo\n";
print OUT "Type\tLength (Bp)\t% in genome\n";
print OUT "DNA\t$denovo_DNA\t"; printf OUT "%.6f\n",$denovo_DNA / $total_len * 100;
print OUT "LINE\t$denovo_LINE\t"; printf OUT "%.6f\n",$denovo_LINE / $total_len * 100;
print OUT "SINE\t$denovo_SINE\t"; printf OUT "%.6f\n",$denovo_SINE / $total_len * 100;
print OUT "LTR\t$denovo_LTR\t"; printf OUT "%.6f\n",$denovo_LTR / $total_len * 100;
print OUT "Other\t$denovo_Other\t"; printf OUT "%.6f\n",$denovo_Other / $total_len * 100;
print OUT "Satellite\t$denovo_Satellite\t"; printf OUT "%.6f\n",$denovo_Satellite / $total_len * 100;
print OUT "Simple_repeat\t$denovo_Simple_repeat\t"; printf OUT "%.6f\n",$denovo_Simple_repeat / $total_len * 100;
print OUT "Unknown\t$denovo_Unknown\t"; printf OUT "%.6f\n",$denovo_Unknown / $total_len * 100;
print OUT "Total\t$denovo_len\t"; printf OUT "%.6f\n",$denovo_len / $total_len * 100;
print OUT "\n";

########## print all result ##########
print OUT "\tRepbase TEs\t\tTE protiens\t\tDe novo\t\tCombined TEs\t\n";
print OUT "Type\tLength (Bp)\t% in genome\tLength (Bp)\t% in genome\tLength (Bp)\t% in genome\tLength (Bp)\t% in genome\n";
### DNA ###
print OUT "DNA\t$rm_DNA\t";printf OUT "%.6f\t",$rm_DNA / $total_len * 100;
print OUT "$pm_DNA\t";printf OUT "%.6f\t",$pm_DNA / $total_len * 100;
print OUT "$denovo_DNA\t";printf OUT "%.6f\t",$denovo_DNA / $total_len * 100;
print OUT "$all_out_trf_DNA\t";printf OUT "%.6f\n",$all_out_trf_DNA / $total_len * 100;
### LINE ###
print OUT "LINE\t$rm_LINE\t";printf OUT "%.6f\t",$rm_LINE / $total_len * 100;
print OUT "$pm_LINE\t";printf OUT "%.6f\t",$pm_LINE / $total_len * 100;
print OUT "$denovo_LINE\t";printf OUT "%.6f\t",$denovo_LINE / $total_len * 100;
print OUT "$all_out_trf_LINE\t";printf OUT "%.6f\n",$all_out_trf_LINE / $total_len * 100;
### SINE ###
print OUT "SINE\t$rm_SINE\t";printf OUT "%.6f\t",$rm_SINE / $total_len * 100;
print OUT "$pm_SINE\t";printf OUT "%.6f\t",$pm_SINE / $total_len * 100;
print OUT "$denovo_SINE\t";printf OUT "%.6f\t",$denovo_SINE / $total_len * 100;
print OUT "$all_out_trf_SINE\t";printf OUT "%.6f\n",$all_out_trf_SINE / $total_len * 100;
### LTR ###
print OUT "LTR\t$rm_LTR\t";printf OUT "%.6f\t",$rm_LTR / $total_len * 100;
print OUT "$pm_LTR\t";printf OUT "%.6f\t",$pm_LTR / $total_len * 100;
print OUT "$denovo_LTR\t";printf OUT "%.6f\t",$denovo_LTR / $total_len * 100;
print OUT "$all_out_trf_LTR\t";printf OUT "%.6f\n",$all_out_trf_LTR / $total_len * 100;
### Other ###
print OUT "Other\t$rm_Other\t";printf OUT "%.6f\t",$rm_Other / $total_len * 100;
print OUT "$pm_Other\t";printf OUT "%.6f\t",$pm_Other / $total_len * 100;
print OUT "$denovo_Other\t";printf OUT "%.6f\t",$denovo_Other / $total_len * 100;
print OUT "$all_out_trf_Other\t";printf OUT "%.6f\n",$all_out_trf_Other / $total_len * 100;
### Unknown ###
print OUT "Unknown\t$rm_Unknown\t";printf OUT "%.6f\t",$rm_Unknown / $total_len * 100;
print OUT "$pm_Unknown\t";printf OUT "%.6f\t",$pm_Unknown / $total_len * 100;
print OUT "$denovo_Unknown\t";printf OUT "%.6f\t",$denovo_Unknown / $total_len * 100;
print OUT "$all_out_trf_Unknown\t";printf OUT "%.6f\n",$all_out_trf_Unknown / $total_len * 100;
### Total ###
print OUT "Total\t$repeatmasker_len\t";printf OUT "%.6f\t",$repeatmasker_len / $total_len * 100;
print OUT "$proteinmask_len\t";printf OUT "%.6f\t",$proteinmask_len / $total_len * 100;

$tmp = $denovo_len - $denovo_Satellite - $denovo_Simple_repeat;
print OUT "$tmp\t";printf OUT "%.6f\t",$tmp / $total_len * 100;

$tmp = $all_out_trf_len - $denovo_Satellite - $denovo_Simple_repeat;
print OUT "$tmp \t";printf OUT "%.6f\n",$tmp / $total_len * 100;
########## Finish ##########
close(OUT);

########## draw ##########
`rm *.tmp`;
if (defined $denovo)
{
	`perl $svg_program denovo.out $total_len -X_step 10 -X_end 40 -Y_step 0.5 -Y_end 2`;
	`mv denovo.out.TEdivergence.svg $species_name.repeat.denovo.sequence.divergence.svg`;
	`perl $png_program -type png $species_name.repeat.denovo.sequence.divergence.svg`;
}
if (defined $repeatmasker)
{
	`perl $svg_program repeatmasker.out $total_len -X_step 10 -X_end 40 -Y_step 0.5 -Y_end 2`;
	`mv repeatmasker.out.TEdivergence.svg $species_name.repeat.known.sequence.divergence.svg`;
	`perl $png_program -type png $species_name.repeat.known.sequence.divergence.svg`;
}
