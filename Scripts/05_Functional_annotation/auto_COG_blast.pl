#!/usr/bin/perl

#Author: liuyabin@genomics.cn
#Data: 2016-05-06

use warnings;
use strict;
use File::Basename qw(fileparse basename dirname);
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use lib "$Bin/../../common_bin";
use GACP qw(parse_config);

die "perl $0 <gene.fasta> <queue> <pro_code> <pep|cds> <cuts> <evalue> <workdir> " if (@ARGV != 7);
my $fasta = shift;
my $Queue = shift;
my $Pro_code = shift;
my $seqType = shift;
my $Cuts = shift;
my $Evalue = shift;
my $workdir = shift;

$fasta = abs_path($fasta);
$workdir = abs_path($workdir);
my $filename = basename($fasta);

my $QP_para;
(defined $Queue) ? $QP_para.="--queue $Queue " : die "Please set the queue for qsub!";
(defined $Pro_code) ? $QP_para.="--pro_code $Pro_code " : die "Please set the project code for qsub!";

my $blast_p = ($seqType =~ /cds/i) ? "blastx" : "blastp";

my $config_file = "$Bin/../../config.txt";
my $RNAdenovo_anno = parse_config($config_file,"RNAdenovo_anno");
my $fastaDeal = parse_config($config_file,"fastaDeal.pl");
my $qsub_sge = parse_config($config_file,"qsub_sge.pl");
my $blast = parse_config($config_file,"blastall");
my $cog = parse_config($config_file,"cog");

my $cog_path = dirname($cog);

system("perl $fastaDeal -cuts $Cuts $fasta -outdir $workdir") == 0 || die $!;

my @subfiles = glob("$workdir/$filename.cut/*.*");
open SH, ">$workdir/$filename.blast.sh" || die $!;
foreach my $subfile (@subfiles){	
	print SH "$blast  -p $blast_p -e $Evalue -a 5 -m 8 -F F -d $cog -i $subfile -o $subfile.blast.cog\n";
}
close SH;

system("perl $qsub_sge $QP_para  --lines 1 --maxjob 30 --resource vf=2G  $filename.blast.sh") == 0 || die $!;

open SH, ">$workdir/blast.cog.process.sh" || die $!;
print SH "#!/bin/bash\n";
print SH "cat $workdir/$filename.cut/*blast.cog > $workdir/$filename.blast.cog && \\\n";
print SH "perl $RNAdenovo_anno/get_annot_info.pl -tophit 5 -topmatch 1 -id $cog_path/cog_clean.fa.id -input $workdir/$filename.blast.cog -out $workdir/$filename.blast.cog.xls && \\\n";
print SH "perl $RNAdenovo_anno/cog_parser.pl $cog_path/whog $cog_path/fun.txt $workdir/$filename.blast.cog.xls && \\\n";
print SH "perl $RNAdenovo_anno/cog_R.pl -catalog $workdir/$filename.COG2Gene.xls -sample $filename -outdir  $workdir\n";
close SH;


system("cat $workdir/$filename.cut/*blast.cog > $workdir/$filename.blast.cog") == 0 || die $!;
system("perl $RNAdenovo_anno/get_annot_info.pl -tophit 5 -topmatch 1 -id $cog_path/cog_clean.fa.id -input $workdir/$filename.blast.cog -out $workdir/$filename.blast.cog.xls") == 0 || die $!;
system("perl $RNAdenovo_anno/cog_parser.pl $cog_path/whog $cog_path/fun.txt $workdir/$filename.blast.cog.xls") == 0 || die $!;
system("perl $RNAdenovo_anno/cog_R.pl -catalog $workdir/$filename.COG2Gene.xls -sample $filename -outdir  $workdir") == 0 || die $!;

