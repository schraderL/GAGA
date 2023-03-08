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

die "perl $0 <gene.fasta> <queue> <pro_code> <animal|plant|all> <pep|cds> <cuts> <evalue> <workdir> " if (@ARGV != 8);
my $fasta = shift;
my $Queue = shift;
my $Pro_code = shift;
my $dbClass = shift;
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

my $blast_p = ($seqType =~ /cds/i) ? "blastn" : "tblastn";

my $config_file = "$Bin/../../config.txt";
my $RNAdenovo_anno = parse_config($config_file,"RNAdenovo_anno");
my $fastaDeal = parse_config($config_file,"fastaDeal.pl");
my $qsub_sge = parse_config($config_file,"qsub_sge.pl");
my $blast = parse_config($config_file,"blastall");
my $nt = parse_config($config_file,"nt");

my ($nt_db, $nt_vf);
if ($dbClass eq "animal"){ $nt_db = "$nt/animal.fa"; $nt_vf = "vf=15G";}
if ($dbClass eq "plant"){ $nt_db = "$nt/Plants.and.Fungi.fa"; $nt_vf = "vf=7G";}
if ($dbClass eq "all"){ $nt_db = "$nt/nt"; $nt_vf = "vf=20G";}

system("perl $fastaDeal -cuts $Cuts $fasta -outdir $workdir") == 0 || die $!;

my @subfiles = glob("$workdir/$filename.cut/*.*");
open SH, ">$workdir/$filename.blast.sh" || die $!;
foreach my $subfile (@subfiles){	
	print SH "$blast  -p $blast_p -e $Evalue -a 8 -F F -d $nt_db -i $subfile -o $subfile.blast.nt\n";
}
close SH;

system("perl $qsub_sge $QP_para  --lines 1 --maxjob 30 --resource $nt_vf  $filename.blast.sh") == 0 || die $!;

open SH, ">$workdir/blast.nt.process.sh" || die $!;
print SH "#!/bin/bash\n";
print SH "cat $workdir/$filename.cut/*blast.nt > $workdir/$filename.blast.nt && \\ \n";
print SH "perl $RNAdenovo_anno/blast_parser.pl -tophit 5 -topmatch 1 $workdir/$filename.blast.nt >  $workdir/$filename.blast.nt.xls && \\ \n";
print SH "rm $workdir/$filename.blast.nt\n";
close SH;

system("cat $workdir/$filename.cut/*blast.nt > $workdir/$filename.blast.nt") == 0 || die $!;
system("perl $RNAdenovo_anno/blast_parser.pl -tophit 5 -topmatch 1 $workdir/$filename.blast.nt >  $workdir/$filename.blast.nt.xls") == 0 || die $!;
system("rm $workdir/$filename.blast.nt") == 0 || die $!;

