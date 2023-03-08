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

my $blast_p = ($seqType =~ /cds/i) ? "blastx" : "blastp";

my $config_file = "$Bin/../../config.txt";
my $RNAdenovo_anno = parse_config($config_file,"RNAdenovo_anno");
my $fastaDeal = parse_config($config_file,"fastaDeal.pl");
my $qsub_sge = parse_config($config_file,"qsub_sge.pl");
my $blast = parse_config($config_file,"blastall");
my $kegg_all = parse_config($config_file,"kegg_all");
my $kegg_animal = parse_config($config_file,"kegg_animal");
my $kegg_plant = parse_config($config_file,"kegg_plant");

my ($kegg_db, $kegg_path, $kegg_id, $kegg_komap);
if ($dbClass eq "animal"){ $kegg_db = $kegg_animal; $kegg_path = dirname($kegg_db); $kegg_id = "$kegg_path/animal.id.annot.xls"; $kegg_komap = "$kegg_path/komap/animal_ko_map.tab";}
if ($dbClass eq "plant"){ $kegg_db = $kegg_plant; $kegg_path = dirname($kegg_db); $kegg_id = "$kegg_path/plant.id.annot.xls"; $kegg_komap = "$kegg_path/komap/plant_ko_map.tab";}
if ($dbClass eq "all"){ $kegg_db = $kegg_all; $kegg_path = dirname($kegg_db); $kegg_id = "$kegg_path/kegg_all_clean.id.annot.xls"; $kegg_komap= "$kegg_path/komap/ko_map.tab";}

system("perl $fastaDeal -cuts $Cuts $fasta -outdir $workdir") == 0 || die $!;

my @subfiles = glob("$workdir/$filename.cut/*.*");
open SH, ">$workdir/$filename.blast.sh" || die $!;
foreach my $subfile (@subfiles){	
	print SH "$blast  -p $blast_p -e $Evalue -m 8 -F F -d $kegg_db -i $subfile -o $subfile.blast.kegg\n";
}
close SH;

system("perl $qsub_sge $QP_para --lines 1 --maxjob 80 --resource vf=2G  $filename.blast.sh") == 0 || die $!; #--reqsub

open SH, ">$workdir/blast.kegg.process.sh" || die $!;
print SH "#!/bin/bash\n";
print SH "cat $workdir/$filename.cut/*blast.kegg > $workdir/$filename.blast.kegg && \\\n";
print SH "perl $RNAdenovo_anno/get_annot_info.pl -tophit 5 -topmatch 1 -id $kegg_id -input $workdir/$filename.blast.kegg -out $workdir/$filename.blast.kegg.xls && \\\n";
print SH "perl $RNAdenovo_anno/blast2ko.pl -input $fasta -output $workdir/$filename.ko -blastout $workdir/$filename.blast.kegg -kegg $kegg_db && \\\n";
print SH "perl $RNAdenovo_anno/pathfind.pl -kegg $kegg_db -maptitle $kegg_path/map_title.tab -komap $kegg_komap -fg $workdir/$filename.ko -output $workdir/$filename.path && \\\n";
print SH "if [ -d $workdir/$filename\_map ]; then rm -rf $workdir/$filename\_map && \\\n";
print SH "mkdir -p $workdir/$filename\_map && \\\n";
print SH "perl $RNAdenovo_anno/keggMap_nodiff.pl -komap $kegg_komap -mapdir $kegg_path/map -ko $workdir/$filename.ko -outdir $workdir/$filename\_map && \\\n";
print SH "perl $RNAdenovo_anno/genPathHTML.pl -indir $workdir\n";
close SH;

system("cat $workdir/$filename.cut/*blast.kegg > $workdir/$filename.blast.kegg") == 0 || die $!;
system("perl $RNAdenovo_anno/get_annot_info.pl -tophit 5 -topmatch 1 -id $kegg_id -input $workdir/$filename.blast.kegg -out $workdir/$filename.blast.kegg.xls") ==0 || die $!;
system("perl $RNAdenovo_anno/blast2ko.pl -input $fasta -output $workdir/$filename.ko -blastout $workdir/$filename.blast.kegg -kegg $kegg_db") ==0 || die $!;
system("perl $RNAdenovo_anno/pathfind.pl -kegg $kegg_db -maptitle $kegg_path/map_title.tab -komap $kegg_komap -fg $workdir/$filename.ko -output $workdir/$filename.path") ==0 || die $!;
system("rm -rf $workdir/$filename\_map") ==0 || die $! if (-d "$workdir/$filename\_map");
system("mkdir -p $workdir/$filename\_map") == 0 || die $!;
system("perl $RNAdenovo_anno/keggMap_nodiff.pl -komap $kegg_komap -mapdir $kegg_path/map -ko $workdir/$filename.ko -outdir $workdir/$filename\_map") ==0 || die $!;
system("perl $RNAdenovo_anno/genPathHTML.pl -indir $workdir") ==0 || die $!;
system("perl $RNAdenovo_anno/drawKEGG.pl -path $workdir/$filename.path -outprefix $workdir/$filename -idCol 3 -level1Col 4 -level2Col 5 -geneCol 6") ==0 || die $!;



