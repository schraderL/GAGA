#!/usr/bin/perl

#Author: liuyabin@genomics.cn
#GO/data: 2016-05-06

use warnings;
use strict;
use File::Basename qw(fileparse basename dirname);
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use lib "$Bin/../../common_bin";
use GACP qw(parse_config);

die "perl $0 <gene.fasta> <queue> <pro_code> <animal|plant|all> <pep|cds> <cuts> <evalue> <workdir> <GO|NO>" if (@ARGV != 9);
my $fasta = shift;
my $Queue = shift;
my $Pro_code = shift;
my $dbClass = shift;
my $seqType = shift;
my $Cuts = shift;
my $Evalue = shift;
my $workdir = shift;
my $GO = shift;

$fasta = abs_path($fasta);
$workdir = abs_path($workdir);
my $filename = basename($fasta);

&mkdir_chdir("$workdir/NR");
my $time = `date`; chomp($time);
print "Start NR function annotation at ... $time!\nMake directory  $workdir/NR && cd  $workdir/NR ...\n";

my $QP_para;
(defined $Queue) ? $QP_para.="--queue $Queue " : die "Please set the queue for qsub!";
(defined $Pro_code) ? $QP_para.="--pro_code $Pro_code " : die "Please set the project code for qsub!";

my $blast_p = ($seqType =~ /cds/i) ? "blastx" : "blastp";

my $config_file = "$Bin/../../config.txt";
my $RNAdenovo_anno = parse_config($config_file,"RNAdenovo_anno");
my $fastaDeal = parse_config($config_file,"fastaDeal.pl");
my $qsub_sge = parse_config($config_file,"qsub_sge.pl");
my $blast = parse_config($config_file,"blastall");
my $nr = parse_config($config_file,"nr");
my $go = parse_config($config_file,"go");

my ($nr_db, $nr_vf);
if ($dbClass eq "animal"){ $nr_db = "$nr/animal.fa"; $nr_vf = "vf=10G";}
if ($dbClass eq "plant"){ $nr_db = "$nr/Plants.fa"; $nr_vf = "vf=5G";}
if ($dbClass eq "all"){ $nr_db = "$nr/nr"; $nr_vf = "vf=15G";}

system("perl $fastaDeal -cuts $Cuts $fasta -outdir $workdir/NR") == 0 || die $!;

my @subfiles = glob("$workdir/NR/$filename.cut/*.*");
open SH, ">$workdir/NR/$filename.blast.sh" || die $!;
foreach my $subfile (@subfiles){	
	print SH "$blast  -p $blast_p -e $Evalue -b 5 -v 5 -a 8 -m 7 -F F -d $nr_db -i $subfile -o $subfile.blast.nr\n";
}
close SH;

system("perl $qsub_sge $QP_para --lines 1 --maxjob 30 --resource $nr_vf  $filename.blast.sh") == 0 || die $!;

open SH, ">$workdir/NR/blast.nr.process.sh" || die $!;
print SH "#!/bin/bash\n";
print SH "cat $workdir/NR/$filename.cut/*blast.nr > $workdir/NR/$filename.blast.nr && \\ \n";
print SH "perl $RNAdenovo_anno/blast_m7_parser.pl $workdir/NR/$filename.blast.nr $workdir/NR/$filename.blast.nr.xls && \\ \n";
print SH "rm $workdir/NR/$filename.blast.nr && \\ \n";
print SH "perl $RNAdenovo_anno/blast_nr_class.pl -nr $workdir/NR/$filename.blast.nr.xls -outdir $workdir/NR\n";
close SH;

system("cat $workdir/NR/$filename.cut/*blast.nr > $workdir/NR/$filename.blast.nr") == 0 || die $!;
system("perl $RNAdenovo_anno/blast_m7_parser.pl $workdir/NR/$filename.blast.nr $workdir/NR/$filename.blast.nr.xls") == 0 || die $!;
system("rm $workdir/NR/$filename.blast.nr") == 0 || die $!;
system("perl $RNAdenovo_anno/blast_nr_class.pl -nr $workdir/NR/$filename.blast.nr.xls -outdir $workdir/NR") == 0 || die $!;

$time = `date`; chomp($time);
if ($GO =~ /GO/i){
	&mkdir_chdir("$workdir/GO");
	print "Start GO function annotation at $time!\nMake directory  $workdir/GO && cd  $workdir/GO ...\n";
	my @nr_blast = glob("$workdir/NR/$filename.cut/*.blast.nr");
	open SH, ">$workdir/GO/$filename.blast2go.sh" || die $!;
	system("mkdir $workdir/GO/blast2go") == 0 || die $!;
	foreach my $nr_blast (@nr_blast){
		my $blast_name = basename($nr_blast);
		print SH "$RNAdenovo_anno/../software/java -Xms1000m -Xmx5000m -cp $RNAdenovo_anno/../software/b2g4pipe/*:$RNAdenovo_anno/../software/b2g4pipe/ext/*: es.blast2go.prog.B2GAnnotPipe -prop $RNAdenovo_anno/../software/b2g4pipe/b2gPipe.properties -in $nr_blast -out $workdir/GO/blast2go/$blast_name -v -annot\n";
	}
	system("perl $qsub_sge $QP_para --lines 1 --maxjob 30 --resource vf=15G  $filename.blast2go.sh") == 0 || die $!;

	open SH, ">$workdir/GO/blast2go.process.sh" || die $!;
	print SH "#!/bin/bash\n";
	print SH "cat $workdir/GO/blast2go/*blast.nr.annot > $workdir/GO/$filename.blast2go.annot && \\\n";
	print SH "if [ -d $workdir/GO/data ];then rm -rf $workdir/GO/data;fi && \\\n";
	print SH "mkdir -p $workdir/GO/data && \\\n";
	print SH "perl $RNAdenovo_anno/annot2goa.pl $go/gene_ontology.1_2.obo $workdir/GO/$filename.blast2go.annot $workdir/GO/data/species && \\\n";
	print SH "perl $RNAdenovo_anno/blast_m7_m8.pl -input \"$workdir/NR/$filename.cut/*.blast.nr\" -output $workdir/GO/data/species.nr.m8 && \\\n";
	print SH "perl $RNAdenovo_anno/getNrDesc.pl -input $workdir/GO/data/species.nr.m8 -rank 1 -nr $nr_db -output $workdir/GO/data/species.nr.desc && \\\n";
	print SH "rm $workdir/GO/data/species.nr.m8 && \\\n";
	print SH "cat $workdir/GO/data/species.[CFP] | cut -f2 | sort -u > $workdir/GO/gene.list && \\\n";
	print SH "perl $RNAdenovo_anno/drawGO.pl -list $workdir/GO/gene.list -goprefix  $workdir/GO/data/species -goclass $go/go.class -outprefix $workdir/GO/$filename\n";

	system("cat $workdir/GO/blast2go/*blast.nr.annot > $workdir/GO/$filename.blast2go.annot") == 0 || die $!;
	system("rm -rf $workdir/GO/data") == 0 || die $! if (-d "$workdir/GO/data");
	system("mkdir -p $workdir/GO/data") == 0 || die $!;
	system("perl $RNAdenovo_anno/annot2goa.pl $go/gene_ontology.1_2.obo $workdir/GO/$filename.blast2go.annot $workdir/GO/data/species") == 0 || die $!;
	system("perl $RNAdenovo_anno/blast_m7_m8.pl -input \"$workdir/NR/$filename.cut/*.blast.nr\" -output $workdir/GO/data/species.nr.m8") == 0 || die $!;
	system("perl $RNAdenovo_anno/getNrDesc.pl -input $workdir/GO/data/species.nr.m8 -rank 1 -nr $nr_db -output $workdir/GO/data/species.nr.desc") == 0 || die $!;
	system("rm $workdir/GO/data/species.nr.m8") == 0 || die $!;
	system("cat $workdir/GO/data/species.[CFP] | cut -f2 | sort -u > $workdir/GO/gene.list") == 0 || die $!;
	system("perl $RNAdenovo_anno/drawGO.pl -list $workdir/GO/gene.list -goprefix  $workdir/GO/data/species -goclass $go/go.class -outprefix $workdir/GO/$filename") == 0 || die $!;
	$time = `date`; chomp($time);
	print "NR && Go annotation DONE ... $time!\n";

}else{
	print "NR annotation DONE ... $time!\n";
}

#####################################sub#################################
sub mkdir_chdir{
	my $dir_name = shift;
	system("rm -rf $dir_name") == 0 || die $! if (-d $dir_name);
	system("mkdir $dir_name") == 0 || die $!;
	chdir $dir_name or die ("ERROR: Cannot change directory to $dir_name: $!\n");
}
