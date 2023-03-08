#!/usr/bin/perl

=head1 Name

	auto_fun_ann.pl  -- the pipeline of function annotation.

=head1 Description

	Base on run_iprscan.pl and blast_database.pl.This program invoke Interpro, KEGG, Swissprot and Trembl.	
	The path of the softwares invoked will be put in the "config.txt".

	Use multiple -appl flags to specify multiple applications. The possible applications
	are listed below:
            ProDom
            PRINTS
            Hamap
            Pfam
            PIRSF
            PANTHER
            TIGRFAM
            SMART
            SUPERFAMILY
            Gene3D
            ProSitePatterns
            ProSiteProfiles
            Coils
            SignalP_EUK
            SignalP_GRAM_POSITIVE
            SignalP_GRAM_NEGATIVE
            Phobius
            TMHMM

=head1 Version
        Author: zhouheling (zhouheling@genomics.org.cn)
        Mender: luchangxin (luchangxin@genomics.org.cn)
        Mender: helijuan (Lyndi.He@genomics.cn)
		Menfer: liuyabin (liuyabin@genomics.cn)
        Version: 1.0    Date: 2010-08-01
        Updata: 2.0     Data: 2015-05-21
        Updata: 2016a     Data: 2016-01-21
		Updata: 2016a     Data: 2016-05-06
        Note:
        (1)Added kegg mapping programm, and could search against kegg database for animal or plant separately.
        (2)Added lines option.

=head1 Usage

	perl auto_fun_ann.pl [options] *.pep

	note: the pep_file must be input by absolute path

	-Interpro		run Interpro
	  -appl <str>		set applications, this option can be used multiple times
	  -iresource <str>	set the required resource used in qsub -l option, default vf=10G #16G
	-KEGG			run KEGG
	-Swissprot		run Swissprot
	-Trembl			run Trembl
	-COG			run COG
	-NT				run NT
	-NR				run NR
	-GO				run GO
	-species <str>		 set the species type, animal or plant; if not set, will compare to all type.
	-seqtype <str>  set the sequence type, pep or cds, default=pep.
	-cuts <int>		cut file with the specified number of seqeunces in each sub file for Interpro,default=50
	-cpu <int>		set the cpu number to use in ipr(former: set th cpu number to use in parallel), default=5
	-Evalue <str>   set the blast evalue, default=1e-5.
	-queue <str>            set the queue in qsub, default no
	-pro_code <str>         set the code for project,default no
	-help			output help information to screen

=head1 Exmple

nohup perl auto_fun_ann.pl -queue bc.q -pro_code APDEnre -Interpro -KEGG -Swissprot -Trembl -COG -NT -NR -GO -species animal -seqtype pep pep.fasta 

=cut

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use lib "$Bin/../../common_bin";
use GACP qw(parse_config);

my ($Interpro,$KEGG,$Swissprot,$Trembl,$COG, $NT, $NR, $GO, $appl_parameter,@Appl,$iResource,$Queue,$Pro_code);
my ($species,$seqtype, $Cutf,$Cuts,$Cpu, $evalue);
my ($help);

GetOptions(
	"Interpro"=>\$Interpro,
		"appl:s"=>\@Appl,
		"iresource:s"=>\$iResource,
	"KEGG"=>\$KEGG,
	"Swissprot"=>\$Swissprot,
	"Trembl"=>\$Trembl,
	"COG" =>\$COG,
	"NT" =>\$NT,
	"NR" =>\$NR,
	"GO" =>\$GO,
	"species:s"=>\$species,
	"seqtype:s"=>\$seqtype,
	"cuts:i"=>\$Cuts,
	"cpu:i"=>\$Cpu,
	"Evalue"=>\$evalue,
	"queue:s"=>\$Queue,
	"pro_code:s"=>\$Pro_code,	
	"Help"=>\$help,
);

die `pod2text $0` if(@ARGV==0 || $help);

$Cuts ||= 100;
$Cpu ||= 5;
$iResource ||= "vf=10G"; #16G
$evalue ||= 1e-5;
$species ||="all";
$seqtype ||="pep";

if (@Appl == 0)
{
	$appl_parameter = "--appl PRINTS,Pfam,SMART,PANTHER,ProSiteProfiles,ProSitePatterns,CDD,SFLD,Gene3D,SUPERFAMILY,TMHMM ";
}
else
{
    $appl_parameter = "-appl ";
	foreach  (@Appl) 
	{ $appl_parameter .= "$_,"; }
    $appl_parameter=~s/,$/ /;
}

my $dir = `pwd`; chomp ($dir);
my $fasta = shift;
my $fasta_name = basename($fasta);
my $species_name = $fasta_name;
#if ($species_name=~/.gene.pep$/) {
#	$species_name=~s/.gene.pep$//;
#}
my $config_file = "$Bin/../../config.txt";
my $run_iprscan = parse_config($config_file,"run_iprscan");
my $auto_blast = parse_config($config_file,"auto_blast");
my $stat_fun_ann = parse_config($config_file,"stat_fun_ann");
my $venny = parse_config($config_file,"venny");
my $addDesc = parse_config($config_file,"addDesc");
my $go =parse_config($config_file,"go");


##add by luchx
my $QP_para;
$QP_para.=" --queue $Queue " if (defined $Queue);
$QP_para.=" --pro_code $Pro_code " if (defined $Pro_code);

open(OUT , ">" . "STEP01_fun_ann_work.sh") or die $!;
	if (defined $Interpro)
	{
		print OUT "##########InterPro##########\n";
		print OUT "mkdir Interpro\ncd Interpro\n";
		print OUT "nohup perl $run_iprscan --cpu $Cpu --cuts $Cuts $QP_para --seqtype $seqtype --resource $iResource $appl_parameter $fasta &\n";
		print OUT "cd ..\n";
	}
	if (defined $KEGG)
	{
		print OUT "##########KEGG##########\n";
		print OUT "mkdir KEGG\ncd KEGG\n";
		print OUT "nohup perl $auto_blast/auto_kegg_blast.pl $fasta $Queue $Pro_code $species $seqtype $Cuts $evalue ./ &\n";
		print OUT "cd ..\n";
	}
	if (defined $Swissprot)
	{
		print OUT "##########SwissProt##########\n";
		print OUT "mkdir Swissprot\ncd Swissprot\n";
		print OUT "nohup perl $auto_blast/auto_swissprot_blast.pl $fasta $Queue $Pro_code $seqtype $Cuts $evalue ./ &\n";
		print OUT "cd ..\n";
	}
	if (defined $Trembl)
	{
		print OUT "##########Trembl##########\n";
		print OUT "mkdir Trembl\ncd Trembl\n";
		print OUT "nohup perl $auto_blast/auto_trembl_blast.pl $fasta $Queue $Pro_code $seqtype $Cuts $evalue ./ &\n";
		print OUT "cd ..\n";
	}
	if (defined $COG)
	{
		print OUT "##########COG##########\n";
		print OUT "mkdir COG\ncd COG\n";
		print OUT "nohup perl $auto_blast/auto_COG_blast.pl $fasta $Queue $Pro_code $seqtype $Cuts $evalue ./ &\n";
		print OUT "cd ..\n";
	}
	if (defined $NT)
	{
		print OUT "##########NT##########\n";
		print OUT "mkdir NT\ncd NT\n";
		print OUT "nohup perl $auto_blast/auto_nt_blast.pl $fasta $Queue $Pro_code $species $seqtype $Cuts $evalue ./ &\n";
		print OUT "cd ..\n";
	}
	if ( $NR && ! $GO )
	{
		print OUT "##########NR##########\n";
		print OUT "nohup perl $auto_blast/auto_nr_blast2GO.pl $fasta $Queue $Pro_code $species $seqtype $Cuts $evalue ./ NO &\n"
	}
	if ( $NR && $GO )
	{
		print OUT "##########NR && GO##########\n";
		print OUT "nohup perl $auto_blast/auto_nr_blast2GO.pl $fasta $Queue $Pro_code $species $seqtype $Cuts $evalue ./ GO &\n"
	}
close(OUT);

open(OUT , ">" . "STEP02_fun_ann_stat.sh") or die $!;
	my $stat_parameter = "";
	my @venn_file;
	my @venn_name;
	print OUT "mkdir fun_ann_stat\n";
	print OUT "cd fun_ann_stat\n";
	print OUT "ln -s $fasta $fasta_name && \\\n";
	print OUT "grep \'^>\' $fasta_name | sed \'s/>//\' | awk \'{print \$1}\' > all_gene.id && \\\n";
	if (defined $Interpro)
	{
		print OUT "ln -s $dir/Interpro/$species_name.iprscan.xls && \\\n";
		$stat_parameter .= "-Interpro $species_name.iprscan.xls ";
		push @venn_file, "$species_name.iprscan.xls";
		push @venn_name, "InterPro";
	}
	if (defined $KEGG)
	{
		print OUT "ln -s $dir/KEGG/$species_name.blast.kegg.xls && \\\n";
		$stat_parameter .= "-kegg $species_name.blast.kegg.xls ";
		push @venn_file, "$species_name.blast.kegg.xls";
		push @venn_name, "KEGG";
	}
	if (defined $Swissprot)
	{
		print OUT "ln -s $dir/Swissprot/$species_name.blast.swissprot.xls && \\\n";
		$stat_parameter .= "-swissprot $species_name.blast.swissprot.xls ";
		push @venn_file, "$species_name.blast.swissprot.xls";
		push @venn_name, "SwissProt";
	}
	if (defined $Trembl)
	{
		print OUT "ln -s $dir/Trembl/$species_name.blast.trembl.xls && \\\n";
		$stat_parameter .= "-trembl $species_name.blast.trembl.xls ";
	}
	if (defined $COG){
		print OUT "ln -s $dir/COG/$species_name.blast.cog.xls && \\\n";
		$stat_parameter .= "-cog $species_name.blast.cog.xls ";
		push @venn_file, "$species_name.blast.cog.xls";
		push @venn_name, "COG";
	}
	if (defined $NT){
		print OUT "ln -s $dir/NT/$species_name.blast.nt.xls && \\\n";
		$stat_parameter .= "-nt $species_name.blast.nt.xls ";
	}
	if (defined $NR){
		print OUT "ln -s $dir/NR/$species_name.blast.nr.xls && \\\n";
		$stat_parameter .= "-nr $species_name.blast.nr.xls ";
		push @venn_file, "$species_name.blast.nr.xls";
		push @venn_name, "NR";
	}
	if ($NR && $GO){
		print OUT "ln -s  $dir/GO/$species_name.Gene2GO.xls && \\\n";
		$stat_parameter .= "-go $species_name.Gene2GO.xls -obo $go/gene_ontology.1_2.obo ";
	}
	print OUT "perl $stat_fun_ann -list all_gene.id $stat_parameter -outxls ./annotation.xls --outstat ./annotation_stat.xls &&\\\n";
	if (@venn_file >= 2 && @venn_file <= 5){
		print OUT "perl $venny -infile " . join(",",@venn_file) . " -name " . join(",",@venn_name) . " -header -color -outdir ./Venny -imgname Annotation.venn && \\\n";
		print OUT "perl $addDesc -input \"\./Venny/*.xls\" -annot annotation.xls && \\\n";
	}
	print OUT "cd ..\n";
close(OUT);

open(OUT , ">" . "STEP03_delete_tmp_files.sh") or die $!;
#	print OUT "cd function_annotation\n";
	if (defined $Interpro)
	{
		print OUT "cd Interpro\n";
		print OUT "rm -r *.cut *.sh.* nohup.out\n";	
		print OUT "cd ..\n";
	}
	if (defined $KEGG)
	{
		print OUT "cd KEGG\n";
		print OUT "rm -r *.cut *.sh.* nohup.out\n";
		print OUT "cd ..\n";
	}
	if (defined $Swissprot)
	{
		print OUT "cd Swissprot\n";
		print OUT "rm -r *.cut *.sh.* nohup.out\n";
		print OUT "cd ..\n";
	}
	if (defined $Trembl)
	{
		print OUT "cd Trembl\n";
		print OUT "rm -r *.cut *.sh.* nohup.out\n";
		print OUT "cd ..\n";
	}
	if (defined $COG)
	{
		print OUT "cd COG\n";
		print OUT "rm -r *.cut *.sh.* nohup.out\n";
		print OUT "cd ..\n";
	}
	if (defined $NT){
		print OUT "cd NT\n";
		print OUT "rm -r *.cut *.sh.* nohup.out\n";
		print OUT "cd ..\n";
	}
	if (defined $NR){
		print OUT "cd NR\n";
		print OUT "rm -r *.cut *.sh.* nohup.out\n";
		print OUT "cd ..\n";
	}
	if ($NR && $GO){
		print OUT "cd GO\n";
		print OUT "rm -r *.cut *.sh.* nohup.out\n";
		print OUT "cd ..\n";
	}
close(OUT);
