#!/usr/bin/perl
=head1 Name

	auto_repeat.pl  -- the pipeline of denovo predict repeat and find repeat.

=head1 Description

	Base on denovo_repeat_find.pl and find_repeat.pl.This program invoke Piler, RepeatScout, LTR_FINDER, RepeatModeler,
	RepeatMasker and  ProteinMask.The path of the softwares invoked will be put in the "config.txt". 

=head1 Version

v1.0

=head1 Usage

	perl auto_repeat.pl [options] seq_file

	note: the seq_file must be input by absolute path

	### denovo predict parameter ###	
	-dcpu <int>		set th cpu number to use in parallel, default=3
	-dcutf <int>		set the number of cutted files
	-drun <str>		set the parallel type, qsub or multi, default=qsub
	-Piler			run Piler
	  -piler_cuts <int>	cut file with the specified number of seqeunces in each sub file
	  -piler_resource <str>	set the required resource used in qsub -l option, default vf=3G
	-RepeatScout		run RepeatScout
	-LTR_FINDER		run LTR_FINDER
	  -tRNA <str>		set the lib file for LTR_FINDER, default Athal-tRNAs.fa 
	-RepeatModeler		run RepeatModeler
	### filter denovo library parameter
	-f_Lcutf <int>		set the number of cutted files for LTRfinder result
	-f_Pcutf <int>		set the number of cutted files for piler result
	-f_Rcutf <int>		set the number of cutted files for repeatscout result
	### find repeat parameter ###
	-kcpu <int>		set th cpu number to use in parallel, default=3
	-kcutf <int>		set the number of cutted files
	-node <str>       set the compute node, eg h=compute-0-151,default is not set
	-krun <str>		set the parallel type, qsub or multi, default=qsub
	-kresource <str>	set the required resource used in qsub -l option, default vf=3G
	-Trf             	run trf
	  -period_size <int>    set the maximum period size for trf, default=2000
	-RepeatMasker    	run RepeatMasker
	  -lib <file>      	set the lib file for RepeatMasker, default Repbase
	-ProteinMask     	run RepeatProteinMask
	  -pvalue <str>    	set the pvalue for RepeatProteinMask, default=1e-4
	-queue <str>            set the queue in qsub, default no
	-pro_code <str>         set the code for project,default no

	-Verbose         	output verbose information to screen
	-Help            	output help information to screen.

=head1 Exmple

	nohup perl auto_repeat.pl -Piler -RepeatScout -LTR_FINDER -Trf -RepeatMasker -ProteinMask /ifs1/GAG/annotation/zhouheling/bin/auto_repeat/RHOziaD.seq &
	
	nohup perl auto_repeat.pl -Piler -piler_cuts 2000 -piler_resource vf=5g /ifs1/GAG/annotation/zhouheling/bin/auto_repeat/RHOziaD.seq &
	
	nohup perl auto_repeat.pl -dcpu 50 -dcutf 100 -Piler -RepeatScout -LTR_FINDER -RepeatModeler /ifs1/GAG/annotation/zhouheling/test/seq_file &

	nohup perl auto_repeat.pl -kcpu 50 -kcutf 50 -kresource vf=1g -Trf -RepeatMasker -ProteinMask /ifs1/GAG/annotation/zhouheling/test/seq_file &

	nohup perl auto_repeat.pl -dcpu 50 -dcutf 100 -LTR_FINDER -kcpu 50 -kcutf 50 -Trf /ifs1/GAG/annotation/zhouheling/test/seq_file &

=cut

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use lib "$Bin/common_bin";
use GACP qw(parse_config);


my ($piler,$piler_cuts,$piler_resource,$repeatscout,$ltr_finder,$tRNA);
my ($repeatmodeler,$denovo_cpu,$denovo_cutf,$denovo_run);
my ($f_Lcutf,$f_Pcutf ,$f_Rcutf);
my ($trf,$trf_period_size,$repeatmasker,$lib,$proteinmask,$pvalue);
my ($known_cpu,$known_cutf,$Node,$known_run,$known_resource,$Queue,$Pro_code);
my ($verbose,$help);

GetOptions(
        "Piler"=>\$piler,
	    "piler_cuts:i"=>\$piler_cuts,
	    "piler_resource:s"=>\$piler_resource,
        "RepeatScout"=>\$repeatscout,
        "LTR_FINDER"=>\$ltr_finder,
	    "tRNA:s"=>\$tRNA,
        "RepeatModeler"=>\$repeatmodeler,
	"dcpu:i"=>\$denovo_cpu,
	"dcutf:i"=>\$denovo_cutf,
	"drun:s"=>\$denovo_run,
	"f_Lcutf:i"=>\$f_Lcutf,
	"f_Pcutf:i"=>\$f_Pcutf,
	"f_Rcutf:i"=>\$f_Rcutf,
        "Trf"=>\$trf,
	    "trf_period_size:i"=>\$trf_period_size,
        "RepeatMasker"=>\$repeatmasker,
	    "lib:s"=>\$lib,
        "ProteinMask"=>\$proteinmask,
	    "pvalue:s"=>\$pvalue,
	"kcpu:i"=>\$known_cpu,
	"kcutf:i"=>\$known_cutf,
	"node:s"=>\$Node,
	"krun:s"=>\$known_run,
	"kresource:s"=>\$known_resource,
	"queue:s"=>\$Queue,
	"pro_code:s"=>\$Pro_code,

        "Verbose"=>\$verbose,
        "Help"=>\$help,
);

die `pod2text $0` if(@ARGV==0 || $help);

#denovo parameter
$denovo_cpu ||= 6; $denovo_cutf ||= 20; $denovo_run ||= "qsub"; 
#piler parameter
$piler_cuts ||= 100; $piler_resource ||= "vf=3G";
#filter parameter
$f_Lcutf||=200;$f_Pcutf||=50 ;$f_Rcutf||=50;
#ltr_finder parameter
$tRNA ||= "/home/xiongzj/GACP/01.repeat_finding/software/LTR_Finder/Dm-tRNAs.fa";

#known parameter
$known_cpu ||= 6; $known_cutf ||= 20; $known_run ||= "qsub"; $known_resource ||= "vf=3G";
#trf parameter
$trf_period_size ||= 2000;
#repeatmasker parameter

my $Lib=(defined $lib)?"-lib $lib":"";
#proteinmask parameter
$pvalue ||= 1e-4;

my $dir = `pwd`; chomp ($dir);
my $seq_file = shift;
my $seq_name = basename($seq_file);
my $species_name=$1 if ($seq_name=~/^(\S+)\.fasta/);
my $config_file = "$Bin/config.txt";
my $denovo_program = parse_config($config_file,"denovo_program");
my $known_program = parse_config($config_file,"known_program");
my $filter_program = parse_config($config_file,"filter_program");
my $stat_prgram = parse_config($config_file,"stat_program");
my $change_ID = parse_config($config_file,"change_ID_path"); 

##add by luchx
my $QP_para;
$QP_para.="--queue $Queue " if (defined $Queue);
$QP_para.="--pro_code $Pro_code " if (defined $Pro_code);

open(OUT , ">" . "STEP01_repeat_work.sh") or die $!;
if ((defined $piler) or (defined $repeatscout) or (defined $ltr_finder) or (defined $repeatmodeler))
{
	print OUT "mkdir denovo\n";
	print OUT "cd denovo\n";
	if (defined $piler)
	{
		print OUT "mkdir piler\n";
		print OUT "cd piler\n";
		if(defined $Node){
			print OUT "nohup perl $denovo_program $QP_para -Piler -cpu $denovo_cpu -cuts $piler_cuts -run $denovo_run -node $Node -resource $piler_resource $seq_file &\n";
		}
		else { print OUT "nohup perl $denovo_program $QP_para -Piler -cpu $denovo_cpu -cuts $piler_cuts -run $denovo_run  -resource $piler_resource $seq_file &\n";}
		print OUT "cd ..\n";
	}
	if (defined $repeatscout)
	{
		print OUT "mkdir repeatscout\n";
		print OUT "cd repeatscout\n";
		if(defined $Node){
			print OUT "nohup perl $denovo_program $QP_para -RepeatScout -cpu $denovo_cpu -cutf $denovo_cutf -run $denovo_run -node $Node $seq_file &\n";
		}
		else {print OUT "nohup perl $denovo_program $QP_para -RepeatScout -cpu $denovo_cpu -cutf $denovo_cutf -run $denovo_run  $seq_file &\n";}
		print OUT "cd ..\n";
	}
	if (defined $ltr_finder)
	{
		print OUT "mkdir ltr_finder\n";
		print OUT "cd ltr_finder\n";
		if(defined $Node){
			print OUT "nohup perl $denovo_program $QP_para -LTR_FINDER -cpu $denovo_cpu -cutf $denovo_cutf -run $denovo_run -node $Node -tRNA $tRNA $seq_file &\n";
		}
		else {print OUT "nohup perl $denovo_program $QP_para -LTR_FINDER -cpu $denovo_cpu -cutf $denovo_cutf -run $denovo_run  -tRNA $tRNA $seq_file &\n";}
		print OUT "cd ..\n";
	}
	if (defined $repeatmodeler)
	{
		print OUT "mkdir repeatmodeler\n";
		print OUT "cd repeatmodeler\n";
		if(defined $Node){
			print OUT "nohup perl $denovo_program $QP_para -RepeatModeler $seq_file -node $Node &\n";
		}
		else {print OUT "nohup perl $denovo_program $QP_para -RepeatModeler $seq_file  &\n";}
		print OUT "cd ..\n";
	}
	print OUT "cd ..\n";
}

if ((defined $trf) or (defined $repeatmasker) or (defined $proteinmask))
{
	print OUT "mkdir known\n";
	print OUT "cd known\n";
	if (defined $trf)
	{
		print OUT "mkdir trf\n";
		print OUT "cd trf\n";
		if(defined $Node){
			print OUT "nohup perl $known_program $QP_para -trf -cpu $known_cpu -cutf $known_cutf -run $known_run -node $Node -resource $known_resource -period_size $trf_period_size $seq_file &\n";
		}
		else{print OUT "nohup perl $known_program $QP_para -trf -cpu $known_cpu -cutf $known_cutf -run $known_run  -resource $known_resource -period_size $trf_period_size $seq_file &\n";}
		print OUT "cd ..\n";
	}
	if (defined $repeatmasker)
	{
		print OUT "mkdir repeatmasker\n";
		print OUT "cd repeatmasker\n";
		if(defined $Node){
			print OUT "nohup perl $known_program $QP_para -repeatmasker -cpu $known_cpu -cutf $known_cutf -run $known_run -node $Node -resource $known_resource -lib $lib $seq_file &\n";
		}
		else{print OUT "nohup perl $known_program $QP_para -repeatmasker -cpu $known_cpu -cutf $known_cutf -run $known_run  -resource $known_resource $Lib $seq_file &\n";}
		print OUT "cd ..\n";
	}
	if (defined $proteinmask)
	{
		print OUT "mkdir proteinmask\n";
		print OUT "cd proteinmask\n";
		if(defined $Node){
			print OUT "nohup perl $known_program $QP_para -proteinmask -cpu $known_cpu -cutf $known_cutf -run $known_run -node $Node -resource $known_resource -pvalue $pvalue $seq_file &\n";
		}
		else{print OUT "nohup perl $known_program $QP_para -proteinmask -cpu $known_cpu -cutf $known_cutf -run $known_run  -resource $known_resource -pvalue $pvalue $seq_file &\n";}
		print OUT "cd ..\n";
	}
	print OUT "cd ..\n";
}	
close(OUT);

open(OUT , ">" . "STEP02_filter_work.sh") or die $!;
open(OUT2, ">" . "STEP03_denovo_repeatmasker.sh") or die $!;
if ((defined $piler) or (defined $repeatscout) or (defined $ltr_finder) or (defined $repeatmodeler))
{
	print OUT "cd denovo\n";
	print OUT2 "cd denovo\nmkdir repeatmasker\ncd repeatmasker\n";
	my $cat_file = "";
	if (defined $piler)
	{
		print OUT "cd piler\nmkdir filter\ncd filter\n";
		print OUT "ln -s $dir/denovo/piler/Piler_Result/$seq_name.piler_library.fa $seq_name.piler_library.fa\n";
		if(defined $Node){
			print OUT "nohup perl $filter_program $QP_para --cpu $f_Pcutf $seq_name.piler_library.fa Piler $Node &\n";
		}
		else{print OUT "nohup perl $filter_program $QP_para --cpu $f_Pcutf $seq_name.piler_library.fa Piler &\n";}
		print OUT "cd ..\ncd ..\n";
		$cat_file = $cat_file . "$dir/denovo/piler/filter/$seq_name.piler_library.fa.final.library" . " ";
	}
	if (defined $repeatscout)
	{
		print OUT "cd repeatscout\nmkdir filter\ncd filter\n";
		print OUT "ln -s $dir/denovo/repeatscout/RepeatScout/$seq_name.17mer.repeatscout.filter.library.fa $seq_name.17mer.repeatscout.filter.library.fa\n";
		if(defined $Node){
			print OUT "nohup perl $filter_program $QP_para --cpu $f_Rcutf $seq_name.17mer.repeatscout.filter.library.fa RepeatScout $Node &\n";
		}
		else{print OUT "nohup perl $filter_program $QP_para --cpu $f_Rcutf $seq_name.17mer.repeatscout.filter.library.fa RepeatScout  &\n";}
		print OUT "cd ..\ncd ..\n";
		$cat_file = $cat_file . "$dir/denovo/repeatscout/filter/$seq_name.17mer.repeatscout.filter.library.fa.final.library" . " ";
	}
	if (defined $ltr_finder)
	{
		print OUT "cd ltr_finder\nmkdir filter\ncd filter\n";
		print OUT "ln -s $dir/denovo/ltr_finder/$seq_name.LTR.fa $seq_name.LTR.fa\n";
		print OUT "nohup perl $filter_program $QP_para --cpu  $f_Lcutf $seq_name.LTR.fa LTR_FINDER &\n";
		#print OUT "ln -s $dir/denovo/ltr_finder/$seq_name.LTR.fa $seq_name.LTR.fa.final.library\n";
		print OUT "cd ..\ncd ..\n";
		$cat_file = $cat_file . "$dir/denovo/ltr_finder/filter/$seq_name.LTR.fa.final.library" . " ";
	}
	if (defined $repeatmodeler)
	{
		$cat_file = $cat_file . "$dir/denovo/repeatmodeler/result/consensi.fa.classified";
	}
	print OUT "cd ..\n";
	print OUT2 "cat $cat_file > final.library\n";
	if(defined $Node){
		print OUT2 "nohup perl $known_program $QP_para -repeatmasker -cpu $known_cpu -cutf $known_cutf -run $known_run -node $Node -resource $known_resource -lib $dir/denovo/repeatmasker/final.library $seq_file &\n";
	}
	else{print OUT2 "nohup perl $known_program $QP_para -repeatmasker -cpu $known_cpu -cutf $known_cutf -run $known_run  -resource $known_resource -lib $dir/denovo/repeatmasker/final.library $seq_file &\n";}
	print OUT2 "cd ..\ncd ..\n";
}
close(OUT);
close(OUT2);

open(OUT , ">" . "STEP04_repeat_statistics.sh") or die $!;
	print OUT "mkdir statistics\ncd statistics\n";
	print OUT "ln -s $seq_file $seq_name\n";
	my $cat_file_with_trf = "";
	my $cat_file_without_trf = "";
	my $stat_parameter = "";
	if ((defined $piler) or (defined $repeatscout) or (defined $ltr_finder) or (defined $repeatmodeler))
	{
		print OUT "ln -s $dir/denovo/repeatmasker/$species_name.denovo.RepeatMasker.gff denovo_tmp.gff\n";
		print OUT "perl $change_ID denovo_tmp.gff\n";
		print OUT "mv denovo_tmp.gff.denovo.gff denovo.gff\n";
		print OUT "rm denovo_tmp.gff\n";
		print OUT "ln -s $dir/denovo/repeatmasker/$species_name.denovo.RepeatMasker.out denovo.out\n";
		$cat_file_with_trf = $cat_file_with_trf . "denovo.gff" . " ";
		$cat_file_without_trf = $cat_file_without_trf . "denovo.gff" . " ";
		$stat_parameter = $stat_parameter . "-denovo" . " ";
	}
	if (defined $trf)
	{
		print OUT "ln -s $dir/known/trf/$species_name.TRF.gff trf.gff\n";
		$cat_file_with_trf = $cat_file_with_trf . "trf.gff" . " ";
		$stat_parameter = $stat_parameter . "-trf" . " ";
	}
	if (defined $repeatmasker)
	{
		print OUT "ln -s $dir/known/repeatmasker/$species_name.known.RepeatMasker.gff repeatmasker.gff\n";
		print OUT "ln -s $dir/known/repeatmasker/$species_name.known.RepeatMasker.out repeatmasker.out\n";
		$cat_file_with_trf = $cat_file_with_trf . "repeatmasker.gff" . " ";
		$cat_file_without_trf = $cat_file_without_trf . "repeatmasker.gff" . " ";
		$stat_parameter = $stat_parameter . "-repeatmasker" . " ";
	}
	if (defined $proteinmask)
	{
		print OUT "ln -s $dir/known/proteinmask/$species_name.RepeatProteinMask.gff proteinmask.gff\n";
		$cat_file_with_trf = $cat_file_with_trf . "proteinmask.gff" . " ";
		$cat_file_without_trf = $cat_file_without_trf . "proteinmask.gff" . " ";
		$stat_parameter = $stat_parameter . "-proteinmask" . " ";
	}
	print OUT "cat $cat_file_with_trf > all.gff\n";
	print OUT "cat $cat_file_without_trf > all_without_trf.gff\n";
	print OUT "nohup perl $stat_prgram $stat_parameter $seq_name &\n";
	print OUT "cd ..\n";
close(OUT);

open(OUT , ">" . "STEP05_delete_tmp_files.sh") or die $!;
	if (defined $piler)
	{
		print OUT "mv $dir/denovo/piler/filter/$seq_name.piler_library.fa.final.library $dir/denovo/piler/$seq_name.piler_library.fa.final.library\n";
		print OUT "mv $dir/denovo/piler/Piler_Result/*.gff $dir/denovo/piler/\n";
		print OUT "rm -r ./denovo/piler/Piler_Result\n";
	}
	if (defined $repeatscout)
	{
		print OUT "mv $dir/denovo/repeatscout/filter/$seq_name.17mer.repeatscout.filter.library.fa.final.library $dir/denovo/repeatscout/$seq_name.17mer.repeatscout.filter.library.fa.final.library\n";
		print OUT "rm -r ./denovo/repeatscout/RepeatScout\n\n";
	}
	if (defined $ltr_finder)
	{
		print OUT "mv $dir/denovo/ltr_finder/filter/$seq_name.LTR.fa.final.library $dir/denovo/ltr_finder/$seq_name.LTR.fa.final.library\n";
		print OUT "rm -r ./denovo/ltr_finder/LTR_Result\n";
	}
	my $tmp_dir = "./*/*";
	print OUT "rm -r $tmp_dir/*.cut $tmp_dir/nohup.out $tmp_dir/*.shell $tmp_dir/*.shell.* $tmp_dir/*.sh $tmp_dir/*.sh.* \n";
	if ((defined $piler) or (defined $repeatscout) or (defined $ltr_finder))
	{
		print OUT "rm -r $tmp_dir/filter\n";
	}
close(OUT);
