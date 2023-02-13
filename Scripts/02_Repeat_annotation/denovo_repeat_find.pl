#!/usr/bin/perl
=head1  Name

	denovo_repeat_find.pl -- the pipeline of denovo predict repeat in sequences.

=head1 Description

	This programm is to denovo predict the repeat sequences (mainly the transponse elements),which contains the softwares such as PILER,RepeatScout,LTR_FINDER and RepeatModeler.You can choose one of these program or all of them according to the vary genome.
	
=head1 Usage:
	
	For small genome sequence(<=600M),you'd better choose the piler and RepeatScout for denovo repeat library construction;and for other genome sequences,you'd better choose the RepeatModeler.The LTR_FINDER is selectable.

=head1 example

	perl denovo_repeat_find.pl -Piler -RepeatScout -LTR_FINDER -RepeatModeler ./test.fa

=head1 Information
	
    Modified:   1. add "-clear -l num_proc=n -binding linear:n", and use Galgal as tRNAdb, remove the reqsub for qsub-sge.pl
                2. add "next, if($ eq "." or $ eq "..");" for each readdir function, so the directory error will not happened
                3. add export TRF_COMMAND and NSEG_COMMAND in repeatscout.sh for RepeatScout
                4. vf=3G to vf=6G

=head1 Usage

        perl denovo_repeat_find.pl [options] input_file
        -Piler                  run Piler
        -RepeatScout            run RepeatScout
        -LTR_FINDER             run LTR_FINDER
        -tRNA <str>             set tRNA lib
        -cpu <int>              set the cpu number to use in parallel, default=3
        -queue <str>            specify the queue to use, default no
        -pro_code <str>         set the code for project,default no
        -run <str>              set the parallel type, qsub or multi, default=qsub
        -node <str>             set the compute node, eg h=compute-0-151,default is not set
        -cutf <int>             cut file with the specified number of subfile in total
        -cuts <int>             cut file with the specified number of seqeunces in each sub file
        -RepeatModeler          run RepeatModeler
        -All_method             run all method
        -resource <str>         set the required resource used in qsub -l option, default vf=6G
        -verbose                output verbose information to screen
        -help                   output help information to screen.

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use lib "$Bin/";
use GACP qw(parse_config);

my ($piler,$repeatscout,$Queue,$Run,$Node,$tRNA,$ltr_finder,$repeatmodeler,$all,$Cpu,$Resource,$cutf,$cuts,$Verbose,$Pro_code,$help);
GetOptions(
	"Piler"=>\$piler,
	"RepeatScout"=>\$repeatscout,
	"LTR_FINDER"=>\$ltr_finder,
	"tRNA:s"=>\$tRNA,
	"cpu:i"=>\$Cpu,
	"queue:s"=>\$Queue,
        "pro_code:s"=>\$Pro_code,
	"run:s"=>\$Run,
	"node:s"=>\$Node,
	"cutf:i"=>\$cutf,
	"cuts:i"=>\$cuts,
	"RepeatModeler"=>\$repeatmodeler,
	"All_method"=>\$all,
	"resource:s"=>\$Resource,
	"verbose"=>\$Verbose,
	"help"=>\$help
);
die `pod2text $0` if(@ARGV==0 || $help);

$Cpu ||= 2;
$Run ||= "qsub";
$cuts ||=100;
$cutf||=20;
$Resource ||= "vf=6G";

my $seq_file=shift;

my $seq_name=basename($seq_file);
my $species_name=$1 if ($seq_name=~/^(\S+)\.fasta/);
my $config_file="$Bin/../config.txt";
my $common_bin = "$Bin/";
my $fastaDeal = "$common_bin/fastaDeal.pl";
my $find_repeat=parse_config($config_file,"known_program");

##############	PILER PATH	#############################################
my $Pals=parse_config($config_file,"pals");
my $Piler=parse_config($config_file,"piler");
my $Piler_path=parse_config($config_file,"piler_path");
my $muscle_path=parse_config($config_file,"muscle");

##############	RepeatScout	PATH	#####################################
my $repeatscout_path=parse_config($config_file,"repeatscout_path");
my $repeatmasker=parse_config($config_file,"repeatmasker");
my $trf=parse_config($config_file,"trf");
my $nseg=parse_config($config_file,"nseg");

##############	RepeatModeler	PATH	#####################################
my $repeatmodeler_path=parse_config($config_file,"repeatmodeler_path");

##############	LTR_FINDER	PATH	#####################################
my $ltr_path=parse_config($config_file,"ltr_path");
$tRNA ||= "$ltr_path/tRNA/Ggal2-tRNAs-all.fa";
my $dir=`pwd`;
chomp $dir;
#print "$dir";

$Queue ||= "st.q";
$Pro_code ||= "P18Z10200N0160";
##add by helijuan
my $QP_para;
$QP_para.="--queue $Queue " if (defined $Queue);
$QP_para.="--pro_code $Pro_code " if (defined $Pro_code);

##################################################
##################################################
##################################################
if(defined $piler){
		my $outdir=$dir."/Piler_Result";
		mkdir $outdir unless (-d $outdir);
		my $hit=$outdir."/".$seq_name.".hit.gff";
		my $trs=$outdir."/".$seq_name.".hit.trs.gff";
		my $fams=$outdir."/".$seq_name."_families";
		my $cons=$outdir."/".$seq_name."_cons";
		my $aligned_fams=$outdir."/".$seq_name."_aligned_fams";
		my $library=$outdir."/".$seq_name.".piler_library.fa";
		my $log=$outdir."/".$seq_name.".piler.log";
		my $piler_shell=$outdir."/"."Run_piler_for_$seq_name.sh";
		my $cut= $dir."/".$seq_name.".cut/";
		`ln -s $seq_file` unless (-e $seq_name);
		my $num=`grep -c ">"  $seq_name`; 
		`perl $fastaDeal --cutf $cutf $seq_name` if($num<=4000);
		`perl $fastaDeal --cuts $cuts $seq_name` if($num>4000);
		my $pipeline="$Piler_path/pipeline.pl";
		`perl $pipeline $cut`;
		`perl $common_bin/multi-process.pl  -cpu $Cpu  find_pals.sh` if($Run eq "multi");
		if(defined $Node){
			`perl $common_bin/qsub-sge.pl $QP_para  -maxjob $Cpu --resource $Resource --node $Node --convert no find_pals.sh` if($Run eq "qsub");
		}
		else{
		`perl $common_bin/qsub-sge.pl $QP_para -maxjob $Cpu --resource $Resource --convert no find_pals.sh` if($Run eq "qsub");}
		`cat $cut/pals*.log >pals.log`;
        `cat $cut/pals*.log2 >pals.log2`;
        `cat $cut/*.hit.gff > final.all.hit.gff`;
		`mv final.all.hit.gff $hit` if ($Run eq "multi");
		`rm *.gff` if ($Run eq "multi");
		`mv final.all.hit.gff $hit` if ($Run eq "qsub");
		`$Piler -trs $hit -out $trs 1>>$log`;
		if(-d $fams){`rm -fr $fams`;}
		unless(-d $fams){`mkdir $fams`;}
		`$Piler -trs2fasta $trs -seq $seq_file -path $fams >>$log`;

		if(-d $aligned_fams){`rm -fr $aligned_fams`;}
		unless(-d $aligned_fams){`mkdir $aligned_fams`;}
		opendir THIS,$fams or die "$!";
		my @family=readdir THIS;
		closedir THIS;
		
		foreach my $fam(@family)
		{
            next, if($fam eq "." or $fam eq ".."); ## add by fangqi
			my $full_name="$fams/$fam";
			my $name=basename($fam);
			my $out="$aligned_fams/$name";
			`$muscle_path -in $full_name -out $out -maxiters 1 -diags 1>>$log ;`
		}
		
		if(-d $cons){`rm -fr $cons`;} 
		unless(-d $cons){print OUT `mkdir $cons`;}
		
		opendir THIS,$aligned_fams or die "$!"; 
		my @consensus=readdir THIS;
		closedir THIS;
		
		foreach my $con(@consensus)
		{
                next, if($con eq "." or $con eq ".."); ## add by fangqi
				my $full_name="$aligned_fams/$con";
				my $name=basename($con);
				my $out="$cons/$name"; 
				`$Piler -cons $full_name -out $out -label $con 1>>$log`;
				`cat $out >> $library`;
		}
}

##################################################
##################################################
##################################################
if(defined $repeatscout){
		my $outdir=$dir."/RepeatScout";
		mkdir $outdir unless(-d $outdir);
		`ln -s $seq_file` unless (-e $seq_name);
		my $repscout=$outdir."/".$seq_name.".17mer.repeatscout";
		my $freq=$outdir."/".$seq_name.".17mer.tbl";
		my $filter=$repscout.".filter";
		my $out=$outdir."/".$seq_name.".out";
		my $library=$filter.".library.fa";
		my $log=$outdir."/".$seq_name.".17mer.repeatscout.log";
		open OUT1,">build_lmer_table.sh" or die "$!";
		open OUT,"> repeatscout.sh" or die "$!";
        print OUT "export NSEG_COMMAND=$nseg\nexport TRF_COMMAND=$trf\n";
		print OUT1 "$repeatscout_path/build_lmer_table -l 17 -sequence $seq_file -freq $freq 1 >> $log\n";
		open OUT2, ">repeatscout_out.sh"or die "$!";
		print OUT2 "$repeatscout_path/RepeatScout -sequence $seq_file -output $repscout -freq $freq -l 17 1 >> $log\n";
		if(defined $Node){
			print OUT "perl $common_bin/qsub-sge.pl $QP_para --lines 1 --resource  $Resource --node $Node  build_lmer_table.sh \n";
		}
		else{print OUT "perl $common_bin/qsub-sge.pl $QP_para  --lines 1 --resource  $Resource   build_lmer_table.sh \n";}	
		if(defined $Node){
			print OUT "perl $common_bin/qsub-sge.pl $QP_para  --lines 1 --resource  $Resource --node $Node   repeatscout_out.sh  \n";
		}
		else{ print OUT "perl $common_bin/qsub-sge.pl $QP_para --lines 1 --resource  $Resource   repeatscout_out.sh  \n";}	
		print OUT "cat $repscout | $repeatscout_path/filter-stage-1.prl 1>$filter 2>>$log\n";
		print OUT "$repeatmasker -nolow -no_is -norna -parallel 1 -lib $filter $seq_file 1>>$log\n  mv $seq_name.RepeatMasker.out $out\n" if($Run eq "multi");
		if(defined $Node){
			print OUT "$find_repeat $QP_para -repeatmasker -cutf $cutf -cpu $Cpu -node $Node -lib $filter $seq_name\n mv $species_name.denovo.RepeatMasker.out $out\n" if($Run eq "qsub");
		}
		else{ print OUT "$find_repeat $QP_para -repeatmasker -cutf $cutf -cpu $Cpu  -lib $filter $seq_name\n mv $species_name.denovo.RepeatMasker.out $out\n" if($Run eq "qsub");}
		print OUT "cat $filter |$repeatscout_path/filter-stage-2.prl --cat=$out --thresh=20 1>$library  2>>$log";
		close OUT;
		close OUT1;
		`sh repeatscout.sh`;
}

##################################################
##################################################
##################################################
if(defined $repeatmodeler)
{
		my $log=$seq_name.".repeatmodeler.log";
		my $rmd_shell="Run_repeatModeler_for_$seq_name.sh";
		open OUT , ">$rmd_shell" or die "$!";
		print OUT "#!/bin/bash\n";
#		print OUT "$repeatmodeler_path/BuildDatabase -name mydb $seq_file > $log \n"; changed by hlj
        print OUT "$repeatmodeler_path/BuildDatabase -engine ncbi -name mydb $seq_file > $log \n";
#		print OUT "$repeatmodeler_path/RepeatModeler  -database mydb > run.out\n"; changed by hlj
		print OUT "$repeatmodeler_path/RepeatModeler -engine ncbi -database mydb -pa 9  > run.out\n";
		print OUT "rm -r RM_*/round-*\n";
		print OUT "mv RM_* result\n";
		close OUT;
#		my $temp_node;
		if(defined $Node){
			$Node= "-l $Node";
			`nohup  qsub -clear -cwd -l num_proc=9,vf=5g -binding linear:9 -l $Node -q $Queue -P $Pro_code $rmd_shell\n`;##coordinated to "-pa 9"
		}
		else { `nohup  qsub -clear -cwd -l num_proc=9,vf=5g -binding linear:9 -q $Queue -P $Pro_code  $rmd_shell\n`;}##coordinated to "-pa 9"

}	

##################################################
##################################################
##################################################
if(defined $ltr_finder)
{
	my $outdir=$dir."/LTR_Result";
	mkdir $outdir unless(-d $outdir);
	`ln -s $seq_file` unless (-e $seq_name);
	my $cut= $dir."/".$seq_name.".cut/";
	`perl $fastaDeal --cutf  $cutf  $seq_name`;
	opendir DIR,$cut or die "$!";
	my @fa= grep {/$seq_name/} readdir DIR;
	closedir DIR;
#	my $ltr_shell="Run_LTR_FINDER_for_$seq_name.shell";
	my $ltr_shell="Run_LTR_FINDER_for_$seq_name.sh";
	open OUT,">$ltr_shell" or die "$!";
	foreach my $subfa(@fa)
	{
        next, if($subfa eq "." or $subfa eq ".."); ## add by fangqi
		my $ltr=$outdir."/".$subfa."ltr_finder";
		my $gff=$ltr.".gff";
		my $log="$outdir/$subfa.ltr_finder.log";
		my $cut_file="$cut$subfa";
#		print OUT "$ltr_path/ltr_finder  -w 2 -s $tRNA  $cut_file 1>$ltr 2>>$log; ";
		print OUT "$ltr_path/ltr_finder  -w 2 -s $tRNA  $cut_file 1>$ltr 2>>$log; ";
		print OUT "perl $ltr_path/Ltr2GFF.pl $ltr 1>$gff 2>>$log; ";
		print OUT "$common_bin/getTE.pl $gff $cut_file > $outdir/$subfa.LTR.fa\n";
	}
    close OUT;

    my $ltr_shell_combin="combin_LTR.sh";
    open  OUT1, ">$ltr_shell_combin" or die "$!";
	print OUT1 "cat $outdir/*.LTR.fa> $seq_name.LTR.fa\;";
	print OUT1 "cat $outdir/*ltr_finder.gff >$seq_name.ltr_finder.gff\;\n";
	close OUT1;

    open OUT2,">LTR_finder_qsub.sh" if (-e $ltr_shell && -e $ltr_shell_combin); 
    if(defined $Node){
#        $Node= "-l $Node";
        print OUT2 "perl $common_bin/qsub-sge.pl $QP_para --lines 1 --resource  $Resource --convert no --maxjob 100 -node $Node  $ltr_shell\n";
        print OUT2 "qsub -clear -cwd -l num_proc=1,$Resource -binding linear:1 -l $Node -q $Queue -P $Pro_code $ltr_shell_combin\n";
    }else {
        print OUT2 "perl $common_bin/qsub-sge.pl $QP_para --lines 1 --resource  $Resource --convert no --maxjob 100 $ltr_shell\n";
        print OUT2 "qsub -clear -cwd -l num_proc=1,$Resource -binding linear:1 -q $Queue -P $Pro_code $ltr_shell_combin\n";
    }
    close OUT2;

    `sh LTR_finder_qsub.sh`;
}
     
