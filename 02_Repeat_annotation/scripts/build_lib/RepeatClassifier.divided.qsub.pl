#!/usr/bin/perl -w

## Modified: add space to $QP_para; add clean/num_proc/-binding linear to qsub
## Author: fangqi@genomics.cn
## Date: 20180807

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use lib "$Bin/../../common_bin";
use GACP qw(parse_config);
use Data::Dumper;
use Cwd qw(abs_path);

my ($Queue,$Pro_code);
GetOptions(
	"queue:s"=>\$Queue,
        "pro_code:s"=>\$Pro_code,
);

my $config_file="$Bin/../../config.txt";
my $repeatclassifer=parse_config($config_file,"repeatclassifer");
my $fastaDeal=parse_config($config_file,"fastaDeal.pl");

&usage if @ARGV<1;

sub usage {
        my $usage = << "USAGE";

        This script run RepeatClassifer by divide db file into sub file.
	NOTE: RepeatClassifer Path: $repeatclassifer	
	Author: zhoujiajian\@genomics.org.cn
	Mender: zhongxiao\@genomics.org.cn
	Mender: luchangxin\@genomics.org.cn

	-queue <str>            set the queue in qsub, default no
	-pro_code <str>         set the code for project,default no

        Usage: $0 -queue bc.q -pro_code gagtest <raw_repeat_db_file> <output dir>
        Example:perl $0 db_file output_dir

USAGE
print "$usage";
exit(1);
};

##add by luchx
my $QP_para=" "; ##add a space, by fangqi 20180807
$QP_para.="-q $Queue " if (defined $Queue);
$QP_para.="-P $Pro_code " if (defined $Pro_code);

my $fasta_long = shift;
my $fasta = basename($fasta_long);
my $outdir = shift;
$outdir = abs_path("$outdir");

## Cut db file to pieces ##
`perl $fastaDeal -cutf 100 --outdir $outdir $fasta`;

my @subfiles = glob("$outdir/$fasta.cut/*");
my @command;
foreach my $subfile (@subfiles){
	push @command,"$repeatclassifer -consensi $subfile";	
}

my $i = "001";
my $Job_mark = "00001";
my $Job_prefix = "work";
foreach(@command){
	`mkdir $outdir/RepeatClassifier.qsub` unless(-e "$outdir/RepeatClassifier.qsub");
	`mkdir $outdir/RepeatClassifier.qsub/$i` unless(-e "$outdir/RepeatClassifier.qsub/$i");	
	open OUT,">","$outdir/RepeatClassifier.qsub/$i/$Job_prefix\_$Job_mark.sh" || die "$!\n";
	print OUT $_."; echo This-Work-is-Completed!\n";
	close OUT;
	$Job_mark++;
	$i++;
}

my @All_subsh = glob("$outdir/RepeatClassifier.qsub/*/work*.sh");

my %qid_h;

open OUT,">","$outdir/RepeatClassifier.sh" || die "$!\n";
foreach my $subsh (@All_subsh){
	#print $subsh,"\n";
	my $sub_dir = dirname($subsh);
	print OUT "cd $sub_dir;/opt/gridengine/bin/lx-amd64/qsub -clear -cwd -l num_proc=1,vf=1g -binding linear:1 $QP_para $subsh;cd $outdir;\n";
	chdir "$sub_dir";
	my $qid_r = `/opt/gridengine/bin/lx-amd64/qsub -clear -cwd -l num_proc=1,vf=1g -binding linear:1 $QP_para $subsh`;
	chomp ($qid_r);
	my ($qid,$qsh_r) = ($1,$2) if($qid_r =~ /Your job (\d+) \("([^"]+)"\) has been submitted/);
	my $qsh = $subsh;
	$qid_h{$qid} = $qsh;
	chdir "$outdir";
}
close OUT;

#Your job 7795328 ("work_00045.sh") has been submitted
#$VAR50 = [
#           '7795837',
#           '/ifs2/BC_GAG/Group/zhoujj/libary_build/test/Gnetum.scafSeq.FG1.LTR.fa.cut100bp.filter_redunance.fa.cut.qsub/050/work_00050.sh'
#         ];

#print STDERR Dumper(@qid);

## Detect job stat ######
open OUT,">","$outdir/RepeatClassifier.log" || die "Can't open file:$!\n";
while(1){
	my $jobid = get_jobid();
	my $flag = 0;
	foreach my $qid (keys %qid_h){
		my $sh = $qid_h{$qid};
		if(exists $jobid->{$qid}){
			$flag = 1;
			next;
		}elsif(-e "$sh.o.$qid"){
			open IN,"$sh.o.$qid" || die "Can't open file:$!\n";
			my @str = <IN>;
			chomp(@str);
			my $str = join " ",@str;
			if($str =~ /This-Work-is-Complieted!/){
				print OUT "Job: $qid is finished!\n";
			}
			close IN;
		}	
	}
	if($flag == 0){
		print OUT "All jobs finished!\n";
		last;
	}
	sleep(10);	
}
close OUT;

my @result = glob("$outdir/$fasta.cut/*.classified");
my $result =  join " ",@result;

system("cat $result > $outdir/$fasta.classified") == 0 || die "cat result failed!\n";

sub get_jobid{
	my $user = `whoami`;
	chomp ($user);
	my $qstat = `/opt/gridengine/bin/lx-amd64/qstat -u $user | awk '{if(\$0 !~ /^job/ && \$0 !~ /^-/) print \$1}'`;
	my @qstat = split /\n/,$qstat;
	my %hash;
	foreach(@qstat){
		$hash{$_} = 1;
	}
	return (\%hash);
}


=cut
### waiting for all jobs finished! ###
foreach $subsh (@All_subsh){
	chdir($dir);

	### check whether the job has been completed! ###
	while (1){
		#	if (-f "*.o*"){
		if (-f "work*.sh.o*"){
			open IN, "work*.sh.o*" ||die "Can't open the file:$!\n";
			my $content;
			$content = join("",<IN>);
			while (1){
				last if ($content =~ /This-Work-is-Completed!/);
				sleep 10;
			}
			close IN;
			last;
		}else{
			sleep 5;
	    }
	}
}
=cut		
