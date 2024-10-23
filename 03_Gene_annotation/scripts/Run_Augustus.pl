#!/usr/bin/perl
use strict;

=head1 Name

	Run_Augustus.pl -- a small pipeline to run augustus.

=head2 Usage
	
	perl run_augustus.pl [option] <seq.fa>
	--maxjob <int>   set cpu [12]
	--prog <str>     set Augustus path
	--species <str>  set specises of parameter to run Augustus  [arabidopsis]
	--run  T/F       T: run pipelie ;F :just generate shell file [T]
	--queue <str>    set the queue, default no
	--pro_code <str> set the project code,default no
	--outdir <str>   set result dir.default ./
	--vf <str>       set the vf
        --node <str>     set the compute node
	--help           print this information

=cut


use strict;
use File::Basename qw(basename);
use FindBin qw($Bin);
use Getopt::Long;

my ($Outdir,$Help);
my ($Maxjob,$Prog,$Species,$Job_q,$Node,$Pro_code);
my $Vf;
my $Run;
GetOptions(
	"outdir:s"=>\$Outdir,
	"help"=>\$Help,
	"maxjob:i"=>\$Maxjob,
	"prog:s"=>\$Prog,
	"species:s"=>\$Species,
	"queue:s"=>\$Job_q,
	"pro_code:s"=>\$Pro_code,
	"run:s"=>\$Run,
        "node:s"=>\$Node,
	"vf:s"=>\$Vf,
);
die `pod2text $0` if (@ARGV==0 ||$Help );

$Job_q ||= "st.q";
$Pro_code ||= "P18Z10200N0102";
$Outdir=~s/\/$//;
$Outdir=$ENV{PWD}."/".$Outdir if ($Outdir !~/^\//);
mkdir($Outdir) unless (-d $Outdir);
$Run||='T';
#$Node ||="h=compute-0-191";
my $Node_para=(defined $Node)?"-node $Node":'';
my %config;
parse_config("/jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/02.augustus/config/config.txt",\%config);

$Prog ||= $config{"augustus"};
$Species ||="GAGA-0014";
$Maxjob||=12;
my $Para="--species=$Species --AUGUSTUS_CONFIG_PATH=$config{augustus_config} --uniqueGeneId=true --noInFrameStop=true --gff3=on --strand=both";
#$Job_q||="all.q";
$Vf||="0.5G";

##add by luchx
my $QP_para;
$QP_para.=" --queue $Job_q" if (defined $Job_q);
$QP_para.=" --pro_code $Pro_code" if (defined $Pro_code);

my $Fastafile=shift;
my $Fastafile_name=basename($Fastafile);

my %sub_seq;
my $Total_len=0;
open IN,$Fastafile or die "can not open $Fastafile:$!";
$/=">";<IN>;$/="\n";
while(<IN>) {
	chomp;
	my $head=$_;
	$/=">";
	my $seq=<IN>;
	chomp $seq;
	$/="\n";
	$seq=~s/\s//g;
	$Total_len+=length($seq);
}
close IN;

my $Sub_len=int ($Total_len/$Maxjob);
my $Cur_len=0;
open IN,$Fastafile or die "can not open $Fastafile:$!";
$/=">";<IN>;$/="\n";
while(<IN>) {
        chomp;
        my $head=$_;
        $/=">";
        my $seq=<IN>;
        chomp $seq;
        $/="\n";
        $seq=~s/\s//g;
	$Cur_len+=length($seq);
	my $sub_id=int($Cur_len/$Sub_len);
	$sub_seq{"$sub_id"}.=">$head\n$seq\n";
}


my $i=0;
my $cut_dir="$Outdir/$Fastafile_name.ag.cut";
`rm -r $cut_dir` unless  ( ! -e $cut_dir);
mkdir($cut_dir) or die "Fail to mkdir $cut_dir:$!";
foreach (values %sub_seq){
	open OUT,">$cut_dir/$Fastafile_name.$i" or die"$!";
	print OUT $_;
	close OUT;
	$i++;	
}

my @subfile=glob("$Outdir/$Fastafile_name.ag.cut/*");
open SH,">$Outdir/$Fastafile_name.ag.sh" or die "can not open $Outdir/$Fastafile_name.ag.sh:$!";
foreach (@subfile) {
	print SH "$Prog $Para $_ > $_.augustus; \n";
}
close SH;

if ($Run eq 'T' ){
	`perl $config{"qsub_sge.pl"} $QP_para --reqsub  --resource vf=$Vf  -maxjob $Maxjob $Node_para --convert no $Outdir/$Fastafile_name.ag.sh`;
	`cat $Outdir/$Fastafile_name.ag.cut/*.augustus >$Outdir/$Fastafile_name.augustus`;
	`perl $Bin/Convert_Augustus.pl $Outdir/$Fastafile_name.augustus >$Outdir/$Fastafile_name.augustus.gff`;
	`perl $config{"getGene.pl"} $Outdir/$Fastafile_name.augustus.gff  $Fastafile >$Outdir/$Fastafile_name.augustus.cds`;
	`perl $config{"cds2aa.pl"} $Outdir/$Fastafile_name.augustus.cds >$Outdir/$Fastafile_name.augustus.pep`;
}

############################################
sub parse_config{
        my $conifg_file = shift;
        my $config_p = shift;
        my $error_status = 0;
        open IN,$conifg_file || die "fail open: $conifg_file";
        while (<IN>) {
                next if(/^#/);
                if (/(\S+)\s*=\s*(\S+)/) {
                        my ($software_name,$software_address) = ($1,$2);
                        $config_p->{$software_name} = $software_address;
                        if (! -e $software_address){
                                warn "Non-exist:  $software_name  $software_address\n";
                                $error_status = 1;
                        }
                }
        }
        close IN;
        die "\nExit due to error of software configuration\n" if($error_status);
}
