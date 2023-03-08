#!/usr/bin/perl

=head1 Name

run_iprscan51-55.pl  --  the pipeline to run iprscan

=head1 Description

this program is the pipeline to run iprscan, note that this program runs 
very slowly.

For qsub_sge way, one job take up 7.6G memory, i.e, one compute-node only hold two 
iprscan jobs at the same time.

For multi_process way, only one iprscan job can run on the host compute-noed at the same time.

Use multiple -appl flags to specify multiple applications. The possible applications
are listed below: ProDom; PRINTS; Hamap; Pfam; PIRSF; PANTHER; TIGRFAM; SMART; SUPERFAMILY; Gene3D; ProSitePatterns; ProSiteProfiles; Coils; SignalP_EUK; SignalP_GRAM_POSITIVE; SignalP_GRAM_NEGATIVE; Phobius; TMHMM 


=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 2.0,  Date: 2008-5-7
  Mender: He Lijuan Lyndi.He@genomics.cn
  Updata: 2016a, Data: 2016-1-21
  Mender: liuyabin liuyabin@genomics.cn
  Updata: 2016a, Data: 2016-05-06
  Note:

=head1 Usage
  
  perl run_iprscan.pl [options] <proteins.fa>
  --appl <str>    set applications, this option can be used multiple times
  --cuts <int>   set the number of sequences in each cutted file, default=100
  --cpu <int>	 set the cpu number to use in parallel, default=3   
  --run <str>    set the parallel type, qsub, or multi, default=qsub
  --node <str>   set the compute node, default h=compute-0-151
  --seqtype <str> set the sequece type, cds or pep, default=pep
  --queue <str>  set the queue, default no
  --pro_code <str>  set the code for project,default no
  --lines <int>   set number of lines to form a job, default 1
  --outdir <str>  set the result directory, default="."
  --resource <str> set the "vf", default=16.9G
  --help      output help information to screen  

=head1 Exmple

   perl ../bin/run_iprscan51-55.pl -cuts 10 -cpu 20 ../input/rice_prot100.fa 

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory
use Cwd qw(abs_path);

my (@Appl,$Cuts,$Cpu,$Run,$Outdir,$resource,$lines,$Node,$Pro_code, $seqtype);
my ($Help);
my $Queue;
GetOptions(
    "cuts:i"=>\$Cuts,
	"appl:s"=>\@Appl,
	"cpu:i"=>\$Cpu,
	"queue:s"=>\$Queue,
	"pro_code:s"=>\$Pro_code,
	"run:s"=>\$Run,
    "node:s"=>\$Node,
	"seqtype:s"=>\$seqtype,
    "lines:i"=>\$lines,
	"outdir:s"=>\$Outdir,
    "resource:s"=>\$resource,
	"help"=>\$Help
);
$Cuts ||= 100;
$Cpu ||= 30;
#$Node ||= "h=compute-0-43";
$seqtype ||= "pep";
my $Node_para=(defined $Node)?"-node $Node":"";
$Run ||= "qsub";
$Outdir ||= ".";
$Outdir = abs_path($Outdir);
$resource ||= "vf=16.9G";
die `pod2text $0` if (@ARGV == 0 || $Help);
$lines ||= 1;

my $seq_file = shift;
my $seq_file_name = basename($seq_file);
my $species_name = $seq_file_name;
#if ($species_name=~/.gene.pep$/) {
#        $species_name=~s/.gene.pep$//;
#}

my %config;
parse_config("$Bin/../../config.txt",\%config);

my $iprscan = "source $config{iprscan_setup}; $config{iprscan} -goterms -f tsv  ";
my $iprscan_parser = $config{"iprscan_parser"};
my $iprscan_parser_xls = $config{"iprscan_parser_xls"}; 
my $fastaDeal = $config{"fastaDeal.pl"};
my $qsub_sge = $config{"qsub_sge.pl"};
my $multi_process = $config{"multi-process.pl"};
my $cds2aa = $config{"cds2aa.pl"};

##add iprscan applications
if (@Appl == 0){
    $iprscan .= "--appl ProDom --appl PRINTS --appl Pfam --appl SMART --appl PANTHER --appl ProSiteProfiles --appl ProSitePatterns ";
}else{
    foreach  (@Appl) {
        $iprscan .= " -appl $_"
    }
}

my $temp="$Outdir/temp";
$iprscan .= " -T  $temp ";
my $iprscan_shell_file = "$Outdir/$seq_file_name.iprscan.sh";
my @subfiles;

$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

if ($seqtype =~ /cds/i){
	`perl $cds2aa $seq_file > $Outdir/$seq_file_name.pep`;
	$seq_file = "$Outdir/$seq_file_name.pep";
	$seq_file_name = basename($seq_file);
}
`perl $fastaDeal -cuts $Cuts $seq_file -outdir $Outdir`;
@subfiles = glob("$Outdir/$seq_file_name.cut/*.*");

##creat shell file
open OUT,">$iprscan_shell_file" || die "fail $iprscan_shell_file";
foreach my $subfile (@subfiles) {
	print OUT "$iprscan -i $subfile -o $subfile.iprscan; \n";
}
close OUT;

##run the shell file
my $opt;
$opt=" -queue $Queue " if defined $Queue;
$opt.=" -pro_code $Pro_code " if (defined $Pro_code);
$opt.=" -lines $lines " if(defined $lines);
`perl $qsub_sge $opt --maxjob $Cpu --reqsub --resource $resource $Node_para $iprscan_shell_file` if ($Run eq "qsub");
`perl $multi_process -cpu 1 $iprscan_shell_file` if ($Run eq "multi");

##cat together the result
`cat $Outdir/$seq_file_name.cut/*.iprscan > $Outdir/$species_name.iprscan`;

##annotations for the iprscan result
`perl $iprscan_parser_xls $Outdir/$species_name.iprscan $Outdir/$species_name.iprscan.xls`;

##extract the IPR and GO annotations for genes
`perl $iprscan_parser $Outdir/$species_name.iprscan -outdir $Outdir`;

#`rm -r $iprscan_shell_file.*.qsub/temp/*`;
`rm -r $temp/*`;

####################################################
################### Sub Routines ###################
####################################################


##parse the config.txt file, and check the existence of each software
####################################################
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


