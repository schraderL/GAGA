#! /usr/local/bin/perl  -w
=head1 Name
build_library.pl  --the pipeline to filter the repeat library predicted by denovo repeat finding software
Version: 2016a
Mender: He Lijuan  Lyndi.He@genomics.cn
=head1 Usage
perl  build_library.pl  <cpu_num> <denovo_predicted_library.fa> <Piler|RepeatScout|LTR_FINDER>
=cut

use strict;
use FindBin qw($Bin $Script);
use lib "$Bin/../";
use GACP qw(parse_config);
use Getopt::Long;
use Cwd;

my($Cpu,$Run,$Outdir,$Queue,$Pro_code);
GetOptions(
	"cpu:s"=>\$Cpu,
	"run:s"=>\$Run,
	"outdir:s"=>\$Outdir,
	"queue:s"=>\$Queue,
	"pro_code:s"=>\$Pro_code
);

$Cpu ||= 3;
$Run ||= "qsub";
$Outdir ||= ".";
$Outdir =~s/\/$//;

if (@ARGV != 2) {
	warn "#Usage: perl $0 <lib> <type> \n";
	exit;
}

my $lib = shift;
my $type = shift;
my $curr_path = getcwd();

my $config_file="$Bin/../../config.txt";
my $formatdb=parse_config($config_file,"formatdb");
my $fastaDeal=parse_config($config_file,"fastaDeal.pl");
my $blast=parse_config($config_file,"blastall");
my $qsub_sge=parse_config($config_file,"qsub_sge.pl");
my $multi=parse_config($config_file,"multi-process.pl");

`perl $Bin/filter_short_seq.pl $lib $type`;
#open OUT, ">cdhit.sh" || die "failed to open :$!\n";
#print OUT "/share/app/cdhit/4.8.1/cd-hit -i $lib.100bp -o $curr_path/$lib.100bp.unredundance\n";
#print OUT "/hwfssz4/BC_COM_FP/bc_paea/Annotation_2018/software/repeat/RepeatModeler-2.0.1/RepeatClassifier -consensi $lib.100bp.unredundance\n";
#close OUT;

my $QP_para = " ";
$QP_para.="--queue $Queue " if (defined $Queue);
$QP_para.="--pro_code $Pro_code " if (defined $Pro_code);

#`perl $qsub_sge $QP_para --lines 1 cdhit.sh`;

my $path = ".";
#chomp($path);
#`perl $Bin/RepeatClassifier.divided.qsub.pl $QP_para $lib.100bp.unredundance $path`;

`perl $Bin/extract_unknown.pl $lib.100bp.unredundance.classified`;

`perl $Bin/blast_database.pl $QP_para --program blastx --swissprot  --cpu $Cpu $lib.100bp.unredundance.classified.unknown`;
 
`perl $Bin/filter_gene.pl $lib.100bp.unredundance.classified.unknown.swissprot.blast.tab $lib.100bp.unredundance.classified.unknown`;
`cat $lib.100bp.unredundance.classified.unknown.filter_gene $lib.100bp.unredundance.classified.known > $lib.final.library`;

`perl $Bin/filter_Nrich.pl  $lib.final.library $type`;
