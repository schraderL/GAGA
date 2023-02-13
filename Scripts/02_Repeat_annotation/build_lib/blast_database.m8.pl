#!/usr/bin/perl

=head1 Name

blast_database.pl  --  the pipeline to run blast against several databases

=head1 Description

This program is the pipeline to run blast. Note that the databases except user defined 
database must be formated previously, better to use formatdb in the same version with
blastall invoked in this program. The "--program" option must be set corrected corresponding
to the databases. 

Nucleotide database: nt,
Protein database: nr, swissprot, trembl, cog, kegg, 

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 3.0,  Date: 2008-12-11
  Note:

=head1 Usage
  
  perl blast_database.pl [options] <sequences.fa>
  --program      set the program of blastall, default blastp
  --evalue       set the E-val for blast, default 1e-5
  --nr           search against Nr database
  --nt           search against Nt database
  --swissprot    search against SwissProt database	
  --trembl       search against TrEMBL database
  --cog          search against Cog database
  --kegg         search against	Kegg database
  --tedna        search against Repbase's RepeatMasker DNA database
  --tepep        search against RepeatProteinMask protein database
  --db <str>     search against user defined database
  --cutf <int>   set the number of subfiles to generate, default=100
  --cpu <int>	 set the cpu number to use in parallel, default=30  
  --node <str>       set the compute node, eg h=compute-0-151,default is not set 
  --run <str>    set the parallel type, qsub, or multi, default=qsub
  --outdir <str>  set the result directory, default="."
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  nohup perl ../bin/blast_database.pl --program blastn --nt --tedna ../input/cucumber.reas.fa -run multi -cpu 2 -outdir cucumber_reas &
  nohup perl ../bin/blast_database.pl --program blastp --nr  --swissprot --trembl --cog  --kegg --tepep ../input/rice_prot3000.fa -outdir ./rice_prot3000 -run multi -cpu 4 &
  

=cut
use FindBin qw($Bin $Script);
use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory
use lib "$Bin/../../common_bin";
use GACP qw(parse_config);

my ($Program,$Evalue,$Nr,$Nt,$Swissprot,$TrEMBL,$Cog,$Kegg,$TEdna,$TEpep,$Userdb,$Cutf,$Cpu,$Node,$Run,$Outdir);
my ($Verbose,$Help);
GetOptions(
	"program:s"=>\$Program,
	"evalue:s"=>\$Evalue,
	"nr"=>\$Nr,
	"nt"=>\$Nt,
	"swissprot"=>\$Swissprot,
	"trembl"=>\$TrEMBL,
	"cog"=>\$Cog,
	"kegg"=>\$Kegg,
	"tedna"=>\$TEdna,
	"tepep"=>\$TEpep,
	"db:s"=>\$Userdb,
	"cutf:i"=>\$Cutf,
	"cpu:i"=>\$Cpu,
	"node:s"=>\$Node,
	"run:s"=>\$Run,
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Program ||= "blastp";
$Evalue ||= 1e-5;
$Cutf ||= 100;
$Cpu ||= 8;
$Run ||= "qsub";
$Outdir ||= `pwd`;
#$Outdir ||= ".";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $seq_file = shift;
my $seq_file_name = basename($seq_file);

my $userdb;
$userdb = basename($Userdb) if(defined $Userdb);
my $config_file = "$Bin/../../config.txt";
my $common_bin = "$Bin/../../common_bin";

my $blast = parse_config($config_file,"blastall")." -F F";;
my $formatdb = parse_config($config_file,"formatdb");

my $nr = parse_config($config_file,"nr")."/nr";
my $nt = parse_config($config_file,"nt")."/nt";
my $swissprot = parse_config($config_file,"swissprot");
my $trembl = parse_config($config_file,"trembl");
my $cog = parse_config($config_file,"cog");
my $kegg = parse_config($config_file,"kegg");
my $tedna = parse_config($config_file,"tedna");
my $tepep = parse_config($config_file,"tepep");

my $fastaDeal = "$common_bin/fastaDeal.pl";
my $blast_parser = "$common_bin/blast_parser.pl";
my $qsub_sge = "$common_bin/qsub-sge.pl";
my $multi_process = "$common_bin/multi-process.pl";

my $blast_shell_file = "./$seq_file_name.blast.$$.sh"; ##只能放在当前目录下
my @subfiles;

$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

`perl $fastaDeal -cutf $Cutf $seq_file -outdir $Outdir`;
@subfiles = glob("$Outdir/$seq_file_name.cut/*.*");

###format the user defined database, judge the type automatically.
if (defined $Userdb) {
		my $seq_type = judge_type($Userdb);
		`$formatdb -p F -o T -i $Userdb` if ($seq_type eq "DNA");
		`$formatdb -p T -o T -i $Userdb` if ($seq_type eq "Protein");
}

##create shell file
open OUT,">$blast_shell_file" || die "fail $blast_shell_file";
foreach my $subfile (@subfiles) {
	print OUT "$blast  -p $Program -e $Evalue   -m 8  -d $nr -i $subfile -o $subfile.nr.blast.m8; \n" if(defined $Nr);
	print OUT "$blast  -p $Program -e $Evalue    -m 8 -d $nt -i $subfile -o $subfile.nt.blast.m8; \n" if(defined $Nt);
	print OUT "$blast  -p $Program -e $Evalue   -m 8 -d $swissprot -i $subfile -o $subfile.swissprot.blast.m8; \n" if(defined $Swissprot);
	print OUT "$blast  -p $Program -e $Evalue    -m 8 -d $trembl -i $subfile -o $subfile.trembl.blast.m8; \n" if(defined $TrEMBL);
	print OUT "$blast  -p $Program -e $Evalue   -m 8 -d $cog -i $subfile -o $subfile.cog.blast.m8; \n" if(defined $Cog);
	print OUT "$blast  -p $Program -e $Evalue    -m 8 -d $kegg -i $subfile -o $subfile.kegg.blast.m8; \n" if(defined $Kegg);
	print OUT "$blast  -p $Program -e $Evalue   -m 8  -d $Userdb -i $subfile -o $subfile.$userdb.blast.m8; \n" if(defined $Userdb);
	print OUT "$blast  -p $Program -e $Evalue   -m 8  -d $tedna -i $subfile -o $subfile.tedna.blast.m8; \n" if(defined $TEdna);
	print OUT "$blast  -p $Program -e $Evalue   -m 8  -d $tepep -i $subfile -o $subfile.tepep.blast.m8; \n" if(defined $TEpep);
}
close OUT;
#run the shell file
if(defined $Node){
`$qsub_sge  --maxjob $Cpu --reqsub --node $Node $blast_shell_file` if ($Run eq "qsub");
}
else{`$qsub_sge  --maxjob $Cpu --reqsub  $blast_shell_file` if ($Run eq "qsub");}
`$multi_process -cpu $Cpu $blast_shell_file` if ($Run eq "multi");


##cat together the result
`cat $Outdir/$seq_file_name.cut/*.nr.blast.m8 > $Outdir/$seq_file_name.nr.blast.m8;` if(defined $Nr);
`cat $Outdir/$seq_file_name.cut/*.nt.blast.m8 > $Outdir/$seq_file_name.nt.blast.m8; ` if(defined $Nt);
`cat $Outdir/$seq_file_name.cut/*.swissprot.blast.m8 > $Outdir/$seq_file_name.swissprot.blast.m8; ` if(defined $Swissprot);
`cat $Outdir/$seq_file_name.cut/*.trembl.blast.m8 > $Outdir/$seq_file_name.trembl.blast.m8; ` if(defined $TrEMBL);
`cat $Outdir/$seq_file_name.cut/*.cog.blast.m8 > $Outdir/$seq_file_name.cog.blast.m8; ` if(defined $Cog);
`cat $Outdir/$seq_file_name.cut/*.kegg.blast.m8 > $Outdir/$seq_file_name.kegg.blast.m8; ` if(defined $Kegg);
`cat $Outdir/$seq_file_name.cut/*.$userdb.blast.m8 > $Outdir/$seq_file_name.$userdb.blast.m8; ` if(defined $Userdb);
`cat $Outdir/$seq_file_name.cut/*.tedna.blast.m8 > $Outdir/$seq_file_name.tedna.blast.m8; ` if(defined $TEdna);
`cat $Outdir/$seq_file_name.cut/*.tepep.blast.m8 > $Outdir/$seq_file_name.tepep.blast.m8; ` if(defined $TEpep);

=head
##covert to table format, and keep all the hits
`$blast_parser   $Outdir/$seq_file_name.nr.blast > $Outdir/$seq_file_name.nr.blast.tab;` if(defined $Nr);
`$blast_parser   $Outdir/$seq_file_name.nt.blast > $Outdir/$seq_file_name.nt.blast.tab; ` if(defined $Nt);
`$blast_parser   $Outdir/$seq_file_name.swissprot.blast > $Outdir/$seq_file_name.swissprot.blast.tab; ` if(defined $Swissprot);
`$blast_parser   $Outdir/$seq_file_name.trembl.blast > $Outdir/$seq_file_name.trembl.blast.tab; ` if(defined $TrEMBL);
`$blast_parser   $Outdir/$seq_file_name.cog.blast > $Outdir/$seq_file_name.cog.blast.tab; ` if(defined $Cog);
`$blast_parser   $Outdir/$seq_file_name.kegg.blast > $Outdir/$seq_file_name.kegg.blast.tab; ` if(defined $Kegg);
`$blast_parser   $Outdir/$seq_file_name.$userdb.blast > $Outdir/$seq_file_name.$userdb.blast.tab; ` if(defined $Userdb);
`$blast_parser   $Outdir/$seq_file_name.tedna.blast  >  $Outdir/$seq_file_name.tedna.blast.tab; ` if(defined $TEdna);
`$blast_parser   $Outdir/$seq_file_name.tepep.blast  >  $Outdir/$seq_file_name.tepep.blast.tab; ` if(defined $TEpep);


##covert to table format, and get the the best hit
`$blast_parser  -tophit 1 -topmatch 1   $Outdir/$seq_file_name.nr.blast > $Outdir/$seq_file_name.nr.blast.tab.best;` if(defined $Nr);
`$blast_parser  -tophit 1 -topmatch 1   $Outdir/$seq_file_name.nt.blast > $Outdir/$seq_file_name.nt.blast.tab.best; ` if(defined $Nt);
`$blast_parser  -tophit 1 -topmatch 1   $Outdir/$seq_file_name.swissprot.blast > $Outdir/$seq_file_name.swissprot.blast.tab.best; ` if(defined $Swissprot);
`$blast_parser  -tophit 1 -topmatch 1   $Outdir/$seq_file_name.trembl.blast > $Outdir/$seq_file_name.trembl.blast.tab.best; ` if(defined $TrEMBL);
`$blast_parser  -tophit 1 -topmatch 1   $Outdir/$seq_file_name.cog.blast > $Outdir/$seq_file_name.cog.blast.tab.best; ` if(defined $Cog);
`$blast_parser  -tophit 1 -topmatch 1   $Outdir/$seq_file_name.kegg.blast > $Outdir/$seq_file_name.kegg.blast.tab.best; ` if(defined $Kegg);
`$blast_parser  -tophit 1 -topmatch 1   $Outdir/$seq_file_name.$userdb.blast > $Outdir/$seq_file_name.$userdb.blast.tab.best; ` if(defined $Userdb);
`$blast_parser  -tophit 1 -topmatch 1   $Outdir/$seq_file_name.tedna.blast  >  $Outdir/$seq_file_name.tedna.blast.tab.best; ` if(defined $TEdna);
`$blast_parser  -tophit 1 -topmatch 1   $Outdir/$seq_file_name.tepep.blast  >  $Outdir/$seq_file_name.tepep.blast.tab.best; ` if(defined $TEpep);
=cut
####################################################
################### Sub Routines ###################
####################################################


##judge DNA or protein automatically
sub judge_type {
	my $file = shift;
	my $sequence;
	open(IN, $file) || die ("can not open $file\n");
	while (<IN>) {
		next if(/^>/);
		$sequence .= $_;
		$sequence =~ s/\s//sg;
		$sequence =~ s/-//sg;
		last if(length($sequence) >= 100000);
	}
	close(IN);
	my $base_num = $sequence=~tr/ACGTNacgtn//;
	my $type = ($base_num / length($sequence) > 0.9) ? "DNA" : "Protein";
	return $type;
}

