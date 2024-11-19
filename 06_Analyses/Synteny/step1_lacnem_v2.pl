#!/usr/bin/perl 
=head1 Description

 lastz, chain, net & maf pipeline;

=head1 Version

 Zijun Xiong, xiongzijun1989@gmail.com
 original developed by Yongli Zeng (zengyongli@genomics.cn)

=head1 Options

 --direction <str>    direction for output files,
                        default "./output";
 --mode <str>         mode selection, "multi" or "single",
                        default "multi";
 --num <int>          split the task into <int> files,
                        default 20;
 --parasuit <str>     easily set parameters suit to define --lpara and --apara.
                        "chimp": for human vs chimp, gorilla, rhesus, marmoset and so on; (near)
                        "chick": for human vs chicken, zebra finch and so on; (far)
                        *** chimp ********************************************************************************
                        * --lpara:
                        *   --hspthresh=4500 --gap=600,150 --ydrop=15000 --notransition
                        *   --scores=chimpMatrix --format=axt
                        * --apara:
                        *   -minScore=5000 -linearGap=medium
                        ******************************************************************************************

                        *** chick ********************************************************************************
                        * --lpara:
                        *   --step=19 --hspthresh=2200 --inner=2000 --ydrop=3400 --gappedthresh=10000
                        *   --scores=birdMatrix --format=axt
                        * --apara:
                        *   -minScore=5000 -linearGap=loose
                        ******************************************************************************************
 --lpara <str>        parameters for lastz,
                        default "--format=axt";
 --apara <str>        parameters for axtChain,
                        default "-linearGap=medium";
 --tn <str>           name for target in maf,
                        default as filename of target.fa;
 --qn <str>           name for query in maf,
                        default as filename of query.fa;
 --step <str>         1:initial; 2:split; 3:lastz; 4:chain; plus; merge; 5:net; 6:maf;
                        default: "123456plusmerge";
 --qsub               use qsub-sge.pl to run lastz & chain;

 --qpara <str>        parameters for qsub-sge.pl,
                        default "--maxjob 50 --resource vf=1.5g --reqsub --convert no";
 --help               show this page;

=head1 Usage

 nohup perl lacnem.pl target.fa query.fa &

=cut

# 2012-08-13	fix path problems;
use strict;
#use FindBin qw($RealBin);
my $RealBin = "/run/media/dell/storage1/User/xiongzj/GAGA_project/manuscript/02.scripts/lastz/bin/";
use File::Basename;
use Getopt::Long;
use Cwd;

#my $pathway = $RealBin;
my $pathway = "/run/media/dell/storage1/User/xiongzj/GAGA_project/manuscript/02.scripts/lastz/lib/";

my ($direction, $mode, $num, $parasuit, $lpara, $apara, $tn, $qn, $step, $qsub, $qpara, $help);
GetOptions(
	"direction:s"	=> \$direction,
	"mode:s"		=> \$mode,
	"num:i"			=> \$num,
	"parasuit:s"    => \$parasuit,
	"lpara:s"    	=> \$lpara,
	"apara:s"       => \$apara,
	"tn:s"			=> \$tn,
	"qn:s"			=> \$qn,
	"step:s"		=> \$step,
	"qsub"			=> \$qsub,
	"qpara:s"	    => \$qpara,
	"help"			=> \$help,
)or die "Unknown option!\n";

my $fasta_target =shift;
my $fasta_query =shift;

die `pod2text $0` if (!($fasta_target && $fasta_query) || $help);
die "$fasta_target not found!\n" unless (-e $fasta_target);
die "$fasta_query not found!\n" unless (-e $fasta_query);

$direction ||= "./output";
$direction =~ s/\/$//;
$mode ||= "multi";
$num ||= 20;
if ($parasuit eq "chimp"){
	$lpara = "--hspthresh=4500 --gap=600,150 --ydrop=15000 --notransition --scores=$pathway/chimpMatrix --format=axt";
	$apara = "-minScore=5000 -linearGap=$pathway/medium";
}elsif ($parasuit eq "chick"){
	$lpara = "--step=19 --hspthresh=2200 --inner=2000 --ydrop=3400 --gappedthresh=10000 --scores=$pathway/birdMatrix --format=axt";
	$apara = "-minScore=5000 -linearGap=$pathway/loose";
}elsif ($parasuit ne ""){
	die "error --parasuit option value: $parasuit\n";
}
$lpara ||= "M=0 K=2200 L=6000 Y=3400 Q=$pathway/HoxD55 E=30 H=2000 O=400 T=1 --format=axt";
$apara ||= "-linearGap=$pathway/loose";

my $n1 = basename($fasta_target);
my $n2 = basename($fasta_query);
$n1 =~ s/\..*$//;
$n2 =~ s/\..*$//;
$tn ||= $n1;
$qn ||= $n2;
$step ||= "123456";
$qpara ||= "--maxjob 50 --resource vf=1.5g --convert no";
$tn .= ".";
$qn .= ".";

#die "$qsub_para\n";
# step1: initial
my $time = time();
if ($step =~ /1/){
	testmkdir("$direction/1.target");

	`$pathway/faToTwoBit $fasta_target $direction/target.2bit`;
	`$pathway/faToTwoBit $fasta_query $direction/query.2bit`;

	`$pathway/faSize $fasta_target -detailed > $direction/target.sizes`;
	`$pathway/faSize $fasta_query -detailed > $direction/query.sizes`;
	
}

#step2: split
if ($step =~ /2/){
	if ($mode eq "single"){
		`$RealBin/split_fasta.pl $fasta_target $direction/1.target 999999999 avg yes`;
	}else{
		`$RealBin/split_fasta.pl $fasta_target $direction/1.target $num avg yes`;
        my @tdir = `ls $direction/1.target/*.fasta`;
        foreach my $tt(@tdir) {
            chomp $tt;
            my $fa = $tt;
            chomp $fa;
            $fa =~ s/\.fasta$//;
            `$pathway/faToTwoBit $tt $fa.2bit`;
        }
	}
}

# step3: run lastz
if ($step =~ /3/){
	testmkdir("$direction/2.lastz");
	testmkdir("$direction/3.chain");
	my $path =getcwd();
	if ($direction =~ /^\//){
		$path = "";
	}
	$path .= "/";

	open SH, ">$direction/lastzshell.sh" or die "can't open lastzshell.sh\n";

	my @tdir = `ls $direction/1.target/*.2bit`;

	my $du = 0;
	foreach my $tt(@tdir){
		my $fa = $tt;
		chomp($fa);
		$fa =~ s/2bit$/fasta/;
		my $size = -s "$fa";
		$du += $size;
	}
	my $seg = int($du / $num);
	my $count = 0;
	my $i = 1;

	foreach my $tt(@tdir){
		chomp($tt);
		my $tinput = $path . "$tt";
		my $tinput_raw = $tinput;
		if ($mode eq "multi"){
			my $fa = $tt;
			chomp($fa);
			$fa =~ s/2bit$/fasta/;
			my $faline = `grep -c ">" $fa`;
			chomp($faline);
			$tinput .= "[multi]" if ($faline != 1);
		}
		my $qinput = $path . "$direction/query.2bit";
		my $nameaxt = basename($tt);
		$nameaxt =~ s/2bit$/axt/;
		$nameaxt = $path . "$direction/2.lastz/$nameaxt";
		my $namechain = basename($nameaxt);
		$namechain =~ s/axt$/chain/;
		$namechain = $path. "$direction/3.chain/$namechain";
		print SH "/run/media/dell/data/public/software/lastz/src/lastz $tinput $qinput $lpara > $nameaxt;\n";
	}
	close SH;
	`perl /run/media/dell/data/User/xiongzj/bin/multi-process.pl $direction/lastzshell.sh --cpu 30`;
}

# step4: chain
if ($step =~ /4/){
	testmkdir("$direction/3.chain");
	my @chr_lastz = `ls $direction/2.lastz`;
	foreach (@chr_lastz){
		chomp;
		my $name = basename($_);
		my $tname = $name;
		$tname =~ s/axt$/2bit/;
		`$pathway/axtChain $apara $direction/2.lastz/$_ $direction/1.target/$tname $direction/query.2bit $direction/3.chain/$name.chain`;
	}
}

# step_plus: plus
if ($step =~ /plus/) {
    testmkdir("$direction/3.chain_plus");
    my @chain = `ls $direction/3.chain/*.chain`;
    foreach (@chain) {
        chomp;
        my $name = basename($_);
        my $tname = $name;
        $tname = s/\.axt\.chain//;
        `/run/media/dell/data/public/software/GenomeAlignmentTools/GenomeAlignmentTools-master/src/patchChain.perl $direction/3.chain/$name $direction/target.2bit $direction/query.2bit $direction/target.sizes $direction/query.sizes -outputDir $direction/3.chain_plus/$name\_pslOutput -numJobs 1 -jobDir $direction/3.chain_plus/$name\_jobs -lastzParameters "--format=axt M=0 K=2400 L=3000 Y=3400 Q=/run/media/dell/storage1/User/xiongzj/GAGA_project/manuscript/02.scripts/lastz/lib//HoxD55 E=30 H=2000 O=400 T=1 "`;
#        `nohup sh $direction/3.chain_plus/$name\_jobs/job0.csh`;
#        `cat $direction/3.chain_plus/$name\_jobs/job0.csh >> $direction/run_patchChain.sh`;
    }
    close OUT;
}

if ($step=~/merge/) {
    testmkdir("$direction/3.chain_plus_merge");
    my @psl = `ls $direction/3.chain_plus/*pslOutput/*.psl`;
    foreach (@psl) {
        chomp;
        `/run/media/dell/storage1/User/xiongzj/GAGA_project/manuscript/02.scripts/lastz/lib/pslToChain $_ $_.chain`;
    }
    `$pathway/chainMergeSort $direction/3.chain/*.chain $direction/3.chain_plus/*pslOutput/*.chain > $direction/3.chain_plus_merge/all.chain`;
    `/run/media/dell/data/public/software/GenomeAlignmentTools/GenomeAlignmentTools-master/src/RepeatFiller.py -c $direction/3.chain_plus_merge/all.chain -T2 $direction/target.2bit -Q2 $direction/query.2bit > $direction/3.chain_plus_merge/all.chain.filler.chain`
}


# step5: net
if ($step =~ /5/){
	testmkdir("$direction/4.prenet");
	testmkdir("$direction/5.net");
    `$pathway/chainSort $direction/3.chain_plus_merge/all.chain.filler.chain $direction/4.prenet/all.chain.filler.sort.chain`;
    `$pathway/chainCleaner $direction/4.prenet/all.chain.filler.sort.chain -tSizes=$direction/target.sizes -qSizes=$direction/query.sizes $direction/target.2bit $direction/query.2bit $direction/4.prenet/all.chain.filler.sort.chain.cleaned.chain $direction/4.prenet/removedSuspects.bed -linearGap=loose`;
	`$pathway/chainNet -tNibDir=$direction/target.2bit -qNibDir=$direction/query.2bit  -rescore $direction/4.prenet/all.chain.filler.sort.chain.cleaned.chain $direction/target.sizes $direction/query.sizes $direction/5.net/UCSC.target.net $direction/5.net/UCSC.query.net -linearGap=loose`;
    `$pathway/netSyntenic $direction/5.net/UCSC.target.net $direction/5.net/UCSC.target.net.temp`;
    `$pathway/NetFilterNonNested.perl -doUCSCSynFilter -keepSynNetsWithScore 5000 -keepInvNetsWithScore 5000 $direction/5.net/UCSC.target.net.temp > $direction/5.net/UCSC.target.filter.net`;
}

# step6: maf
if ($step =~ /6/){
	testmkdir("$direction/6.net_to_axt");
	testmkdir("$direction/7.maf");
	`$pathway/netToAxt $direction/5.net/target.net $direction/4.prenet/all_sort.chain $direction/target.2bit $direction/query.2bit $direction/6.net_to_axt/all.axt`;
	`$pathway/axtSort $direction/6.net_to_axt/all.axt $direction/6.net_to_axt/all_sort.axt`;
	`$pathway/axtToMaf -tPrefix=$tn -qPrefix=$qn $direction/6.net_to_axt/all_sort.axt $direction/target.sizes $direction/query.sizes $direction/7.maf/all.maf`;
}

$time = time() - $time;
my $hour = int($time / 3600);
my $minute = int(($time - $hour * 3600) / 60);
my $second = int($time % 60);
print "\nTotal time cost: $hour h $minute m $second s.\n";

#################
sub testmkdir(){
	my $dir = shift;
	if (-e $dir){
		warn "Warning: Folder ($dir) exists! all files in it will be deleted!\n";
#		`rm -r $dir`;
	}else{
	    `mkdir -p $dir`;
    }
}

