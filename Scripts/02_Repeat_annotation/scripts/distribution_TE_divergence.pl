#!/usr/bin/perl -w

=head1 Name

	distribution_TE_distribution_new.pl

=head1 Programm Description

	This Program is designed to be draw the distribution of repeat diveragence.

=head1 Contact & Version

	Author: XiongZijun, xiongzijun@genomics.cn
	Version: 1.0

=head1 Command-line Option

	perl distribution_TE_distribution_new.pl <.out> <length>| [options]
	--x_left_border 	set the left border of X axis, default=60
	--x_right_border 	set the right border of X axis, default=540
	--y_up_border 		set the up border of Y axis, default=50
	--y_down_border 	set the down border of Y axis, default=410
	--x_start 		set the start number of X, default=0
	--y_start 		set the start number of Y, default=0
	--x_end 		set the end number of X, default=40
	--y_end 		set the end number of Y, defaul=2
	--x_step 		set the step number of X, defaul=10
	--y_step 		set the step number of Y, defaul=0.5
	--help			output help information to screen

=head1 Usage Exmple

	perl distribution_TE_distribution_new.pl repeatmask.out total_len  --x_end 60 --x_step 10 --y_end 2  --y_step  0.5

=cut



use strict;
use FindBin qw($Bin $Script);
use lib "$Bin/SVG";
use SVG;
use Getopt::Long;

BEGIN {
    my $script_path = "/opt/perl/lib/site_perl/5.14.2/" if(-d "/opt/perl/lib/site_perl/5.14.2/");
    push @INC, $script_path;
}
use SVG;


my ($x1_border,$x2_border,$y1_border,$y2_border,$x_start,$x_end,$y_start,$y_end,$x_step,$y_step,$Help);
GetOptions(
	"--x_left_border:s"=>\$x1_border,
	"--x_right_border:s"=>\$x2_border,
	"--y_up_border:s"=>\$y1_border,
	"--y_down_border:s"=>\$y2_border,
	"--x_start:s"=>\$x_start,
	"--y_start:s"=>\$y_start,
	"--x_end:s"=>\$x_end,
	"--y_end:s"=>\$y_end,
	"--x_step:s"=>\$x_step,
	"--y_step:s"=>\$y_step,
	"--help"=>\$Help
);

die `pod2text $0` if @ARGV < 2 or $Help;


my $svg = SVG->new('width',680,'height',540);
#my ($x1_border,$x2_border,$y1_border,$y2_border) = (60,540,50,410);

$x1_border ||= 100;
$x2_border ||= 580;
$y1_border ||= 50;
$y2_border ||= 450;

#my ($x_start,$x_end,$y_start,$y_end) = (0,40,0,0.4);
$x_start ||= 0;
$x_end ||= 40;
$y_start ||= 0;
$y_end ||= 2;
#my ($x_step,$y_step) = (10,0.1);

$x_step ||= 10;
$y_step ||= 0.5;



my $file = shift;
my $total_len = shift;

open IN,"<$file" or die "fail to open $file:$!\n";
my %hash;
while(<IN>) {

#	    score  div               scaffold    begin      end                       repeat class
	if(/\d+\s+(\S+)\s+\S+\s+\S+\s+ \S+   \s+ (\d+) \s+ (\d+)  \s+\S+\s+\S+\s+\S+\s+   (\S+)/x) {
		my $div = int $1;
		next unless($div >0);
		my $len = $3-$2+1;
		my $class = $4;
		$hash{$div}{"DNA"} += $len if($class =~ /DNA/);
		$hash{$div}{"LINE"} += $len if($class =~ /LINE/);
		$hash{$div}{"LTR"} += $len if($class =~ /LTR/);
		$hash{$div}{"SINE"} += $len if($class =~ /SINE/);
	}
}
close IN;

my $x_unit = ($x2_border-$x1_border)/(($x_end-$x_start)*2);
my $y_unit = ($y2_border-$y1_border)/20;

my $x_covert = $x_unit*2;
my $y_covert = ($y2_border-$y1_border)/($y_end-$y_start);

my @colors = qw/green yellow orange red blue /;
my @elements = qw/DNA LINE LTR SINE/;

################ draw x y border ###################
$svg->rect('x',$x1_border,'y',$y1_border,'width',($x2_border-$x1_border),'height',($y2_border-$y1_border+5),'stroke','black','fill','white','stroke-width',2);

################ draw rect description #############
$svg->rect('x',$x2_border-$x_unit*16-10,'y',$y1_border+5,'width',$x_unit*16+5,'height',$y_unit*4.8+10,'stroke','black','fill','white','stroke-width',2);
foreach (1..4) {
	$svg->rect('x',$x2_border-$x_unit*16-10+2,'y',$y1_border+5+2*$_+($_-1)*1.2*$y_unit,'width',3*$x_unit,'height',1.2*$y_unit,'stroke','black','fill',$colors[4-$_]);
	$svg->text('x',$x2_border-$x_unit*16-10+2+3*$x_unit+5,'y',$y1_border+5+2*$_+$_*1.2*$y_unit,'font-size','20','font-family','ArialNarrow-Bold','fill','black','-cdata',$elements[$_-1]);
}
############### draw x y coordinate ##############
my $x_number = ($x_end-$x_start)/$x_step;
my $y_number = ($y_end-$y_start)/$y_step;
foreach (0..$x_number) {
	$svg->line('x1',$x1_border+$_*$x_step*$x_covert,'y1',$y2_border+5,'x2',$x1_border+$_*$x_step*$x_covert,'y2',$y2_border+10,'stroke','black','fill','black','stroke-width',2);
	$svg->text('x',$x1_border+$_*$x_step*$x_covert,'y',$y2_border+45,'font-size','24','font-family','ArialNarrow-Bold','text-anchor','middle','fill','black','-cdata',$_*$x_step);
}
foreach (0..$y_number) {
	$svg->line('x1',$x1_border,'y1',$y2_border-$_*$y_step*$y_covert,'x2',$x1_border-5,'y2',$y2_border-$_*$y_step*$y_covert,'stroke','black','fill','black','stroke-width',2);
	$svg->text('x',$x1_border-8,'y',$y2_border-$_*$y_step*$y_covert+16,'font-size','24','font-family','ArialNarrow-Bold','text-anchor','end','fill','black','-cdata',$_*$y_step);	
}
############# draw x y description ###############
$svg->text('x',$x1_border+($x2_border-$x1_border)/2,'y',$y2_border+80,'font-size','28','font-family','ArialNarrow-Bold','text-anchor','middle','fill','black','-cdata',"Sequence divergence rate(%)");
my $g = $svg->group("transform"=>"rotate(-90)");
my ($x,$y) = (-($y1_border-20+($y2_border-$y1_border)/2),$x1_border-70);
$g->text('x',$x,'y',$y,'font-size','28','font-family','ArialNarrow-Bold','text-anchor','middle','fill','black','-cdata','Percentage of genome(%)');


############## draw rect distribution ############
#my $total_len = 2308415131;
foreach(1..$x_end) {
	my $y_sine = exists $hash{$_}{"SINE"} ? $hash{$_}{"SINE"} : 0;
	my $y_ltr = exists $hash{$_}{"LTR"} ? $hash{$_}{"LTR"} : 0;
	my $y_line = exists $hash{$_}{"LINE"} ? $hash{$_}{"LINE"} : 0;
	my $y_dna = exists $hash{$_}{"DNA"} ? $hash{$_}{"DNA"} : 0;

	my $y_sine_covert = ($y_sine*$y_covert*100)/$total_len;
	my $y_ltr_covert = ($y_ltr*$y_covert*100)/$total_len;
	my $y_line_covert = ($y_line*$y_covert*100)/$total_len;
	my $y_dna_covert = ($y_dna*$y_covert*100)/$total_len;

	$svg->rect('x',$x1_border+(2*$_-1)*$x_unit,'y',$y2_border-$y_sine_covert,'width',$x_unit,'height',$y_sine_covert,'stroke',$colors[0],'fill',$colors[0]);	
	$svg->rect('x',$x1_border+(2*$_-1)*$x_unit,'y',$y2_border-$y_sine_covert-$y_ltr_covert,'width',$x_unit,'height',$y_ltr_covert,'stroke',$colors[1],'fill',$colors[1]);
	$svg->rect('x',$x1_border+(2*$_-1)*$x_unit,'y',$y2_border-$y_sine_covert-$y_ltr_covert-$y_line_covert,'width',$x_unit,'height',$y_line_covert,'stroke',$colors[2],'fill',$colors[2]);
	$svg->rect('x',$x1_border+(2*$_-1)*$x_unit,'y',$y2_border-$y_sine_covert-$y_ltr_covert-$y_line_covert-$y_dna_covert,'width',$x_unit,'height',$y_dna_covert,'stroke',$colors[3],'fill',$colors[3]);
} 
chomp $file;
 open OUT,">$file\.TEdivergence\.svg" or die "$!";
        print OUT $svg->xmlify();
        close OUT;
#print $svg->xmlify();



