#!/usr/bin/perl

=head1 Name

  go_svg.pl -- draw WEGO-like Gene Ontology Classification figure

=head1 Version

  Author:Hanc,Zheng  zhenghch@genomics.org.cn
  Version:2.1  ,Date: 5/26/2008

=head1 Description

 This stcript would read from go class file and *.go.id file and then call discribute_svg.pl(author:List) to draw .
 It also has the use to call parse_GO.pl(author:Fanw) to read GO flatfiles, and get GO iterms: (GO id, GO description) and generate a  class file;

=head1 Usage
 
 perl go_svg.pl [options] name1.go.id [name2.go.id.....]
 --outdir       output directory(default current directory)
 --name         output name(default out.*)
 --class        change class file(only calculate go_id in the file;default go.class)
 --alias        change alias file (only calculate go_id in the file;default go.alias)
 --displayturn  change the display turn of  class.if you don't set this parameter it will be sorted by ASCII
                (biological_process   cellular_component   molecular_function)
 --color        change color of the rectangle
 --mark	        change mark text
 --note         change the title of the picture
 --verbose      output verbose information to screen
 --help         output help information to screen

=head1 Example
 
 ##use default to draw
 perl go_svg.pl name1.go.id [name2.go.id.....]			
 
 ##use special defined options, mainly for file
 perl go_svg.pl name1.go.id [name2.go.id.....] -outdir  special_outdir -class class_file -alias alias_file -name output_prefix
 
 ##use special defined options, mainly for figure
 perl go_svg.pl name1.go.id [name2.go.id .....] -displayturn first_class:second_class:... -color  color1:color2:... -mark mark1:mark2... -note note -y_name name

=cut

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname);
use FindBin qw ($Bin $Script);
my ($class_file,$alias_file);
my($turn,$mark,$note,$color,$outdir,$name);
my ($verbose,$help);

my $error_out;				##print go_id which can't find in the class file			
my $txt_out;					##the statistics of the data
my $svg_out;

my @cmp;					##$cmp[$i]{$go_id}=[$ref_id...];
my @total;					##$total[$i] store total ref_gene in the file;
my @num;					## $num{$class_st}{$class_se}{$gene_id}=\d+;

my $xend;					##x end
my $is_goid;					##check if go_id in the class_file;

my @default_color=("#99330C","#9B97C8","#D8A51E","green","blue","red");
my @default_mark=qw(cmp1 cmp2 cmp3 cmp4 cmp5 cmp6);

my %mod=(
	"general"=>"Type:Rect\nWidth:1500\nHeight:500\nBothYAxis:1\nMultiRY:1\n#MultiY:1\nScaleLen:8\nWholeScale:1\nMovePer:0.125\nXStep:1\nRYStep:6\nXScalePos:0.5\nXScaleRoate:75\nXStart:0\nRYStart:0.01\nRYEnd:100\nYNeedLog:10i\nYStep:25\nYEnd:100\nMarkStyle:horizonal\nMarkPos:left\nMarkNoBorder:1\nFontFamily:Arial\nFontSize:30\nXUnit:1",
	"YStart"=>"0.01",
	"Note"=>"Comparison of GO classification",
	"y"=>"Percent of genes",
	"RY"=>"Number of genes",
);

my $dir_up=$1 if $Bin=~/(.*)\//;
GetOptions(
	"outdir:s"=>\$outdir,
	"name:s"=>\$name,
	"class:s"=>\$class_file,
	"alias:s"=>\$alias_file,
	"displayturn:s"=>\$turn,
	"color:s"=>\$color,
	"mark:s"=>\$mark,
	"note:s"=>\$note,
	"verbose"=>\$verbose,
	"help"=>\$help
);

die `pod2text $0` if (@ARGV == 0 || $help);

my %config;
parse_config("$Bin/config.txt",\%config);

$outdir ||= ".";
$outdir =~ s/\/$//;
$name ||= "out";
$class_file ||= $config{go_class};
$alias_file ||= $config{go_alias};

my $distribute_svg = $config{distribute_svg};

my $svg_file="$outdir/$name.svg";
my $lst_file="$outdir/$name.lst";
my $error_file="$outdir/$name.error";
my $txt_file="$outdir/$name.txt";

my @color=split /:/,$color if $color;
@color=(@color >= @ARGV)?@color:@default_color;
my @mark=split /:/,$mark if $mark;
@mark=(@mark>=@ARGV)?@mark:@default_mark;		##if parameter is not enough then set it to default

$mod{"Note"}=$note if $note;

my ($class,$is_goid)=read_class($class_file);		##read class file

my @turn=split /:/,$turn if $turn;
@turn=(sort keys %$class) unless (@turn>0 && @turn eq scalar keys %$class);

my $alias=read_alias($alias_file);				##read alias in order to transform the synonymy go_id

my $temp;
for (my $i=0;$i<@ARGV;$i++){
	($total[$i],$cmp[$i],$temp)=read_go($ARGV[$i],$is_goid,$alias);	##read go.id file
	$error_out.="$temp";
	$num[$i]=cal_class($class,$cmp[$i]);		##calculate  every class' gene_id
}


###generate lst file
$svg_out="$mod{'general'}\n";
$svg_out.="YStart:$mod{'YStart'}\n";
$svg_out.="UnitPer:@{[0.6/@ARGV]}\nOffsetPer:@{[0.6/@ARGV]}\n";	##auto ajust the width of rectangle by @ARGV(number of input file)
$svg_out.="Note:$mod{'Note'}\ny:$mod{'y'}\nRY:$mod{'RY'}\n";

($xend,$temp)=set_Xscale($class);			##add the x-axis mark to the lst file,and count x-end to fit the data
$svg_out.="XEnd:$xend\n$temp\n";

for (my $i=0;$i<@num;$i++){					##this seg is tend to set RY 	
	$svg_out.="Color:$color[$i]\nMark:$mark[$i]\nYMark:r\nStart:0\nEND:2\nStep:1\n";
	my $half=$total[$i]/100;				##
	$svg_out.="Scale:\n0\n$half\n$total[$i]\n:END\n";
	$temp=set_value($class,$num[$i],$total[$i],$mod{'YStart'});
	$svg_out.="$temp\n";
}


open OUT,">$lst_file" or die "fail $lst_file:$!";
print OUT $svg_out;
close OUT;

####
if($error_out){
	open OUT,">$error_file" or die "fail $error_file:$!";
	print OUT $error_out;
	close OUT;
}

$txt_out=set_text($class,\@num,\@total);
open OUT,">$txt_file" or die "fail $txt_file:$!";
print OUT $txt_out;
close OUT;

$svg_out="";
$txt_out="";
$error_out="";
system ("perl $distribute_svg $lst_file $svg_file");

############################
#calculate gene_id of every class
 sub cal_class{
	my ($class,$cmp)=@_;
	my %cmp_num;
	for  my $class_st(sort keys %$class){
		for my $class_se(sort keys %{$class->{$class_st}}){
			for my $go_id(sort @{$class->{$class_st}{$class_se}}){
				if (exists $cmp->{$go_id}){
					for my $gene_id(@{$cmp->{$go_id}}){
						$cmp_num{$class_st}{$class_se}{$gene_id}++;
					}
				}
			}
		}
	}
return \%cmp_num;
}

#############################
##read from class file
sub read_class{
	my $class_file=shift;
	my %class;
	my %is_goid;
	open IN,$class_file or die "fail $class_file:$!";
	while(<IN>){
		chomp;
		my ($class_st,$class_se,$go_id)=(split /\t+/)[0..2];
		$class_st=~s/^\s+//;$class_st=~s/\s+$//;		##delete blank begin and end  
		$class_se=~s/^\s+//;$class_se=~s/\s+$//;		##
		$class_se=~s/activity$//g;
		push @{$class{$class_st}{$class_se}},$go_id;
		$is_goid{$go_id}++;
	}
	close IN;
return (\%class,\%is_goid);
}	

############################
#read from gene-go_id file 
sub read_go{
	my %GO;
	my ($go_file,$is_goid,$alias)=@_;
	my $num;
	my @error;
	my $error;
	my $error_out;
	open IN,$go_file or die "fail $go_file:$!";
	while(<IN>){
		chomp;
		my ($ref_id,@go_id)=split/\t+/;
		$num++;
		for (my $i=0;$i<@go_id;$i++){
			if (exists $is_goid->{$go_id[$i]}){			##go_id is in class_file;
				push @{$GO{$go_id[$i]}},$ref_id;
			}else{		
					if (exists $alias->{$go_id[$i]}){	
						push @{$GO{$alias->{$go_id[$i]}}},$ref_id; ##chang it to go_id used now
					}else{
					
						push @error,$go_id[$i];	##can't find it in class_file and it has not synonymy go_id in alias file;
						$error++;
					}
				}
		}
	}
	close IN;
	
	if($error ){
		$error_out.="Error in file $go_file:can't find go id. You have chosen part of the go_id to compare or the class file is  different version of $go_file.Please check it.\n";
		$error_out.="$_\n" foreach ( @error);
	}

return ($num,\%GO,$error_out);
}

############################
#read alias file 
sub read_alias{
	my $alias_file=shift;
	my %alias;
	open IN,$alias_file or die "fail $alias_file:$!";
	while(<IN>){
		chomp;
		my ($go_id,@alias_id)=split /\s+/;
		$alias{$_}=$go_id foreach (@alias_id);
	}
return \%alias;
}

#############################################
#set X-axis info
sub set_Xscale{
	my $class=shift;
	my $out_scale="Scale:\n";
	my $out_group="Group:\n";
	my $out;	
	my $xend;	
	for  my $class_st(@turn){
		my $n=keys %{$class->{$class_st}};
		$out_group.="$n:$class_st\n";
		for my $class_se(sort keys %{$class->{$class_st}}){
			$out_scale.="$class_se\n";
			$xend++;
		}
	}
	$out="$out_scale:End\n$out_group:End\n";
return ($xend,$out);
}
			
	
#############################################
#calculate value of every rectangle
sub set_value{
	my ($class,$cmp_num,$max,$min)=@_;
	my $out;
	my $i=0;
	for  my $class_st(@turn){
		for my $class_se(sort keys %{$class->{$class_st}}){
			my $num_cmp=(exists $cmp_num->{$class_st}{$class_se})? (scalar keys %{$cmp_num->{$class_st}{$class_se}})/$max*100:$min;
			$num_cmp=($num_cmp<$min)?$min:$num_cmp;
			$out.="$i:$num_cmp\n";
			$i++;
		}
	}
return $out;
}
############################################
#generate text file
sub set_text{
	my ($class,$num,$total)=@_;
	my $txt_out;
	for  my $class_st(sort keys %$class){
		for my $class_se(sort keys %{$class->{$class_st}}){
			$txt_out.="$class_st\t$class_se\t";
			for (my $i=0;$i<@$num;$i++){
				my $num_cmp=(exists $num->[$i]{$class_st}{$class_se})?scalar (keys %{$num->[$i]{$class_st}{$class_se}}):0;
				$txt_out.="$num_cmp\t@{[$num_cmp/$total->[$i]]}\t";
			}
			$txt_out.="\n";
		}
	}
return $txt_out;
}
##########################################


##parse the software.config file, and check the existence of each software
####################################################
sub parse_config{
	my $conifg_file = shift;
	my $config_p = shift;
	
	my $error_status = 0;
	open IN,$conifg_file || die "fail open: $conifg_file";
	while (<IN>) {
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
