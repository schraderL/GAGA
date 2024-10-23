#!/usr/bin/perl

=head1 Name

kegg_parser.pl  -- the program to do kegg analysis

=head1 Description

This program extract the ko number and its class, also color pink up on the corresponding ko map.
The input is the table format blast result against KEGG database.

There is one kegg class type that has not been understand:
Protein Families; Cellular Processes and Signaling; Cytoskeleton proteins [BR:ko04812]

=head1 Version

  Author: Sun Juan, sunjuan@genomics.org.cn
  Author: Fan wei, fanw@genomics.org.cn
  Version: 2.0,  Date: 2008-5-21

=head1 Usage

  perl kegg_parser.pl <blast.tab> 
  --outdir <str>  set the result directory, default="Kegg_map"
  --verbose      output running progress information to screen  
  --help         output help information to screen  

=head1 Exmple

  perl kegg_parser.pl mmbao.fasta.allgene.pep.blast.kegg.tab

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use GD;

my ($Outdir);
my ($Verbose,$Help);
GetOptions(
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);

$Outdir ||= "./";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $blast_tab = shift;
my $blast_tab_base = basename($blast_tab);
use Data::Dumper;
my %config;
parse_config("$Bin/../../config.txt",\%config);

my $map_path = $config{"map"};
my $ko_path = $config{"ko"};

my %KO_Class; ##存储全部的KO与其他属性的对应关系
my %KO_EC;
my %KO_Name;
my %KO_Defi;
my %Map_Gene; ##存储实际要画的map图和所包含的基因

$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

##read ko file for useful information
open KO,$ko_path || "fail $ko_path";
$/="///";
while (<KO>) {
	my ($ko,$name,$defi,$class,$ec);
	$ko = $1 if (/ENTRY\s+(K\d+)/s);  $ko =~ s/\s+/ /sg;
	$name = $1 if(/\nNAME\s+(.+?)\n\w/s);  $name =~ s/\s+/ /sg;
	##$defi = $1 if(/\nDEFINITION\s+(.+?)\n\w/s);  $defi =~ s/\s+/ /sg;
	$defi=$1 if (/\nDEFINITION\s+([^\[\]\n]+)\n\w/s);
	($defi,$ec)=($1,$2) if (/\nDEFINITION\s+(.+) \[EC:(\S+)\]\n\w/s);
	$defi =~ s/\s+/ /sg;$ec =~ s/\s+/ /sg;
	$class = $1 if(/\nCLASS\s+(.+?)\n\w/s);  $class =~ s/\s+/ /sg;
	#$ec = $1 if(/\nDBLINKS.+?\s+EC:\s+(\S+)/s);  $ec =~ s/\s+/ /sg;
	#print "$ko\t$name\t$defi\t$class\t$ec\n";
	$KO_Name{$ko} = ($name) ? $name : "--";
	$KO_Defi{$ko} = ($defi) ? $defi : "--";
	$KO_EC{$ko} = ($ec) ? $ec : "--";
	$KO_Class{$ko} = ($class) ? $class : "--";
}
$/="\n";
close KO;

##read the table format blast result file
open TAB,$blast_tab || die "fail $blast_tab";
open OUT,">$Outdir/$blast_tab_base.ko.class" || die "fail $Outdir/$blast_tab_base.ko.class";
while (<TAB>) {
	chomp;
	my ($gene_id,$kegg_gene,$Eval,$kegg_info) = (split /\t/)[0,4,13,15];
	my $ko = $1 if($kegg_info =~ /\s+(K\d+)\s*/);
	print OUT "$gene_id\t$kegg_gene\t$Eval\t$ko\t$KO_Name{$ko}\t$KO_Defi{$ko}\t$KO_EC{$ko}\t$KO_Class{$ko}\n";
	
	if ($KO_EC{$ko} ne '--') { ##根据酶EC号,将基因定到map图上。
		while ($KO_Class{$ko} =~ /\[PATH:ko(\d+)\]/g) {
			push @{$Map_Gene{"map$1"}}, [$gene_id,$ko,$KO_EC{$ko}];
			#print "$1\t$gene_id\t$ko\t$KO_EC{$ko}\n";
		}
	}
}
close OUT;
close TAB;

##generate a record file
open OUT,">$Outdir/$blast_tab_base.map.gene" || die "fail $Outdir/$blast_tab_base.map.gene";
foreach my $map_id (sort keys %Map_Gene) {
	my $map_p = $Map_Gene{$map_id};
	my $gene_num = @$map_p;
	print OUT "$map_id\t$gene_num";
	foreach my $gene_p (@$map_p) {
		my $gene_str = join(",",@$gene_p);
		print OUT "\t$gene_str";
	}
	print OUT "\n";
}
close OUT;

##draw the figures
my $map_dir = "$Outdir/$blast_tab_base.map.fig";
#my $map_dir = "$Outdir/$blast_tab_base.map.png";
mkdir($map_dir);
foreach my $map_id (sort keys %Map_Gene) {
	my $map_p = $Map_Gene{$map_id};
	
	my %conf_ec;
	if ( -e "$map_path/$map_id.html"){
	open (CONF,"$map_path/$map_id.html") || die "fail $map_path/$map_id.html";
	while (<CONF>) {
		chomp;
		s///g;
		#$conf_ec{$5} = [$1,$2,$3,$4] if (/rect \((\d+),(\d+)\) \((\d+),(\d+)\).+\?enzyme\+(\S+)/);
		if (/^\<area shape=rect.+coords\=(\d+),(\d+),(\d+),(\d+).+(\d+\.\d+\.\d+\.[\d\-]+).+/){
			push @{$conf_ec{$5}}, [$1,$2,$3,$4];
			#print "$5\t$1\t$2\t$3\t$4\n";
			#$conf_ec{$5} = [$1,$2,$3,$4];

		}
		#$conf_ec{$5} = [$1,$2,$3,$4] if (/^\<area shape=rect.+coords\=(\d+),(\d+),(\d+),(\d+).+(\d+\.\d+\.\d+\.\d+).+/);
		#print "$5\t$1\t$2\t$3\t$4\n";
	}
	}
	close CONF;
	warn "$map_id not exist",next unless(-f "$map_path/$map_id.png" && -f "$map_path/$map_id.html");
	#open (GIF,"$map_path/$map_id.gif") || die "fail $map_path/$map_id.gif";
	open (PNG,"$map_path/$map_id.png")|| die "fail $map_path/$map_id.png";
	open (RES,">$map_dir/$map_id.png") || die "fail $map_dir/$map_id.png";
	my $im = GD::Image->new(*PNG);
	my $pink = $im->colorAllocate(255,204,255);
	my $black = $im->colorAllocate(0,0,0);
	foreach my $gene_p (@$map_p) {
		my $ec = $gene_p->[2];
		foreach my $temp (@{$conf_ec{$ec}}){
			#$im->filledRectangle($conf_ec{$ec}[0],$conf_ec{$ec}[1],$conf_ec{$ec}[2],$conf_ec{$ec}[3],$pink); ##画填充的矩形
			#$im->string(gdSmallFont,$conf_ec{$ec}[0],$conf_ec{$ec}[1],$ec,$black); 
			$im->filledRectangle($temp->[0],$temp->[1],$temp->[2],$temp->[3],$pink);
			$im->string(gdSmallFont,$temp->[0],$temp->[1],$ec,$black);
		}
	}
	binmode RES;
	print RES $im->png;
	close GIF;
	close RES;
}



##parse the software.config file, and check the existence of each software
####################################################
sub parse_config{
	my $conifg_file = shift;
	my $config_p = shift;
	
	my $error_status = 0;
	open IN,$conifg_file || die "fail open: $conifg_file";
	while (<IN>) {
		next if(/#/);
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













