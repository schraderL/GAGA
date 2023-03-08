#!/usr/bin/perl

=head1 Name

iprscan_parser51-55.pl  --  parse and convert the iprscan result

=head1 Description

This program is designed to be a general iprscan raw result parser.
It will do both IPR and GO analysis.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Mender: He lijuan, Lyndi.He@genomics.cn
  Version: 2.0,  Date: 2008-5-22
  Version: 3.0,  Data: 2015-5-21

=head1 Usage
 
  --outdir <str>  set the result directory, default="."
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl ../bin/iprscan_parser51-55.pl rice_prot100.fa.iprscan

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;

my ($Verbose,$Help,$Outdir);
GetOptions(
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Outdir ||= ".";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $iprscan_file = shift;
my $iprscan_file_base = basename($iprscan_file);

my %config;
parse_config("$Bin/../../config.txt",\%config);
my $ipr2go=$config{"ipr2go"};
my %go;
open IPR2GO,"<$ipr2go"||die;
#print "$ipr2go\n";
while(<IPR2GO>){
   chomp;
   my @tmp=split /\|/,$_;
   my $go_id=$tmp[0]."_".$tmp[1];
   $go{$go_id}{desc}=$tmp[2];
   $go{$go_id}{class}=$tmp[3];
}


my %Gene_IPR;
my %Gene_GO;
my %IPR_Gene;
my %IPR;
my %GO;

$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);


##parse the iprscan raw result file
open IN,$iprscan_file || die "fail open $iprscan_file";
while (<IN>) {
	chomp;
	my @t = split /\t/;
	my $gene = $t[0];
	my $ipr_id = $t[11];
	my $ipr_desc = $t[12];
	my $eval = $t[8];
	my $method = $t[3];

	if($ipr_id =~ /IPR\d+/){
		$Gene_IPR{$gene}{$ipr_id} = 1;
		$IPR_Gene{$ipr_id}{$gene} = 1;
		$IPR{$ipr_id} = $ipr_desc;
	}

	my $go_str = $t[13];
#    print "1: $go_str\n";
    my @go_str=split /\|/,$go_str;
#    print "@go_str\n";
	foreach my $tmp(@go_str) {
        my $go_id=$ipr_id."_".$tmp;
        $Gene_GO{$gene}{$tmp} = 1;
        if(exists $go{$go_id}){
            $GO{$tmp}{desc} = $go{$go_id}{desc};
            $GO{$tmp}{class} = $go{$go_id}{class};
#            print "$tmp\t$GO{$tmp}{desc}\t$GO{$tmp}{class}\n";
        }
    }

}
close IN;


##output gene.IPR result
open OUT, ">$Outdir/$iprscan_file_base.gene.ipr" || die $!;
foreach my $gene (sort keys %Gene_IPR) {
	my $IPR_p = $Gene_IPR{$gene};
	my $IPR_num = keys %$IPR_p;
	print OUT "$gene\t$IPR_num";
	foreach my $ipr_id (sort keys %$IPR_p) {
		print OUT "\t$ipr_id; $IPR{$ipr_id}";
	}
	print OUT "\n";
}
close OUT;

##output IPR.gene result
open OUT, ">$Outdir/$iprscan_file_base.ipr.gene" || die $!;
foreach my $ipr_id (sort keys %IPR_Gene) {
	my $gene_p = $IPR_Gene{$ipr_id};
	my $gene_num = keys %$gene_p;
	print OUT "$ipr_id\t$IPR{$ipr_id}\t$gene_num";
	foreach my $gene (sort keys %$gene_p) {
		print OUT "\t$gene";
	}
	print OUT "\n";
}
close OUT;

##output gene.GO result
open OUT, ">$Outdir/$iprscan_file_base.gene.GO" || die $!;
open WEGO, ">$Outdir/$iprscan_file_base.gene.wego" || die $!;
foreach my $gene (sort keys %Gene_GO) {
	my $GO_p = $Gene_GO{$gene};
	my $GO_num = keys %$GO_p;
	print OUT "$gene\t$GO_num";
	print WEGO "$gene";
	foreach my $go_id (sort keys %$GO_p) {
		print OUT "\t$go_id; $GO{$go_id}{desc}; $GO{$go_id}{class}";
		print WEGO "\t$go_id";
	}
	print OUT "\n";
	print WEGO "\n";
}
close OUT;
close WEGO;


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

