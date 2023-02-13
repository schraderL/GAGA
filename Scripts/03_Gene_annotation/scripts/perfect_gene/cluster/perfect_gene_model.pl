#!/usr/bin/perl

=head1 Name

perfect_gene_model.pl  -- the pipeline to get perfect gene models

=head1 Description

This pipeline takes the gene gff file and genome fasta file as input. It will
check the cds models(start and stop codons, frameshift, and premature stop codon ), 
check splicing site(only GT/AG is allowed for eukaryote now),
and filter the redundance on genomic level. At last generate a perfect representative
gene set, which do not have any above drawback.

It will also remove the small cds genes and less exon genes, based on the min_cds_len
and min_exon_num.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-4-19

=head1 Usage
  perl perfect_gene_model.pl <gene.gff> <genome.fa>
  --min_cds_len <num>   set the minimum cds length, default 150bp
  --min_exon_num <num>  set the minimum exon number, default 1
  --outdir <str>   set the result directory, default="./" 
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

 

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

my ($Verbose,$Help,$Outdir,$Min_cds_len,$Min_exon_num);
GetOptions(
	"min_cds_len:i"=>\$Min_cds_len,
	"min_exon_num:i"=>\$Min_exon_num,
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Min_cds_len ||= 150;
$Min_exon_num ||= 1;
$Outdir ||= "./";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $gff_file = shift;
my $genome_file = shift;

my $gff_file_basename = basename($gff_file);

$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

my %config;
parse_config("$Bin/config.txt",\%config);

`perl $config{"getGene.pl"} $gff_file $genome_file > $Outdir/$gff_file_basename.cds`;
`perl $config{"cds2aa.pl"} -check $Outdir/$gff_file_basename.cds > $Outdir/$gff_file_basename.cds.wrong`;
`perl $config{"fishInWinter.pl"} -except -ff gff $Outdir/$gff_file_basename.cds.wrong $gff_file > $Outdir/$gff_file_basename.cds.correct.gff`;

`perl $config{"getGene.pl"} -type splice $gff_file $genome_file > $Outdir/$gff_file_basename.splice`;
`more $Outdir/$gff_file_basename.splice | awk '\$3!="GT/AG"' > $Outdir/$gff_file_basename.splice.wrong`;
`perl $config{"fishInWinter.pl"}  -except -ff gff $Outdir/$gff_file_basename.splice.wrong $Outdir/$gff_file_basename.cds.correct.gff > $Outdir/$gff_file_basename.cds.correct.splice.correct.gff`;

`perl $Bin/filter_on_genomic_level.pl $Outdir/$gff_file_basename.cds.correct.splice.correct.gff -outdir $Outdir`;
`mv $Outdir/$gff_file_basename.cds.correct.splice.correct.gff.nonredundant $Outdir/$gff_file_basename.cds.correct.splice.correct.nr.gff`;

`perl $Bin/filter_gene_gff.pl --min_cds_len $Min_cds_len --min_exon_num $Min_exon_num $Outdir/$gff_file_basename.cds.correct.splice.correct.nr.gff > $Outdir/$gff_file_basename.cds.correct.splice.correct.nr.final.gff`;



####################################################
################### Sub Routines ###################
####################################################



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

