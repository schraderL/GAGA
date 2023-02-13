#!/usr/bin/perl

=head1 Name

findOverlap.pl  --  find overlap relations between two block sets on chromosome

=head1 Description

This program is designed to find overlap relations of blocks, the input file
can be in psl gff or table format. For gene, it considers the whole gene locus
instead of cds or exon. The strand on chromosome is also considered.

There are 3 types of distance:
(1) real distance, is defined as: 1 - overlap_size * 2 / (ref_size + pre_size);
(2) length distance, is 0 if the overlap size is equal or bigger than the cutoff, is 1 for otherwise. 
(3) percent distance, is 0 if cutoff% of one gene is overlapped with another gene, is 1 for otherwise.
    If percent is set to 100, this means that the two genes are same, or one is included inside another.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 2.0,  Date: 2008-6-16

=head1 Usage
  --real          use the real distance, range from 0 to 1
  --length <num>  use the 0/1 distance decided by overlapped size (bp)
  --percent <num> use the 0/1 distance decided by overlapped percent (both)
  --format <str>     set the format of query file
  --verbose   output verbose information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl findOverlap_CalcuDist.pl  mRNA.gff >  mRNA.dist

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

my ($Real,$Length,$Percent);
my ($QueryFormat,$SubjectFormat);
my ($Verbose,$Help);
GetOptions(
	"real"=>\$Real,
	"length:i"=>\$Length,
	"percent:i"=>\$Percent,
	"format:s"=>\$QueryFormat,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV < 1 || $Help);

print STDERR "length is $Length bp\n" if(defined $Length && defined $Verbose);
print STDERR "percent is $Percent %\n" if(defined $Percent && defined $Verbose);

my $ref_file=shift;
my $pre_file=$ref_file;
$SubjectFormat = $SubjectFormat;
my ( %Ref, %Pre );

if($QueryFormat =~ /psl/ || $ref_file=~/\.psl$/){
	read_psl($ref_file,\%Ref);
}elsif($QueryFormat =~ /gff/ || $ref_file=~/\.gff$/){
	read_gff($ref_file,\%Ref);
}else{
	read_table($ref_file,\%Ref);
}


print STDERR "read ref_file done\n" if($Verbose);


if($SubjectFormat =~ /psl/ || $pre_file=~/\.psl$/){
	read_psl($pre_file,\%Pre);
}elsif($SubjectFormat =~ /gff/ || $pre_file=~/\.gff$/){
	read_gff($pre_file,\%Pre);
}else{
	read_table($pre_file,\%Pre);
}

print STDERR "read pre_file done\n" if($Verbose);


find_overlap(\%Ref,\%Pre);

print STDERR "find overlap done\n" if($Verbose);


####################################################
################### Sub Routines ###################
####################################################

##利用两组元素之间找overlap的机制，找到overlap,并计算用于聚类的距离
####################################################
sub find_overlap{
	my $Ref_hp=shift;
	my $Pre_hp=shift;
	
	foreach  my $chr (sort keys %$Ref_hp) {
		
		my $output;
		my @ref_chr = (exists $Ref_hp->{$chr})  ? (sort {$a->[0] <=> $b->[0]} @{$Ref_hp->{$chr}}) : ();
		my @pre_chr = (exists $Pre_hp->{$chr})  ? (sort {$a->[0] <=> $b->[0]} @{$Pre_hp->{$chr}}) : ();
		
		print STDERR "find overlap on $chr\n" if($Verbose);
		
		my $pre_pos = 0;
		for (my $i=0; $i<@ref_chr; $i++) {
			my $ref_gene = $ref_chr[$i][2];
			my $ref_size = $ref_chr[$i][1] - $ref_chr[$i][0] + 1;
			my $ref_strand = $ref_chr[$i][3];
			my @overlap;
			
			for (my $j=$pre_pos; $j<@pre_chr; $j++) {
				if ($pre_chr[$j][1] < $ref_chr[$i][0]) {
					$pre_pos++;
					next;
				}
				if ($pre_chr[$j][0] > $ref_chr[$i][1]) {
					last;
				}
				my $pre_strand = $pre_chr[$j][3];
				my $pre_size = $pre_chr[$j][1] - $pre_chr[$j][0] + 1;
				my $overlap_size = overlap_size($pre_chr[$j],$ref_chr[$i]);
				my $distance_rate;
				if(defined $Length){
					$distance_rate = ($overlap_size >= $Length) ? 0 : 1;
				}elsif(defined $Percent){
					$distance_rate = ($overlap_size/$ref_size*100 >= $Percent || $overlap_size/$pre_size*100 >= $Percent) ? 0 : 1;
	
				}else{
					$distance_rate = 1 - $overlap_size * 2 / ($ref_size + $pre_size);
				}
				
				##push @overlap,"$pre_chr[$j][2],$pre_size,$overlap_size";
				next if($ref_gene eq $pre_chr[$j][2]);
				$output .= "$ref_gene\t$pre_chr[$j][2]\t$distance_rate\t$overlap_size\t$chr";
				$output .= "\t$ref_strand\t$ref_size\t$ref_chr[$i][0]\t$ref_chr[$i][1]";
				$output .= "\t$pre_strand\t$pre_size\t$pre_chr[$j][0]\t$pre_chr[$j][1]\n";
			}
			
		}

		print $output;
	}

}


sub overlap_size {
	my $block1_p = shift;
	my $block2_p = shift;
	
	my $combine_start = ($block1_p->[0] < $block2_p->[0]) ?  $block1_p->[0] : $block2_p->[0];
	my $combine_end   = ($block1_p->[1] > $block2_p->[1]) ?  $block1_p->[1] : $block2_p->[1];
	
	my $overlap_size = ($block1_p->[1]-$block1_p->[0]+1) + ($block2_p->[1]-$block2_p->[0]+1) - ($combine_end-$combine_start+1);

	return $overlap_size;
}



##read gff3 format
sub read_gff{
	my $file=shift;
	my $ref=shift;

	open (IN,$file) || die ("fail open $file\n");
	while (<IN>) {
		chomp;
		s/^\s+//;
		my @t = split(/\t/);
		my $tname = $t[0];

		if ($t[2] eq 'mRNA' || $t[2] eq 'match') {
			my $qname = $1 if($t[8] =~ /ID=([^;]+);*/);
			my $start = $t[3];
			my $end = $t[4];
			my $strand = $t[6];
			push @{$ref->{$tname}},[$start,$end,$qname,$strand];
		}
	}
	close(IN);
	
}

sub read_psl{
	my $file=shift;
	my $ref=shift;
	open(REF,$file)||die("fail to open $file\n");
	while (<REF>) {
		chomp;
		my @temp=split(/\s+/,$_);
		my $chr=$temp[13];
		my $gene=$temp[9];
		my $start=$temp[15]+1;
		my $end=$temp[16];
		my $strand = $temp[8];

		push @{$ref->{$chr}},[$start,$end,$gene,$strand];
	}
	close(REF);
}

sub read_table{
	my $file=shift;
	my $ref=shift;
	open(REF,$file)||die("fail to open $file\n");
	while (<REF>) {
		chomp;
		my @temp=split(/\t/,$_);
		my ($gene,$chr,$strand,$start,$end) = @temp;
		push @{$ref->{$chr}},[$start,$end,$gene,$strand];
	}
	close(REF);
}

