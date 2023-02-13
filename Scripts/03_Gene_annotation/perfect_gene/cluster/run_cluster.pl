#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($Bin);

# Run the culster program

die "Usage:$0 <all.gff > <all.gff.cds.len> [cutoff, default=1; cpu,default:10] \n\n" if @ARGV<2;

my $gff=shift;
my $cds_len=shift;
my $Cutoff=shift;
my $Cpu=shift;

$Cutoff ||=1;
$Cpu ||=10;

my $name=$gff;
$name=~s/.+\/([^\/]+)$/$1/;

my %Gene;
read_gff($gff,\%Gene);

my %Len;
read_len($cds_len,\%Len);

my @New_gff;
my @New_len;

my $num=0;
foreach my $chr (sort { scalar(keys %{$Gene{$b}}) <=> scalar(keys %{$Gene{$a}}) } keys %Gene ){
	my $idx = $num % $Cpu;
	$num++;
	foreach my $id (sort keys %{$Gene{$chr}} ){
		$New_gff[$idx].=join("\t",@{$Gene{$chr}{$id}{mRNA}})."\n";
		foreach my $cds ( @{$Gene{$chr}{$id}{CDS}} ){
			$New_gff[$idx].=join("\t",@$cds)."\n";
		}
		$New_len[$idx].=join("\t",$id,$Len{$id})."\n";
	}
}

my $Dir="./cluster".time();
`rm -r $Dir ` unless ( ! -e $Dir );
mkdir ($Dir) unless ( -d $Dir );

open SH,">cluster.all.sh" or die "$!";
for(my $i=0;$i<@New_gff;$i++){
	my $sub_gff="$Dir/$i.gff";
	my $sub_len="$Dir/$i.cds.len";
	open OUT,">$sub_gff" or die "$!";
	print OUT $New_gff[$i];
	close OUT;
	open OUT,">$sub_len" or die "$!";
	print OUT $New_len[$i];
	close OUT;
	print SH "perl $Bin/cluster_genes_on_genomic_locus.pl -cutoff $Cutoff -outdir $Dir $sub_gff $sub_len \n";
}
close SH;

#`qsub-sge.pl -resource vf=0.1G -maxjob $Cpu cluster.all.sh`;
`multi-process.pl -cpu $Cpu cluster.all.sh`;

`cat $Dir/*.dist.hcluster.reference > $name.dist.hcluster.reference`;
`cat $Dir/*.uncluster.gff >$name.uncluster.gff`;
`cat $Dir/*.nr.gff >$name.nr.gff`;

sub read_gff{
	my ($file,$p)=@_;
	open IN,$file or die "$!";
	while(<IN>){
		next if /^#/;
		chomp;
		my @c=split(/\t/);
		if ($c[2] eq 'mRNA' && $c[8]=~/ID=([^;]+)/){
			@{$p->{$c[0]}{$1}{mRNA}}=@c;
		}elsif (($c[2] eq 'CDS' || $c[2] eq 'exon') && $c[8]=~/Parent=([^;]+)/){
			push @{$p->{$c[0]}{$1}{CDS}},[@c];
		}
	}
	close IN;
}

sub read_len{
	my ($file,$p)=@_;
	open IN,$file or die "$!";
	while(<IN>){
		if (/^(\S+)\s+(\d+)/){
			$p->{$1}=$2;
		}
	}
	close IN;
}
