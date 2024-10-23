#!/usr/bin/perl
#Author: liaoxinhui@genomics.org.cn
#Date:  Sun Jul  5 16:39:28 CST 2015


use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use Cwd 'abs_path';
use File::Basename;

my ($list, $nr, $nt, $swissprot, $cog, $kegg, $go, $obo, $interpro, $trembl, $outxls, $outstat, $help);
GetOptions (
	"list:s" => \$list,
	"nr:s" => \$nr,
	"nt:s" => \$nt,
	"go:s" => \$go,
	"obo:s" => \$obo,
	"cog:s" => \$cog,
	"kegg:s" => \$kegg,
	"swissprot:s" => \$swissprot,
	"interpro:s" => \$interpro,
	"trembl:s" => \$trembl,
	"outxls:s" => \$outxls,
	"outstat:s" => \$outstat,
	"help|?" => \$help
);

if ( !$list || (!$nr && !$nt && !$go && !$cog && !$kegg && !$swissprot && !$interpro && !$trembl) || $help || ($go && !$obo)) {
	die <<USAGE;
=========================================================================================
Description: Merge all annotation files
Usage: perl $0 [options]
Options:
  * -list           input all_gene_id.list, no header, format: <geneID  somethingOrnot>
	-nr				input nr annotation file
	-nt				input nt annotation file
	-go				input go annotation file
	-obo			input obo database file, set with -go
	-cog			input cog annotation file
	-kegg			input kegg annotation file
	-interpro		input interpro annotation file
	-swissprot		input swissprot annotation file
	-trembl			input trembl annotation file
	-outxls			output xls, default: [./all_annotation.xls]
	-outstat        output stat, default: [./annotation_stat.xls]
	-help			print this help information
E.g.:
	perl $0 -list all_gene_id.list -nr nr.xls -nt nt.xls -kegg kegg.xls
=========================================================================================
USAGE
}

$outxls ||= "./all_annotation.xls";
$outstat ||= "./annotation_stat.xls";

my %info = ();
my @db = ();
if ($nr) {
	phrase_m8(\%info,$nr,"Nr");
	push @db, "Nr";
}
if ($nt) {
	phrase_m8(\%info,$nt,"Nt");
	push @db, "Nt";
}
if ($swissprot) {
	phrase_m8(\%info,$swissprot,"Swissprot");
	push @db,"Swissprot";
}
if ($kegg) {
	phrase_m8(\%info,$kegg,"KEGG");
	push @db,"KEGG";
}
if ($cog) {
	phrase_m8(\%info,$cog,"COG");
	push @db,"COG";
}
if ($trembl) {
	phrase_m8(\%info,$trembl,"TrEMBL");
	push @db,"TrEMBL";
}
if ($interpro) {
	phrase_ipr(\%info,$interpro,"Interpro");
	push @db,"Interpro";
}
if ($go) {
	phrase_go(\%info,$go,$obo,"GO");
	push @db,"GO";
}

my %stat = ();

open OXLS,">$outxls" or die $!;
open STAT,">$outstat" or die $!;
print OXLS "Unigene\t" . join("\t", @db) . "\n";
print STAT "Values\tTotal\t" . join("\t", map{"$_-Annotated"} @db) . "\tOverall\n";

my $total_unigene = 0;

open LIST,$list or die $!;
while (<LIST>) {
	next if (/^#/ || /^\s*$/);
	chomp; my @a = split /\s+/;
	$total_unigene++;
	my $out = $a[0];
	my $flag = 0;
	foreach my $db (@db) {
		if ($info{$a[0]}{$db}) {
			$out .= "\t$info{$a[0]}{$db}";
			$stat{$db}++;
			$flag = 1;
		}else {
			$out .= "\tNA";
		}
	}
	$stat{'Overall'}++ if ($flag == 1);
	print OXLS "$out\n";
}
close LIST;
close OXLS;

push @db,"Overall";
my %ptg = ();
foreach my $k (@db) {
	$ptg{$k} = $stat{$k} ? sprintf("%.2f%%",$stat{$k}/$total_unigene*100) : "0.00%";
}

while ($total_unigene =~  s/(\d+)(\d{3}(,\d{3})*)/$1,$2/) {};
my $stat_num = "Number\t$total_unigene";
my $stat_ptg = "Percentage\t100%";
foreach my $db (@db) {
	while($stat{$db} =~ s/(\d+)(\d{3}(,\d{3})*)/$1,$2/) {};
	$stat_num .= ($stat{$db}) ? "\t$stat{$db}" : "\t0";
	$stat_ptg .= ($ptg{$db}) ? "\t$ptg{$db}" : "\t0.00%";
}
print STAT "$stat_num\n";
print STAT "$stat_ptg\n";
close STAT;

sub phrase_m8 {
	my ($hash, $file, $key) = @_;
	open M8, $file or die $!;<M8>;
	while (<M8>) {
		next if (/^#/ || /^\s*$/);
		chomp; my @a = split /\t+/,$_,13;
		$$hash{$a[0]}{$key} ||= "$a[1]/$a[-3]/$a[-1]"; # subject_id / evalue / description
	}
	close M8;
}

sub phrase_ipr {
	my ($hash, $file, $key) = @_;
	open IPR, $file or die $!;<IPR>;
	while (<IPR>) {
		next if (/^#/ || /^\s*$/);
		chomp; my @a = split /\t+/;
		$$hash{$a[0]}{$key} ||= "$a[1]/$a[-2]/$a[-1]";
	}
	close IPR;
}

sub phrase_go {
#	use PerlIO::gzip;
	my ($hash, $file, $db, $key) = @_;
	my %annot = ();
	my $mod = ($db =~ /\.gz$/) ? "<:gzip " : "<";
	open DB,$mod,$db or die $!;
	my ($id,$name,$namespace);
	while (<DB>) {
		chomp;
		if (/^id:\s+(\S+)/) {
	       $id = $1;
   		}elsif (/^name:\s+(.*)/) {
     		$name = $1;
    	}elsif (/^namespace:\s+(\S+)/) { # molecular_function, biological_process, cellular_component
        	$namespace = $1;
        	$annot{$id} = [$name, $namespace];
		}
	}
	close DB;
	open GO,$file or die $!;
	while (<GO>) {
		chomp; my @a = split /\s+/;
		my $gene = shift @a;
		my %tmp = ();
		foreach my $go (@a) {
			next unless ($annot{$go}); # oh fuck! juck ignore it!
			push @{$tmp{$annot{$go}[1]}}, "$go//$annot{$go}[0]";
		}
		foreach my $k (sort keys %tmp) {
			$$hash{$gene}{$key} .= "$k:".join(";",@{$tmp{$k}}).";";
		}
	}
	close GO;
}
