#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# Quick way: grep ' p-value =' absrel/*output | less -S


# usage: perl analyze_absrel.pl

# INPUT
my $cdsdir = "/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Cleaned_alignments/run_all_withgard/hmmcleaner_50_aln/";
my $absreldir = "/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Hyphy_absrel/absrel_withgardpartitions_hmmclean50aln/absrel_allhog/";

my $outdir = "absrel_output_tables";
system ("mkdir -p $outdir");

my $pvalue = "0.001"; # p-value to filter significant hits 
my $fdr = "0.01"; # FDR to filter significant hits

my ($line, $name);

# Reading each alignment file for which absrel was run

#my %protfasta;
system ("ls $cdsdir\/\*\.nonstop.aln > tmp_protlist.txt");
open(File, "<", "tmp_protlist.txt");
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	

	my $id = "";
	if ($line =~ /.*\/(\S+)\.cds/){
		$id = $1;
	} else {
		die "Can't find OG id in $line\n";
	}  

	my $absrelout = "$absreldir\/$id\_ABSREL.json";
	my $tableout = "$outdir\/$id\_absrel_wdnds.tsv";
	my $tablesignout = "$outdir\/$id\_absrel_significant\_pval$pvalue\.tsv";
	my $tablesignfdrout = "$outdir\/$id\_absrel_significant\_fdr$fdr\.tsv";


	system ("python3 /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/significant_aBSREL_2_table.py $absrelout $tablesignout 0 $pvalue \n");
	system ("python3 /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/significant_aBSREL_2_table.py $absrelout $tablesignfdrout $fdr 0\n");
	system ("python3 /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/aBSREL_2_table_withdnds.py $absrelout $tableout \n");


}
close File;
system ("rm tmp_protlist.txt");





