#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# usage: perl get_gemoma_gff.pl $chem\/gemoma_outdir/final_annotation.gff $chem\/gemoma_outdir/comparison.tabular $chem\/gemoma_outdir/ $chem $genome

my ($line, $name, $nameout);

my $ingffrep = $ARGV[0];
my $ingffok = $ARGV[1];
my $intable = $ARGV[2];
my $outgff = $ARGV[3];


my $excludegenes = "";

open (File , "<", $intable); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	next if ($line =~ /^#/);
	my @subl = split (/\t/, $line);

	my $gene = uc($subl[3]); ### uc added to avoid gemoma upper case or not in different versions
	unless ($subl[11] =~ /^NA$/){
		$excludegenes .= " $gene ";
	}
}

open (Resultsgff, ">", "$outgff");
print Resultsgff "##gff-version 3\n";


open (GFFfile , "<", $ingffrep); 
while (<GFFfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);

	if ($line =~ /^#/){
#		print Resultsgff "$line\n";
		next;
	}
	
	my @subline = split (/\t/, $line);

	if ($subline[2] =~ /CDS/){
		my $genename = "";
		my $rest = "";
		if ($subline[8] =~ /Parent=([^;]+)(\S+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
			$rest = $2;
		}
		else {die "ERROR in run_OR_classification.pl: It fails detecting Parent ID in $line\n";}

		my $nnamef = "Cb"."$genename";

		my $ucgenename = uc($genename);

		unless ($excludegenes =~ / $ucgenename /){
			print Resultsgff "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tParent=$nnamef"."$rest\;\n";

		}
		
	}
	elsif ($subline[2] =~ /mRNA/){
		my $genename = "";
		my $rest = "";
		if ($subline[8] =~ /ID=([^;]+)\;Parent=[^;]+(\;\S+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
			$rest = $2;
		}
		else {print "ERROR in run_OR_classification.pl: It fails detecting ID in $line\n";}

		my $nnamef = "Cb"."$genename";

		my $ucgenename = uc($genename);

		unless ($excludegenes =~ / $ucgenename /){
			print Resultsgff "$subline[0]\t$subline[1]\tgene\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=g"."$nnamef"."$rest;\n";		
			print Resultsgff "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=$nnamef\;Parent=g"."$nnamef"."$rest;\n";		

		}

	}
}
close GFFfile;

open (GFFfile , "<", $ingffok); 
while (<GFFfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);

	if ($line =~ /^#/){
#		print Resultsgff "$line\n";
		next;
	}
	
	my @subline = split (/\t/, $line);

	if ($subline[2] =~ /CDS/){
		my $genename = "";
		my $rest = "";
		if ($subline[8] =~ /Parent=([^;]+)(\S+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
			$rest = $2;
		}
		else {die "ERROR in run_OR_classification.pl: It fails detecting Parent ID in $line\n";}

		my $nnamef = "Cb"."$genename";

		my $ucgenename = uc($genename);

#		unless ($excludegenes =~ / $ucgenename /){ # No filtering here, we want to keep all genes from this GFF
			print Resultsgff "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tParent=$nnamef"."$rest\;\n";

#		}
		
	}
	elsif ($subline[2] =~ /mRNA/){
		my $genename = "";
		my $rest = "";
		if ($subline[8] =~ /ID=([^;]+)\;Parent=[^;]+(\;\S+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
			$rest = $2;
		}
		else {print "ERROR in run_OR_classification.pl: It fails detecting ID in $line\n";}

		my $nnamef = "Cb"."$genename";

		my $ucgenename = uc($genename);

#		unless ($excludegenes =~ / $ucgenename /){
			print Resultsgff "$subline[0]\t$subline[1]\tgene\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=g"."$nnamef"."$rest;\n";		
			print Resultsgff "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=$nnamef\;Parent=g"."$nnamef"."$rest;\n";		

#		}

	}
}
close GFFfile;


close Resultsgff;



