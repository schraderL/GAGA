#!/usr/bin/perl
use strict;
use warnings;

#########################################################
#####
#
# Filter bacterial scaffolds from the assembly
#
#####

# Usage
# perl filter_bacterial_scaffolds.pl folder_with_final_assemblies_infasta Directory_with_output_from_bacterial_pipeline

## Input variables

my $genomedir = $ARGV[0];
my $contdir= $ARGV[1];
#my $contdir = "/home/projects/ku_00039/people/dinghe/working_dr/metagenome_lgt/GAGA";

## Start

open (Outtable, ">", "GAGA_bacterial_andhuman_scaffolds_stats_20210420.tsv");
print Outtable "GAGA-ID\tBacterial scaffolds\tTotal length of bacterial scaffolds\tHuman scaffolds\tTotal length of human scaffolds\t List of discarded scaffolds\tBacterial taxa\n";

system("ls $genomedir\/*fasta > Input_genomes.txt");
open(File, "<", "Input_genomes.txt");
while(<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	my $id = "";
    my $fullid = "";
	if ($line =~ /\/+(\S{3,4}\-\d{4})\_(\S+).fasta/){
		$id = $1;
        $fullid = "$1\_$2";
	} else {
		die "Cannot find GAGA-ID in $line\n";
	}

    my $contfile = "$contdir\/$id\/results/contaminantscaffolds.csv";
    
    if (-f $contfile) {
        #Analysis done
    } else {
        print "No contaminants file found for $id\n";
        next;
    }

    my $contnum = 0;
    my $contlen = 0;
    my $contlist = "";
    my $conttaxa = "";
    open(Filelog, "<", $contfile);
    while(<Filelog>){
        chomp;
        my $nline = $_;
        next if ($nline !~ /\S+/);
        next if ($nline =~ /LongestContProWindowSize/);
        my @subl = split (/\t/, $nline);
        $contnum++;
        $contlen+=$subl[17];
        $contlist .= "$subl[0] ";
        $conttaxa .= "$subl[1] "
    }
    close Filelog;


    my $humancontfile = "$contdir\/$id\/results/contaminants.human.csv";
    
    if (-f $humancontfile) {
        #Analysis done
    } else {
        print "No human contaminants file found for $id\n";
        next;
    }

    my $humancontnum = 0;
    my $humancontlen = 0;
    open(Filelog, "<", $humancontfile);
    while(<Filelog>){
        chomp;
        my $nline = $_;
        next if ($nline !~ /\S+/);
        next if ($nline =~ /LongestContProWindowSize/);
        my @subl = split (/\t/, $nline);
        $humancontnum++;
        $humancontlen+=$subl[1];
        $contlist .= "$subl[0] ";
    }
    close Filelog;

    
    print Outtable "$id\t$contnum\t$contlen\t$humancontnum\t$humancontlen\t$contlist\t$conttaxa\n";
    
    
    open (Outfasta, ">", "$fullid\_filt.fasta");
    
    my %fasta;
    my $name = "";
    open(Fasta, "<", $line);
    while(<Fasta>){
        chomp;
        my $nline = $_;
        next if ($nline !~ /\S+/);
        if ($nline =~ />(\S+)/){
            $name = $1;
        } else {
            $fasta{$name} .= uc($nline); # Change all sequences to upper case nucleotides
        }
    }
    close Filelog;
    
#    foreach my $key (sort {$a <=> $b} (keys %fasta)){
    foreach my $key (sort (keys %fasta)){
        if ($contlist =~ /$key /){
            next;
        } else {
            print Outfasta ">$key\n$fasta{$key}\n";
        }
    
    }
    
	
}
close File;

close Outtable;





