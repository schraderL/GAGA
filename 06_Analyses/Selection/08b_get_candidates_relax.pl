#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

## Changed line 32 with file extension
# usage: perl get_candidates_relax.pl Relax_output_folder(output folder from hyphy RELAX containing the JSON files) Relax_results(folder to store the results from this script)

# module load ngs tools gcc/7.4.0 intel/perflibs R/4.2.0

my $relaxdir = "$ARGV[0]";
my $outprefix = $ARGV[1];


## Read dir with all json files - write them in a table

my $outrelaxtable = "$outprefix\_summary_table.txt";
open (Results, ">", "$outrelaxtable");
print Results "HOG\tpvalue\trelaxation or intensification parameter\n";

my $ogrelaxfail = "0";
system ("ls $relaxdir\/*RELAX.json > tmp_relaxfiles.txt");
open(Filef, "<", "tmp_relaxfiles.txt");
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );

    my $og = "";
#    if ($line =~ /(HOG\d\d\d\d\d\d\d)/){
    if ($line =~ /(HOG\S+)\_RELAX\.json/){   # Here do it for each partition
        $og = $1;
    } else {
        die "Cannot find OG in $line\n";
    }

    my $ogok = 0;
    system ("grep -A 1 \"p-value\" $line".' | sed \'s/\s//g\' | sed \'s/\"p-value\"//g\' | sed \'s/\"relaxationorintensificationparameter\"//g\' | sed \'s/\://g\' | sed \'s/\,//g\' | sed \'N;s/\n/\t/\''." > tmp_relaxvals.txt");
    open(Filetmp, "<", "tmp_relaxvals.txt");
    while(<Filetmp>){
        chomp;
        my $line2 = $_;
        if ($line2 !~ /\S+/ ){
            $ogrelaxfail++;
            print "$og does not contain any character, filling with 1\n";
            print Results "$og\t1\t1\n";
#            print "$og does not contain any character, skipping\n";
            next;
        } else {
            print Results "$og\t$line2\n";
            $ogok++;
        }
    }
    close Filetmp;
    system ("rm tmp_relaxvals.txt");

    if ($ogok == 0){
        $ogrelaxfail++;
        print "$og does not contain any character from relax, filling with 1\n";
        print Results "$og\t1\t1\n";
#       print "$og does not contain any character, skipping\n";
    }

}
close Filef;
close Results;
system ("rm tmp_relaxfiles.txt");

print "\nHOGs with failed RELAX: $ogrelaxfail\n\n";


## Then, merge partitions into single HOGs, check if relax/intensification occurs, and check if congruent, if not just reduce the p-val to 1. 
open (Results, ">", "$outprefix\_summary_table_mergeparts.txt");
print Results "HOG\tpvalue\trelaxation or intensification parameter\tRELAX\n";

my %hogtable;
my $hogincongruent = 0;
open(Filef, "<", "$outprefix\_summary_table.txt");
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    next if ($line =~ /relaxation or intensification parameter/ );    
    my $og = "";
    if ($line =~ /(HOG\d\d\d\d\d\d\d)/){
        $og = $1;
    }

    my @cursubl = split (/\t/, $line);

    if (exists $hogtable{$og}){

        my @prevsubl = split (/\t/, $hogtable{$og});

        if ($hogtable{$og} =~ /Incongruent/){
            next; # Already incongruent with some partitions
        } elsif ($cursubl[1] < $prevsubl[1]){
            if ($cursubl[1] <= 0.05 && $prevsubl[1] <= 0.05){
                if ($cursubl[2] > 1 && $prevsubl[2] > 1){ # Ok Intensifition in test
                    $hogtable{$og} = "$og\t$cursubl[1]\t$cursubl[2]";
                } elsif ($cursubl[2] < 1 && $prevsubl[2] < 1){ # Ok Relaxation in test
                    $hogtable{$og} = "$og\t$cursubl[1]\t$cursubl[2]";
                } else {
                    $hogincongruent++;
                    $hogtable{$og} = "$og\t1\tIncongruent";
                }
            }

        } elsif ($cursubl[1] >= $prevsubl[1]){
            if ($cursubl[1] <= 0.05 && $prevsubl[1] <= 0.05){
                if ($cursubl[2] > 1 && $prevsubl[2] > 1){ # Ok Intensifition in test
                    #$hogtable{$og} = $line;
                    next; # Previous better p-value
                } elsif ($cursubl[2] < 1 && $prevsubl[2] < 1){ # Ok Relaxation in test
                    #$hogtable{$og} = $line;
                    next; # Previous better p-value
                } else {
                    $hogincongruent++;
                    $hogtable{$og} = "$og\t1\tIncongruent";
                }
            }

        }
        
    } else {
        $hogtable{$og} = "$og\t$cursubl[1]\t$cursubl[2]";
    }

}
close Filef;

foreach my $og (sort keys %hogtable){
    my @subl = split (/\t/, $hogtable{$og});
    my $rpar = "";
    if ($subl[2] > 1){
        $rpar = "Intensification";
    } else {
        $rpar = "Relaxation";
    }

    print Results "$hogtable{$og}\t$rpar\n";
}
close Results;

print "\nHOGs with incongruent partitions in RELAX (i.e.: one partition relaxed, other intensified): $hogincongruent\n\n";


## Correct the table and get fdr - R script in the table

open (Resultsr, ">", "$outprefix\_relax_fdr.R");

print Resultsr '
table <- read.csv("'."$outprefix\_summary_table_mergeparts.txt".'", header = T, sep = "\t")

table$FDR <- p.adjust(table$pvalue, "fdr")
table$Bonferroni <- p.adjust(table$pvalue, "bonferroni")

write.table(table, file="'."$outprefix\_summary_table_mergeparts_fdr.txt".'", quote=FALSE, sep="\t", col.names = NA)

';

system("Rscript $outprefix\_relax_fdr.R");


## Get candidate gene list, and reference set filtering pval and fdr

open (Results, ">", "$outprefix\_candidates_relax_pvalue001_table.txt");
print Results "HOG\tpvalue\trelaxation or intensification parameter\tRELAX\tFDR\tBonferroni\n";
open (Resultsr, ">", "$outprefix\_candidates_relax_pvalue001_table_relaxedtest.txt");
print Resultsr "HOG\tpvalue\trelaxation or intensification parameter\tRELAX\tFDR\tBonferroni\n";
open (Resultsi, ">", "$outprefix\_candidates_relax_pvalue001_table_intensifiedtest.txt");
print Resultsi "HOG\tpvalue\trelaxation or intensification parameter\tRELAX\tFDR\tBonferroni\n";

open (Resultsf, ">", "$outprefix\_candidates_relax_fdr001_table.txt");
print Resultsf "HOG\tpvalue\trelaxation or intensification parameter\tRELAX\tFDR\tBonferroni\n";
open (Resultsfr, ">", "$outprefix\_candidates_relax_fdr001_table_relaxedtest.txt");
print Resultsfr "HOG\tpvalue\trelaxation or intensification parameter\tRELAX\tFDR\tBonferroni\n";
open (Resultsfi, ">", "$outprefix\_candidates_relax_fdr001_table_intensifiedtest.txt");
print Resultsfi "HOG\tpvalue\trelaxation or intensification parameter\tRELAX\tFDR\tBonferroni\n";

open(Filef, "<", "$outprefix\_summary_table_mergeparts_fdr.txt");
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    next if ($line =~ /Bonferroni/ );    
    my $og = "";
    if ($line =~ /(HOG\d\d\d\d\d\d\d)/){
        $og = $1;
    }

    my @subl = split(/\t/, $line);

    if ($subl[2] <= 0.01){ # P-value 0.01
        print Results "$subl[1]\t$subl[2]\t$subl[3]\t$subl[4]\t$subl[5]\t$subl[6]\n";
        if ($subl[4] =~ /Relaxation/){
            print Resultsr "$subl[1]\t$subl[2]\t$subl[3]\t$subl[4]\t$subl[5]\t$subl[6]\n";
        } elsif ($subl[4] =~ /Intensification/){
            print Resultsi "$subl[1]\t$subl[2]\t$subl[3]\t$subl[4]\t$subl[5]\t$subl[6]\n";
        } else {
            die "Cannot find Relaxation or Intensification in $line\n";
        }
    }

    if ($subl[5] <= 0.01){ # FDR 0.01
        print Resultsf "$subl[1]\t$subl[2]\t$subl[3]\t$subl[4]\t$subl[5]\t$subl[6]\n";
        if ($subl[4] =~ /Relaxation/){
            print Resultsfr "$subl[1]\t$subl[2]\t$subl[3]\t$subl[4]\t$subl[5]\t$subl[6]\n";
        } elsif ($subl[4] =~ /Intensification/){
            print Resultsfi "$subl[1]\t$subl[2]\t$subl[3]\t$subl[4]\t$subl[5]\t$subl[6]\n";
        } else {
            die "Cannot find Relaxation or Intensification in $line\n";
        }
    }

}
close Filef;
close Results;

