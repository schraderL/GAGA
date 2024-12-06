#!/usr/bin/perl
use strict; use warnings;

# Script to filter codons with more than the stablished percent of gaps
# Filter species with less than the number of codons and remove then those genes from the gene tree
# Root the gene trees in Leptanilla (GAGA-0392)

# Load required modules: module load ngs tools newick-utils/1.6

# Usage: perl clean_aln_fromhmmcleaner.pl inputaln.fasta(HOGID.cds.nonstop_hmmclean.ali) folderwithgenetree(named as HOGID.treefile) outputdir percentgapstofilter(i.e. 0.50 will filter codons with gaps in more than 50% genes)
#perl ../clean_aln_fromhmmcleaner.pl input/HOG0000492_parts_SPAN_1.cds.nonstop_clean.ali ../testfiles hmmcleaner_50_aln 0.50

# It will print a curated output from hmmclean, only filtering genes with low number of codons, in the same dir with extension cds.nonstop_hmmclean_curated.aln, and also the genetree filtering these sequences with name hmmclean.treefile
# In the outputdir, it will print the alignment filtering also codons with high percentage of gaps, and the corresponding gene trees

# Parameters
my $percgap = "$ARGV[3]"; # Percentage of gaps in a codon to filter (i.e. 0.50 will filter codons with gaps in more than 50% genes, that is keeping codons in at least 50% of the genes); Use 0.25 to keep codons that are at least in 75% of the genes.
my $minperclength = "0.2"; # Minimum percentage of nucleotides respect the average hog length to retain a gene. Total number is calculated from all genes (not counting gaps). For example, if a hog has average length of 100 nt, and a gene has only 19, it is discarded
my $generoot = "GAGA-0392"; # Pattern in gene name for the species to root the gene trees | It will only root when there is one gene, otherwise the gene tree is rooted in the mid point

# Input files
my $alnfile = "$ARGV[0]";
my $treefolder = "$ARGV[1]";
my $outdir = "$ARGV[2]";

my $alndir = ""; my $filename = "";
if ($alnfile =~ /(.*)\/(\S+)\.cds/){
    $alndir = $1;
    $filename = $2;
}

# Controls in the script
# Check length and add gaps to have same lenght across all seqs in the alignment (hmmclean trims the last codons if they are cleaned)
# Check codons, if only one/two nt (i.e.: one or two gaps only in the codon), print warning and correct codon with three gaps
# Check that the alignment is multiple of three

# Read alignment
my $line = ""; my $name = "";
my %fasta;
my $ognumseqs = "0"; my $ogalnlength = "0"; my $maxogalnlength = "0";
my $ogunalignedlength = "0";
my $outgroupgene = ""; my $outgroupgenecount = "0";
open(Filef, "<", $alnfile);
while(<Filef>){
    chomp;
    my $line2 = $_;
    next if ($line2 !~ /\S+/);
    next if ($line2 =~ /^\#/);
    if ($line2 =~ />(\S+)/){
        $name = $1;
        $fasta{$name} .= ""; # Initialize fasta to avoid skipping empty sequences
        $ognumseqs++;
        if ($ogalnlength > $maxogalnlength){ # Get maximum length to correct shorter sequences from hmmclean
            $maxogalnlength = $ogalnlength;
        }
        $ogalnlength = 0;
        if ($name =~ /$generoot/){
            $outgroupgene = $name;
            $outgroupgenecount++;
        }
    } else {
        $line2 =~ s/\ /-/g; ### HMMcleaner introduced spaces instead of - at the end of the sequences if there is trimming
        $fasta{$name} .= "$line2";
        #if ($ognumseqs == 1){ # only count aln length for the first sequence 
            $ogalnlength += length($line2);
        #}
        my $nogapseq = $line2;
        $nogapseq =~ s/\-//g;
        $ogunalignedlength += length($nogapseq);
    }
}
close Filef;

my %newfasta;
my %countgapcodon; my %countseqcodon;
my $genestodiscard = ""; my $genestodiscardcount = "0";
my $firstseq = "0";

foreach my $key (sort keys (%fasta)){
    my $seq = "";
    if (exists $fasta{$key}){
        $seq = $fasta{$key};
    }
    my $newseq = "";
    my $lengthprev = length($seq);
    if ($lengthprev < $maxogalnlength){
        print "Warning: $alnfile sequence $key is shorter than the alignment, filling with gaps...\n";
        my $diflen = $maxogalnlength - $lengthprev;
        my $addgap = "-" x $diflen;
        $newseq = "$seq"."$addgap";
    } else {
        $newseq = $seq;
    }

    my $lengthcheck = length($newseq)/3;
    my $nseq = "";

    if ($lengthcheck =~ /\.33/){
        #print Results "$line2"."NN\n"; 
        $nseq = "$newseq"."--";
        print "Warning: $alnfile sequence $seq is not multiple of 3!\n";   # Print seq id to check
    } elsif ($lengthcheck =~ /\.66/) {
        #print Results "$line2"."N\n";
        $nseq = "$newseq"."-";
        print "Warning: $alnfile sequence $seq is not multiple of 3!\n";   # Print seq id to check
    } else {
        #print Results "$line2\n";
        $nseq = "$newseq";
    }

    # Discard gene if it contains less than the percentage of nucleotides
    my $aveoglength = $ogunalignedlength/$ognumseqs;
    my $nogapseq = $nseq;
    $nogapseq =~ s/\-//g;
    my $nogapseqlength = length($nogapseq);
    my $filtval = $aveoglength * $minperclength;
    if ($nogapseqlength < $filtval){
        $genestodiscard .= "$key ";
        $genestodiscardcount++;
        next; # Filter this gene
    }

    # Go through the codon and count how many are gaps

    my $nseqcodoncorrected = "";

    my $length = length($nseq);
    for (my $n = 0; $n < ($length - 2); $n += 3) {

        if ($firstseq == 0){ # Initialize the values for all codons
            $countgapcodon{$n} = 0;
            $countseqcodon{$n} = 0;
        }

        my $codon = substr ($nseq, $n, 3);
        if ($codon =~ /TAA/ || $codon =~ /TAG/ || $codon =~ /TGA/ ){ # Stop codons
#           print Nonstop "---";
#           if ($n >= ($length - 3)){ # Only discard stop if it is the last one
#               $nseqcodoncorrected .= "";    
#           } else {
                $nseqcodoncorrected .= "---";  # Discard stop codons from the alignment
#           }
            $countgapcodon{$n}++;
        } elsif ($codon =~ /\w+\-/){ # Incomplete codon with gap
#           print Nonstop "---";
            print "Warning: Incomplete codon with gaps $codon in sequence $seq in file $alnfile \n";
            $nseqcodoncorrected .= "---";
            $countgapcodon{$n}++;  
        } elsif ($codon =~ /\-/){ # Gap regions
#           print Nonstop "---";
            $nseqcodoncorrected .= "---";
            $countgapcodon{$n}++;  
        }
        else {
#           print Nonstop "$codon";
            $nseqcodoncorrected .= "$codon";
            $countseqcodon{$n}++;
        }
    }
    $firstseq++;

    # Save final sequence
    $newfasta{$key} = $nseqcodoncorrected;

}

# Exclude filtered genes from gene tree
my $genetree = "$treefolder\/$filename\.treefile";
my $genetreeout = "$alndir\/$filename\.hmmclean.treefile";
if ($genestodiscardcount > 0){
    print "Excluding $genestodiscardcount genes $genestodiscard in $alndir\/$filename\.cds.nonstop_hmmclean_curated.aln\n";
    system ("nw_prune $genetree $genestodiscard > $genetreeout");
} else {
    system ("cp $genetree $genetreeout");
}


# Go now through the codons again to filter blocks with the percentage of gaps
my %newfastanogap; my %newfastacurated;
my $numbergenesretained = "0"; my $numbercodonsretained = "0";
#open (Resultscur, ">", "$alndir\/$filename\.cds.nonstop_hmmclean_curated.aln"); # Also filtering blocks with no alignment

foreach my $key (sort keys (%newfasta)){
    my $seq = $newfasta{$key};
#    print Resultscur ">$key\n$seq\n"; # Print the curated hmmalign file

   # Go through the codons and filter those with less codon sequences than the required percentage

    my $nseqcodonfilt = "";
    my $nseqcodoncurated = "";

    my $length = length($seq);
    for (my $n = 0; $n < ($length - 2); $n += 3) {
        my $codon = substr ($seq, $n, 3);

        unless (exists $countgapcodon{$n}){
            die "Error it cannot find the gap count for codon $n in $alnfile\n";
        }
        my $gappercent = $countgapcodon{$n}/$ognumseqs;

        # Save here the codon alignment, filtering only blocks with no sequence
#       if ($gappercent < 1 ){ ## 1 means 100% of the genes contain gaps in this codon
        if ($countseqcodon{$n} > 0 ){ ## At least one gene that has a codon sequence in this region
            $nseqcodoncurated .= "$codon";
        }

        # Now save the codon alignment filtering gap blocks with the required percentage
        if ($gappercent > $percgap ){
            # Skip as it contains more gaps that the required percentage
            next;
        } else {
            $nseqcodonfilt .= "$codon";
            $numbercodonsretained++;
        }

    }

    # Save curated hmmcleaner output
    $newfastacurated{$key} = $nseqcodoncurated;

    # Discard gene if it contains less than 15 codons after filtering gaps
    my $nogapseq = $nseqcodonfilt;
    $nogapseq =~ s/\-//g;
    my $nogapseqlength = length($nogapseq);
    if ($nogapseqlength < 15){
        $genestodiscard .= "$key ";
        $genestodiscardcount++;
        next; # Filter this gene
    }


    # Save final sequence
    $newfastanogap{$key} = $nseqcodonfilt;
    $numbergenesretained++;

}
#close Resultscur;


# Go now through the final seq with filtered codons to print the alignment

# Skip this hog if the number of final kept sequences is below 20, or the codon length is below 20.

if ($numbergenesretained < 16 || $numbercodonsretained < 30){
    # Skip this HOG
    print "Skipping $alnfile as it contains $numbergenesretained genes with $numbercodonsretained aligned codons\n";
    system ("mv $alndir\/$filename\.cds.nonstop_hmmclean_curated.aln $alndir\/$filename\.cds.nonstop_hmmclean_curated_oldbadqual.aln");
} else {

    open (Resultsc, ">", "$outdir\/$filename\.cds.nonstop_hmmclean_gapfilt.aln");
    open (Resultscur, ">", "$alndir\/$filename\.cds.nonstop_hmmclean_curated.aln"); # Also filtering blocks with no alignment


    foreach my $key (sort keys (%newfastanogap)){
        my $seq = $newfastanogap{$key};
        print Resultsc ">$key\n$seq\n"; # Print the curated hmmalign file
    }
    foreach my $key (sort keys (%newfastacurated)){
        my $seq = $newfastacurated{$key};
        print Resultscur ">$key\n$seq\n"; # Print the curated hmmalign file
    }    
    close Resultsc;
    close Resultscur;


    # Remove the genes from the gene tree in the final filtering file
    my $originalgenetree = "$treefolder\/$filename\.treefile";
    my $finalgenetreeout = "$outdir\/$filename\.hmmclean_gapfilt.treefile";
    if ($genestodiscardcount > 0){
        print "Excluding $genestodiscardcount genes $genestodiscard in $outdir\/$filename\.cds.nonstop_hmmclean_gapfilt.aln\n";
        system ("nw_prune $originalgenetree $genestodiscard > $finalgenetreeout");
    } else {
        system ("cp $originalgenetree $finalgenetreeout");
    }

    # Finally, root trees with leptanilla
    #my $finalgenetreeout = "$outdir\/$filename\.hmmclean_gapfilt.treefile";
    my $finalgenetreeoutroot = "$outdir\/$filename\.hmmclean_gapfilt.rooted.treefile";

    if ($outgroupgenecount == 1 && $genestodiscard !~ /$outgroupgene/){
        system ("nw_reroot $finalgenetreeout $outgroupgene > $finalgenetreeoutroot");
    } else {
        system ("nw_reroot $finalgenetreeout > $finalgenetreeoutroot");
    }


    # Check that gene tree and sequences printed match
    #nw_labels -I ../HOG0000308_parts_SPAN_1.hmmclean_gapfilt.treefile | wc -l

    system("nw_labels -I $genetreeout | wc -l > tmp_countreetips_$filename\.txt");
    my $tipcount = "0";
    my $seqcount = scalar keys %newfastacurated;
    open (Tmptree, "<", "tmp_countreetips_$filename\.txt");
    while (<Tmptree>){
        chomp;
        my $nline = $_;
        next if ($nline !~ /\S+/);
        $tipcount = $nline;
    }
    close Tmptree;
    if ($tipcount == $seqcount){
#        print "OK tree\n";
    } else {
        print "Error in $filename, number of tips in tree do not match with sequences in fasta: TREE $genetreeout ALN $alndir\/$filename\.cds.nonstop_hmmclean_curated.aln\n";
    }
    system("rm tmp_countreetips_$filename\.txt");

    system("nw_labels -I $finalgenetreeout | wc -l > tmp_countreetips_$filename\.txt");
    my $tipcountnogap = "0";
    my $seqcountnogap = scalar keys %newfastanogap;
    open (Tmptree, "<", "tmp_countreetips_$filename\.txt");
    while (<Tmptree>){
        chomp;
        my $nline = $_;
        next if ($nline !~ /\S+/);
        $tipcountnogap = $nline;
    }
    close Tmptree;
    if ($tipcountnogap == $seqcountnogap){
#        print "OK tree\n";
    } else {
        print "Error in $filename, number of tips in tree do not match with sequences in fasta: TREE $finalgenetreeout ALN $outdir\/$filename\.cds.nonstop_hmmclean_gapfilt.aln\n";
    }
    system("rm tmp_countreetips_$filename\.txt");

}





