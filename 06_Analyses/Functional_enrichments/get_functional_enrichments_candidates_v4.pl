#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

### Script to get the functional tables from a list of candidates, and conduct functional enrichments | NOTE the orthogroups are for the ants, if using other orthogroups (i.e. CAFE), change lines 55,56,57 and comment some parts, such as getting GAGA DE as it is based in ant orthogroups
## Specific steps
# Add the columns with expression data in the output functional annotation table - lines 472-473
# Add the columns with GAGA caste biased expression from Yue's analyses - lines 476 and 477 | NOTE change the path in the corresponding perl scripts
# GO enrichments: uncomment lines 703-731  | lines for simplyfy enrichment plots 1757 - 1761
# GO figure NOT USEFUL PLOTS - indicate script and uncomment lines 734 - 772
# KEGG enrichments: uncomment line 1676 that runs the R script 
# Mpha development plots: unzip files with plots and uncomment lines 1680 - 1740
# GAGA expression plots: uncomment lines 1744-1745 | for plots separated by groups, add species file and uncomment lines 1750-1752

## Specify the required files in lines 

# Usage: perl get_functional_enrichments_candidates_v4.pl Absrel_candidates_mergepart_single_test.txt Absrel_candidates_mergepart_single_fullstats.txt
# Example usage: perl ~/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/get_functional_enrichments_candidates_v4.pl HOG_summary_DEG_ants_75pc_species.tsv ../../Queen_worker_DE_summaries/HOG_summary_DEG_ants.tsv

# Computerome - Load modules: module load ngs tools gcc/7.4.0 intel/perflibs R/4.2.0 anaconda3/2023.03
# Necessary installed R modules: library(GOstats) library(GSEABase) library(GOstatsPlus) library(topGO) library(GO.db) library(ggplot2) library(scales)
# GOfigure only works with module load anaconda3/4.4.0 !!! Not with any 2023 versions

# Indicate if there are headers, it can be added as an input varialble by uncomment lines 26 and 27
my $clist = "$ARGV[0]"; # Candidate list file
my $flist = "$ARGV[1]"; # File with all tested HOGs
#my $headerclist = "$ARGV[2]"; # Add if the candidate list contains a first-line header (T) or not (F)
#my $headerflist = "$ARGV[3]"; # Add if the reference list contains a first-line header (T) or not (F)
my $headerclist = "T"; # Add if the candidate list contains a first-line header (T) or not (F)
my $headerflist = "T"; # Add if the reference list contains a first-line header (T) or not (F)
my $cafecandidates = "F"; # Indicate if the candidate files are from CAFE, as the HOG ids are in the second column

# Output prefix name
my $outname = "";
if ($clist =~ /(\S+)\.txt/){
    $outname = $1;
} else {
    $outname = $clist;
}

# Remove partitions from candidate list, and just get the HOG
my $outfile = "$outname\_nopart_ids.txt";
if ($headerclist =~ /T/){
    system ("tail -n+2 $clist | awk '{print \$1}' | sort | uniq > $outfile");
} else {
    system ("cat $clist | awk '{print \$1}' | sort | uniq > $outfile"); # No header in candidate list
}
if ($cafecandidates =~ /T/){ # The HOG ID in the CAFE candidate files is in the second column
    system ("cat $clist | awk '{print \$2}' | sort | uniq > $outfile"); # No header in candidate list
}

## Needed files for the script
my $goannotfile = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/Orthogroup_functional_annotations/Orthogroups/Orthogroups_functional_annotation_GO.annot";
my $fannotfile = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/Orthogroup_functional_annotations/Orthogroups/Orthogroups_functional_annotation_summary.tsv";
my $keggfile = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/Orthogroup_functional_annotations/Orthogroups/Orthogroups_functional_annotation_KEGG.annot";
#my $goannotfile = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/Orthogroup_functional_annotations/Orthogroups_withoutgroups/Orthogroups_functional_annotation_GO.annot"; # CAFE candidates
#my $fannotfile = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/Orthogroup_functional_annotations/Orthogroups_withoutgroups/Orthogroups_functional_annotation_summary.tsv"; # CAFE candidates
#my $keggfile = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/Orthogroup_functional_annotations/Orthogroups_withoutgroups/Orthogroups_functional_annotation_KEGG.annot"; # CAFE candidates
my $keggRdata = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/KEGG_levs.pathway.module.Rdata";
my $keggdescfile = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/KEGG_KO_description_fromR.txt";
my $obogofile = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Candidate_genes/GO_19Apr2023.obo"; # Read ovo file with GO descriptions to include in output table
my $mphafile = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/Mpha_gene_association_GAGA2NCBI.txt";
#my $mphafolder = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/Bitao_NEE/pharaoh_developmantal_exp_plot_v2/";
my $mphafolder = "/Users/joel/Documents/Results/Analyses_traits/Analyses_folder/Bitao_NEE/pharaoh_developmantal_exp_plot_v2/";
my $bitaocanalized = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/Bitao_NEE/NEE_supp_file_caste_bias_canalized_expression.tsv";
my $bitaode = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/Bitao_NEE/NEE_supp_file_caste_differential_expression.tsv";
my $bitaodempha = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/Bitao_NEE/NEE_supp_file_caste_differential_expression_larvaeMpha.tsv";
my $dmelfile = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/Dmel_FlyAtlas2_gene_data_Mpha.txt";


#my $gofigurescript = "/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/gofigure/GO-Figure/gofigure.py"; # GO figure script path
my $simplifyenrichmentscript = "";
my $scriptspath = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/"; # Path to the scripts, like get_expression topGO ...

## Start script

# Read GO descriptions
my %ovodesc; 
my $obogo = ""; my $oboname = "";
open(Filef, "<", $obogofile);
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    #my @subl = split (/\t/, $line);
    if ($line =~ /^id: (GO\S+)/){
        $obogo = $1;
    } elsif ($line =~ /^name: (.+)/) {
        $oboname = "$1";
    } elsif ($line =~ /^namespace: molecular_function/) {
        $oboname .= " [MF]";
        $ovodesc{$obogo} = $oboname;
    } elsif ($line =~ /^namespace: biological_process/) {
        $oboname .= " [BP]";
        $ovodesc{$obogo} = $oboname;
    } elsif ($line =~ /^namespace: cellular_component/) {
        $oboname .= " [CC]";
        $ovodesc{$obogo} = $oboname;
    }  elsif ($line =~ /^alt_id: (.+)/) {
        $obogo = "$1";
        $ovodesc{$obogo} = $oboname;
    }
}
close Filef;

# Read Mpha association table
my %mphatable; 
open(Filef, "<", $mphafile);
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    my @subl = split (/\t/, $line);
    $mphatable{$subl[0]} = $subl[2];
}
close Filef;

# Read Bitao candidates
my %bitaotablecanalized; my $bitaotablecanalizedheader = "NA"; my %bitaotablecanalizedfull;
open(Filef, "<", $bitaocanalized);
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    if ($line =~ /Canalized expression pattern/){ # Header
        $bitaotablecanalizedheader = $line;
        next;
    }
    my @subl = split (/\t/, $line);
    if (exists $bitaotablecanalized{$subl[0]}){
        next; # Already saved, might be repeated those DE in Mpha and Aech, that appear later only in Mpha
    } else {
        $bitaotablecanalized{$subl[0]} = $subl[1];
        $bitaotablecanalizedfull{$subl[0]} = $line;
    }
}
close Filef;

# Read Bitao candidates
my %bitaotablede; my $bitaotabledeheader = "NA"; my %bitaotabledefull;
open(Filef, "<", $bitaode);
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    if ($line =~ /Stage for caste DE/){ # Header
        $bitaotabledeheader = $line;
        next;
    }
    my @subl = split (/\t/, $line);
    if (exists $bitaotablede{$subl[0]}){ # Add the new DE stage into the same gene
        $bitaotablede{$subl[0]} .= ", $subl[1]";
        $bitaotabledefull{$subl[0]} = "$subl[0]\t$bitaotablede{$subl[0]}\t$subl[2]\t$subl[3]\t$subl[4]\t$subl[5]\t$subl[6]\t$subl[7]\t$subl[8]\n";        
    } else {
        $bitaotablede{$subl[0]} = $subl[1];
        $bitaotabledefull{$subl[0]} = $line;
    }
}
close Filef;

# Read Bitao candidates
my %bitaotabledempha; my $bitaotabledemphaheader = "NA"; my %bitaotabledemphafull;
open(Filef, "<", $bitaodempha);
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    if ($line =~ /Caste DE/){ # Header
        $bitaotabledemphaheader = $line;
        next;
    }
    my @subl = split (/\t/, $line);
    if (exists $bitaotabledempha{$subl[0]}){
        next; # Already saved, might be repeated those DE in Mpha and Aech, that appear later only in Mpha
    } else {
        $bitaotabledempha{$subl[0]} = $subl[1];
        $bitaotabledemphafull{$subl[0]} = $line;
    }
}
close Filef;

# Read Dmel ortholog and expression table
my %dmeltable; my %dmelexprtable; my $dmeltableheader = "NA"; my %dmeltablefull;
open(Filef, "<", $dmelfile);
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    if ($line =~ /Top5 Expression Level/){ # Header
        $dmeltableheader = $line;
        next;
    }
    $line =~ s/ +/ /g; # Replace multiple spaces with just one
    my @subl = split (/\t/, $line);
    if (exists $dmeltable{$subl[0]}){
        next; # Already saved, might be repeated those DE in Mpha and Aech, that appear later only in Mpha
    } else {
        $dmeltable{$subl[1]} = "$subl[3]\t$subl[5]"; # old line before editing:      $dmeltable{$subl[1]} = "$subl[5] $subl[3]";
        $dmelexprtable{$subl[1]} = "$subl[6]";
        $dmeltablefull{$subl[1]} = $line;
    }
}
close Filef;

#########

# Create reference set
open (Results, ">", "$outname\_population_GO.annot");
open (Resultstab, ">", "$outname\_population_GO.tsv");
open (Resultskegg, ">", "$outname\_population_KEGG.annot");
open (Resultsmpha, ">", "$outname\_population_Mphaids.annot");

my $referencegolist = "";
my $headercount = "0";
open(Filef, "<", $flist);
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    #next if ($line =~ /Branches under selection/ || $line =~ /relaxation/ || $line =~ /Average dN/); # Header for absrel, or relax tables
    if ($headerflist =~ /T/){ # Skip first line if header = T in the reference set of genes
        if ($headercount == 0){
            $headercount++;
            next;
        }
    }
    #my @subl = split (/\t/, $line);
    my $og = "";
    if ($line =~ /(HOG\d\d\d\d\d\d\d)/){
        $og = $1;
    } else {
        die "Cannot find OG in line 208 file: $flist line $line\n";
    }
    if ($referencegolist =~ / $og /){
        next; # OG already saved
    } else {
        $referencegolist .= " $og ";
    }
}
close Filef;

my %fullgotable;
my %fullgotablewdesc;
my %popgotable;

open(Filef, "<", $goannotfile);
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    #my @subl = split (/\t/, $line);
    my $og = "";
    if ($line =~ /(HOG\d\d\d\d\d\d\d)/){
        $og = $1;
    } else {
        die "Cannot find OG in line 232, file: $goannotfile line $line\n";
    }

    my @subl = split(/\t/, $line);
    $fullgotable{$subl[0]} .= "$subl[1] ";

    if ($referencegolist =~ / $og /){
        print Results "$line\n"; # HOG in reference, save.   
        if (exists $popgotable{$og}){
            $popgotable{$og} .= ",$subl[1]"; 
        } else {
            $popgotable{$og} = "$subl[1]"; # First GO, no comma needed
        }
    } else {
        next; # HOG not tested for selection
    }

    # Get the description for each GO
    if (exists $ovodesc{$subl[1]}){
        $fullgotablewdesc{$subl[0]} .= "$subl[1]: $ovodesc{$subl[1]}; ";
    } else {
        print "Warning: no description in obo file for $subl[1] in $subl[0]\n";
        $fullgotablewdesc{$subl[0]} .= "$subl[1]; ";
    }

}
close Filef;
close Results;

foreach my $hog (keys %popgotable){
    print Resultstab "$hog\t$popgotable{$hog}\n";
}
close Resultstab;

# Create candidate GO set
open (Results, ">", "$outname\_GO.annot");
my $candidategolist = "";
open(Filef, "<", $outfile);
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    next if ($line =~ /Branches under selection/ || $line =~ /relaxation/);
    #my @subl = split (/\t/, $line);
    my $og = "";
    if ($line =~ /(HOG\d\d\d\d\d\d\d)/){
        $og = $1;
    } else {
        die "Cannot find OG in line 280, file: $clist line $line\n";
    }

    if ($candidategolist =~ / $og /){
        next; # OG saved
    } else {
        $candidategolist .= " $og ";
    }
}
close Filef;

open(Filef, "<", $goannotfile);
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    #my @subl = split (/\t/, $line);
    my $og = "";
    if ($line =~ /(HOG\d\d\d\d\d\d\d)/){
        $og = $1;
    } else {
        die "Cannot find OG in line 301, file: $goannotfile line $line\n";
    }

    if ($candidategolist =~ / $og /){
        print Results "$line\n"; # HOG in candidate list, save.    
    } else {
        next; # HOG not in candidate
    }
}
close Filef;
close Results;

# KEGG
my %hogkegg;
open(Filef, "<", $keggfile);
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    #my @subl = split (/\t/, $line);
    my $og = "";
    if ($line =~ /(HOG\d\d\d\d\d\d\d)/){
        $og = $1;
    } else {
        die "Cannot find OG in line 324, file: $keggfile line $line\n";
    }

    my @subl = split(/\t/, $line);

    $hogkegg{$og} .= "$subl[1] ";

    if ($referencegolist =~ / $og /){
        print Resultskegg "$line\n"; # HOG in reference, save.   
    } else {
        next; # HOG not tested for selection
    }

}
close Filef;
close Resultskegg;

# Create list of Dmel genes in test and population
open (Resultskeggdmelpop, ">", "$outname\_population_KEGG_Dmel.txt");
open (Resultskeggdmel, ">", "$outname\_KEGG_Dmel.txt");

# Create candidate table with gene functions
open (Results, ">", "$outname\_funct_annot.tsv");
print Results "Orthogroup\tNumber of sequences\tNumber of species\tSwissprot\tNumber of seqs with that Swissprot\tDrosophila melanogaster swissprotbest hit\tNumber of seqs with that Dmel seq as best hit\tEggNog\tNumber of seqs with that EggNog\tTrEMBL\tNumber of seqs with that TrEML\tInterPro\tNumber of seqs with that InterPro\tPfam\tNumber of seqs with that Pfam\tSignalP\tNumber of seqs with SignalP detected\tNumber of GOs (at least 33\% sequences)\tTotal number of GOs\tNumber of KEGGs (at least 33\% sequences)\tTotal number of KEGGs\tMonomorium pharaonis gene names\tKEGG id\tMonomorium pharaonis ncbi gene names\tDrosophila melanogaster best repiprocal hit\tExpression in Dmel top5 tissues\tCaste bias canalized expression NEE\tCaste differential expression NEE\tCaste differential expression Mpha NEE\tGOs\tGOs with description\n";

open(Filef, "<", $fannotfile);
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    next if ($line =~ /Number of seqs with that Swissprot/);
    my @subl = split (/\t/, $line);
    my $og = "";
    if ($line =~ /(HOG\d\d\d\d\d\d\d)/){
        $og = $1;
    } else {
        die "Cannot find OG in line 354, file: $fannotfile line $line\n";
    }

    if ($candidategolist =~ / $og /){ # HOG in candidate list, save. 
        my $mphancbi = ""; # Getting the Mpha ncbi ids to add them in the output table
        my $dmelhit = ""; # Getting Dmel best reciprocal hit
        my $dmelexpr = "";
        my $neecan = ""; # NEE candidates from canalized expression
        my $neede = "";
        my $needempha = ""; 
        if (exists ($subl[21])){
            my $mphaids = $subl[21];
            my @sepmphaids = split(/\,/, $mphaids);
            foreach my $sepid (@sepmphaids){
                $sepid =~ s/\s//g;
                if (exists ($mphatable{$sepid})){
                    $mphancbi .= "$mphatable{$sepid} ";
                    print Resultsmpha "$og\t$mphatable{$sepid}\n";

                    if (exists $bitaotablecanalized{$mphatable{$sepid}}){
                        $neecan = "$bitaotablecanalized{$mphatable{$sepid}} ";
                    }
                    if (exists $bitaotablede{$mphatable{$sepid}}){
                        $neede = "$bitaotablede{$mphatable{$sepid}} ";
                    }
                    if (exists $bitaotabledempha{$mphatable{$sepid}}){
                        $needempha = "$bitaotabledempha{$mphatable{$sepid}} ";
                    }

                }
                if (exists $dmeltable{$sepid}){
                    $dmelhit .= "$dmeltable{$sepid} ";
                    $dmelexpr .= "$dmelexprtable{$sepid} ";
                    print Resultskeggdmelpop "$dmeltable{$sepid}\n";
                    print Resultskeggdmel "$dmeltable{$sepid}\n";
                }

            }
        }

        my $keggid = "NA";
        if (exists $hogkegg{$og}){
            $keggid = "$hogkegg{$og}";
        }

        if (exists $fullgotable{$og}){
            print Results "$line\t$keggid\t$mphancbi\t$dmelhit\t$dmelexpr\t$neecan\t$neede\t$needempha\t$fullgotable{$og}\t$fullgotablewdesc{$og}\n";  
        } else {
            print Results "$line\t$keggid\t$mphancbi\t$dmelhit\t$dmelexpr\t$neecan\t$neede\t$needempha\tNA\tNA\n";  
        }
    } elsif ($referencegolist =~ / $og /){ # HOG in the tested list, save the Mpha id for kegg enrichments as population
        my $dmelhit = ""; # Getting Dmel best reciprocal hit
        my $mphancbi = ""; # Getting the Mpha ncbi ids to add them in the output table
        if (exists ($subl[21])){
            my $mphaids = $subl[21];
            my @sepmphaids = split(/\,/, $mphaids);
            foreach my $sepid (@sepmphaids){
                $sepid =~ s/\s//g;

                if (exists ($mphatable{$sepid})){
                    $mphancbi .= "$mphatable{$sepid} ";
                    print Resultsmpha "$og\t$mphatable{$sepid}\n";
                }

                if (exists $dmeltable{$sepid}){

                    #$dmelhit .= "$dmeltable{$sepid}\n";
                    print Resultskeggdmelpop "$dmeltable{$sepid}\n";
                }

            }
        }

    } else {
        next; # HOG not in candidate or population
    }
}
close Filef;
close Results;
close Resultsmpha;
close Resultskeggdmel;
close Resultskeggdmelpop;

# Add additional expression data 
#system ("perl /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/get_expression_tables_additionaldata.pl $outname\_funct_annot.tsv");
system ("perl $scriptspath\/get_expression_tables_additionaldata_wkeggidcol_v2.pl $outname\_funct_annot.tsv");
system ("sed \'s/\\\./,/g\' $outname\_funct_annot_add_expression.tsv > $outname\_funct_annot_add_expression_commadec.tsv");

# Add GAGA expression results
system ("perl $scriptspath\/get_DE_results_for_hog_list.pl $outname\_funct_annot_add_expression.tsv ");
system ("perl $scriptspath\/get_DE_workers_results_for_hog_list_v2.pl $outname\_funct_annot_add_expression.tsv_addGAGADE.tsv ");

## GO enrichments

# GO stats package

# Create R script
open (ResultsR, ">", "$outname\_GOenrichment_script.R");
print ResultsR 'library(GOstats)
library(GSEABase)
library(GOstatsPlus)

# load function b2g_to_gsc
b2g_to_gsc = function(file, organism="my-organism"){
  blast2GOannot = read.table(file=file, sep="\t", fill=TRUE)
  colnames(blast2GOannot)<-c("Seq", "GO")
  
  # Create GO package
  GOdata<-data.frame(GO.id = blast2GOannot$GO,
                     evidence.code = rep("ISA", dim(blast2GOannot)[1]),
                     gene.id = as.character(blast2GOannot$Seq))
  goFrame = GOFrame(GOdata, organism = organism)
  goAllFrame = GOAllFrame(goFrame)
  
  GO_gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
  #  GO_mappings<-getGOFrameData(goAllFrame)
  
  return(GO_gsc)
}

# Here it is the reference set (i.e.; all orthologs)
GO_gsc = b2g_to_gsc(file="'."$outname\_population_GO.annot".'",organism="Anthog") # file Annot consist in "Gene GO"


# Gene list to run the GO enrichment
file = read.csv2(file = "'."$outfile".'", header = F) # Gene list

# enrichment (p-value = 0.05)
gene_IDs_of_interest = as.character(file$V1)

GO_BPpval = test_GO(gene_IDs_of_interest, ontology = "BP", gsc=GO_gsc, pval = 0.05)
GO_MFpval = test_GO(gene_IDs_of_interest, ontology = "MF", gsc=GO_gsc, pval = 0.05)
GO_CCpval = test_GO(gene_IDs_of_interest, ontology = "CC", gsc=GO_gsc, pval = 0.05)

# Print tables
write.table(summary(GO_BPpval), file = "'."$outname\_GOenrichment_BP_pval005.txt".'",quote=FALSE, sep = "\t")
write.table(summary(GO_BPpval, categorySize = 10), file = "'."$outname\_GOenrichment_BP_pval005_catsize10.txt".'",quote=FALSE, sep = "\t")
write.table(summary(GO_MFpval), file = "'."$outname\_GOenrichment_MF_pval005.txt".'",quote=FALSE, sep = "\t")

go_termspval = unlist(lapply(c(GO_BPpval, GO_MFpval, GO_CCpval),
                         function(d){significant_terms(d, cutoff = 0.05)}
))
#go_terms
go_termspval_df <- data.frame(GOterm = names(go_termspval), enrichment_P_value = go_termspval)

# Table to visualize in REVIGO (http://revigo.irb.hr/)
write.table(go_termspval_df, file = "'."$outname\_GOenrichment_GOterms_pval005.txt".'", quote = FALSE, sep = "\t", col.names = c("% GOterm", "enrichment_P-value"), row.names = FALSE)


# enrichment (p-value = 0.01)
gene_IDs_of_interest = as.character(file$V1)

GO_BP = test_GO(gene_IDs_of_interest, ontology = "BP", gsc=GO_gsc, pval = 0.01)
GO_MF = test_GO(gene_IDs_of_interest, ontology = "MF", gsc=GO_gsc, pval = 0.01)
GO_CC = test_GO(gene_IDs_of_interest, ontology = "CC", gsc=GO_gsc, pval = 0.01)

# result summary
#head(summary(GO_BP))
#summary(GO_MF)
#summary(GO_BP)
#summary(GO_BP, categorySize = 10) # Because significant results for GO terms with fewer than 10 genes may reflect mere chance rather than evidence for enrichment

# Print tables
write.table(summary(GO_MF), file = "'."$outname\_GOenrichment_MF_pval001.txt".'",quote=FALSE, sep = "\t")
write.table(summary(GO_BP), file = "'."$outname\_GOenrichment_BP_pval001.txt".'",quote=FALSE, sep = "\t")
write.table(summary(GO_MF, categorySize = 10), file = "'."$outname\_GOenrichment_MF_pval001_catsize10.txt".'",quote=FALSE, sep = "\t")
write.table(summary(GO_BP, categorySize = 10), file = "'."$outname\_GOenrichment_BP_pval001_catsize10.txt".'",quote=FALSE, sep = "\t")

go_terms = unlist(lapply(c(GO_BP, GO_MF, GO_CC),
                         function(d){significant_terms(d, cutoff = 0.01)}
))
#go_terms
go_terms_df <- data.frame(GOterm = names(go_terms), enrichment_P_value = go_terms)

# Table to visualize in REVIGO (http://revigo.irb.hr/)
#write.table(go_terms, file = "'."$outname\_GOenrichment_GOterms_pval001.txt".'",quote=FALSE)
write.table(go_terms_df, file = "'."$outname\_GOenrichment_GOterms_pval001.txt".'", quote = FALSE, sep = "\t", col.names = c("% GOterm", "enrichment_P-value"), row.names = FALSE)

#### Create some plots with the BP and MF enrichment results
library(ggplot2)
library(scales)

goEnrichment_result <- "'."$outname\_GOenrichment_BP_pval001.txt".'"
topNum <- 20 ## Number of GOs to show
goEnrichment<-read.table(goEnrichment_result,sep="\t",header=T)
goEnrichment$Pvalue <- as.numeric(goEnrichment$Pvalue)
goEnrichment <- goEnrichment[goEnrichment$Pvalue < 0.05,] 
goEnrichment <- goEnrichment[,c("GOBPID","Term","Pvalue")]

ntop <- min(topNum, nrow(goEnrichment)) # If there are less than the specified number of representative GOs, just use all from the table
ggdata <- goEnrichment[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
pdf(paste(goEnrichment_result,".PvalPlot.pdf",sep=""),height=16,width=12)
ggplot(ggdata,
       aes(x = Term, y = -log10(Pvalue), size = -log10(Pvalue), fill = -log10(Pvalue))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = \'royalblue\', high = \'red4\') +
  
  xlab(\'\') + ylab(\'Enrichment score\') +
  labs(
    title = \'GO enrichment, Top\',
    subtitle = paste(topNum," terms ordered by p-value"),
    caption = \'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001\') +
  
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "dotted", "dotted"),
             colour = c("red", "red", "red"),
             size = c(1, 1, 1)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = \'right\',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = \'bold\', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = \'bold\', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = \'bold\', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = \'bold\', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = \'bold\', vjust = 0.5),
    axis.title = element_text(size = 12, face = \'bold\'),
    axis.title.x = element_text(size = 12, face = \'bold\'),
    axis.title.y = element_text(size = 12, face = \'bold\'),
    axis.line = element_line(colour = \'black\'),
    
    legend.key = element_blank(), 
    legend.key.size = unit(1, "cm"), 
    legend.text = element_text(size = 14, face = "bold"),
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()

dev.off()

goEnrichment_result <- "'."$outname\_GOenrichment_MF_pval001.txt".'"
topNum <- 20
goEnrichment<-read.table(goEnrichment_result,sep="\t",header=T)
goEnrichment$Pvalue <- as.numeric(goEnrichment$Pvalue)
goEnrichment <- goEnrichment[goEnrichment$Pvalue < 0.05,] 
goEnrichment <- goEnrichment[,c("GOMFID","Term","Pvalue")]

ntop <- min(topNum, nrow(goEnrichment)) # If there are less than the specified number of representative GOs, just use all from the table
ggdata <- goEnrichment[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
pdf(paste(goEnrichment_result,".PvalPlot.pdf",sep=""),height=16,width=12)
ggplot(ggdata,
       aes(x = Term, y = -log10(Pvalue), size = -log10(Pvalue), fill = -log10(Pvalue))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = \'royalblue\', high = \'red4\') +
  
  xlab(\'\') + ylab(\'Enrichment score\') +
  labs(
    title = \'GO enrichment, Top\',
    subtitle = paste(topNum," terms ordered by p-value"),
    caption = \'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001\') +
  
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "dotted", "dotted"),
             colour = c("red", "red", "red"),
             size = c(1, 1, 1)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = \'right\',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = \'bold\', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = \'bold\', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = \'bold\', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = \'bold\', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = \'bold\', vjust = 0.5),
    axis.title = element_text(size = 12, face = \'bold\'),
    axis.title.x = element_text(size = 12, face = \'bold\'),
    axis.title.y = element_text(size = 12, face = \'bold\'),
    axis.line = element_line(colour = \'black\'),
    
    legend.key = element_blank(), 
    legend.key.size = unit(1, "cm"), 
    legend.text = element_text(size = 14, face = "bold"),
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()

dev.off()

# Uglier bar plot
#pdf(paste(goEnrichment_result,".PvalbarPlot.pdf",sep=""),height=16,width=12)
#ggplot(ggdata, aes(x=Term, y=-log10(Pvalue))) +
#  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
#  xlab("Biological process") +
#  ylab("Enrichment") +
#  ggtitle("GO enrichment") +
#  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$Pvalue)), by = 2), 1)) +
#  theme_bw(base_size=24) +
#  theme(
#    legend.position=\'none\',
#    legend.background=element_rect(),
#    plot.title=element_text(angle=0, size=18, face="bold", vjust=1),
#    axis.text.x=element_text(angle=0, size=12, face="bold", hjust=1.10),
#    axis.text.y=element_text(angle=0, size=12, face="bold", vjust=0.5),
#    axis.title=element_text(size=16, face="bold"),
#    legend.key=element_blank(),     #removes the border
#    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
#    legend.text=element_text(size=12),  #Text size
#    title=element_text(size=16)) +
#  guides(colour=guide_legend(override.aes=list(size=2.5))) +
#  coord_flip()
#dev.off()

';
close ResultsR;
#=head
# Run Rscript
#print "Rscript $outname\_GOenrichment_script.R\n";
system ("Rscript $outname\_GOenrichment_script.R");

## TopGO
#print "Rscript /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/topGO_enrich.R $outname\_population_GO.tsv $outfile\n";
system ("Rscript $scriptspath\/topGO_enrich.R $outname\_population_GO.tsv $outfile");

# TopGO figures
system ("Rscript $scriptspath\/topGO_showTopTerms.R $outfile\.CC_elim_pval001.txt 20");
system ("Rscript $scriptspath\/topGO_showTopTerms.R $outfile\.CC_classic_pval001.txt 20");
system ("Rscript $scriptspath\/topGO_showTopTerms.R $outfile\.CC_weight01_pval001.txt 20");
system ("Rscript $scriptspath\/topGO_showTopTerms.R $outfile\.MF_elim_pval001.txt 20");
system ("Rscript $scriptspath\/topGO_showTopTerms.R $outfile\.MF_classic_pval001.txt 20");
system ("Rscript $scriptspath\/topGO_showTopTerms.R $outfile\.MF_weight01_pval001.txt 20");
system ("Rscript $scriptspath\/topGO_showTopTerms.R $outfile\.BP_elim_pval001.txt 20");
system ("Rscript $scriptspath\/topGO_showTopTerms.R $outfile\.BP_classic_pval001.txt 20");
system ("Rscript $scriptspath\/topGO_showTopTerms.R $outfile\.BP_weight01_pval001.txt 20");


# Revigo type plots - rrvgo package
system ("Rscript $scriptspath\/analyse_GOenrichment_rrvgo.R $outname\_GOenrichment_BP_pval001.txt");
system ("Rscript $scriptspath\/analyse_GOenrichment_rrvgo.R $outname\_GOenrichment_BP_pval005.txt");
system ("Rscript $scriptspath\/analyse_GOenrichment_rrvgo.R $outfile\.BP_classic_pval001.txt");
system ("Rscript $scriptspath\/analyse_GOenrichment_rrvgo.R $outfile\.BP_classic_pval005.txt");
system ("Rscript $scriptspath\/analyse_GOenrichment_rrvgo.R $outfile\.BP_elim_pval001.txt");
system ("Rscript $scriptspath\/analyse_GOenrichment_rrvgo.R $outfile\.BP_elim_pval005.txt");
#=cut

# GO figures
=h
#GOFigure
system ("python $gofigurescript -i $outname\_GOenrichment_GOterms_pval001.txt -o $outname\_GOfigure_GOenrichment_GOstats_pval001_si02 -f small -m 20 -e 50 -si 0.2");
system ("python $gofigurescript -i $outname\_GOenrichment_GOterms_pval001.txt -o $outname\_GOfigure_GOenrichment_GOstats_pval001_si05 -f small -m 20 -e 50 -si 0.5");
system ("python $gofigurescript -i $outname\_GOenrichment_GOterms_pval001.txt -o $outname\_GOfigure_GOenrichment_GOstats_pval001_si08 -f small -m 20 -e 50 -si 0.8");

#Additional figures
system ("Rscript $scriptspath\/topGO_showTopTerms_GOFigure.R $outname\_GOfigure_GOenrichment_GOstats_pval001_si02 20\n");
system ("Rscript $scriptspath\/topGO_showTopTerms_GOFigure.R $outname\_GOfigure_GOenrichment_GOstats_pval001_si05 20\n");
system ("Rscript $scriptspath\/topGO_showTopTerms_GOFigure.R $outname\_GOfigure_GOenrichment_GOstats_pval001_si08 20\n");

#GOfigure topgo results
#GOFigure
system ("python $gofigurescript -i $outfile\.GOterms_elim_pval001.txt -o $outname\_GOfigure_GOenrichment_TopGO_elim_pval001_si02 -f small -m 20 -e 50 -si 0.2");
system ("python $gofigurescript -i $outfile\.GOterms_elim_pval001.txt -o $outname\_GOfigure_GOenrichment_TopGO_elim_pval001_si05 -f small -m 20 -e 50 -si 0.5");
system ("python $gofigurescript -i $outfile\.GOterms_elim_pval001.txt -o $outname\_GOfigure_GOenrichment_TopGO_elim_pval001_si08 -f small -m 20 -e 50 -si 0.8");

system ("python $gofigurescript -i $outfile\.GOterms_classic_pval001.txt -o $outname\_GOfigure_GOenrichment_TopGO_classic_pval001_si02 -f small -m 20 -e 50 -si 0.2");
system ("python $gofigurescript -i $outfile\.GOterms_classic_pval001.txt -o $outname\_GOfigure_GOenrichment_TopGO_classic_pval001_si05 -f small -m 20 -e 50 -si 0.5");
system ("python $gofigurescript -i $outfile\.GOterms_classic_pval001.txt -o $outname\_GOfigure_GOenrichment_TopGO_classic_pval001_si08 -f small -m 20 -e 50 -si 0.8");

system ("python $gofigurescript -i $outfile\.GOterms_weight01_pval001.txt -o $outname\_GOfigure_GOenrichment_TopGO_weight01_pval001_si02 -f small -m 20 -e 50 -si 0.2");
system ("python $gofigurescript -i $outfile\.GOterms_weight01_pval001.txt -o $outname\_GOfigure_GOenrichment_TopGO_weight01_pval001_si05 -f small -m 20 -e 50 -si 0.5");
system ("python $gofigurescript -i $outfile\.GOterms_weight01_pval001.txt -o $outname\_GOfigure_GOenrichment_TopGO_weight01_pval001_si08 -f small -m 20 -e 50 -si 0.8");

#Additional figures
system ("Rscript $scriptspath\/topGO_showTopTerms_GOFigure.R $outname\_GOfigure_GOenrichment_TopGO_elim_pval001_si02 20\n");
system ("Rscript $scriptspath\/topGO_showTopTerms_GOFigure.R $outname\_GOfigure_GOenrichment_TopGO_elim_pval001_si05 20\n");
system ("Rscript $scriptspath\/topGO_showTopTerms_GOFigure.R $outname\_GOfigure_GOenrichment_TopGO_elim_pval001_si08 20\n");

system ("Rscript $scriptspath\/topGO_showTopTerms_GOFigure.R $outname\_GOfigure_GOenrichment_TopGO_classic_pval001_si02 20\n");
system ("Rscript $scriptspath\/topGO_showTopTerms_GOFigure.R $outname\_GOfigure_GOenrichment_TopGO_classic_pval001_si05 20\n");
system ("Rscript $scriptspath\/topGO_showTopTerms_GOFigure.R $outname\_GOfigure_GOenrichment_TopGO_classic_pval001_si08 20\n");

system ("Rscript $scriptspath\/topGO_showTopTerms_GOFigure.R $outname\_GOfigure_GOenrichment_TopGO_weight01_pval001_si02 20\n");
system ("Rscript $scriptspath\/topGO_showTopTerms_GOFigure.R $outname\_GOfigure_GOenrichment_TopGO_weight01_pval001_si05 20\n");
system ("Rscript $scriptspath\/topGO_showTopTerms_GOFigure.R $outname\_GOfigure_GOenrichment_TopGO_weight01_pval001_si08 20\n");

=cut

## KEGG enrichment - vold
=h
# Create R script
open (ResultsR, ">", "$outname\_KEGGenrichment_script.R");
print ResultsR 'library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(dplyr)

load("'."$keggRdata".'")

# input files
ko2geneList <- "'."$outname\_population_KEGG.annot".'"
mpha2geneList <- "'."$outname\_population_Mphaids.annot".'"
targetGeneList <- "'."$outfile".'"

# read files
#ko2geneID <- read_table(ko2geneList,col_names=c("Ko","GID"))
ko2geneID <- read_table(ko2geneList,col_names=c("GID","Ko"))
ko2geneID <- ko2geneID[c("Ko", "GID")]

targetGene <- read_table(targetGeneList,col_names=c("GID"))

loc2geneID <- read_table(mpha2geneList,col_names=c("GID","LOCID"))

# Enrichment of enzymes, K_ID

kegg_enrich <- enricher(targetGene$GID,TERM2GENE = ko2geneID,TERM2NAME = keg_lvs[,c(1,4)],pAdjustMethod = "BH",pvalueCutoff  = 0.01)
#write_csv(as.data.frame(kegg_enrich),paste0(targetGeneList,".KEGG_enzyme.csv"))
write.table(as.data.frame(kegg_enrich),file = paste(targetGeneList,".KEGG_enzyme.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# Enrichment of pathways

pathway2gene <- inner_join(ko2geneID, ko2pathway, by = "Ko")
pathway2gene <- pathway2gene[c("Pathway", "GID")]
colnames(pathway2gene) <- c("Ko", "GID")

kegg_enrich <- enricher(targetGene$GID,TERM2GENE = pathway2gene,TERM2NAME = keg_lvs[,c(1,4)],pAdjustMethod = "BH",pvalueCutoff  = 0.01)
write.table(as.data.frame(kegg_enrich),file = paste(targetGeneList,".KEGG_pathway.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

topNum <- 20
#goEnrichment<-read.table(goEnrichment_result,sep="\t",header=T)
goEnrichment<-as.data.frame(kegg_enrich)
goEnrichment$pvalue <- as.numeric(goEnrichment$pvalue)
goEnrichment <- goEnrichment[goEnrichment$pvalue < 0.05,] 
goEnrichment <- goEnrichment[,c("ID","Description","pvalue")]

ntop <- min(as.numeric(topNum), nrow(goEnrichment)) # If there are less than the specified number of representative GOs, just use all from the table
ggdata <- goEnrichment[1:ntop,]
ggdata$Description <- factor(ggdata$Description, levels = rev(ggdata$Description)) # fixes order
pdf(paste(targetGeneList,".KEGG_pathway_Plot.pdf",sep=""),height=16,width=12)
ggplot(ggdata,
       aes(x = Description, y = -log10(pvalue), size = -log10(pvalue), fill = -log10(pvalue))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = "royalblue", high = "red4") +
  
  xlab("") + ylab("Enrichment score") +
  labs(
    title = "KEGG enrichment, Top",
    subtitle = paste(ntop," terms ordered by p-value"),
    caption = "Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001") +
  
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "dotted", "dotted"),
             colour = c("red", "red", "red"),
             size = c(1, 1, 1)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = "right",
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = "bold", vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = "bold", vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = "bold", vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = "bold", hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = "bold", vjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.line = element_line(colour = "black"),
    
    legend.key = element_blank(), 
    legend.key.size = unit(1, "cm"), 
    legend.text = element_text(size = 14, face = "bold"),
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()

dev.off()

### KEGG enrichment using enrichKEGG
# Using all KEGG universe

targetGene2Kid <- inner_join(targetGene, ko2geneID, by = "GID")
targetGeneKid <- targetGene2Kid$Ko
universeGeneKid <- ko2geneID$Ko

kegg_path_enrich <- enrichKEGG(targetGeneKid, organism="ko", keyType="kegg")
write.table(as.data.frame(kegg_path_enrich),file = paste(targetGeneList,".KEGG_allko_pathway.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

topNum <- 20
#goEnrichment<-read.table(goEnrichment_result,sep="\t",header=T)
goEnrichment<-as.data.frame(kegg_path_enrich)
goEnrichment$pvalue <- as.numeric(goEnrichment$pvalue)
goEnrichment <- goEnrichment[goEnrichment$pvalue < 0.05,] 
goEnrichment <- goEnrichment[,c("ID","Description","pvalue")]

ntop <- min(as.numeric(topNum), nrow(goEnrichment)) # If there are less than the specified number of representative GOs, just use all from the table
ggdata <- goEnrichment[1:ntop,]
ggdata$Description <- factor(ggdata$Description, levels = rev(ggdata$Description)) # fixes order
pdf(paste(targetGeneList,".KEGG_allko_pathway_Plot.pdf",sep=""),height=16,width=12)
ggplot(ggdata,
       aes(x = Description, y = -log10(pvalue), size = -log10(pvalue), fill = -log10(pvalue))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = "royalblue", high = "red4") +
  
  xlab("") + ylab("Enrichment score") +
  labs(
    title = "KEGG enrichment, Top",
    subtitle = paste(ntop," terms ordered by p-value"),
    caption = "Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001") +
  
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "dotted", "dotted"),
             colour = c("red", "red", "red"),
             size = c(1, 1, 1)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = "right",
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = "bold", vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = "bold", vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = "bold", vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = "bold", hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = "bold", vjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.line = element_line(colour = "black"),
    
    legend.key = element_blank(), 
    legend.key.size = unit(1, "cm"), 
    legend.text = element_text(size = 14, face = "bold"),
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()

dev.off()

# Using only the total set of KEGG universe

kegg_path_enrich <- enrichKEGG(targetGeneKid, organism="ko", keyType="kegg", universe = universeGeneKid)
write.table(as.data.frame(kegg_path_enrich),file = paste(targetGeneList,".KEGG_ko_pathway.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

topNum <- 20
#goEnrichment<-read.table(goEnrichment_result,sep="\t",header=T)
goEnrichment<-as.data.frame(kegg_path_enrich)
goEnrichment$pvalue <- as.numeric(goEnrichment$pvalue)
goEnrichment <- goEnrichment[goEnrichment$pvalue < 0.05,] 
goEnrichment <- goEnrichment[,c("ID","Description","pvalue")]

ntop <- min(as.numeric(topNum), nrow(goEnrichment)) # If there are less than the specified number of representative GOs, just use all from the table
ggdata <- goEnrichment[1:ntop,]
ggdata$Description <- factor(ggdata$Description, levels = rev(ggdata$Description)) # fixes order
pdf(paste(targetGeneList,".KEGG_ko_pathway_Plot.pdf",sep=""),height=16,width=12)
ggplot(ggdata,
       aes(x = Description, y = -log10(pvalue), size = -log10(pvalue), fill = -log10(pvalue))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = "royalblue", high = "red4") +
  
  xlab("") + ylab("Enrichment score") +
  labs(
    title = "KEGG enrichment, Top",
    subtitle = paste(ntop," terms ordered by p-value"),
    caption = "Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001") +
  
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "dotted", "dotted"),
             colour = c("red", "red", "red"),
             size = c(1, 1, 1)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = "right",
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = "bold", vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = "bold", vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = "bold", vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = "bold", hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = "bold", vjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.line = element_line(colour = "black"),
    
    legend.key = element_blank(), 
    legend.key.size = unit(1, "cm"), 
    legend.text = element_text(size = 14, face = "bold"),
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()

dev.off()

# Using Mpha list from KEGG

loc2geneID <- loc2geneID %>% mutate(LOCID = stringr::str_remove(LOCID, "LOC"))
targetGene2loc <- inner_join(targetGene, loc2geneID, by = "GID")
targetGeneloc <- targetGene2loc$LOCID
universeGeneloc <- loc2geneID$LOCID

kegg_path_enrich <- enrichKEGG(targetGeneloc, organism="mpha", keyType="kegg")
#kegg_path_enrich2 <- enrichKEGG(targetGeneloc, organism="mpha", keyType="kegg", universe = universeGeneloc)
write.table(as.data.frame(kegg_path_enrich),file = paste(targetGeneList,".KEGG_allmpha_pathway.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

topNum <- 20
#goEnrichment<-read.table(goEnrichment_result,sep="\t",header=T)
goEnrichment<-as.data.frame(kegg_path_enrich)
goEnrichment$pvalue <- as.numeric(goEnrichment$pvalue)
goEnrichment <- goEnrichment[goEnrichment$pvalue < 0.05,] 
goEnrichment <- goEnrichment[,c("ID","Description","pvalue")]

ntop <- min(as.numeric(topNum), nrow(goEnrichment)) # If there are less than the specified number of representative GOs, just use all from the table
ggdata <- goEnrichment[1:ntop,]
ggdata$Description <- factor(ggdata$Description, levels = rev(ggdata$Description)) # fixes order
pdf(paste(targetGeneList,".KEGG_allmpha_pathway_Plot.pdf",sep=""),height=16,width=12)
ggplot(ggdata,
       aes(x = Description, y = -log10(pvalue), size = -log10(pvalue), fill = -log10(pvalue))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = "royalblue", high = "red4") +
  
  xlab("") + ylab("Enrichment score") +
  labs(
    title = "KEGG enrichment, Top",
    subtitle = paste(ntop," terms ordered by p-value"),
    caption = "Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001") +
  
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "dotted", "dotted"),
             colour = c("red", "red", "red"),
             size = c(1, 1, 1)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = "right",
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = "bold", vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = "bold", vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = "bold", vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = "bold", hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = "bold", vjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.line = element_line(colour = "black"),
    
    legend.key = element_blank(), 
    legend.key.size = unit(1, "cm"), 
    legend.text = element_text(size = 14, face = "bold"),
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()

dev.off()

# Using Mpha list from KEGG, with specific universe

kegg_path_enrich <- enrichKEGG(targetGeneloc, organism="mpha", keyType="kegg", universe = universeGeneloc)
write.table(as.data.frame(kegg_path_enrich),file = paste(targetGeneList,".KEGG_mpha_pathway.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

topNum <- 20
#goEnrichment<-read.table(goEnrichment_result,sep="\t",header=T)
goEnrichment<-as.data.frame(kegg_path_enrich)
goEnrichment$pvalue <- as.numeric(goEnrichment$pvalue)
goEnrichment <- goEnrichment[goEnrichment$pvalue < 0.05,] 
goEnrichment <- goEnrichment[,c("ID","Description","pvalue")]

ntop <- min(as.numeric(topNum), nrow(goEnrichment)) # If there are less than the specified number of representative GOs, just use all from the table
ggdata <- goEnrichment[1:ntop,]
ggdata$Description <- factor(ggdata$Description, levels = rev(ggdata$Description)) # fixes order
pdf(paste(targetGeneList,".KEGG_mpha_pathway_Plot.pdf",sep=""),height=16,width=12)
ggplot(ggdata,
       aes(x = Description, y = -log10(pvalue), size = -log10(pvalue), fill = -log10(pvalue))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = "royalblue", high = "red4") +
  
  xlab("") + ylab("Enrichment score") +
  labs(
    title = "KEGG enrichment, Top",
    subtitle = paste(ntop," terms ordered by p-value"),
    caption = "Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001") +
  
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "dotted", "dotted"),
             colour = c("red", "red", "red"),
             size = c(1, 1, 1)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = "right",
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = "bold", vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = "bold", vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = "bold", vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = "bold", hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = "bold", vjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.line = element_line(colour = "black"),
    
    legend.key = element_blank(), 
    legend.key.size = unit(1, "cm"), 
    legend.text = element_text(size = 14, face = "bold"),
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()

dev.off()

';
close ResultsR;

# Run Rscript
#print "Rscript $outname\_KEGGenrichment_script.R\n";
system ("Rscript $outname\_KEGGenrichment_script.R");
=cut

# KEGG enrichment v3
## KEGG enrichment
# Create R script
open (ResultsR, ">", "$outname\_KEGGenrichment_script.R");
print ResultsR 'library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(dplyr)
library(pathview)

load("'."$keggRdata".'")

# input files
ko2geneList <- "'."$outname\_population_KEGG.annot".'"
mpha2geneList <- "'."$outname\_population_Mphaids.annot".'"
dmel2geneList <- "'."$outname\_population_KEGG_Dmel.txt".'"

targetGeneList <- "'."$outfile".'"
targetdmelGeneList <- "'."$outname\_KEGG_Dmel.txt".'"


# read files
#ko2geneID <- read_table(ko2geneList,col_names=c("Ko","GID"))
ko2geneID <- read_table(ko2geneList,col_names=c("GID","Ko"))
ko2geneID <- ko2geneID[c("Ko", "GID")]

targetGene <- read_table(targetGeneList,col_names=c("GID"))

loc2geneID <- read_table(mpha2geneList,col_names=c("GID","LOCID"))

dmelpopulation <- read.table(dmel2geneList, sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names = c("GID", "Function"))
# Remove spaces in columns
dmelpopulation$GID <- gsub(" ", "", dmelpopulation$GID)
dmelpopulation$Function <- gsub(" ", "_", dmelpopulation$Function)

dmeltargetGene <- read.table(targetdmelGeneList, sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names = c("GID", "Function"))
dmeltargetGene$GID <- gsub(" ", "", dmeltargetGene$GID)
dmeltargetGene$Function <- gsub(" ", "_", dmeltargetGene$Function)


# Enrichment of enzymes, K_ID

KOdesc <- read.csv("'."$keggdescfile".'", sep = "\t", header = F)

kegg_enrich <- enricher(targetGene$GID,TERM2GENE = ko2geneID,TERM2NAME = KOdesc[,c(1,2)],pAdjustMethod = "BH",pvalueCutoff  = 0.05)
#write_csv(as.data.frame(kegg_enrich),paste0(targetGeneList,".KEGG_enzyme.csv"))
write.table(as.data.frame(kegg_enrich),file = paste(targetGeneList,".KEGG_enzyme.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


# Enrichment of pathways

pathway2gene <- inner_join(ko2geneID, ko2pathway, by = "Ko")
pathway2gene <- pathway2gene[c("Pathway", "GID")]
colnames(pathway2gene) <- c("Ko", "GID")

kegg_enrich <- enricher(targetGene$GID,TERM2GENE = pathway2gene,TERM2NAME = keg_lvs[,c(1,4)],pAdjustMethod = "BH",pvalueCutoff  = 0.05)
write.table(as.data.frame(kegg_enrich),file = paste(targetGeneList,".KEGG_pathway.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# Loop through each pathway in the enrichment result and generate pathview plots
for (i in 1:nrow(kegg_enrich)) {
  pathway_id <- kegg_enrich$ID[i]
  #print(pathway_id)
  # Generate pathview plot for the pathway
  pathview_plot <- try(pathview(
    gene.data  = targetGene$GID,
    pathway.id = pathway_id,
    species    = "ko"
  ))
}

topNum <- 20
#goEnrichment<-read.table(goEnrichment_result,sep="\t",header=T)
goEnrichment<-as.data.frame(kegg_enrich)
goEnrichment$pvalue <- as.numeric(goEnrichment$pvalue)
goEnrichment <- goEnrichment[goEnrichment$pvalue < 0.05,] 
goEnrichment <- goEnrichment[,c("ID","Description","pvalue")]

ntop <- min(as.numeric(topNum), nrow(goEnrichment)) # If there are less than the specified number of representative GOs, just use all from the table
ggdata <- goEnrichment[1:ntop,]
ggdata$Description <- factor(ggdata$Description, levels = rev(ggdata$Description)) # fixes order
pdf(paste(targetGeneList,".KEGG_pathway_Plot.pdf",sep=""),height=16,width=12)
ggplot(ggdata,
       aes(x = Description, y = -log10(pvalue), size = -log10(pvalue), fill = -log10(pvalue))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = "royalblue", high = "red4") +
  
  xlab("") + ylab("Enrichment score") +
  labs(
    title = "KEGG enrichment, Top",
    subtitle = paste(ntop," terms ordered by p-value"),
    caption = "Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001") +
  
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "dotted", "dotted"),
             colour = c("red", "red", "red"),
             size = c(1, 1, 1)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = "right",
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = "bold", vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = "bold", vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = "bold", vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = "bold", hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = "bold", vjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.line = element_line(colour = "black"),
    
    legend.key = element_blank(), 
    legend.key.size = unit(1, "cm"), 
    legend.text = element_text(size = 14, face = "bold"),
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()

dev.off()

### KEGG enrichment using enrichKEGG
# Using all KEGG universe

targetGene2Kid <- inner_join(targetGene, ko2geneID, by = "GID")
targetGeneKid <- targetGene2Kid$Ko
universeGeneKid <- ko2geneID$Ko


## Using Dmel
# Using Dmel list from KEGG, with general universe
library(org.Dm.eg.db)
# See all available keytypes
#available_keytypes <- columns(org.Dm.eg.db)
#print(available_keytypes)

# Get ENTREZ IDs 
idsdmeltarget<-bitr(dmeltargetGene$GID, fromType = "FLYBASECG", toType = "ENTREZID", OrgDb="org.Dm.eg.db")
#dedup_ids = ids[!duplicated(ids[c("FLYBASECG")]),] # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
idsdmeluniverse<-bitr(dmelpopulation$GID, fromType = "FLYBASECG", toType = "ENTREZID", OrgDb="org.Dm.eg.db")

write.table(idsdmeltarget$ENTREZID, file = paste(targetGeneList,"_KEGG_dmel_entrezid_list.txt",sep=""), col.names = FALSE, row.names = FALSE,  quote = FALSE)


kegg_organism = "dme"
kegg_path_enrich <- enrichKEGG(gene=idsdmeltarget$ENTREZID, universe=idsdmeluniverse$ENTREZID,organism=kegg_organism, pvalueCutoff = 0.1, keyType = "ncbi-geneid")
write.table(as.data.frame(kegg_path_enrich),file = paste(targetGeneList,".KEGG_dmel_pathway.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# Loop through each pathway in the enrichment result and generate pathview plots
for (i in 1:nrow(kegg_path_enrich)) {
  pathway_id <- kegg_path_enrich$ID[i]
  #print(pathway_id)
  # Generate pathview plot for the pathway
  pathview_plot <- try(pathview(
    gene.data  = idsdmeltarget$ENTREZID,
    pathway.id = pathway_id,
    species    = "dme"
  ))
}

kegg_organism = "dme"
kegg_path_enrich <- enrichKEGG(gene=idsdmeltarget$ENTREZID, organism=kegg_organism, pvalueCutoff = 0.1, keyType = "ncbi-geneid")
write.table(as.data.frame(kegg_path_enrich),file = paste(targetGeneList,".KEGG_dmel_allpathway.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# Loop through each pathway in the enrichment result and generate pathview plots
for (i in 1:nrow(kegg_path_enrich)) {
  pathway_id <- kegg_path_enrich$ID[i]
  #print(pathway_id)
  # Generate pathview plot for the pathway
  pathview_plot <- try(pathview(
    gene.data  = idsdmeltarget$ENTREZID,
    pathway.id = pathway_id,
    species    = "dme"
  ))
}

topNum <- 20
#goEnrichment<-read.table(goEnrichment_result,sep="\t",header=T)
goEnrichment<-as.data.frame(kegg_path_enrich)
goEnrichment$pvalue <- as.numeric(goEnrichment$pvalue)
goEnrichment <- goEnrichment[goEnrichment$pvalue < 0.05,] 
goEnrichment <- goEnrichment[,c("ID","Description","pvalue")]

ntop <- min(as.numeric(topNum), nrow(goEnrichment)) # If there are less than the specified number of representative GOs, just use all from the table
ggdata <- goEnrichment[1:ntop,]
ggdata$Description <- factor(ggdata$Description, levels = rev(ggdata$Description)) # fixes order
pdf(paste(targetGeneList,".KEGG_dmel_pathway_Plot.pdf",sep=""),height=16,width=12)
ggplot(ggdata,
       aes(x = Description, y = -log10(pvalue), size = -log10(pvalue), fill = -log10(pvalue))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = "royalblue", high = "red4") +
  
  xlab("") + ylab("Enrichment score") +
  labs(
    title = "KEGG enrichment, Top",
    subtitle = paste(ntop," terms ordered by p-value"),
    caption = "Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001") +
  
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "dotted", "dotted"),
             colour = c("red", "red", "red"),
             size = c(1, 1, 1)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = "right",
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = "bold", vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = "bold", vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = "bold", vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = "bold", hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = "bold", vjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.line = element_line(colour = "black"),
    
    legend.key = element_blank(), 
    legend.key.size = unit(1, "cm"), 
    legend.text = element_text(size = 14, face = "bold"),
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()

dev.off()

# Using Dmel as organism but with kegg ids - nothing it does not associate the kegg id with pathway

#kegg_path_enrich <- enrichKEGG(targetGeneKid, organism="dme", keyType="kegg")
#write.table(as.data.frame(kegg_path_enrich),file = paste(targetGeneList,".KEGG_dmel_bykegg_allpathway.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#kegg_path_enrich <- enrichKEGG(targetGeneKid, organism="dme", keyType="kegg", universe = universeGeneKid )
#write.table(as.data.frame(kegg_path_enrich),file = paste(targetGeneList,".KEGG_dmel_bykegg_pathway.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


# Using Mpha list from KEGG

loc2geneID <- loc2geneID %>% mutate(LOCID = stringr::str_remove(LOCID, "LOC"))
targetGene2loc <- inner_join(targetGene, loc2geneID, by = "GID")
targetGeneloc <- targetGene2loc$LOCID
universeGeneloc <- loc2geneID$LOCID

kegg_path_enrich <- enrichKEGG(targetGeneloc, organism="mpha", keyType="kegg")
#kegg_path_enrich2 <- enrichKEGG(targetGeneloc, organism="mpha", keyType="kegg", universe = universeGeneloc)
write.table(as.data.frame(kegg_path_enrich),file = paste(targetGeneList,".KEGG_allmpha_pathway.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# Loop through each pathway in the enrichment result and generate pathview plots
for (i in 1:nrow(kegg_path_enrich)) {
  pathway_id <- kegg_path_enrich$ID[i]
  #print(pathway_id)
  # Generate pathview plot for the pathway
  pathview_plot <- try(pathview(
    gene.data  = targetGeneloc,
    pathway.id = pathway_id,
    species    = "mpha"
  ))
}

topNum <- 20
#goEnrichment<-read.table(goEnrichment_result,sep="\t",header=T)
goEnrichment<-as.data.frame(kegg_path_enrich)
goEnrichment$pvalue <- as.numeric(goEnrichment$pvalue)
goEnrichment <- goEnrichment[goEnrichment$pvalue < 0.05,] 
goEnrichment <- goEnrichment[,c("ID","Description","pvalue")]

ntop <- min(as.numeric(topNum), nrow(goEnrichment)) # If there are less than the specified number of representative GOs, just use all from the table
ggdata <- goEnrichment[1:ntop,]
ggdata$Description <- factor(ggdata$Description, levels = rev(ggdata$Description)) # fixes order
pdf(paste(targetGeneList,".KEGG_allmpha_pathway_Plot.pdf",sep=""),height=16,width=12)
ggplot(ggdata,
       aes(x = Description, y = -log10(pvalue), size = -log10(pvalue), fill = -log10(pvalue))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = "royalblue", high = "red4") +
  
  xlab("") + ylab("Enrichment score") +
  labs(
    title = "KEGG enrichment, Top",
    subtitle = paste(ntop," terms ordered by p-value"),
    caption = "Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001") +
  
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "dotted", "dotted"),
             colour = c("red", "red", "red"),
             size = c(1, 1, 1)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = "right",
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = "bold", vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = "bold", vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = "bold", vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = "bold", hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = "bold", vjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.line = element_line(colour = "black"),
    
    legend.key = element_blank(), 
    legend.key.size = unit(1, "cm"), 
    legend.text = element_text(size = 14, face = "bold"),
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()

dev.off()

# Using Mpha list from KEGG, with specific universe

kegg_path_enrich <- enrichKEGG(targetGeneloc, organism="mpha", keyType="kegg", universe = universeGeneloc)
write.table(as.data.frame(kegg_path_enrich),file = paste(targetGeneList,".KEGG_mpha_pathway.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# Loop through each pathway in the enrichment result and generate pathview plots
for (i in 1:nrow(kegg_path_enrich)) {
  pathway_id <- kegg_path_enrich$ID[i]
  #print(pathway_id)
  # Generate pathview plot for the pathway
  pathview_plot <- try(pathview(
    gene.data  = targetGeneloc,
    pathway.id = pathway_id,
    species    = "mpha"
  ))
}

topNum <- 20
#goEnrichment<-read.table(goEnrichment_result,sep="\t",header=T)
goEnrichment<-as.data.frame(kegg_path_enrich)
goEnrichment$pvalue <- as.numeric(goEnrichment$pvalue)
goEnrichment <- goEnrichment[goEnrichment$pvalue < 0.05,] 
goEnrichment <- goEnrichment[,c("ID","Description","pvalue")]

ntop <- min(as.numeric(topNum), nrow(goEnrichment)) # If there are less than the specified number of representative GOs, just use all from the table
ggdata <- goEnrichment[1:ntop,]
ggdata$Description <- factor(ggdata$Description, levels = rev(ggdata$Description)) # fixes order
pdf(paste(targetGeneList,".KEGG_mpha_pathway_Plot.pdf",sep=""),height=16,width=12)
ggplot(ggdata,
       aes(x = Description, y = -log10(pvalue), size = -log10(pvalue), fill = -log10(pvalue))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = "royalblue", high = "red4") +
  
  xlab("") + ylab("Enrichment score") +
  labs(
    title = "KEGG enrichment, Top",
    subtitle = paste(ntop," terms ordered by p-value"),
    caption = "Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001") +
  
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "dotted", "dotted"),
             colour = c("red", "red", "red"),
             size = c(1, 1, 1)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = "right",
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = "bold", vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = "bold", vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = "bold", vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = "bold", hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = "bold", vjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.line = element_line(colour = "black"),
    
    legend.key = element_blank(), 
    legend.key.size = unit(1, "cm"), 
    legend.text = element_text(size = 14, face = "bold"),
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()

dev.off()


# Using all KEGG universe

kegg_path_enrich <- enrichKEGG(targetGeneKid, organism="ko", keyType="kegg")
write.table(as.data.frame(kegg_path_enrich),file = paste(targetGeneList,".KEGG_allko_pathway.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# Visualize pathways | see https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
#browseKEGG(kegg_path_enrich, "ko04310") #ko04140   Autophagy - animal
#library("pathview")
##hsa04110 <- pathview(gene.data  = geneList,
##                     pathway.id = "hsa04110",
##                     species    = "hsa",
##                     limit      = list(gene=max(abs(geneList)), cpd=1))
#ko04310 <- pathview(gene.data  = targetGeneKid,
#                    pathway.id = "ko04310",
#                    species    = "ko")

# Loop through each pathway in the enrichment result and generate pathview plots
for (i in 1:nrow(kegg_path_enrich)) {
  pathway_id <- kegg_path_enrich$ID[i]
  #print(pathway_id)
  # Generate pathview plot for the pathway
  pathview_plot <- try(pathview(
    gene.data  = targetGeneKid,
    pathway.id = pathway_id,
    species    = "ko"
  ))
}


# Plot
topNum <- 20
#goEnrichment<-read.table(goEnrichment_result,sep="\t",header=T)
goEnrichment<-as.data.frame(kegg_path_enrich)
goEnrichment$pvalue <- as.numeric(goEnrichment$pvalue)
goEnrichment <- goEnrichment[goEnrichment$pvalue < 0.05,] 
goEnrichment <- goEnrichment[,c("ID","Description","pvalue")]

ntop <- min(as.numeric(topNum), nrow(goEnrichment)) # If there are less than the specified number of representative GOs, just use all from the table
ggdata <- goEnrichment[1:ntop,]
ggdata$Description <- factor(ggdata$Description, levels = rev(ggdata$Description)) # fixes order
pdf(paste(targetGeneList,".KEGG_allko_pathway_Plot.pdf",sep=""),height=16,width=12)
ggplot(ggdata,
       aes(x = Description, y = -log10(pvalue), size = -log10(pvalue), fill = -log10(pvalue))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = "royalblue", high = "red4") +
  
  xlab("") + ylab("Enrichment score") +
  labs(
    title = "KEGG enrichment, Top",
    subtitle = paste(ntop," terms ordered by p-value"),
    caption = "Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001") +
  
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "dotted", "dotted"),
             colour = c("red", "red", "red"),
             size = c(1, 1, 1)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = "right",
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = "bold", vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = "bold", vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = "bold", vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = "bold", hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = "bold", vjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.line = element_line(colour = "black"),
    
    legend.key = element_blank(), 
    legend.key.size = unit(1, "cm"), 
    legend.text = element_text(size = 14, face = "bold"),
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()

dev.off()

# Using only the total set of KEGG universe

kegg_path_enrich <- enrichKEGG(targetGeneKid, organism="ko", keyType="kegg", universe = universeGeneKid)
write.table(as.data.frame(kegg_path_enrich),file = paste(targetGeneList,".KEGG_ko_pathway.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


topNum <- 20
#goEnrichment<-read.table(goEnrichment_result,sep="\t",header=T)
goEnrichment<-as.data.frame(kegg_path_enrich)
goEnrichment$pvalue <- as.numeric(goEnrichment$pvalue)
goEnrichment <- goEnrichment[goEnrichment$pvalue < 0.05,] 
goEnrichment <- goEnrichment[,c("ID","Description","pvalue")]

ntop <- min(as.numeric(topNum), nrow(goEnrichment)) # If there are less than the specified number of representative GOs, just use all from the table
ggdata <- goEnrichment[1:ntop,]
ggdata$Description <- factor(ggdata$Description, levels = rev(ggdata$Description)) # fixes order
pdf(paste(targetGeneList,".KEGG_ko_pathway_Plot.pdf",sep=""),height=16,width=12)
ggplot(ggdata,
       aes(x = Description, y = -log10(pvalue), size = -log10(pvalue), fill = -log10(pvalue))) +
  
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = "royalblue", high = "red4") +
  
  xlab("") + ylab("Enrichment score") +
  labs(
    title = "KEGG enrichment, Top",
    subtitle = paste(ntop," terms ordered by p-value"),
    caption = "Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001") +
  
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "dotted", "dotted"),
             colour = c("red", "red", "red"),
             size = c(1, 1, 1)) +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = "right",
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = "bold", vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = "bold", vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = "bold", vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = "bold", hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = "bold", vjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.line = element_line(colour = "black"),
    
    legend.key = element_blank(), 
    legend.key.size = unit(1, "cm"), 
    legend.text = element_text(size = 14, face = "bold"),
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()

dev.off()

# Loop through each pathway in the enrichment result and generate pathview plots
for (i in 1:nrow(kegg_path_enrich)) {
  pathway_id <- kegg_path_enrich$ID[i]
  #print(pathway_id)
  # Generate pathview plot for the pathway
  pathview_plot <- try(pathview(
    gene.data  = targetGeneKid,
    pathway.id = pathway_id,
    species    = "ko"
  ))
}


';
close ResultsR;

# Run Rscript
#print "Rscript $outname\_KEGGenrichment_script.R\n";
system ("Rscript $outname\_KEGGenrichment_script.R");


## Bitao NEE expression and get tables
#=head
system("mkdir -p $outname\_Mpha_expression_NEE");

open (Resultscan, ">", "$outname\_NEE_caste_bias_canalized.tsv");
print Resultscan "$bitaotablecanalizedheader\n";

open (Resultsde, ">", "$outname\_NEE_caste_DE.tsv");
print Resultsde "$bitaotabledeheader\n";

open (Resultsdempha, ">", "$outname\_NEE_caste_DE_Mpha.tsv");
print Resultsdempha "$bitaotabledemphaheader\n";

open(Filef, "<", $fannotfile);
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    next if ($line =~ /Number of seqs with that Swissprot/);
    my @subl = split (/\t/, $line);
    my $og = "";
    if ($line =~ /(HOG\d\d\d\d\d\d\d)/){
        $og = $1;
    } else {
        die "Cannot find OG in line 1054, file: $fannotfile and line $line\n";
    }

    if ($candidategolist =~ / $og /){ # HOG in candidate list, save. 
        my $mphancbi = ""; # Getting the Mpha ncbi ids to add them in the output table
        if (exists ($subl[21])){
            my $mphaids = $subl[21];
            my @sepmphaids = split(/\,/, $mphaids);
            foreach my $sepid (@sepmphaids){
                $sepid =~ s/\s//g;
                if (exists ($mphatable{$sepid})){
 
                    system ("cp $mphafolder\/$mphatable{$sepid}\.pdf $outname\_Mpha_expression_NEE/$og\_$sepid\_$mphatable{$sepid}\.pdf");

                    if (exists $bitaotablecanalizedfull{$mphatable{$sepid}}){
                        print Resultscan "$bitaotablecanalizedfull{$mphatable{$sepid}}\n";
                    }
                    if (exists $bitaotabledefull{$mphatable{$sepid}}){
                        print Resultsde "$bitaotabledefull{$mphatable{$sepid}}\n";
                    }
                    if (exists $bitaotabledemphafull{$mphatable{$sepid}}){
                        print Resultsdempha "$bitaotabledemphafull{$mphatable{$sepid}}\n";
                    }

                }
            }
        }


    } else {
        next; # HOG not in candidate
    }
}
close Filef;
close Resultscan;
close Resultsde;
close Resultsdempha;
#=cut

# Ploting the GAGA data
#=h
system("perl $scriptspath\/get_expression_plots.pl $outfile");
#system("perl $scriptspath\/get_expression_plots.pl $outfile");
#=cut

# Plotting the GAGA data separated by groups
=head
my $specieslistfile = ""; # Define here the species comparison file used in the trait selection analyses
system("perl $scriptspath\/get_expression_plots_traitcomps.pl $outfile $specieslistfile");
system("perl $scriptspath\/get_expression_plots_traitcomps.pl $outfile $specieslistfile");
=cut

# jokergoo simplifyenrichment - MOVE AT THE END THE PVAL005, it takes very long
#=h
system ("Rscript $scriptspath\/simplifyEnrichment_script.R $outname\_GOenrichment_BP_pval001.txt");
system ("Rscript $scriptspath\/simplifyEnrichment_script.R $outname\_GOenrichment_BP_pval005.txt");
#system ("Rscript $scriptspath\/simplifyEnrichment_script.R $outfile\.BP_classic_pval001.txt");
system ("Rscript $scriptspath\/simplifyEnrichment_script.R $outfile\.BP_classic_pval005.txt");
#system ("Rscript $scriptspath\/simplifyEnrichment_script.R $outfile\.BP_elim_pval005.txt");
#=cut

## Organize all files into folders
system("mv $outname\_KEGG_enrichment_v3 $outname\_KEGG_enrichment_v3old 2>/dev/null"); # Rename if there are previous folders
system("mkdir -p $outname\_KEGG_enrichment_v3");
system("mv $outname*KEGG* $outname*pulation_Mphaids.annot $outname\_KEGG_enrichment_v3/ 2>/dev/null");
system("mkdir -p $outname\_KEGG_enrichment_v3/pathway_plots");
system("mv *.png *.xml $outname\_KEGG_enrichment_v3/pathway_plots/ 2>/dev/null");


# Remove population files that are quite big in size
#system ("rm $outname\_population_*");

## Organize all files into folders
system("rm -rf $outname\_GO_enrichment_v3 2>/dev/null"); # Remove if there is a previous folder
system("mkdir -p $outname\_GO_enrichment_v3");
system("rm tGO.* 2>/dev/null");
system("mv $outname*GOenrich* *tGO* $outname\_GOfigure* $outname*classic* $outname*weight01* $outname*elim* $outname\_GO_enrichment_v3/ 2>/dev/null");



# rerun some steps manually or clean folder interrupted as simplifyenrichment gets stacked
=head
rm Rplots.pdf tGO* *population*
PREFIX=Absrel_candidates_reference_enrichedPS
mkdir -p "$PREFIX"_GO_enrichment_v3
mv "$PREFIX"_GOenrichment* "$PREFIX"_nopart_ids.txt.* "$PREFIX"_GO_enrichment_v3
rm -rf "$PREFIX"_GOfigure*

PREFIX=Relax_candidates_relax_fdr001_table_relaxedtest
PREFIX=Wings_hogtable
PREFIX=All_candidates_intensified_PS_queenless_gamergate_ergatoid.tsv
----

FILER=Absrel_candidates_reference_enrichedPS_nopart_ids.txt.BP_classic_pval001.txt
FILER=All_candidates_relaxed_PS_qwdimorphism_polymorphism_colonysize.tsv_nopart_ids.txt.BP_elim_pval005.txt
Rscript /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/simplifyEnrichment_script_onlyreduced.R $FILER > xsimplifyenrichment_onlyreduced.out 2> xsimplifyenrichment_onlyreduced.err &

Rscript /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/simplifyEnrichment_script_onlydmel.R $FILER > xsimplifyenrichment_onlyfly.out 2> xsimplifyenrichment_onlyfly.err &

Rscript /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/simplifyEnrichment_script_onlydefault.R $FILER > xsimplifyenrichment_onlydef.out 2> xsimplifyenrichment_onlydef.err &

Rscript /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/simplifyEnrichment_script_onlydmel_Rupdate.R $FILER > xsimplifyenrichment_onlyfly.out 2> xsimplifyenrichment_onlyfly.err &




FILER=Absrel_candidates_reference_enrichedPS_GOenrichment_BP_pval001.txt
FILER=All_candidates_intensified_PS_queenless_gamergate_ergatoid.tsv_GOenrichment_BP_pval001.txt

Rscript /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/simplifyEnrichment_script_onlyreduced.R $FILER > xsimplifyenrichment_onlyreducedpv1.out 2> xsimplifyenrichment_onlyreducedpv1.err &

#Rscript /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/simplifyEnrichment_script_onlydmel.R $FILER > xsimplifyenrichment_onlyflypv1.out 2> xsimplifyenrichment_onlyflypv1.err &

Rscript /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/simplifyEnrichment_script_onlydefault.R $FILER > xsimplifyenrichment_onlydefpv1.out 2> xsimplifyenrichment_onlydefpv1.err &

Rscript /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/simplifyEnrichment_script_onlydmel_Rupdate.R $FILER > xsimplifyenrichment_onlyflydefpv1.out 2> xsimplifyenrichment_onlyflydefpv1.err &


FILER=Absrel_candidates_reference_enrichedPS_GOenrichment_BP_pval005.txt
FILER=All_candidates_intensified_PS_queenless_gamergate_ergatoid.tsv_GOenrichment_BP_pval005.txt

Rscript /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/simplifyEnrichment_script_onlyreduced.R $FILER > xsimplifyenrichment_onlyreducedpv5.out 2> xsimplifyenrichment_onlyreducedpv5.err &

#Rscript /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/simplifyEnrichment_script_onlydmel.R $FILER > xsimplifyenrichment_onlyflypv5.out 2> xsimplifyenrichment_onlyflypv5.err &

Rscript /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/simplifyEnrichment_script_onlydefault.R $FILER > xsimplifyenrichment_onlydefpv5.out 2> xsimplifyenrichment_onlydefpv5.err &

Rscript /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/simplifyEnrichment_script_onlydmel_Rupdate.R $FILER > xsimplifyenrichment_onlyflydefpv5.out 2> xsimplifyenrichment_onlyflydefpv5.err &


=cut



