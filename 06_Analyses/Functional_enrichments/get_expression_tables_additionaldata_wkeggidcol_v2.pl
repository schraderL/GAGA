#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# Script to add additional columns in the functional annotation table with the expression patterns in ovary, mpha gyne tissues..
# v2 includes updated information from Mpha developmental transcriptomes across caste (Bitao 2022)

# usage: perl get_expression_tables_additionaldata.pl file_funcannot.tsv 


my $clist = "$ARGV[0]"; # Candidate list file functional annotation .tsv


# Output prefix name
my $outname = "";
if ($clist =~ /(\S+)\.tsv/){
    $outname = $1;
} else {
    $outname = $clist;
}

# Needed files for the script
#my $goannotfile = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_functional_annotation_wKEGG_v2/Orthogroups/Orthogroups_functional_annotation_GO.annot";
#my $fannotfile = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_functional_annotation_wKEGG_v2/Orthogroups/Orthogroups_functional_annotation_summary.tsv";
#my $keggfile = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_functional_annotation_wKEGG_v2/Orthogroups/Orthogroups_functional_annotation_KEGG.annot";
#my $keggRdata = "/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/KEGG_levs.pathway.module.Rdata";
#my $obogofile = "/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/GO_19Apr2023.obo"; # Read ovo file with GO descriptions to include in output table

my $mphafilezijun = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/Mpha_gene_association_GAGA2NCBI.txt";
my $mphafile = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/aMpha_gaga_ncbi_correspondance_table.tsv";

#my $mphafolder = "/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/Bitao_NEE/pharaoh_developmantal_exp_plot_v2/";
my $bitaocanalized = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/Bitao_NEE/NEE_supp_file_caste_bias_canalized_expression.tsv";
my $bitaode = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/Bitao_NEE/NEE_supp_file_caste_differential_expression.tsv";
my $bitaodemphadev = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/Bitao_NEE/NEE_supp_file_caste_differential_expression_developmentMpha.tsv";
my $bitaodempha = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/Bitao_NEE/NEE_supp_file_caste_differential_expression_larvaeMpha.tsv";
my $dmelfile = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/Dmel_FlyAtlas2_gene_data_Mpha.txt";

my $mphabraindefile = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/GAGA_expression/Bitao_brain/DEG.LOC2Mpha.combined_v2.csv";
my $mphabraindeantsfile = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/GAGA_expression/Bitao_brain/NEE2018qiu_DEG_genes.csv";

my $mphagynefile = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/GAGA_expression/Gyne_tissues_china/abundance.addPvalue3_gyne_tissues_DE_withinfo.tsv";

my $mphaovaryfile = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/GAGA_expression/Ovary_Yue/ov_DEG.res.05.lfc0_5d_Dmel.txt";

my $mphaopticlobefile = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/GAGA_expression/Optic_lobe/Opticlobe_genes.txt";
my $mphaolidfile = "/Users/joel/Dropbox/Joel/trabajo/GAGAproject/Analyses_functional_tables_enrich/GAGA_expression/Optic_lobe/Pharaoh_ant_id_correspondance_GAGA.txt";

#my $gofigurescript = "/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Analisis/gofigure/GO-Figure/gofigure.py"; # GO figure script path

# Start script


# Read Mpha association table from zijun
my %mphatablezijun; 
open(Filef, "<", $mphafilezijun);
while(<Filef>){
    chomp;
    my $line = $_;
    $line =~ s/\r//g; # Replace all shit
    $line =~ s/\n//g; # Replace all shit
    next if ($line !~ /\S+/ );
    my @subl = split (/\t/, $line);
    $mphatablezijun{$subl[0]} = $subl[2];
}
close Filef;

# Read Mpha association table from orthofinder
my %mphatable; my $mphatableheader = ""; my $mphatablehadercount = "0";
open(File, "<", "$mphafile");
while(<File>){
    chomp;
    my $line = $_;
    $line =~ s/\r//g; # Replace all shit
    $line =~ s/\n//g; # Replace all shit
    next if ($line !~ /\S+/);
    my @subl = split (/\t/, $line);

    if ($mphatablehadercount == 0){
        # Header
        $mphatableheader = $line;
        $mphatablehadercount++;
        next;
    }

    #my $gagaid = "$subl[4]";
    my $ncbiid = "$subl[1]"; 

    if ($subl[4] =~ /\,/){ # There is more than one gene
        my @subgagaid = split(/,/, $subl[4]);
        foreach my $gagaid (@subgagaid){
            $gagaid =~ s/\s//g; # Replace any space
             if (exists $mphatable{$gagaid}){
                if ($mphatable{$gagaid} =~ /NA/){ # Replace with new entry
                    $mphatable{$gagaid} = $ncbiid;
                } elsif ($gagaid =~ /NA/){ # Skip this entry with NA
                    next;
                } elsif ($mphatable{$gagaid} =~ /$ncbiid/){ # Same entry
                    next;
                } else {
                    #print "Different Mpha ids for gene $subl[1] protein $subl[0]\n";
                }

            } else { # New entry
                $mphatable{$gagaid} = $ncbiid;
            }            
        } 
    } else {
        my $gagaid = "$subl[4]";
        if (exists $mphatable{$gagaid}){
            if ($mphatable{$gagaid} =~ /NA/){ # Replace with new entry
                $mphatable{$gagaid} = $ncbiid;
            } elsif ($gagaid =~ /NA/){ # Skip this entry with NA
                next;
            } elsif ($mphatable{$gagaid} =~ /$ncbiid/){ # Same entry
                next;
            } else {
                #print "Different Mpha ids for gene $subl[1] protein $subl[0]\n";
            }

        } else { # New entry
            $mphatable{$gagaid} = $ncbiid;
        }        
    } 

}
close File;

#=h
# Read Bitao candidates - canalized
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

# Read Bitao candidates - DE across caste in mpha and aech
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
    my $logfc = "$subl[2]";
    my $castebias = "";
    if ($logfc =~ /\-/){
        $castebias = "Worker";
    } else {
        $castebias = "Gyne";
    } 

    if (exists $bitaotablede{$subl[0]}){ # Add the new DE stage into the same gene
        $bitaotablede{$subl[0]} .= ", $castebias $subl[1]";
        $bitaotabledefull{$subl[0]} = "$subl[0]\t$bitaotablede{$subl[0]}\t$subl[2]\t$subl[3]\t$subl[4]\t$subl[5]\t$subl[6]\t$subl[7]\t$subl[8]\n";        
    } else {
        $bitaotablede{$subl[0]} = "$castebias $subl[1]";
        $bitaotabledefull{$subl[0]} = $line;
    }
}
close Filef;

# Read Bitao candidates - DE across Mpha developmental stages
my %bitaotablededev; my $bitaotablededevheader = "NA"; my %bitaotablededevfull;
open(Filef, "<", $bitaodemphadev);
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    if ($line =~ /Caste DE/){ # Header
        $bitaotablededevheader = $line;
        next;
    }
    my @subl = split (/\t/, $line);
    my $logfc = "$subl[4]";
    my $castebias = "";
    if ($logfc =~ /\-/){
        $castebias = "Worker";
    } else {
        $castebias = "Gyne";
    } 

    if (exists $bitaotablededev{$subl[0]}){ # Add the new DE stage into the same gene
        $bitaotablededev{$subl[0]} .= ", $castebias $subl[1]";
        $bitaotablededevfull{$subl[0]} = "$subl[0]\t$bitaotablededev{$subl[0]}\t$subl[2]\t$subl[3]\t$subl[4]\t$subl[5]\t$subl[6]\t$subl[7]\t$subl[8]\n"; 
    } else {
        $bitaotablededev{$subl[0]} = "$castebias $subl[1]";
        $bitaotablededevfull{$subl[0]} = $line;
    }
}
close Filef;

# Read Bitao candidates - DE across Mpha larval caste transcriptomes
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
    my $logfc = "$subl[4]";
    my $castebias = "";
    if ($logfc =~ /\-/){
        $castebias = "Worker";
    } else {
        $castebias = "Gyne";
    } 

    if (exists $bitaotabledempha{$subl[0]}){
        next; # Already saved, might be repeated those DE in Mpha and Aech, that appear later only in Mpha
    } else {
        $bitaotabledempha{$subl[0]} = "$castebias $subl[1]";
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
        $dmeltable{$subl[1]} = "$subl[5] $subl[3]";
        $dmelexprtable{$subl[1]} = "$subl[6]";
        $dmeltablefull{$subl[1]} = $line;
    }
}
close Filef;
#=cut

# Read monomorium brain expression - 2 files | saved by gaga id

my %mphabraindetable; my $mphabraindetableheader = ""; my $mphabraindetablehadercount = "0"; my %mphabraindetablefull;
my $mphabraindetableheadersaved = "Mpha brain expression\tMpha brain DE worker-gyne\tMpha brain Fold change\tMpha worker brain TPM\tMpha gyne brain TPM";
open(File, "<", "$mphabraindefile");
while(<File>){
    chomp;
    my $line = $_;
    $line =~ s/\r//g; # Replace all shit
    $line =~ s/\n//g; # Replace all shit
    next if ($line !~ /\S+/);
    my @subl = split (/\t/, $line);

    if ($mphabraindetablehadercount == 0){
        # Header
        $mphabraindetableheader = $line;
        $mphabraindetablehadercount++;
        next;
    }

    $mphabraindetablefull{$subl[1]} = $line;

    my $fdr = $subl[16];
    my $foldchange = $subl[14];
    my $workertpm = $subl[12];
    my $queentpm = $subl[13];

    my $totaltpm = $workertpm + $queentpm;
    my $brainexpr = " ";

    if ($totaltpm > 250){
        $brainexpr = "Brain highly expressed TPM: $totaltpm";
    }   

    if ($fdr =~ /NA/){
        $mphabraindetable{$subl[1]} = "$brainexpr\t\t$foldchange\t$workertpm\t$queentpm"; # Saving 5 colums for brain transcriptome
    } elsif ($fdr < 0.01){
        if ($foldchange > 0 ){
            $mphabraindetable{$subl[1]} = "$brainexpr\tWorker brain\t$foldchange\t$workertpm\t$queentpm"; # Saving 5 colums for brain transcriptome
        } else {
            $mphabraindetable{$subl[1]} = "$brainexpr\tGyne brain\t$foldchange\t$workertpm\t$queentpm"; # Saving 5 colums for brain transcriptome
        } 
    } else {
        $mphabraindetable{$subl[1]} = "$brainexpr\t\t$foldchange\t$workertpm\t$queentpm"; # Saving 5 colums for brain transcriptome
    }
 
}
close File;

my %mphabraindeantstable; my $mphabraindeantstableheader = ""; my $mphabraindeantstablehadercount = "0"; my %mphabraindeantstablefull;
my $mphabraindeantstableheadersaved = "Mpha brain DE ants Bitao";
open(File, "<", "$mphabraindeantsfile");
while(<File>){
    chomp;
    my $line = $_;
    $line =~ s/\r//g; # Replace all shit
    $line =~ s/\n//g; # Replace all shit
    next if ($line !~ /\S+/);
    my @subl = split (/\,/, $line);

    if ($mphabraindeantstablehadercount == 0){
        # Header
        $mphabraindeantstableheader = $line;
        $mphabraindeantstablehadercount++;
        next;
    }

    $mphabraindeantstable{$subl[0]} = "Brain DE ants";
 
}
close File;


# Read monomorium gyne tissues | Expression is per gene - so consider it for annotated with gMphaOR, or Mpha_g0001 no iso

my %mphagynetable; my $mphagynetableheader = ""; my $mphagynetablehadercount = "0"; my %mphagynetablefull;
my $mphagynetableheadersaved = "Mpha gyne tissue with differential expression\tMpha gyne tissue most expressed\tMpha gyne tissue most expressed TPM>100\tMpha gyne Brain TPM\tMpha gyne Fatbody TPM\tMpha gyne Ovary TPM\tMpha gyne Thorax TPM\tMpha gyne Intestine TPM";
open(File, "<", "$mphagynefile");
while(<File>){
    chomp;
    my $line = $_;
    $line =~ s/\r//g; # Replace all shit
    $line =~ s/\n//g; # Replace all shit
    next if ($line !~ /\S+/);
    my @subl = split (/\t/, $line);

    if ($mphagynetablehadercount == 0){
        # Header
        $mphagynetableheader = $line;
        $mphagynetablehadercount++;
        next;
    }

    $mphagynetablefull{$subl[0]} = $line;

    my $tissuede = $subl[73];
    my $maintissue = $subl[68];
    my $maintissuetpm = $subl[69];
    my $braintpm = $subl[63];
    my $fatbodytpm = $subl[64];
    my $intestinetpm = $subl[65];
    my $ovarytpm = $subl[66];
    my $thoraxtpm = $subl[67];

    $mphagynetable{$subl[0]} = "$tissuede\t$maintissue\t$maintissuetpm\t$braintpm\t$fatbodytpm\t$ovarytpm\t$thoraxtpm\t$intestinetpm"; # Saving 3 colums for ovary transcriptome

}
close File;



# Read monomorium ovary queen Vs gyne | saved by ncbi id - LOC

my %mphaovarytable; my $mphaovarytableheader = ""; my $mphaovarytablehadercount = "0"; my %mphaovarytablefull;
my $mphaovarytableheadersaved = "Mpha ovary expression gyne Vs queen\tMpha ovary DE Fold Change\tMpha ovary DE FDR";
open(File, "<", "$mphaovaryfile");
while(<File>){
    chomp;
    my $line = $_;
    $line =~ s/\r//g; # Replace all shit
    $line =~ s/\n//g; # Replace all shit
    next if ($line !~ /\S+/);
    my @subl = split (/\t/, $line);

    if ($mphaovarytablehadercount == 0){
        # Header
        $mphaovarytableheader = $line;
        $mphaovarytablehadercount++;
        next;
    }

    $mphaovarytablefull{$subl[0]} = $line;

    my $fdr = $subl[3];
    my $foldchange = $subl[1];
    my $bias = $subl[6];

    if ($fdr < 0.01){
        if ($foldchange > 1 || $foldchange < -1 ){
            $mphaovarytable{$subl[0]} = "Ovary $bias\t$foldchange\t$fdr"; # Saving 3 colums for ovary transcriptome
        } else {
            $mphaovarytable{$subl[0]} = "Ovary $bias low FC\t$foldchange\t$fdr"; # # Saving 3 colums for ovary transcriptome
        } 
    } else {
        $mphaovarytable{$subl[1]} = "\t$foldchange\t$fdr"; # Saving 3 colums for ovary transcriptome
    }

}
close File;



# Read monomorium optic lobe | saved by gaga id

my %mphaolidtable; my $mphaolidtableheader = ""; my $mphaolidtablehadercount = "0";  
open(Filef, "<", $mphaolidfile);
while(<Filef>){
    chomp;
    my $line = $_;
    $line =~ s/\r//g; # Replace all shit
    $line =~ s/\n//g; # Replace all shit
    next if ($line !~ /\S+/ );
    my @subl = split (/\t/, $line);

    if ($mphaolidtablehadercount == 0){
        # Header
        $mphaolidtableheader = $line;
        $mphaolidtablehadercount++;
        next;
    }

    $mphaolidtable{$subl[0]} = $subl[1];
}
close Filef;

my %mphaoltable; my $mphaoltableheader = ""; my $mphaoltablehadercount = "0"; my %mphaoltablefull;
open(Filef, "<", $mphaopticlobefile);
while(<Filef>){
    chomp;
    my $line = $_;
    $line =~ s/\r//g; # Replace all shit
    $line =~ s/\n//g; # Replace all shit
    next if ($line !~ /\S+/ );
    my @subl = split (/\t/, $line);

    if ($mphaoltablehadercount == 0){
        # Header
        $mphaoltableheader = $line;
        $mphaoltablehadercount++;
        next;
    }

    my $mphaid = "$subl[0]";
    if (exists $mphaolidtable{$subl[0]}){
         $mphaoltable{$mphaolidtable{$subl[0]}} = "$subl[2]";
         $mphaoltablefull{$mphaolidtable{$subl[0]}} = "$line";
    }

}
close Filef;




#########

# Create candidate table with gene functions and expression
open (Results, ">", "$outname\_add_expression.tsv");
print Results "Orthogroup\tNumber of sequences\tNumber of species\tSwissprot\tNumber of seqs with that Swissprot\tDrosophila melanogaster swissprotbest hit\tNumber of seqs with that Dmel seq as best hit\tEggNog\tNumber of seqs with that EggNog\tTrEMBL\tNumber of seqs with that TrEML\tInterPro\tNumber of seqs with that InterPro\tPfam\tNumber of seqs with that Pfam\tSignalP\tNumber of seqs with SignalP detected\tNumber of GOs (at least 33\% sequences)\tTotal number of GOs\tNumber of KEGGs (at least 33\% sequences)\tTotal number of KEGGs\tKEGG ids\tMonomorium pharaonis gene names\tMonomorium pharaonis ncbi gene names from Zijun table\tMonomorium pharaonis ncbi gene names newversion\tDrosophila melanogaster best repiprocal hit\tExpression in Dmel top5 tissues\tCaste bias canalized expression NEE\tCaste differential expression Mpha & Aech NEE\tCaste differential expression Mpha development NEE\tCaste differential expression Mpha larvae NEE\t$mphabraindetableheadersaved\tMpha brain DE across ants\t$mphagynetableheadersaved\t$mphaovarytableheadersaved\tMpha Optic lobe\tGOs\n";

open(Filef, "<", $clist);
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
        die "Cannot find OG in line 354, file: $clist line $line\n";
    }

        my $mphancbi = ""; # Getting the Mpha ncbi ids to add them in the output table

        my $brainde = "\t\t\t\t";
        my $braindeants = "";
        my $gynetissues = "\t\t\t\t\t\t\t";
        my $ovaryde = "\t\t"; # 3 columns within
        my $opticlobe = "";

        #my $dmelhit = ""; # Getting Dmel best reciprocal hit
        #my $dmelexpr = "";
        my $neecan = ""; # NEE candidates from canalized expression
        my $neede = "";
        my $neededev = "";
        my $needempha = ""; 
        if (exists ($subl[21])){
            my $mphaids = $subl[21];
            my @sepmphaids = split(/\,/, $mphaids);
            my $mphaidscount = @sepmphaids;
            foreach my $sepid (@sepmphaids){
                $sepid =~ s/\s//g;

                if (exists $mphaoltable{$sepid}){
                    $opticlobe = $mphaoltable{$sepid};
                }

                #if (exists $mphabraindetable{$sepid}){ # Not here, it is saved with gene not mrna id
                #    $brainde = $mphabraindetable{$sepid};
                #}
                if (exists $mphabraindeantstable{$sepid}){
                    $braindeants = $mphabraindeantstable{$sepid};
                }

                if (exists ($mphatable{$sepid})){ # Get id in ncbi
                    if ($line =~ /MphaOR/ || $line =~ /MphaGR/ || $line =~ /Cytochrome_P450_Mpha/ || $mphaidscount > 4){ # For ORs, GRs, CP450, or orthogroups with more than 5 - use Zijun assignment as it is more accurate
                        if (exists ($mphatablezijun{$sepid})){
                            $mphancbi = "$mphatablezijun{$sepid}";
                        } else {
                            $mphancbi = "$mphatable{$sepid}";
                        }
                    } else {
                        $mphancbi = "$mphatable{$sepid}";
                    }
                    #$mphancbi .= "$mphatable{$sepid} ";

                    if (exists $bitaotablecanalized{$mphancbi}){
                        $neecan = "$bitaotablecanalized{$mphancbi} ";
                    }
                    if (exists $bitaotablede{$mphancbi}){
                        $neede = "$bitaotablede{$mphancbi} ";
                    }
                    if (exists $bitaotablededev{$mphancbi}){
                        $neededev = "$bitaotablededev{$mphancbi} ";
                    }
                    if (exists $bitaotabledempha{$mphancbi}){
                        $needempha = "$bitaotabledempha{$mphancbi} ";
                    }

                    if (exists $mphaovarytable{$mphancbi}){
                        $ovaryde = $mphaovarytable{$mphancbi};
                    }

                }

                #if (exists $dmeltable{$sepid}){
                #    $dmelhit .= "$dmeltable{$sepid} ";
                #    $dmelexpr .= "$dmelexprtable{$sepid} ";
                #}

                my $genename = "";
                if ($sepid =~ /(Mpha\S+)\_i/){ # Get gene id
                    $genename = $1;
                    if (exists $mphagynetable{$genename}){
                        $gynetissues = $mphagynetable{$genename}
                    }
                } elsif ($sepid =~ /(\S+)/){
                    $genename = "g$1"; # Chemosensory with gMphaOBP1
                    if (exists $mphagynetable{$genename}){
                        $gynetissues = $mphagynetable{$genename}
                    }
                }

                if ($sepid =~ /(Mpha\S+)\_i/){ # Get gene id
                    $genename = $1;
                    if (exists $mphabraindetable{$genename}){
                        $brainde = $mphabraindetable{$genename}
                    }
                } elsif ($sepid =~ /(\S+)/){
                    $genename = "$1"; # Chemosensory with gMphaOBP1, but here it does not have the initial g!
                    if (exists $mphabraindetable{$genename}){
                        $brainde = $mphabraindetable{$genename}
                    }
                }                

            }
        }


    #print Results "$line\t$mphancbi\t$dmelhit\t$dmelexpr\t$neecan\t$neede\t$needempha\tNA\tNA\n";  
    print Results "$subl[0]\t$subl[1]\t$subl[2]\t$subl[3]\t$subl[4]\t$subl[5]\t$subl[6]\t$subl[7]\t$subl[8]\t$subl[9]\t$subl[10]\t$subl[11]\t$subl[12]\t$subl[13]\t$subl[14]\t$subl[15]\t$subl[16]\t$subl[17]\t$subl[18]\t$subl[19]\t$subl[20]\t$subl[22]\t$subl[21]\t$subl[23]\t";
#    print Results "$mphancbi\t$subl[24]\t$subl[25]\t$subl[26]\t$subl[27]\t$subl[28]\t";
    print Results "$mphancbi\t$subl[24]\t$subl[25]\t$neecan\t$neede\t$neededev\t$needempha\t"; # printing updated Bitao NEE DE
    print Results "$brainde\t$braindeants\t$gynetissues\t$ovaryde\t$opticlobe\t";
    print Results "$subl[29]\n";

}
close Filef;
close Results;




