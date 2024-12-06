#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# Script to generate the gene expression plots from queen, worker and male transcriptome data in GAGA

# usage: perl get_expression_plots.pl Absrel_candidates_mergepart_single_test.txt[Candidate list of orthogroups]

# Load modules: module load ngs tools gcc/7.4.0 intel/perflibs R/4.2.0 anaconda3/2023.03

my $clist = "$ARGV[0]"; # Candidate list file
#my $headerclist = "F"; # Input without header, otherwise uncomment this line, set to T, and uncomment lines 33-40 | Add if the candidate list contains a first-line header (T) or not (F)

# Output prefix name
my $outname = "";
if ($clist =~ /(\S+)\.txt/){
    $outname = $1;
} else {
    $outname = $clist;
}

# Needed files for the script
my $expressionfile = "/Users/joel/Documents/Results/Analyses_traits/Analyses_folder/GAGA_expression/all_wb_abun_tpm_addTree.txt";
my $samplefile = "/Users/joel/Documents/Results/Analyses_traits/Analyses_folder/GAGA_expression/all_samp_sp_beforeDE_newID_group.txt";

# Start script

# Step already done when doing the functional enrichments, using here that input directly || Remove partitions from candidate list, and just get the HOG
my $outfile = "$clist";
#my $outfile = "$outname\_nopart_ids.txt";
#if ($headerclist =~ /T/){
#    system ("tail -n+2 $clist | awk '{print \$1}' | sort | uniq > $outfile");
#} else {
#    system ("cat $clist | awk '{print \$1}' | sort | uniq > $outfile"); # No header in candidate list
#}

# Create output folders
system("mkdir -p $outname\_GAGA_expression_caste");
system("mkdir -p $outname\_GAGA_expression_caste_subgroups");


# Read list of orthogroups with expression (single copy)

my $expressionhoglist = "";
open(Filef, "<", $expressionfile);
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    next if ($line =~ /worker/ || $line =~ /group/ || $line =~ /cluster/); # Headers for the sample table
#    if ($headerflist =~ /T/){ # Skip first line if header = T in the reference set of genes
#        if ($headercount == 0){
#            $headercount++;
#            next;
#        }
#    }
    #my @subl = split (/\t/, $line);
    my $og = "";
    if ($line =~ /(HOG\d\d\d\d\d\d\d)/){
        $og = $1;
    } else {
        die "Cannot find HOG, script get_expression_plots.pl line 64 in $line\n";
    }
    if ($expressionhoglist =~ / $og /){
        next; # OG already saved
    } else {
        $expressionhoglist .= " $og ";
    }
}
close Filef;


# Start R script
open (ResultsR, ">", "$outname\_GAGA_expression_script.R");
print ResultsR 'library(ggplot2)
library(preprocessCore)
library(ggplot2)
library(RColorBrewer)
#library(hrbrthemes)
library(ggsci)
library(dplyr)
library(reshape2)

# Read files
all_tpm <- read.table("'."$expressionfile".'", header = TRUE, row.names = 1, check.names = FALSE)
colnames(all_tpm) <- gsub("-", "_", colnames(all_tpm)) #replace "-" with "_"
all_tpm <- all_tpm[-c(1,2,3),]
all_tpm <- replace(all_tpm, all_tpm == "na", NA)
all_tpm_bk <- all_tpm


all_tpm <- apply(all_tpm, 2, function(x) as.numeric(x))
all_tpm_norm <- normalize.quantiles(as.matrix(all_tpm))
colnames(all_tpm_norm) <- colnames(all_tpm)
rownames(all_tpm_norm) <- gsub("-", "_", rownames(all_tpm_bk))
t1 <- apply(all_tpm_norm, 2, function(x) quantile(x, probs = seq(0, 1, 1/4), na.rm = TRUE)) #To check the quantiles
df_all_tpm_norm <- as.data.frame(all_tpm_norm)

c_tissue <- c("brain", "antenna", "pupa", "gaster", "legs", "larva", "pupa", "head", "thorax", "abdomen")

# Read sample list
sample_newid_wb <- read.table("'."$samplefile".'", header = FALSE, col.names = c("species",  "sample", "group",  "newID", "caste")) 

sample_newid <- sample_newid_wb
sample_newid <- sample_newid[which(sample_newid_wb$caste != "unknown"), ] # remove unknown castes but tissue samples are still included.
sample_newid <- sample_newid[-grep(paste(c_tissue, collapse = "|"), sample_newid$sample),] # remove tissues
sample_newid <- sample_newid[which(sample_newid$species != "GAGA-0350"), ] # There 7 samples from GAGA-0350, which have wrong sample info.

subgroup = F
if (subgroup == T){
  sample_newid$caste[which(sample_newid$caste == "soldier")] <- "Large_worker" 
  sample_newid$caste[which(sample_newid$caste == "small_worker")] <- "Worker"
  sample_newid$caste[which(sample_newid$caste == "worker")] <- "Worker"
  sample_newid$caste[which(sample_newid$caste == "large_worker")] <- "Large_worker"
  sample_newid$caste[which(sample_newid$caste == "queen")] <- "Queen"
  sample_newid$caste[which(sample_newid$caste == "gyne")] <- "Gyne"
  sample_newid$caste[which(sample_newid$caste == "male")] <- "Male"
  mylevel = c("Queen", "Gyne", "Worker", "Large_worker", "Male")
}else{
  sample_newid$caste[which(sample_newid$caste == "gyne")] <- "Queen"
  sample_newid$caste[which(sample_newid$caste == "queen")] <- "Queen"
  sample_newid$caste[which(sample_newid$caste == "worker")] <- "Worker" 
  sample_newid$caste[which(sample_newid$caste == "soldier")] <- "Worker" 
  sample_newid$caste[which(sample_newid$caste == "small_worker")] <- "Worker"
  sample_newid$caste[which(sample_newid$caste == "large_worker")] <- "Worker"
  sample_newid$caste[which(sample_newid$caste == "male")] <- "Male"
  mylevel = c("Queen", "Worker", "Male")
}
table(sample_newid$caste)

sample_newid$sample <- gsub("-", "_", sample_newid$sample) #replace "-" with "_"

';
#close ResultsR;

# Read candidate list, and add them into the script
open(Filef, "<", $clist);
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    #next if ($line =~ /Branches under selection/ || $line =~ /relaxation/);
    #my @subl = split (/\t/, $line);
    my $og = "";
    if ($line =~ /(HOG\d\d\d\d\d\d\d)/){
        $og = $1;
    } else {
        die "Cannot find HOG, script get_expression_plots.pl line 150 in $line\n";
    }

    unless ($expressionhoglist =~ / $og /){ # Only write the code for HOGs with expression data (mostly single copy)
        next;
    }

    # Print R code to create the plot
    print ResultsR '
geneID <- "'."$og".'"

select_genes_tpm <- t(df_all_tpm_norm[geneID,])
select_genes_tpm <- merge(sample_newid, select_genes_tpm, by.x = "sample", by.y = "row.names")
names(select_genes_tpm)[names(select_genes_tpm) == "'."$og".'"] <- "norm_tpm"
select_genes_tpm$caste <- factor(select_genes_tpm$caste, levels = mylevel) 

a <- min(log2(select_genes_tpm$norm_tpm), na.rm=T)
if (a < 0){ a <- ceiling(-a) }

#q <- ggplot(select_genes_tpm, aes(x = factor(caste), y = log2(norm_tpm) + a, color = caste)) +
#  geom_boxplot(outlier.colour="white", outlier.shape = 10, outlier.size = 1, alpha = 0.8, size = 1, position = position_dodge(width=1)) +
#  #scale_color_manual(values = c(pal_nejm("default")(5)), labels = mylevel) +
#  scale_color_manual(values = c("#BC3C29FF", "#0072B5FF", "#000000"), labels = mylevel) +
#  scale_x_discrete(labels = mylevel) +
#  geom_point(aes(fill = caste), size=2, shape=21, alpha = 0.2, position = position_jitterdodge(jitter.width = 0, dodge.width = 1)) +
#  #scale_fill_manual(values = c(pal_nejm("default")(5)), labels = mylevel) +
#  scale_fill_manual(values = c("#BC3C29FF", "#0072B5FF", "#000000"), labels = mylevel) +
#  theme_bw() +
#  theme(legend.position = "none",
#        axis.text = element_text(size = 11, colour = "black"),
#        axis.text.x = element_text(size = 11, colour = "black", angle = 45, hjust = 1),
#        axis.title = element_text(size =11, face = "bold"),
#        legend.text = element_text(size = 11),
#        legend.title = element_text(size = 11, face = "bold"),
#        strip.text = element_text(size = 11, face = "bold")) +
#  labs(title = paste("Expresion abundance ", geneID), x = "Castes", y = "Expression abundance log2(TPM)" ) +
#  theme(plot.title = element_text(face = "bold", color = "black", size = 11, hjust = 0.5), strip.text.x = element_text(size = 12))
  
#print(q)
#ggsave(filename = paste0("'."$outname\_GAGA_expression_caste/".'expression_abundance", geneID, "_whole_body_samples_boxplot.pdf"), height = 4, width = 4, dpi = 300)

# Violin with boxplot
q <- ggplot(select_genes_tpm, aes(x = factor(caste), y = log2(norm_tpm) + a, color = caste)) +
  geom_violin(alpha = 0.8, size = 1, scale = "width") +
  geom_boxplot(width=0.25, outlier.colour="white", outlier.shape = 10, outlier.size = 1, alpha = 0.8, size = 1, position = position_dodge(width=0.25)) +
  # scale_color_manual(values = c(pal_nejm("default")(5)), labels = mylevel) +
  scale_color_manual(values = c("#BC3C29FF", "#0072B5FF", "#000000"), labels = mylevel) +
  scale_x_discrete(labels = mylevel) +
  geom_point(aes(fill = caste), size=2, shape=21, alpha = 0.2, position = position_jitterdodge(jitter.width = 0, dodge.width = 1)) +
  #  scale_fill_manual(values = c(pal_nejm("default")(5)), labels = mylevel) +
  scale_fill_manual(values = c("#BC3C29FF", "#0072B5FF", "#000000"), labels = mylevel) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 11, colour = "black"),
        axis.text.x = element_text(size = 11, colour = "black", angle = 45, hjust = 1),
        axis.title = element_text(size =11, face = "bold"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 11, face = "bold")) +
  labs(title = paste("Expresion abundance", geneID), x = "Castes", y = "Expression abundance log2(TPM)" ) +
  theme(plot.title = element_text(face = "bold", color = "black", size = 11, hjust = 0.5), strip.text.x = element_text(size = 12))

print(q)
ggsave(filename = paste0("'."$outname\_GAGA_expression_caste/".'Expression_abundance", geneID, "_whole_body_samples_violinplot.pdf"), height = 4, width = 4, dpi = 300)

';


}
close Filef;
close ResultsR;

# Run R script
system ("Rscript $outname\_GAGA_expression_script.R");
system ("mv $outname\_GAGA_expression_script.R $outname\_GAGA_expression_caste/");


#### Now do the same but to create the plot with subgroups: gyne, queen, worker, major, male. 
# Start R script
open (ResultsR, ">", "$outname\_GAGA_expression_subgroup_script.R");
print ResultsR 'library(ggplot2)
library(preprocessCore)
library(ggplot2)
library(RColorBrewer)
#library(hrbrthemes)
library(ggsci)
library(dplyr)
library(reshape2)

# Read files
all_tpm <- read.table("'."$expressionfile".'", header = TRUE, row.names = 1, check.names = FALSE)
colnames(all_tpm) <- gsub("-", "_", colnames(all_tpm)) #replace "-" with "_"
all_tpm <- all_tpm[-c(1,2,3),]
all_tpm <- replace(all_tpm, all_tpm == "na", NA)
all_tpm_bk <- all_tpm


all_tpm <- apply(all_tpm, 2, function(x) as.numeric(x))
all_tpm_norm <- normalize.quantiles(as.matrix(all_tpm))
colnames(all_tpm_norm) <- colnames(all_tpm)
rownames(all_tpm_norm) <- gsub("-", "_", rownames(all_tpm_bk))
t1 <- apply(all_tpm_norm, 2, function(x) quantile(x, probs = seq(0, 1, 1/4), na.rm = TRUE)) #To check the quantiles
df_all_tpm_norm <- as.data.frame(all_tpm_norm)

c_tissue <- c("brain", "antenna", "pupa", "gaster", "legs", "larva", "pupa", "head", "thorax", "abdomen")

# Read sample list
sample_newid_wb <- read.table("'."$samplefile".'", header = FALSE, col.names = c("species",  "sample", "group",  "newID", "caste")) 

sample_newid <- sample_newid_wb
sample_newid <- sample_newid[which(sample_newid_wb$caste != "unknown"), ] # remove unknown castes but tissue samples are still included.
sample_newid <- sample_newid[-grep(paste(c_tissue, collapse = "|"), sample_newid$sample),] # remove tissues
sample_newid <- sample_newid[which(sample_newid$species != "GAGA-0350"), ] # There 7 samples from GAGA-0350, which have wrong sample info.

subgroup = T
if (subgroup == T){
  sample_newid$caste[which(sample_newid$caste == "soldier")] <- "Large_worker" 
  sample_newid$caste[which(sample_newid$caste == "small_worker")] <- "Worker"
  sample_newid$caste[which(sample_newid$caste == "worker")] <- "Worker"
  sample_newid$caste[which(sample_newid$caste == "large_worker")] <- "Large_worker"
  sample_newid$caste[which(sample_newid$caste == "queen")] <- "Queen"
  sample_newid$caste[which(sample_newid$caste == "gyne")] <- "Gyne"
  sample_newid$caste[which(sample_newid$caste == "male")] <- "Male"
  mylevel = c("Queen", "Gyne", "Worker", "Large_worker", "Male")
}else{
  sample_newid$caste[which(sample_newid$caste == "gyne")] <- "Queen"
  sample_newid$caste[which(sample_newid$caste == "queen")] <- "Queen"
  sample_newid$caste[which(sample_newid$caste == "worker")] <- "Worker" 
  sample_newid$caste[which(sample_newid$caste == "soldier")] <- "Worker" 
  sample_newid$caste[which(sample_newid$caste == "small_worker")] <- "Worker"
  sample_newid$caste[which(sample_newid$caste == "large_worker")] <- "Worker"
  sample_newid$caste[which(sample_newid$caste == "male")] <- "Male"
  mylevel = c("Queen", "Worker", "Male")
}
table(sample_newid$caste)

sample_newid$sample <- gsub("-", "_", sample_newid$sample) #replace "-" with "_"

';
#close ResultsR;

# Read candidate list, and add them into the script
open(Filef, "<", $clist);
while(<Filef>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/ );
    #next if ($line =~ /Branches under selection/ || $line =~ /relaxation/);
    #my @subl = split (/\t/, $line);
    my $og = "";
    if ($line =~ /(HOG\d\d\d\d\d\d\d)/){
        $og = $1;
    } else {
        die "Cannot find OG in $line\n";
    }

    unless ($expressionhoglist =~ / $og /){ # Only write the code for HOGs with expression data (mostly single copy)
        next;
    }

    # Print R code to create the plot
    print ResultsR '
geneID <- "'."$og".'"

select_genes_tpm <- t(df_all_tpm_norm[geneID,])
select_genes_tpm <- merge(sample_newid, select_genes_tpm, by.x = "sample", by.y = "row.names")
names(select_genes_tpm)[names(select_genes_tpm) == "'."$og".'"] <- "norm_tpm"
select_genes_tpm$caste <- factor(select_genes_tpm$caste, levels = mylevel) 

a <- min(log2(select_genes_tpm$norm_tpm), na.rm=T)
if (a < 0){ a <- ceiling(-a) }

q <- ggplot(select_genes_tpm, aes(x = factor(caste), y = log2(norm_tpm) + a, color = caste)) +
  geom_violin(alpha = 0.8, size = 1, scale = "width") +
  geom_boxplot(width=0.25, outlier.colour="white", outlier.shape = 10, outlier.size = 1, alpha = 0.8, size = 1, position = position_dodge(width=0.25)) +
  # scale_color_manual(values = c(pal_nejm("default")(5)), labels = mylevel) +
  scale_color_manual(values = c("#BC3C29FF", "#e05b48", "#048fe0", "#136494", "#000000"), labels = mylevel) +
  scale_x_discrete(labels = c("Queen", "Gyne", "Worker", "Large\nWorker", "Male")) +  # diagonal labels
  #scale_x_discrete(labels = mylevel) +
  geom_point(aes(fill = caste), size=2, shape=21, alpha = 0.2, position = position_jitterdodge(jitter.width = 0, dodge.width = 1)) +
  #  scale_fill_manual(values = c(pal_nejm("default")(5)), labels = mylevel) +
  scale_fill_manual(values = c("#BC3C29FF", "#e05b48", "#048fe0", "#136494", "#000000"), labels = mylevel) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 11, colour = "black"),
        axis.text.x = element_text(size = 11, colour = "black", angle = 45, hjust = 1),
        axis.title = element_text(size =11, face = "bold"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 11, face = "bold")) +
  labs(title = paste("Expresion abundance", geneID), x = "Castes", y = "Expression abundance log2(TPM)" ) +
  theme(plot.title = element_text(face = "bold", color = "black", size = 11, hjust = 0.5), strip.text.x = element_text(size = 12))

print(q)
ggsave(filename = paste0("'."$outname\_GAGA_expression_caste_subgroups/".'Expression_abundance", geneID, "_whole_body_samples_subgroups_violinplot.pdf"), height = 4, width = 5, dpi = 300)

';


}
close Filef;
close ResultsR;

# Run R script
system ("Rscript $outname\_GAGA_expression_subgroup_script.R");

system ("mv $outname\_GAGA_expression_subgroup_script.R $outname\_GAGA_expression_caste_subgroups/");

system ("rm Rplots.pdf"); # Remove a file that R prints with all plots together
