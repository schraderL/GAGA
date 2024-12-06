#!/usr/bin/perl
use strict;
use warnings;


# Usage: perl get_pgls_genefam_vs_traits.pl


my $intableref = "/home/projects/ku_00039/people/joeviz/PGLS/Traits_table_v6_18Dec_toPGLS_reduced.tsv";
#my $intableref = "/home/projects/ku_00039/people/joeviz/PGLS/Traits_table_v6_18Dec_toPGLS_reduced_nostLFRwLept.tsv";

#my $intable = "/home/projects/ku_00039/people/joeviz/PGLS/GAGA_genefamilies_reannot_withORsubfam_toPGLS_nostLFR.txt"; 
my $intable = "/home/projects/ku_00039/people/joeviz/PGLS/Traits_table_v6_18Dec_toPGLS_reduced.tsv";

my $intree = "/home/projects/ku_00039/people/joeviz/PGLS/GAGA_dated_phylogeny_newick.tre";
#my $intree = "/home/projects/ku_00039/people/joeviz/PGLS/GAGA_dated_phylogeny_newick_nostLFR.tre";


#my $outputfolder = "/home/projects/ku_00039/people/joeviz/PGLS/PGLS_final_traits_Vs_traits";
#my $outputfolder = "/home/projects/ku_00039/people/joeviz/PGLS/PGLS_final_Genefamreannot1to1filt_Vs_traits";

system ("mkdir $outputfolder");
system ("mkdir $outputfolder\/Venn_plots");
system ("mkdir $outputfolder\/PGLS_output");

# Copy submit script
system ("cp submit_Rpgls.sh $outputfolder\/");

my $output = "$outputfolder\/Rscript_PGLS_traitcomp.R";
#if ($sptable =~ /(\S+)\.txt/ ){
#	$output = "$1\_fulltable.txt";
#} else {
#	$output = "sptable\_fulltable.txt";
#}

my @headerref; my $headerrefcheck = 0;
my $headerrefline = ""; my $headerreflinecomma = "";

open (Table, "<", $intableref);
while(<Table>){
	chomp;	
	my $line = $_;
	next if ($line !~ /\S+/);
	$line =~ s/ /\./g; # Remove spaces for "."	
	$line =~ s/\-/\./g; # Remove spaces for "."	Because R transform them to .
	$line =~ s/\r//g; # Remove new lines missed by chomp	

	#$line =~ s/[^a-zA-Z0-9,]+//g; # Remove all non alphanumeric characters
	#$line =~ s/[^\S]+//g; # Remove all non alphanumeric characters

	my @subl = split (/\t/, $line);
	if ($headerrefcheck < 1 ){
		# Header
		$headerrefcheck++;
		foreach my $traithead (@subl){
			#push (@header, $traithead);
			unless ($traithead =~ /GAGA\.ID/ || $traithead =~ /Species\.name/){
#				if ($traithead =~ /Polygyny/ || $traithead =~ /Polymorphism/){
					$headerrefline .= "$traithead\t";
					$headerreflinecomma .= "$traithead\, ";
					push (@headerref, $traithead);
#				}
			}
		}
	} 

}
close Table;

print "Reference trait table: $headerrefline\n"; # Check traits

# Opening second table
my @header; my $headercheck = 0;
my $headerline = ""; my $headerlinecomma = "";

open (Table, "<", $intable);
while(<Table>){
	chomp;	
	my $line = $_;
	next if ($line !~ /\S+/);
	$line =~ s/ /\./g; # Remove spaces for "."	
	$line =~ s/\-/\./g; # Remove spaces for "."	Because R transform them to .
	$line =~ s/\r//g; # Remove new lines missed by chomp	

	#$line =~ s/[^a-zA-Z0-9,]+//g; # Remove all non alphanumeric characters
	#$line =~ s/[^\S]+//g; # Remove all non alphanumeric characters

	my @subl = split (/\t/, $line);
	if ($headercheck < 1 ){
		# Header
		$headercheck++;
		foreach my $traithead (@subl){
			#push (@header, $traithead);
			unless ($traithead =~ /GAGA\.ID/ || $traithead =~ /Species\.name/){
#				if ($traithead =~ /Polygyny/ || $traithead =~ /Polymorphism/){
					$headerline .= "$traithead\t";
					$headerlinecomma .= "$traithead\, ";
					push (@header, $traithead);
#				}
			}
		}
	} 
	
}
close Table;

print "Second trait tables: $headerline\n"; # Check traits


## Open now the script and print all comparisons

open (Results, ">", $output);
print Results 'setwd("'."$outputfolder".'")

library(ape)
library(phytools)
library(tidyverse)
library(caper)
require(MASS)
require(dplyr)

ALL_TRAITSREF <- read.csv("'."$intableref".'",sep="\t",dec = ".", quote = "")
ALL_TRAITS <- read.csv("'."$intable".'",sep="\t",dec = ".", quote = "")

tree <- read.tree("'."$intree".'")

#sink("PGLS_results_all.txt")

sink("PGLS_results_matrix_all_coeffs.txt") # To print a matrix table with regression coefficients
cat("\t")
cat("'."$headerrefline".'")
cat("\n")
sink()

sink("PGLS_results_matrix_all_signtrend.txt") # To print a matrix table with the trend from regression coefficients in significant cases
cat("\t")
cat("'."$headerrefline".'")
cat("\n")
sink()

sink("PGLS_results_matrix_all.txt") # To print a matrix table
cat("\t")
cat("'."$headerrefline".'")
cat("\n")
sink()

';

#print Results 'ALL_TRAITS %>%  dplyr::select(GAGA.ID,'."$headerlinecomma".') %>% as.data.frame() -> mm'."\n"; # To use the mm with all traits and see if there are differences


my $trn = 0;
foreach my $tr (@header){

	#print Results "cat(\"$tr\\t\")\n"; # For matrix output
	
	print Results 'sink("PGLS_results_matrix_all.txt", append = T)'."\n"; # Open matrix to add trait
	print Results 'cat("'."$tr".'\t")'."\n";
	print Results 'sink()'."\n"; 
	print Results 'sink("PGLS_results_matrix_all_coeffs.txt", append = T)'."\n"; # Open matrix to add trait line
	print Results 'cat("'."$tr".'\t")'."\n";
	print Results 'sink()'."\n"; 
	print Results 'sink("PGLS_results_matrix_all_signtrend.txt", append = T)'."\n"; # Open matrix to add trait line
	print Results 'cat("'."$tr".'\t")'."\n";
	print Results 'sink()'."\n"; 

	my $tr2n = 0;
	foreach my $tr2 (@headerref){
#		if ($trn == $tr2n){ # Skip with itself
		if ($tr =~ /^$tr2$/){ # Skip with itself
			$tr2n++;
			
			#print Results "cat(\"\\t\")\n"; # For matrix output
			
			print Results 'sink("PGLS_results_matrix_all.txt", append = T)'."\n"; # 
			print Results 'cat("\t")'."\n";
			print Results 'sink()'."\n"; 
			print Results 'sink("PGLS_results_matrix_all_coeffs.txt", append = T)'."\n"; # 
			print Results 'cat("\t")'."\n";
			print Results 'sink()'."\n"; 
			print Results 'sink("PGLS_results_matrix_all_signtrend.txt", append = T)'."\n"; # 
			print Results 'cat("\t")'."\n";
			print Results 'sink()'."\n"; 

			next;
		}

		print Results "\#Comparing $tr with $tr2\n";
		print Results 'ALL_TRAITS %>%  dplyr::select(GAGA.ID, '."$tr ".') %>% as.data.frame() -> mmtable'."\n";
		print Results 'ALL_TRAITSREF %>%  dplyr::select(GAGA.ID, '."$tr2 ".') %>% as.data.frame() -> mmref'."\n";
		print Results 'mm <- inner_join(x = mmtable, y = mmref, by = "GAGA.ID")'."\n";
		print Results 'comp <- comparative.data(tree, mm, GAGA.ID, vcv=TRUE, vcv.dim=3)'."\n";
		print Results 'Model_1 = try(pgls('."$tr \~ $tr2".', comp, lambda="ML"))'."\n";
#		print Results "summary(Model_1)\n";  # Comment for matrix output
		print Results "pval <- try(anova(Model_1))\n";
#		print Results 'pval$`Pr(>F)`'."\n";
#		print Results 'cat("'."$tr VERSUS $tr2 correlation pvalue ".'")'."\n".'cat(pval$`Pr(>F)`)'."\n".'cat("\n")'."\n"; # Comment for matrix output
#		print Results 'sink()'."\n"; # Exit the sink, so it does not print "[1] 1" in case that there is error 
		print Results 'pvalue <- tryCatch(unlist(pval$`Pr(>F)`), error = function(e) print("1"))'."\n";
		print Results 'sink("PGLS_results_matrix_all.txt", append = T)'."\n"; # Open again the matrix file to print the pvalue		
		print Results 'cat(pvalue[1])'."\n".'cat("\t")'."\n"; # Uncomment for matrix output
		print Results 'sink()'."\n".'sink("'."$outputfolder\/PGLS_output\/$tr\_Vs_$tr2\_summary.txt".'")'."\n".'try(summary(Model_1))'."\n"; # Print here the PGLS summary
		print Results 'sink()'."\n";

		# print coefficients
		print Results 'coef_values <- try(coef(summary(Model_1))[, "Estimate"])
		estimate_value <- tryCatch(unlist(as.numeric(coef_values["'."$tr2".'"])), error = function(e) print("0"))
		sink("PGLS_results_matrix_all_coeffs.txt", append = T)
		cat(estimate_value)
		cat("\t")
		sink()

		sink("PGLS_results_matrix_all_signtrend.txt", append = T)
		try(if (!is.na(estimate_value) && !is.na(pvalue[1])) {
		  if (estimate_value > 0) {
		    if (pvalue[1] < 0.05) {
		      cat("+")
		    }
		    if (pvalue[1] < 0.01) {
		      cat("+")
		    }
		    if (pvalue[1] < 0.001) {
		      cat("+")
		    }
		  } else {
		    if (pvalue[1] < 0.05) {
		      cat("-")
		    }
		    if (pvalue[1] < 0.01) {
		      cat("-")
		    }
		    if (pvalue[1] < 0.0001) {
		      cat("-")
		    }
		  }
		})
		cat("\t")
		sink()


		'."\n"; # print coefficients

		# Print violin plot commands; if pvalue es > 0.05, hacer el violin plot; si hay mas de 10 categories (i.e: numerical traits), hacer simplemente la correlacion
		print Results 'if (pvalue[1] < 0.05){
  numcat <- length(unique(mm$'."$tr2".'))
  if (numcat < 10){			
	  mmplot <- mm
	  mmplot$'."$tr2".' <- as.factor(mmplot$'."$tr2".')
	  p1 <-  ggplot(mmplot ,aes(y='."$tr".', x='."$tr2".', color='."$tr2".')) + 
	    geom_violin() 
	  pdf("'."$outputfolder\/Venn_plots\/$tr\_Vs_$tr2\_violin_boxplot.pdf".'")
	  try(print(p1 + geom_boxplot(width=0.1) +scale_color_brewer(palette="Set1")))
	  dev.off()
	  pdf("'."$outputfolder\/Venn_plots\/$tr\_Vs_$tr2\_violin_meanplot.pdf".'")
	  try(print(p1 + stat_summary(fun.data="mean_sdl", geom="pointrange") +scale_color_brewer(palette="Set1")))
	  dev.off()  
  }	  
  if (numcat >= 10) {
    plot_data = tibble(Y=Model_1$y, X=Model_1$x[,2], Fitted=Model_1$fitted )

    # Correlation plots
    p1 <- ggplot(plot_data, aes(x = X, y = Y)) +
      geom_point(alpha = 0.6) +
      theme_bw() +
      geom_smooth(method = "lm")+
      xlab("'."$tr2".'") + ylab("'."$tr".'")

    pdf("'."$outputfolder\/Venn_plots\/$tr\_Vs_$tr2\_correlation_plot.pdf".'")
    try(print(p1))
    dev.off()

    p2 <- ggplot(plot_data, aes(x = X, y = Y)) +
      geom_point(alpha = 0.6) +
      theme_bw() +
      scale_y_log10() +
      scale_x_log10() +
      geom_smooth(method = "lm")+
      xlab("log10 '."$tr2".'") + ylab("log10 '."$tr".'")  

    pdf("'."$outputfolder\/Venn_plots\/$tr\_Vs_$tr2\_correlation_plot_log10.pdf".'")
    try(print(p2))
    dev.off()

    p3 <- ggplot(plot_data, aes(x = X, y = Y)) +
      geom_point(alpha = 0.6) +
      theme_bw() +
      scale_y_log10() +
      geom_smooth(method = "lm")+
      xlab("'."$tr2".'") + ylab("log10 '."$tr".'")  

    pdf("'."$outputfolder\/Venn_plots\/$tr\_Vs_$tr2\_correlation_plot_log10y.pdf".'")
    try(print(p3))
    dev.off()


  }
}'."\n";

		#print Results 'sink("PGLS_results_matrix_all.txt", append = T)'."\n"; # Open again the matrix file || Commented, it is opened in each step
		print Results "\n\n";

		$tr2n++;

	}

	$trn++;

	print Results 'sink("PGLS_results_matrix_all.txt", append = T)'."\n"; # Open matrix to add new line
	print Results 'cat("\n")'."\n";
	print Results 'sink()'."\n"; 
	print Results 'sink("PGLS_results_matrix_all_coeffs.txt", append = T)'."\n"; # Open matrix to add new line
	print Results 'cat("\n")'."\n";
	print Results 'sink()'."\n"; 
	print Results 'sink("PGLS_results_matrix_all_signtrend.txt", append = T)'."\n"; # Open matrix to add new line
	print Results 'cat("\n")'."\n";
	print Results 'sink()'."\n"; 
}

print Results 'sink()'."\n\n\n";

# Get FDR corrected matrix (per column)
print Results 'pvaltable <- read_tsv(file="PGLS_results_matrix_all.txt")
fdrtable <- as.data.frame(lapply(pvaltable[-1], \(x) p.adjust (x, method = "fdr")))
fdrtablefull <- cbind(pvaltable$...1,fdrtable)
write_tsv(fdrtablefull, file="PGLS_results_matrix_all_fdr.txt")
';

close Results;


# Run these lines when the R script finishes to organize the anayses for specific traits: script_get_trait_results_order.sh
=h
system ("mkdir $outputfolder\/Traits_colony_size_plots");
system ("cp $outputfolder\/Venn_plots/*olony* $outputfolder\/Traits_colony_size_plots");
system ("mkdir $outputfolder\/Traits_colony_size_txt");
system ("cp $outputfolder\/PGLS_output/*olony* $outputfolder\/Traits_colony_size_txt");

system ("mkdir $outputfolder\/Traits_dimorphism_plots");
system ("cp $outputfolder\/Venn_plots/*imorphis* $outputfolder\/Traits_dimorphism_plots");
system ("mkdir $outputfolder\/Traits_dimorphism_txt");
system ("cp $outputfolder\/PGLS_output/*imorphis* $outputfolder\/Traits_dimorphism_txt");

system ("mkdir $outputfolder\/Traits_gamergates_plots");
system ("cp $outputfolder\/Venn_plots/*amergate* $outputfolder\/Traits_gamergates_plots");
system ("mkdir $outputfolder\/Traits_gamergates_txt");
system ("cp $outputfolder\/PGLS_output/*amergate* $outputfolder\/Traits_gamergates_txt");

system ("mkdir $outputfolder\/Traits_polyandry_plots");
system ("cp $outputfolder\/Venn_plots/*olyandry* $outputfolder\/Traits_polyandry_plots");
system ("mkdir $outputfolder\/Traits_polyandry_txt");
system ("cp $outputfolder\/PGLS_output/*olyandry* $outputfolder\/Traits_polyandry_txt");

system ("mkdir $outputfolder\/Traits_social_parasitism_plots");
system ("cp $outputfolder\/Venn_plots/*parasitism* $outputfolder\/Traits_social_parasitism_plots");
system ("mkdir $outputfolder\/Traits_social_parasitism_txt");
system ("cp $outputfolder\/PGLS_output/*parasitism* $outputfolder\/Traits_social_parasitism_txt");

system ("mkdir $outputfolder\/Traits_trophallaxis_plots");
system ("cp $outputfolder\/Venn_plots/*ophallaxi* $outputfolder\/Traits_trophallaxis_plots");
system ("mkdir $outputfolder\/Traits_trophallaxis_txt");
system ("cp $outputfolder\/PGLS_output/*ophallaxi* $outputfolder\/Traits_trophallaxis_txt");

system ("mkdir $outputfolder\/Traits_polygyny_plots");
system ("cp $outputfolder\/Venn_plots/*olygyn* $outputfolder\/Traits_polygyny_plots");
system ("mkdir $outputfolder\/Traits_polygyny_txt");
system ("cp $outputfolder\/PGLS_output/*olygyn* $outputfolder\/Traits_polygyny_txt");

system ("mkdir $outputfolder\/Traits_polymorphism_plots");
system ("cp $outputfolder\/Venn_plots/*olymorphis* $outputfolder\/Traits_polymorphism_plots");
system ("cp $outputfolder\/Venn_plots/*oldier* $outputfolder\/Traits_polymorphism_plots");
system ("mkdir $outputfolder\/Traits_polymorphism_txt");
system ("cp $outputfolder\/PGLS_output/*olymorphis* $outputfolder\/Traits_polymorphism_txt");
system ("cp $outputfolder\/PGLS_output/*oldier* $outputfolder\/Traits_polymorphism_txt");

system ("mkdir $outputfolder\/Traits_trophobiosis_plots");
system ("cp $outputfolder\/Venn_plots/*ophobiosi* $outputfolder\/Traits_trophobiosis_plots");
system ("mkdir $outputfolder\/Traits_trophobiosis_txt");
system ("cp $outputfolder\/PGLS_output/*ophobiosi* $outputfolder\/Traits_trophobiosis_txt");
=cut


