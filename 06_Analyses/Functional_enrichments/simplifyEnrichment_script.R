# Script to plot GO enrichments using jokergoo summarizeGO package

# Usage:  Rscript simplifyEnrichment_script.R Formicoids_genecandidates_positiveselection_GO_enrichment/Formicoids_genecandidates_positiveselection.txt.BP_weight01_pval001.txt
# It will print all output files with the same prefix as the input and in the same directory

# Help page: https://jokergoo.github.io/2023/10/02/simplified-simplifyenrichment-plot/
# No model organism : https://jokergoo.github.io/2023/03/30/use-simplifyenrichment-for-non-model-organisms/
# Although it has much less GO, just generate an additional plot with fly

# install
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("summarizeGO")

# Install through github - latest version including summarizeGO function
#install.packages("remotes")
#remotes::install_github("jokergoo/simplifyGO")

# Load modules
library(simplifyEnrichment)

# here read input file entered in the first argument
file <- commandArgs(trailingOnly = TRUE)[1]
GOtable <- read.csv(file, header = TRUE, sep = "\t")

# Detect if it is GOstats or TopGO input
# Use df as input table
df = GOtable

# Get the list of ids in first column
if("GOBPID" %in% colnames(df)) {
  # If input is GOstats
  go_id = df$GOBPID[df$Pvalue < 0.05]
} else {
  # Else input is TopGO
  go_id = df$GO.ID[df$Fisher < 0.05]
}


# get output name, being file without .txt 
outname <- sub("\\.txt$", "", file)
#and adding _simplifyGO.pdf
outsimpgo <- paste0(outname, "_simplifyGO_default.pdf") 

# Run simplify enrichment
matdef = GO_similarity(go_id)

#pdf("Formicoids_genecandidates_positiveselection_GOenrichment_BP_pval001_catsize10_simplifyGO.pdf") # Use filename depending on file and using full path
pdf(outsimpgo)
simplifyGO(matdef)
dev.off()


### Using Drosophila as model 
library("org.Dm.eg.db")
outsimpgodmel <- paste0(outname, "_simplifyGO_Dmel.pdf") 
matdmel = GO_similarity(go_id, db = "org.Dm.eg.db")
pdf(outsimpgodmel)
simplifyGO(matdmel)
dev.off()


# Simpler plot without heatmap
l = if("GOBPID" %in% colnames(df)) df$Pvalue < 0.05 else df$Fisher < 0.05


pdf(paste0(outname, "_simplifyGO_sum.pdf"))
if("GOBPID" %in% colnames(df)) {
  summarizeGO(df$GOBPID[l], -log10(df$Pvalue)[l], axis_label = "average -log10(p.adjust)")
} else {
  summarizeGO(df$GO.ID[l], -log10(df$Fisher)[l], axis_label = "average -log10(p.adjust)")
}
dev.off()
#system(paste0("pdftk ", outname, "_summarizeGO.pdf cat 2-end output ", outname, "_summarizeGO_ok.pdf"))

# Heatmap with other non model organisms as orgdb : Double comment ## is in lines that should be used
##library(AnnotationHub)
##ah = AnnotationHub()
##ah = ah[ah$rdataclass == "OrgDb"]
#meta = mcols(ah)
#rownames(meta)[ which(meta$species == "Apis mellifera") ]
#rownames(meta)[ which(meta$species == "Monomorium pharaonis") ]
#rownames(meta)[ which(meta$species == "Ooceraea biroi") ]
#rownames(meta)[ which(meta$species == "Solenopsis invicta") ]
#rownames(meta)[ which(meta$species == "Drosophila melanogaster") ]

#orgdb = ah[["AH108737"]] # Mpha
#orgdb = ah[["AH109059"]] # Obir
#orgdb = ah[["AH108746"]] # Sinv

#orgdb = ah[["AH108106"]] # Dolphin example
##orgdb = ah[["AH107058"]] # Fly
#orgdb = ah[["AH109223"]] # honeybee

#orgdb
#columns(orgdb)

# Visualize now using the organism selected
##mat = GO_similarity(go_id, db = orgdb)
#mat = GO_similarity(go_id, db = orgdb, ont = "BP")
##pdf(paste0(outname, "_orgdbfly.pdf"))
##simplifyGO(mat)
##dev.off()


