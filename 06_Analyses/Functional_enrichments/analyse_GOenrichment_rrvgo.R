# Script to plot GO enrichments using rrvgo package

# https://ssayols.github.io/rrvgo/articles/rrvgo.html

# Load modules
library(tidyverse)
library(vroom)
library(ggplot2)
library(ggrepel)
#devtools::install_github("ssayols/rrvgo")
library(rrvgo)
#BiocManager::install("org.Dm.eg.db")


############################################################
# Load the GO enrichment analysis
# read input file entered in the first argument
#file <- commandArgs(trailingOnly = TRUE)[1]
file <- commandArgs(trailingOnly = TRUE)[1]
GOtable <- read.csv(file, header = TRUE, sep = "\t")

#GOtable <- read.csv("/Users/lukas/Downloads/All_pairwise_comp_genes_shared_module_noduphogsawk.tsv_GOenrichment_BP_pval001.txt", header = TRUE, sep = "\t")

# Detect if it is GOstats or TopGO input
# Use df as input table
df = GOtable

# Get the list of ids in first column
if("GOBPID" %in% colnames(df)) {
  # If input is GOstats
  go_id = df$GOBPID[df$Pvalue < 0.05]
  scores <- setNames(-log10(df$Pvalue), df$GOBPID)
} else {
  # Else input is TopGO
  go_id = df$GO.ID[df$Fisher < 0.05]
  scores <- setNames(-log10(df$Fisher), df$GO.ID)
}


simMatrix <- calculateSimMatrix(go_id,
                                orgdb="org.Dm.eg.db",
                                ont="BP",
                                method="Rel")

#scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Dm.eg.db")

############################################################
##Similarity matrix heatmap
############################################################

# get output name, being file without .txt 
outname <- sub("\\.txt$", "", file)
#and adding file extension
outheatmap <- paste0(outname, "_rrvgo_heatmap.pdf") 


raw.plot.heatmap <- heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)

#raw.plot.heatmap
ggsave(file=outheatmap,raw.plot.heatmap,width=15,height=15)

############################################################
## Scatter plot depicting groups and distance between terms
############################################################

outscatter <- paste0(outname, "_rrvgo_scatterplot.pdf") 

raw.plot.scatter <- scatterPlot(simMatrix, reducedTerms)
#raw.plot.scatter
ggsave(file=outscatter,raw.plot.scatter,width=15,height=15)

#outscatterv <- paste0(outname, "_rrvgo_scatterplotv.pdf") 
#pdf(outscatterv) # It losses quality this way
#raw.plot.scatter
#dev.off()

#outscatterlabel <- paste0(outname, "_rrvgo_scatterplot_addlabel.pdf") 
#raw.plot.scatterlabel <- scatterPlot(simMatrix, reducedTerms,addLabel=T)
#raw.plot.scatterlabel
#ggsave(file=outscatterlabel,raw.plot.scatterlabel,width=15,height=15)


# use edited scatterPlot function
outscattermod <- paste0(outname, "_rrvgo_scatterplot_mod.pdf") 

df.n10<-scatterPlot(simMatrix,reducedTerms,addLabel=T,onlyParents=T,size="size")$data
#wordcloudPlot(reducedTerms.n10, min.freq=0)

# Add shifted labels

labelSize <- 1.5
p.n10<-ggplot2::ggplot(df.n10, ggplot2::aes(x = V1, y = V2, color = parentTerm)) +
  ggplot2::geom_point(ggplot2::aes(size = size), alpha = 0.5) +
  ggplot2::scale_color_discrete(guide = "none") + ggplot2::scale_size_continuous(range = c(0, 20)) + ggplot2::scale_x_continuous(name = "") +
  ggplot2::scale_y_continuous(name = "") + ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank()
  ) + geom_text_repel(
  data          = subset(df.n10, V1 > 0 & parent == rownames(df.n10)),
  aes(label         = parentTerm),
  nudge_x       = 1.4 - subset(df.n10, V1 > 0 & parent == rownames(df.n10))$V1,
  segment.size  = 0.2,
  segment.color = "grey50",
  direction     = "y",
  hjust         = 0
) +
  geom_text_repel(
    data          = subset(df.n10, V1 <= 0 & parent == rownames(df.n10)),
    aes(label         = parentTerm),
    nudge_x       = -1.2 - subset(df.n10, V1 <= 0 & parent == rownames(df.n10))$V1,
    segment.size  = 0.2,
    segment.color = "grey50",
    direction     = "y",
    hjust         = 1
  )+coord_cartesian(xlim=c(-10,10))

#p.n10
ggsave(file=outscattermod,p.n10,width=15,height=8)


############################################################
## Treemap plot
############################################################

outtreemap <- paste0(outname, "_rrvgo_treemap.pdf") 

#raw.plot.treemap <- treemapPlot(reducedTerms)
#raw.plot.treemap
#ggsave(file=outtreemap,raw.plot.treemap,width=15,height=15)

#pdf(outtreemap)
pdf(outtreemap, width = 15, height = 15)
treemapPlot(reducedTerms)
dev.off()

############################################################
## Word cloud
############################################################

outword <- paste0(outname, "_rrvgo_wordcloud.pdf") 

raw.plot.word <- wordcloudPlot(reducedTerms, min.freq=1, colors="black")
#raw.plot.word
ggsave(file=outword,raw.plot.word,width=15,height=15)

pdf(outword, width = 15, height = 15)
wordcloudPlot(reducedTerms, min.freq=1, colors="black")
dev.off()


#outword <- paste0(outname, "_rrvgo_wordcloud_red.pdf") 
#pdf(outword, width = 15, height = 15)
#wordcloudPlot(reducedTerms, min.freq=2, colors="black")
#dev.off()

outword <- paste0(outname, "_rrvgo_wordcloud_red_color.pdf") 
pdf(outword, width = 15, height = 15)
#wordcloudPlot(reducedTerms, min.freq=2, colors=c("blue", "orange", "green"))
wordcloudPlot(reducedTerms, min.freq=2, colors=c("#99CCFF", "#FFD699", "#B2DFB2"))
dev.off()

