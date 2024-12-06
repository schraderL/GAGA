#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
    cat("Rscript topGO_showTopTerms.R <goEnrichment_result> <topNum> \n")
    quit('no')
}else if (length(args) > 1) {
    goEnrichment_result <- args[1]
    topNum <- args[2]
}

library(ggplot2)
library(scales)


goEnrichment<-read.table(goEnrichment_result,sep="\t",header=T)
goEnrichment$Fisher <- as.numeric(goEnrichment$Fisher)
goEnrichment <- goEnrichment[goEnrichment$Fisher < 0.05,] 
goEnrichment <- goEnrichment[,c("GO.ID","Term","Fisher")]

ntop <- min(as.numeric(topNum), nrow(goEnrichment)) # If there are less than the specified number of representative GOs, just use all from the table
ggdata <- goEnrichment[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
pdf(paste(goEnrichment_result,".FisherPlot.pdf",sep=""),height=16,width=12)
ggplot(ggdata,
  aes(x = Term, y = -log10(Fisher), size = -log10(Fisher), fill = -log10(Fisher))) +

  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +

  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'GO enrichment, Top',
    subtitle = paste(ntop," terms ordered by Fisher Exact p-value"),
    caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +

  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
    linetype = c("dotted", "dotted", "dotted"),
    colour = c("red", "red", "red"),
    size = c(1, 1, 1)) +

  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),

    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),

    legend.key = element_blank(), 
    legend.key.size = unit(1, "cm"), 
    legend.text = element_text(size = 14, face = "bold"),
    title = element_text(size = 14, face = "bold")) +

  coord_flip()

dev.off()

ntop <- nrow(goEnrichment) # If there are less than the specified number of representative GOs, just use all from the table
ggdata <- goEnrichment
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
pdf(paste(goEnrichment_result,".FisherPlot_all.pdf",sep=""),height=16,width=12)
ggplot(ggdata,
  aes(x = Term, y = -log10(Fisher), size = -log10(Fisher), fill = -log10(Fisher))) +

  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +

  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'GO enrichment, Top',
    subtitle = paste(ntop," terms ordered by Fisher Exact p-value"),
    caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +

  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
    linetype = c("dotted", "dotted", "dotted"),
    colour = c("red", "red", "red"),
    size = c(1, 1, 1)) +

  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),

    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),

    legend.key = element_blank(), 
    legend.key.size = unit(1, "cm"), 
    legend.text = element_text(size = 14, face = "bold"),
    title = element_text(size = 14, face = "bold")) +

  coord_flip()

dev.off()