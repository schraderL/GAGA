#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
    cat("Rscript topGO_enrich.R <gene2GO map_file> <relax genelist> \n")
    quit('no')
}else if (length(args) > 1) {
    mapfile <- args[1]
    relaxInput<- args[2]
}


library(topGO)
library(GO.db)

gene_id = readMappings(file = mapfile)
gene_names = names(gene_id)
gene_list = rep(0,length(gene_id))
names(gene_list) = names(gene_id)
my_genes = read.table(relaxInput)[,1]
gene_list[match(my_genes,names(gene_list))] = 1
top_diff_genes = function(x)(x == 1)

## cellular component
CC_go <- new("topGOdata",nodeSize=5,ontology = "CC",allGenes = gene_list,annot = annFUN.gene2GO,gene2GO = gene_id,geneSel=top_diff_genes)
result_fisher.elim = runTest(CC_go,algorithm = "elim",statistic = "fisher")
allres = GenTable(CC_go,Fisher = result_fisher.elim,ranksOf = "classicFisher",numChar = "9999",topNodes = attributes(result_fisher.elim)$geneData[4])
filtres <- allres[allres$Fisher < 0.05,] 
write.table(filtres,file = paste(relaxInput,".CC_elim_pval005.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
filtrescc <- allres[allres$Fisher < 0.01,] 
write.table(filtrescc,file = paste(relaxInput,".CC_elim_pval001.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
pdf(paste(relaxInput,".CC_elim.pdf",sep=""))
showSigOfNodes(CC_go,score(result_fisher.elim),firstSigNodes = 5,useInfo = "all")
dev.off()
printGraph(CC_go,result_fisher.elim,firstSigNodes = 5,fn.prefix = "tGO.CC",useInfo = "pval",pdfSW = TRUE)

## molecular function
MF_go <- new("topGOdata",nodeSize=5,ontology = "MF",allGenes = gene_list,annot = annFUN.gene2GO,gene2GO = gene_id,geneSel=top_diff_genes)
result_fisher.elim = runTest(MF_go,algorithm = "elim",statistic = "fisher")
allres = GenTable(MF_go,Fisher = result_fisher.elim,ranksOf = "classicFisher",numChar = "9999",topNodes = attributes(result_fisher.elim)$geneData[4])
filtres <- allres[allres$Fisher < 0.05,] 
write.table(filtres,file = paste(relaxInput,".MF_elim_pval005.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
filtresmf <- allres[allres$Fisher < 0.01,] 
write.table(filtresmf,file = paste(relaxInput,".MF_elim_pval001.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
pdf(paste(relaxInput,".MF_elim.pdf",sep=""))
showSigOfNodes(MF_go,score(result_fisher.elim),firstSigNodes = 5,useInfo = "all")
dev.off()
printGraph(MF_go,result_fisher.elim,firstSigNodes = 5,fn.prefix = "tGO.MF",useInfo = "pval",pdfSW = TRUE)

## biological process
BP_go <- new("topGOdata",nodeSize=5,ontology = "BP",allGenes = gene_list,annot = annFUN.gene2GO,gene2GO = gene_id,geneSel=top_diff_genes)
result_fisher.elim = runTest(BP_go,algorithm = "elim",statistic = "fisher")
allres = GenTable(BP_go,Fisher = result_fisher.elim,ranksOf = "classicFisher",numChar = "9999",topNodes = attributes(result_fisher.elim)$geneData[4])
filtres <- allres[allres$Fisher < 0.05,] 
write.table(filtres,file = paste(relaxInput,".BP_elim_pval005.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
filtresbp <- allres[allres$Fisher < 0.01,] 
write.table(filtresbp,file = paste(relaxInput,".BP_elim_pval001.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
pdf(paste(relaxInput,".BP_elim.pdf",sep=""))
showSigOfNodes(BP_go,score(result_fisher.elim),firstSigNodes = 5,useInfo = "all")
dev.off()
printGraph(BP_go,result_fisher.elim,firstSigNodes = 5,fn.prefix = "tGO.BP",useInfo = "pval",pdfSW = TRUE)

#Combine tables with BP, CC and MF
combined_table <- rbind(filtresbp, filtresmf, filtrescc)
colnames(combined_table) <- c("GO.ID", "Term", "Annotated", "Significant", "Expected", "P_value")
subtable <- combined_table[,c("GO.ID","P_value")]
colnames(subtable) <- c("% GOterm", "enrichment_P-value")
subtable <- subtable[order(subtable$`enrichment_P-value`), ]

write.table(subtable,file = paste(relaxInput,".GOterms_elim_pval001.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


# Now with classic method
## cellular component
CC_go <- new("topGOdata",nodeSize=5,ontology = "CC",allGenes = gene_list,annot = annFUN.gene2GO,gene2GO = gene_id,geneSel=top_diff_genes)
result_fisher.elim = runTest(CC_go,algorithm = "classic",statistic = "fisher")
allres = GenTable(CC_go,Fisher = result_fisher.elim,ranksOf = "classicFisher",numChar = "9999",topNodes = attributes(result_fisher.elim)$geneData[4])
filtres <- allres[allres$Fisher < 0.05,] 
write.table(filtres,file = paste(relaxInput,".CC_classic_pval005.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
filtrescc <- allres[allres$Fisher < 0.01,] 
write.table(filtrescc,file = paste(relaxInput,".CC_classic_pval001.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
pdf(paste(relaxInput,".CC_classic.pdf",sep=""))
showSigOfNodes(CC_go,score(result_fisher.elim),firstSigNodes = 5,useInfo = "all")
dev.off()
printGraph(CC_go,result_fisher.elim,firstSigNodes = 5,fn.prefix = "tGO.CC",useInfo = "pval",pdfSW = TRUE)

## molecular function
MF_go <- new("topGOdata",nodeSize=5,ontology = "MF",allGenes = gene_list,annot = annFUN.gene2GO,gene2GO = gene_id,geneSel=top_diff_genes)
result_fisher.elim = runTest(MF_go,algorithm = "classic",statistic = "fisher")
allres = GenTable(MF_go,Fisher = result_fisher.elim,ranksOf = "classicFisher",numChar = "9999",topNodes = attributes(result_fisher.elim)$geneData[4])
filtres <- allres[allres$Fisher < 0.05,] 
write.table(filtres,file = paste(relaxInput,".MF_classic_pval005.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
filtresmf <- allres[allres$Fisher < 0.01,] 
write.table(filtresmf,file = paste(relaxInput,".MF_classic_pval001.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
pdf(paste(relaxInput,".MF_classic.pdf",sep=""))
showSigOfNodes(MF_go,score(result_fisher.elim),firstSigNodes = 5,useInfo = "all")
dev.off()
printGraph(MF_go,result_fisher.elim,firstSigNodes = 5,fn.prefix = "tGO.MF",useInfo = "pval",pdfSW = TRUE)

## biological process
BP_go <- new("topGOdata",nodeSize=5,ontology = "BP",allGenes = gene_list,annot = annFUN.gene2GO,gene2GO = gene_id,geneSel=top_diff_genes)
result_fisher.elim = runTest(BP_go,algorithm = "classic",statistic = "fisher")
allres = GenTable(BP_go,Fisher = result_fisher.elim,ranksOf = "classicFisher",numChar = "9999",topNodes = attributes(result_fisher.elim)$geneData[4])
filtres <- allres[allres$Fisher < 0.05,] 
write.table(filtres,file = paste(relaxInput,".BP_classic_pval005.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
filtresbp <- allres[allres$Fisher < 0.01,] 
write.table(filtresbp,file = paste(relaxInput,".BP_classic_pval001.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
pdf(paste(relaxInput,".BP_classic.pdf",sep=""))
showSigOfNodes(BP_go,score(result_fisher.elim),firstSigNodes = 5,useInfo = "all")
dev.off()
printGraph(BP_go,result_fisher.elim,firstSigNodes = 5,fn.prefix = "tGO.BP",useInfo = "pval",pdfSW = TRUE)

#Combine tables with BP, CC and MF
combined_table <- rbind(filtresbp, filtresmf, filtrescc)
colnames(combined_table) <- c("GO.ID", "Term", "Annotated", "Significant", "Expected", "P_value")
subtable <- combined_table[,c("GO.ID","P_value")]
colnames(subtable) <- c("% GOterm", "enrichment_P-value")
subtable <- subtable[order(subtable$`enrichment_P-value`), ]

write.table(subtable,file = paste(relaxInput,".GOterms_classic_pval001.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


# Default weight01 method
## cellular component
CC_go <- new("topGOdata",nodeSize=5,ontology = "CC",allGenes = gene_list,annot = annFUN.gene2GO,gene2GO = gene_id,geneSel=top_diff_genes)
result_fisher.elim = runTest(CC_go,algorithm = "weight01",statistic = "fisher")
allres = GenTable(CC_go,Fisher = result_fisher.elim,ranksOf = "classicFisher",numChar = "9999",topNodes = attributes(result_fisher.elim)$geneData[4])
filtres <- allres[allres$Fisher < 0.05,] 
write.table(filtres,file = paste(relaxInput,".CC_weight01_pval005.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
filtrescc <- allres[allres$Fisher < 0.01,] 
write.table(filtrescc,file = paste(relaxInput,".CC_weight01_pval001.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
pdf(paste(relaxInput,".CC_weight01.pdf",sep=""))
showSigOfNodes(CC_go,score(result_fisher.elim),firstSigNodes = 5,useInfo = "all")
dev.off()
printGraph(CC_go,result_fisher.elim,firstSigNodes = 5,fn.prefix = "tGO.CC",useInfo = "pval",pdfSW = TRUE)

## molecular function
MF_go <- new("topGOdata",nodeSize=5,ontology = "MF",allGenes = gene_list,annot = annFUN.gene2GO,gene2GO = gene_id,geneSel=top_diff_genes)
result_fisher.elim = runTest(MF_go,algorithm = "weight01",statistic = "fisher")
allres = GenTable(MF_go,Fisher = result_fisher.elim,ranksOf = "classicFisher",numChar = "9999",topNodes = attributes(result_fisher.elim)$geneData[4])
filtres <- allres[allres$Fisher < 0.05,] 
write.table(filtres,file = paste(relaxInput,".MF_weight01_pval005.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
filtresmf <- allres[allres$Fisher < 0.01,] 
write.table(filtresmf,file = paste(relaxInput,".MF_weight01_pval001.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
pdf(paste(relaxInput,".MF_weight01.pdf",sep=""))
showSigOfNodes(MF_go,score(result_fisher.elim),firstSigNodes = 5,useInfo = "all")
dev.off()
printGraph(MF_go,result_fisher.elim,firstSigNodes = 5,fn.prefix = "tGO.MF",useInfo = "pval",pdfSW = TRUE)

## biological process
BP_go <- new("topGOdata",nodeSize=5,ontology = "BP",allGenes = gene_list,annot = annFUN.gene2GO,gene2GO = gene_id,geneSel=top_diff_genes)
result_fisher.elim = runTest(BP_go,algorithm = "weight01",statistic = "fisher")
allres = GenTable(BP_go,Fisher = result_fisher.elim,ranksOf = "classicFisher",numChar = "9999",topNodes = attributes(result_fisher.elim)$geneData[4])
filtres <- allres[allres$Fisher < 0.05,] 
write.table(filtres,file = paste(relaxInput,".BP_weight01_pval005.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
filtresbp <- allres[allres$Fisher < 0.01,] 
write.table(filtresbp,file = paste(relaxInput,".BP_weight01_pval001.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
pdf(paste(relaxInput,".BP_weight01.pdf",sep=""))
showSigOfNodes(BP_go,score(result_fisher.elim),firstSigNodes = 5,useInfo = "all")
dev.off()
printGraph(BP_go,result_fisher.elim,firstSigNodes = 5,fn.prefix = "tGO.BP",useInfo = "pval",pdfSW = TRUE)

#Combine tables with BP, CC and MF
combined_table <- rbind(filtresbp, filtresmf, filtrescc)
colnames(combined_table) <- c("GO.ID", "Term", "Annotated", "Significant", "Expected", "P_value")
subtable <- combined_table[,c("GO.ID","P_value")]
colnames(subtable) <- c("% GOterm", "enrichment_P-value")
subtable <- subtable[order(subtable$`enrichment_P-value`), ]

write.table(subtable,file = paste(relaxInput,".GOterms_weight01_pval001.txt",sep=""),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
