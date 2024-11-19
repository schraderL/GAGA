setwd("D:/Project/项目汇报/20231207/differential_co-expression/")
library(WGCNA)
library(tidyverse)

Gyne_data<-read_table("./Gyne_abundance_200samples.csv", col_types = cols())
Worker_data<-read_table("./Worker_abundance_170samples.csv", col_types = cols())

Metadata <- read_table("./Gyne_worker_SRA_info.1.csv",col_types = cols())

Gyne_data_wide <- Gyne_data %>% pivot_longer(cols = !gene_ID,names_to = "library",values_to = "tpm") %>% mutate(logTPM=log2(tpm+1))
Worker_data_wide <- Worker_data %>% pivot_longer(cols = !gene_ID,names_to = "library",values_to = "tpm") %>% mutate(logTPM=log2(tpm+1))

Gyne_data_wide %>% full_join(Metadata %>% select(Run, dev_stage, Sample), by = c("library"="Run")) -> Gyne_data_log_wide
Worker_data_wide %>% full_join(Metadata %>% select(Run, dev_stage, Sample), by = c("library"="Run")) -> Worker_data_log_wide

Gyne_data_log_wide %>%  group_by(gene_ID, dev_stage) %>% summarise(mean.logTPM = mean(logTPM)) %>% ungroup() -> Gyne_data_wide_averaged
Worker_data_log_wide %>%  group_by(gene_ID, dev_stage) %>% summarise(mean.logTPM = mean(logTPM)) %>% ungroup() -> Worker_data_wide_averaged

Gyne_data_wide_averaged %>% mutate(sample=str_c("Gyne",dev_stage)) -> Gyne_data_wide_averaged_add
Worker_data_wide_averaged %>% mutate(sample=str_c("Worker",dev_stage)) -> Worker_data_wide_averaged_add

Worker_data_wide_averaged %>% mutate(sample=str_c("Worker.",dev_stage)) -> Worker_data_wide_averaged_add
Gyne_data_wide_averaged %>% mutate(sample=str_c("Gyne.",dev_stage)) -> Gyne_data_wide_averaged_add

Worker_data_wide_averaged_add <- Worker_data_wide_averaged_add %>% dplyr::rename(dev_stage_ori=dev_stage,dev_stage=sample) %>% select(gene_ID,dev_stage,mean.logTPM)
Gyne_data_wide_averaged_add <- Gyne_data_wide_averaged_add %>% dplyr::rename(dev_stage_ori=dev_stage,dev_stage=sample) %>% select(gene_ID,dev_stage,mean.logTPM)

all_data_averaged<-rbind(Gyne_data_wide_averaged_add,Worker_data_wide_averaged_add)

all_data_averaged_wide <- all_data_averaged %>% select(gene_ID, dev_stage, mean.logTPM) %>% pivot_wider(names_from = dev_stage, values_from = mean.logTPM) %>% as.data.frame()

all_data_averaged_wide <- all_data_averaged_wide %>% filter(gene_ID!="NA")
row.names(all_data_averaged_wide) <- all_data_averaged_wide$gene_ID
all_data_averaged_wide <- all_data_averaged_wide %>% select(-gene_ID)

all_data_averaged_wide1<-t(all_data_averaged_wide)

powers=18
exp_tpm <- all_data_averaged_wide1
datExpr<-exp_tpm
nGenes=ncol(datExpr)
nSamples=nrow(datExpr)

sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")

plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],labels = powers, cex = 1.2, col = "red")

net <- blockwiseModules(datExpr,power = powers,maxBlockSize = ncol(datExpr),corType = "bicor",networkType = "signed",TOMType = "signed", minModuleSize = 30, reassignThreshold=0, mergeCutHeight = 0.25, numericLabels = TRUE, saveTOMs = F,verbose = 3)

moduleColors1 <- labels2colors(net$colors)
table(moduleColors1)

MES0 <- moduleEigengenes(datExpr,moduleColors1)$eigengenes
MEs = orderMEs(MES0)
cp<-corAndPvalue(MEs,design,alternative ="greater")

