library(tidyverse)
library(igraph)
library(ggraph)


Gyne_data<-read_table("./Gyne_abundance_200samples.csv", col_types = cols())
Worker_data<-read_table("./Worker_abundance_170samples.csv", col_types = cols())
Metadata <- read_table("./Gyne_worker_SRA_info.1.csv",col_types = cols())
module_input <- read_table("./module4.txt.GN.DEG.lst",col_types = cols())

Gyne_data_wide <- Gyne_data %>% pivot_longer(cols = !gene_ID,names_to = "library",values_to = "tpm") %>% mutate(logTPM=log2(tpm+1))
Worker_data_wide <- Worker_data %>% pivot_longer(cols = !gene_ID,names_to = "library",values_to = "tpm") %>% mutate(logTPM=log2(tpm+1))
Gyne_data_wide %>% full_join(Metadata %>% select(Run, dev_stage, Sample), by = c("library"="Run")) -> Gyne_data_log_wide
Worker_data_wide %>% full_join(Metadata %>% select(Run, dev_stage, Sample), by = c("library"="Run")) -> Worker_data_log_wide
Gyne_data_log_wide %>%  group_by(gene_ID, dev_stage) %>% summarise(mean.logTPM = mean(logTPM)) %>% ungroup() -> Gyne_data_wide_averaged
Worker_data_log_wide %>%  group_by(gene_ID, dev_stage) %>% summarise(mean.logTPM = mean(logTPM)) %>% ungroup() -> Worker_data_wide_averaged
Worker_data_wide_averaged %>% mutate(sample=str_c("Worker.",dev_stage)) -> Worker_data_wide_averaged_add
Gyne_data_wide_averaged %>% mutate(sample=str_c("Gyne.",dev_stage)) -> Gyne_data_wide_averaged_add
Worker_data_wide_averaged_add <- Worker_data_wide_averaged_add %>% dplyr::rename(dev_stage_ori=dev_stage,dev_stage=sample) %>% select(gene_ID,dev_stage,mean.logTPM)
Gyne_data_wide_averaged_add <- Gyne_data_wide_averaged_add %>% dplyr::rename(dev_stage_ori=dev_stage,dev_stage=sample) %>% select(gene_ID,dev_stage,mean.logTPM)

all_data_averaged<-rbind(Gyne_data_wide_averaged_add,Worker_data_wide_averaged_add)
all_data_averaged_wide <- all_data_averaged %>% select(gene_ID, dev_stage, mean.logTPM) %>% pivot_wider(names_from = dev_stage, values_from = mean.logTPM) %>% as.data.frame()
all_data_averaged_wide <- all_data_averaged_wide %>% filter(gene_ID!="NA")
row.names(all_data_averaged_wide) <- all_data_averaged_wide$gene_ID
all_data_averaged_wide <- all_data_averaged_wide %>% select(-gene_ID)

cor_matrix <- cor(t(all_data_averaged_wide[, -1]))
cor_matrix_upper_tri <- cor_matrix
cor_matrix_upper_tri[lower.tri(cor_matrix_upper_tri)] <- NA
number_of_tissue_stage = 12 
edge_table <- cor_matrix_upper_tri %>% as.data.frame() %>% mutate(from = row.names(cor_matrix)) %>% pivot_longer(cols = !from, names_to = "to", values_to = "r") %>% filter(is.na(r) == F) %>% filter(from != to) %>% mutate(t = r*sqrt((number_of_tissue_stage-2)/(1-r^2))) %>% mutate(p.value = case_when(t > 0 ~ pt(t, df = number_of_tissue_stage-2, lower.tail = F),t <=0 ~ pt(t, df = number_of_tissue_stage-2, lower.tail = T)))  %>% mutate(FDR = p.adjust(p.value, method = "fdr"))
edge_table_filter <- edge_table %>% filter(FDR < 0.01)

edge_table_input <- edge_table_filter %>% filter(from %in% module_input$Gene & to %in% module_input$Gene)
node_table <- data.frame(gene_ID = c(edge_table_input$from, edge_table_input$to) %>% unique())
node_table <- node_table %>% left_join(module_input,by=c(gene_ID="Gene")) %>% select(gene_ID,Name,DEG)

my_network <- graph_from_data_frame(edge_table_input,vertices = node_table,directed = F)
my_network %>% ggraph(layout = "kk") + geom_edge_arc(strength = 0.2, width = 0.25, alpha = 0.2) + geom_node_point(alpha = 0.8, color = "white", size=6,shape = 21,aes(fill = DEG)) + scale_fill_manual(values=c("#EA8379","#299D8F")) +  geom_node_text(aes(label = Name),size=4, color="black",repel = TRUE) + theme_void() + theme(text = element_text(size = 14),legend.position = "bottom",legend.justification = 1,title = element_text(size = 12))

ggsave("./Module4.network.svg", height = 8, width = 6, bg = "white")



