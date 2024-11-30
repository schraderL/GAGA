library(tidyverse)
library(RColorBrewer)
library(ggh4x)




Gyne_dat<-read_table("./Gyne_abundance_200samples.csv", col_types = cols())
Worker_dat<-read_table("./Worker_abundance_170samples.csv", col_types = cols())
Metadata <- read_table("./Gyne_worker_SRA_info.1.csv",col_types = cols())


Gyne_data_wide <- Gyne_dat %>% pivot_longer(cols = !gene_ID,names_to = "library",values_to = "tpm") %>% mutate(logTPM=log2(tpm+1))
Worker_data_wide <- Worker_dat %>% pivot_longer(cols = !gene_ID,names_to = "library",values_to = "tpm") %>% mutate(logTPM=log2(tpm+1))
Gyne_data_wide %>% left_join(Metadata %>% select(Run, dev_stage, Sample), by = c("library"="Run")) -> Gyne_data_wide
Worker_data_wide %>% left_join(Metadata %>% select(Run, dev_stage, Sample), by = c("library"="Run")) -> Worker_data_wide


Gyne_data_wide_add <- Gyne_data_wide %>% mutate(tissue=case_when(str_detect(Sample,"Gyne")~"Gyne",str_detect(Sample,"Worker")~"Worker"))
Worker_data_wide_add <- Worker_data_wide %>% mutate(tissue=case_when(str_detect(Sample,"Gyne")~"Gyne",str_detect(Sample,"Worker")~"Worker"))


Exp_table_wide <- rbind(Gyne_data_wide_add,Worker_data_wide_add)
my_genes<-data.frame(gene_ID=c("Mpha_g13810","Mpha_g13811","Mpha_g13812","Mpha_g13815"),symbol=c("CG8745","Vitellogenin-2","Vitellogenin-1","CG18063"))
genes_TPM <- Exp_table_wide %>% filter(gene_ID %in% my_genes$gene_ID)
genes_TPM <- genes_TPM %>% left_join(my_genes,by="gene_ID")


df<-genes_TPM %>%
mutate(order_x = case_when(
            str_detect(dev_stage, "2nd") ~ 1,
            str_detect(dev_stage, "3rd") ~ 2,
            str_detect(dev_stage, "Pre.Pupa") ~ 3,
            str_detect(dev_stage, "Pupa.Young") ~ 4,
            str_detect(dev_stage, "Pupa.Old") ~ 5,
            str_detect(dev_stage, "Imago") ~ 6
            )) %>%
mutate(dev_stage = reorder(dev_stage, order_x)) %>% 
mutate(symbol = fct_relevel(symbol, "CG8745", "Vitellogenin-2", "Vitellogenin-1", "CG18063")) %>% 
mutate(tag = symbol)  


scales <- list(
  # Here you have to specify all the scales, one for each facet row in your case
  scale_y_continuous(limits = c(0,6)),
  scale_y_continuous(limits = c(0, 3)),
  scale_y_continuous(limits = c(0, 4)),
  scale_y_continuous(limits = c(0, 3)))
  
p.main <- df %>% ggplot(aes(fill=tissue,x = as.numeric(dev_stage), y = logTPM,group=tissue,color=tissue)) +
facet_grid(symbol~ ., scales = "free_y") +
  geom_rect(xmin=-3,xmax=10,ymin=0,ymax=Inf,fill="grey98",alpha=1,color=NA) +
  geom_point(size=.4,alpha=.3,position=position_dodge(width=0.1))+
  geom_vline(xintercept = 1:6,color="grey90",alpha=.8) +
  geom_smooth(method = 'loess',size=0.8,se = T,show.legend = FALSE,linetype = "solid",alpha=.1) +
  stat_summary(
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median,
    pch=21,
    col="black",
    position=position_dodge(width=0.1))+
scale_fill_manual(values=c("red2","steelblue3")) +
  scale_color_manual(values=c("red2","steelblue3")) +
  labs(x = NULL,y = "log2(TPM)") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.spacing = unit(0, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    strip.background = element_blank()
  ) 

p.supplement <- p.main +
  facetted_pos_scales(y = scales)

p.main
ggsave("./syntenic_gene.Vitellogenin.svg", height = 6, width = 6, bg = "white")
p.supplement
ggsave("./syntenic_gene.Vitellogenin.sup.svg", height = 6, width = 6, bg = "white")
