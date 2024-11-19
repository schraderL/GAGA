remotes::install_github("flr/ggplotFL")
library(ggplot2)
library(tidyverse)
library(ggplotFL)

setwd("/Users/lukas/sciebo/Projects/GAGA/mainpaper/Figure2/data/input/")
fold_data<-read_table("./log2fd.lst",col_types = cols())
module<-read_table("./4module.txt",col_types=cols())

fold_data_select <- fold_data %>% filter(gene_ID %in% module$gene_ID)

fold_data_module <- fold_data_select %>% left_join(module,by="gene_ID")

fold_data_module_mean <- fold_data_module  %>%
     group_by(module, dev_stage) %>%
     summarise(mean.log2fd = mean(log2fd)) %>%
     ungroup()

df <- fold_data_module  %>%
   mutate(order_x = case_when(
   str_detect(dev_stage, "2nd") ~ 1,
   str_detect(dev_stage, "3rd") ~ 2,
   str_detect(dev_stage, "Pre.Pupa") ~ 3,
   str_detect(dev_stage, "Pupa.Young") ~ 4,
   str_detect(dev_stage, "Pupa.Old") ~ 5,
   str_detect(dev_stage, "Imago") ~ 6,
)) %>%
mutate(dev_stage = reorder(dev_stage, order_x)) 




#df<-df %>% group_by(dev_stage,module) %>% mutate(q_stage=quantile(log2fd,0.25), med_stage=median(log2fd))

df %>% ggplot(aes(x = dev_stage, y = log2fd )) +
  coord_cartesian(ylim = c(-3, 3)) +
  geom_rect(xmin=-3,xmax=10,ymin=0,ymax=Inf,fill="ghostwhite",alpha=1) +
  geom_rect(xmin=-3,xmax=10,ymin=0,ymax=-Inf,fill="grey90",alpha=1)+
  facet_wrap(~module,strip.position = "right",ncol=1) +
  geom_line(aes(group = gene_ID), alpha = 0.1, color = "grey60") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_flquantiles(aes(group=module),fill="steelblue1",alpha=0.3,probs = c(0.1, 0.5, 0.9))+
  
  geom_boxplot(aes(group = dev_stage),fill="white", alpha = 1, color = "black",width=.1,outlier.colour = NA,lwd=.5) + 
  labs(x = NULL,
       y = "log2foldchange") +
  theme_light() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(0.5, "line")
  )

ggsave("./4module.dev_stage.log2fd.svg", height = 6, width = 4, bg = "white")
