library(tidyverse)
library(RIdeogram)

colors=data.frame(row.names = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21"),col=c("699ECA","FF8C00","F898CB","4DAF4A","D65190","731A73","FFCB5B","E87B1E","0076B9","3D505A","0098B2","1F78B4","33A02C","E31A1C","FFFF99","B15928","D95F02","7570B3","E7298A","66A61E","6A3D9A"))

k_target_ord=c("18","1","17","3","2","4","5","6","15","11","12","16","10","14","13","19","7","8","9") ### chromosomes of target genomes
k_query_ord=c("1","5","6","2","9","8","7","10","4","11","3")  ### chromosomes of query genomes

lgall=read.table("input.lst.all",h=T)
lg=read.table("input.lst",h=T)

lgall %>%  group_by(s1chr) %>%   summarise(End=max(End_1)) %>% arrange(-End) %>% mutate(Start=1,species='GAGA-0003',fill='GAGA-0003',size=12,color=252525) %>% rename(Chr=s1chr) %>% relocate(Start,.after=Chr)-> k_target
lgall %>%  group_by(s2chr) %>%   summarise(End=max(End_2)) %>% arrange(-End) %>% mutate(Start=1,species='GAGA-0527',fill='GAGA-0527',size=12,color=252525) %>% rename(Chr=s2chr) %>% relocate(Start,.after=Chr)-> k_query
k_target=k_target[order(match(k_target$Chr,k_target_ord)),]
k_query=k_query[order(match(k_query$Chr,k_query_ord)),]

lg$Species_1=match(lg$s1chr,k_target$Chr)
lg$Species_2=match(lg$s2chr,k_query$Chr)
rbind(k_target,k_query) -> chroms

lg %>% select(s1chr,Species_1,Start_1,End_1,Species_2,Start_2,End_2) %>% mutate(fill=colors[match(s1chr,rownames(colors)),]) %>%  relocate(End_1,.after=Start_1) %>% relocate(End_2,.after=Start_2)  %>% select(-s1chr) -> lgt

ideogram(karyotype = chroms, synteny = lgt,output="GAGA-0003_GAGA-0527.all.svg")
convertSVG("GAGA-0003_GAGA-0527.all.svg", file="GAGA-0003_GAGA-0527.all",device = "png")
