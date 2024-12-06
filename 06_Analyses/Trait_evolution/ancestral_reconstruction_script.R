library(treeio)
library(ggtree)
library(phytools)
library(ggplot2)
library(ggnewscale)
library(wesanderson)
library(ggimage)
library(geiger)

setwd("/Users/joel/GAGAproject/Ancestral_reconstruction/")

### READ TREE
treeu <- read.tree("GAGA_dated_phylogeny_newick.tre")

## READ TRAITS TABLE
traits <- read.csv("Traits_table_v6_18Dec_toPGLS_reduced.tsv",sep="\t",dec = ".", quote = "")

# READ SUBFAMILIES & NODES TABLE: to plot the subfamily ranges with the corresponding colors from ITOL
subfamilies <- read.csv("subfamilies_nodes.csv",sep="\t",dec = ",")
subfamilies_colors <- c("#04aa08", "#c80a0a","#020bf3","#e77a01","#1cc2c4","#60c304","#03783f","#d2ae04","#4505d4","#6d3c12","#f41e23","#069ffe")
names(subfamilies_colors) <- subfamilies$Subfamily

### Phylogenetic signal
traittotest <- traits$Colony.size.cat
names(traittotest) <- traits$GAGA.ID
phylosig(treeu, traittotest,  method = "lambda", test = TRUE)

lambda.traittotest<-fitContinuous(phy = treeu, dat = traittotest, model = "lambda")
lambda.traittotest

#################### DISCRETE CHARACTERS: Ancestral Character Estimation ##################################################
#The model of trait evolution used for each trait is determined by comparing the log likelihoods of alternative models

TRAIT <- traits$Colony.size.cat # Select here the column with trait of interest to conduct the ancestral reconstruction
names(TRAIT) <- paste(traits$GAGA.ID, traits$Species.name) # Retrieve GAGA ID for each trait and the species name
TRAIT<- replace(TRAIT, TRAIT=="", NA)

unique(TRAIT) # See the categories for the trait
numtraits <- 4 # Write here the number of traits (do not count NA)
traitname <- "Colony size" # Name for the trait

TRAIT_colors <- as.vector(wes_palette("Zissou1", numtraits, type = "continuous")) # colors for the plot

# Visualize the color chosen
plot(NULL, xlim=c(0,length(TRAIT_colors)), ylim=c(0,1), 
     xlab="", ylab="", xaxt="n", yaxt="n")
rect(0:(length(TRAIT_colors)-1), 0, 1:length(TRAIT_colors), 1, col=TRAIT_colors)


names(TRAIT_colors) <- c("0", "1", "2", "3") # Name of the categories


# rename tree labels as "GAGA-ID species name"
for (name in 1:length(treeu$tip.label)) {
  treeu$tip.label[name] <- paste(treeu$tip.label[name],
                                      traits$Species.name[which(traits$GAGA.ID == treeu$tip.label[name])])
}


# Ancestral character estimation, testing models
fitER<-ace(TRAIT, treeu, model="ER", type="discrete") # Ancestral Character Estimation: ER model
fitARD<-ace(TRAIT, treeu, model="ARD", type="discrete") # Ancestral Character Estimation: ARD model
fitSYM<-ace(TRAIT, treeu, model="SYM", type="discrete") # Ancestral Character Estimation: SYM model

# Using other alternative models to test more relevant biological hypotheses for specific trait catgorizations
fit2single<-ace(TRAIT, treeu, type="discrete", model=matrix(c(0, 0, 1, 0), 2)) # Ancestral Character Estimation
fit2singlerev<-ace(TRAIT, treeu, type="discrete", model=matrix(c(0, 1, 0, 0), 2), marginal = TRUE, use.expm = FALSE) # Ancestral Character Estimation
fitpp<-ace(TRAIT, treeu, type="discrete", model=matrix(c(0, 4, 0, 0,
                                                         4, 0, 5, 0,
                                                         0, 3, 0, 0,
                                                         0, 2, 1, 0), 4), marginal = TRUE, use.expm = FALSE) # Ancestral Character Estimation
fitpp<-ace(TRAIT, treeu, type="discrete", model=matrix(c(0, 1, 2, 3, 4, 5, 0,
                                                         1, 0, 7, 8, 9, 10, 0,
                                                         2, 7, 0, 12, 13, 14, 0,
                                                         3, 8, 12, 0, 16, 17, 0,
                                                         4, 9, 13, 16, 0, 19, 0,
                                                         5, 10, 14, 17, 19, 0, 0, 
                                                         6, 11, 15, 18, 20, 21, 0), 7), marginal = TRUE, use.expm = FALSE) # Ancestral Character Estimation
fitpp<-ace(TRAIT, treeu, type="discrete", model=matrix(c(0, 1, 2,
                                                         1, 0, 3,
                                                         0, 0, 0), 3), marginal = TRUE, use.expm = FALSE) # Ancestral Character Estimation
fitpp2<-ace(TRAIT, treeu, type="discrete", model=matrix(c(0, 1, 2,
                                                         0, 0, 3,
                                                         0, 0, 0), 3), marginal = TRUE, use.expm = FALSE) # Ancestral Character Estimation

fitmanc<-ace(TRAIT, treeu, type="discrete", model=matrix(c(0, 0, 0,
                                                           1, 0, 3,
                                                           2, 3, 0), 3), marginal = TRUE, use.expm = FALSE) # Ancestral Character Estimation
fitmancard<-ace(TRAIT, treeu, type="discrete", model=matrix(c(0, 4, 0,
                                                              1, 0, 5,
                                                              2, 3, 0), 3), marginal = TRUE, use.expm = FALSE) # Ancestral Character Estimation
fitmancard2<-ace(TRAIT, treeu, type="discrete", model=matrix(c(0, 4, 0,
                                                              1, 0, 0,
                                                              2, 3, 0), 3), marginal = TRUE, use.expm = FALSE) # Ancestral Character Estimation
fitmancsym<-ace(TRAIT, treeu, type="discrete", model=matrix(c(0, 1, 0,
                                                              1, 0, 0,
                                                              2, 3, 0), 3), marginal = TRUE, use.expm = FALSE) # Ancestral Character Estimation
fitmanc2<-ace(TRAIT, treeu, type="discrete", model=matrix(c(0, 0, 0,
                                                            1, 0, 0,
                                                            2, 3, 0), 3), marginal = TRUE, use.expm = FALSE) # Ancestral Character Estimation
fitpp6<-ace(TRAIT, treeu, type="discrete", model=matrix(c(0, 0, 0, 0, 0, 0,
                                                         1, 0, 6, 7, 8, 9,
                                                         2, 6, 0, 10, 11, 12,
                                                         3, 7, 10, 0, 13, 14,
                                                         4, 8, 11, 13, 0, 15,
                                                         5, 9, 12, 14, 15, 0), 6), marginal = TRUE, use.expm = FALSE) # Ancestral Character Estimation

# Checking likelihoods and AIC's
#fitER$loglik
AIC(fitER)
#fitARD$loglik
AIC(fitARD)
#fitSYM$loglik
AIC(fitSYM)
#fit2single$loglik
AIC(fit2single)
#fit2ard$loglik
AIC(fit2ard)
#fitpp$loglik
AIC(fitpp)
AIC(fitmancsym)
AIC(fitmancard2)

# Type here the best retrieved model (lowest AIC)
best_model <- fitER # type here the best model

# Creating data frame with the reconstructed states for the best model
mas <- as.data.frame(round(best_model$lik.anc,2))
#mas <- as.data.frame(round(best_model$lik.anc, digits = 0))
mas$node <- row.names(mas)

# Plotting the tree
p1 <-ggtree(treeu, size = 0.2) +
  xlim(0,300) +
  geom_tiplab(offset=10,align= TRUE, size=0.9, linesize = 0.2)+
  #  geom_hilight(data=subfamilies, mapping = aes(node=node,fill=Subfamily, ), align = "right", type = "gradient", alpha=0.3)+
  geom_hilight(data=subfamilies, mapping = aes(node=node,fill=Subfamily, ), align = "right", type = "rect", alpha=0.2)+
  scale_fill_manual(values=subfamilies_colors) + guides(fill = guide_legend(override.aes = list(alpha = 1)))+  
  scale_y_reverse()+ new_scale_fill()
p1

# add heatmap to the tree (column at the right of the tree with the trait values)
p2 <- gheatmap(p1, data.frame(factor(TRAIT)), width=0.02, low="#D95F02",
               high="#7570B3", colnames = FALSE, font.size=2, color="white", offset=4) +
  scale_fill_manual(values=TRAIT_colors, labels = names(TRAIT_colors), name = traitname, na.value = NA)+
  scale_y_reverse()
p2

# Add the pie charts with the ancestral reconstructions
pies <- nodepie(mas, cols=1:numtraits, color=TRAIT_colors, alpha=0.8)
#p12 <- inset(p2, pies, width = 0.06, height = 0.06, reverse_y = TRUE, hjust = 0.005, vjust = -0.2)
#p12 <- inset(p2, pies, width = 0.06, height = 0.06, reverse_y = TRUE, hjust = 1, vjust = -0.4)
p12 <- inset(p2, pies, width = 0.09, height = 0.09, reverse_y = TRUE, hjust = 1, vjust = -0.4) # Mac updated to R 4.2.2- change the scale

p12

# Save the tree into a file
ggsave(
  "TraitNNN_ult_ancestral_fitER_discrete.pdf",
  device = "pdf",
  plot = last_plot(),
  height = 6.33,
  width = 5.14,
  units = "in"
)

dev.off()



#################### ANCESTRAL RECONSTRUCTION FOR CONTINUOUS CHARACTERS ####################################################################
### Note: the continuous trait ancestral reconstruction does not accepted missing values or NA
### Then, if the dataset contains missing data, those species need to be excluded from the table and tree, here the commands to do it:
## Subset tree to only contain species without NA (here as example the trait column MATmean)
#tree.subset<-treeio::drop.tip(treeu@phylo,tip=traits$GAGA.ID[is.na(traits$MATmean)])
## trim trait data to same set of species
#TRAIT<-traits$MATmean[!is.na(traits$MATmean)]
#names(TRAIT)<-traits$GAGA.ID[!is.na(traits$MATmean)]

setwd("/Users/joel/GAGAproject/Ancestral_reconstruction/")
treeu <- read.tree("GAGA_dated_phylogeny_newick.tre")
traits <- read.csv("Traits_table_v6_18Dec_toPGLS_reduced.tsv",sep="\t",dec = ".", quote = "")
subfamilies <- read.csv("subfamilies_nodes.csv",sep="\t",dec = ",")
subfamilies_colors <- c("#04aa08", "#c80a0a","#020bf3","#e77a01","#1cc2c4","#60c304","#03783f","#d2ae04","#4505d4","#6d3c12","#f41e23","#069ffe")
names(subfamilies_colors) <- subfamilies$Subfamily


TRAIT_colors <- as.vector(wes_palette("Zissou1", 5, type = "continuous")) # change by number of categories/colours to plot the trait

TRAIT <- traits$Colony.size.log.max # select column with CONTINUOUS trait of interest
#TRAIT <- as.numeric(gsub(",","",traits$Genome.size)) # select column with CONTINUOUS trait of interest
traitname <- "Colony"
names(TRAIT) <-traits$GAGA.ID # obtain GAGA ID for each trait


## Ancestral reconstruction for the continuous trait
fit<-fastAnc(treeu, TRAIT) # fast estimation of ML ancestral states
# Parsing the results
traitvector <- data.frame(node = ggtree::nodeid(treeu, names(TRAIT)), trait = TRAIT)
traitfit <- data.frame(node = names(fit), trait = as.vector(fit))
traitfit2 <- rbind(traitvector, traitfit)
traitfit2$node<-as.numeric(traitfit2$node)
tree2 <- as.treedata(tidytree::full_join(as_tibble(treeu), traitfit2, by = 'node'))

# rename tree labels in format "GAGAID species name"
for (name in 1:length(tree2@phylo$tip.label)) {
  tree2@phylo$tip.label[name] <- paste(tree2@phylo$tip.label[name],
                                       traits$Species.name[which(traits$GAGA.ID == tree2@phylo$tip.label[name])])
}

# Plot the tree
p1 <-  ggtree(tree2,ladderize = T, aes(color=trait), continuous = "color", size=0.3,
                                     lineend='square')+  scale_y_reverse() +
  geom_tiplab(offset=1,align = TRUE ,color = "black", size = 1, linesize = 0.2)+
  scale_color_gradientn(colours=TRAIT_colors, name = traitname)+ 
#  scale_color_gradientn(colours=TRAIT_colors, labels = names(TRAIT_colors), name = "Colony size")+ 
  coord_cartesian(clip = 'off') + 
  xlim(0, 300)+ theme(legend.position = "left")

# Add subfamily to the right of the tree
#for (i in 1:length(subfamilies_colors)) { # Changed offset here to use in my home mac
#  p1 <- p1 + geom_cladelabel(align = TRUE,subfamilies$node[i], names(subfamilies_colors[i]), offset=80, barsize=0.5,
#                             angle=0, offset.text=1, hjust=0, fontsize=1, color = subfamilies_colors[i])
#}
for (i in 1:length(subfamilies_colors)) {
  p1 <- p1 + geom_cladelabel(align = TRUE,subfamilies$node[i], names(subfamilies_colors[i]), offset=100, barsize=0.5,
                             angle=0, offset.text=1, hjust=0, fontsize=1, color = subfamilies_colors[i])
}
p1

ggsave(
  "Colony_size_log_max_ancestral_continuous.pdf",
  device = "pdf",
  plot = last_plot(),
  height = 6.33,
  width = 5.14,
  units = "in"
)




