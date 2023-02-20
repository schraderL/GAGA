# GAGA Gene family evolution

We used [CAFE](https://github.com/hahnlab/CAFE5) to analyze changes in gene family size across the ant phylogeny. This program uses a birth and death process to model gene gain and loss across the phylogenetic tree. Briefly, we used CAFE to estimate gene family evolutionary rates, both global (the whole phylogeny share the same _birth-and-death_ rate) and local (we tested different scenarios in which we model different _birth-and-death_ rates in the phylogeny, i.e.: one rate for poneroids and other for formicoids). In addition, we infered gene family counts for each orthologous group in all internal nodes of the phylogeny, allowing us to retreive gene families that have expanded or contracted in specific lineages. 


### CAFE input
First, we use as input in CAFE the orthogroup table that we retrieved in the steps detailed in [Orthology](../Orthology). This table contains the gene count of each infered orthogroup for each species, including our 163 ant genomes and the 8 outgroups from Apoidea. The format of the table should be (example with three genomes):
```
Desc    HOG     GAGA-0001   GAGA-0003   GAGA-0004   
HOG0000754  HOG0000754  5   6   4  
HOG0000866  HOG0000866  3   4   4  
```

In addition, we also use the calibrated phylogeny that contains the phylogenetic relationships and dates retrieved for the 163 ant species sequenced in GAGA, see [Phylogeny](../Phylogeny)for more details on these steps. For the 8 Apoidea species used as outgroups, we reconstructed the phylogeny using these outgroups and single copy ortholog groups, including the 163 ant genomes, and used the dates from [Peters et al. 2017](http://dx.doi.org/10.1016/j.cub.2017.01.027). 

Note that we used the error model implemented in CAFE, which accounts for non-biological factors (e.g., genome sequencing and coverage differences, gene family clustering errors, etc.) leading to incorrect gene family counts in input files. However, we opted to conduct the CAFE analyses both keeping and discarding all species with low genome contiguity (Scaffold N50 < 500 Kb), as the gene family numbers might be biased due to the fragmentation of the assembly, resulting for instance in multiple partial genes or missed copies. And therefore we evaluated both approaches. 
- We removed the columns for the species with low contiguity genomes in the orthogroup gene count table. 
- We pruned the dated species tree to generate a tree without the low contiguity genome species using ***nw_prune*** from [Newick Utilities](https://bio.tools/newick_utilities). 
```
nw_prune GAGA_dated_phylogeny_newick_woutgroup.tre GAGA-0002 GAGA-0055 GAGA-0220 GAGA-0234 GAGA-0331 GAGA-0366 GAGA-0386 GAGA-0387 GAGA-0389 GAGA-0392 GAGA-0406 GAGA-0408 GAGA-0495 GAGA-0535 GAGA-0536 GAGA-0537 GAGA-0553 GAGA-0554 GAGA-0577 > GAGA_dated_phylogeny_newick_woutgroup_nostLFR.tre
```

### Running CAFE
We used the script [run_cafe_analyses.sh](run_cafe_analyses.sh) to run the CAFE analyses. The script generates a number of independent runs in order to ensure convergence. The script contains the code to run all the following steps that we conducted, just comment/uncomment the lines to run each step. 

- First, we used


### Generating CAFE results


Folder with the scripts used in GAGA
 We use CAFE, using the N0.tsv gene count table retrieved in orthofinder. But we excluded all stLFR short read assemblies given the low accuracy and incompleteness in the annotiation .
- Add scripts to filter out species from the species tree (for instance the stLFR)

Script blabla

1 - run all to see
2 - filter HOG with high variation
3 - run cafe 
4 - run error model
5 fix lambda an re run with all orthogroups


