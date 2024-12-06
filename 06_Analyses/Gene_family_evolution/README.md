# Gene family evolution

We used [CAFE](https://github.com/hahnlab/CAFE5) to analyze changes in gene family size across the ant phylogeny. This program uses a birth and death process to model gene gain and loss across the phylogenetic tree. Briefly, we used CAFE to estimate gene family evolutionary rates, both global (the whole phylogeny share the same _birth-and-death_ rate) and local (we tested different scenarios in which we model different _birth-and-death_ rates in the phylogeny, i.e.: one rate for the ants and another for the outgroups). In addition, we infered gene family counts for each orthologous group in all internal nodes of the phylogeny, allowing us to retrieve gene families that have expanded or contracted in specific lineages. 


### CAFE input

We use as input in CAFE the orthogroup table that we detail in the [Orthology assessment](../Orthology). This table contains the gene count of each infered orthogroup for each species, including our 163 ant genomes and the 8 outgroups from the Apoidea. The format of the table should be (example with three genomes):
```
Desc    HOG     GAGA-0001   GAGA-0003   GAGA-0004   
HOG0000754  HOG0000754  5   6   4  
HOG0000866  HOG0000866  3   4   4  
```

Note that we used the error model implemented in CAFE, which accounts for non-biological factors (e.g., genome sequencing and coverage differences, gene family clustering errors, etc.) leading to incorrect gene family counts. However, we opted to conduct the CAFE analyses removing all species with low genome contiguity (Scaffold N50 < 500 Kb), as the gene family numbers might be biased due to the fragmentation of the assembly, resulting for instance in multiple partial genes or missed copies.

- We removed the columns for the species with low contiguity genomes in the orthogroup gene count table: [N0_GeneCounts_nostLFR.tsv file](N0_GeneCounts_nostLFR.tsv.zip). 
- We pruned the dated species tree to generate a tree without the low contiguity genome species using ***nw_prune*** from [Newick Utilities](https://bio.tools/newick_utilities): [GAGA_dated_phylogeny_newick_woutgroup_nostLFR.tre file](GAGA_dated_phylogeny_newick_woutgroup_nostLFR.tre).

```
nw_prune GAGA_dated_phylogeny_newick_woutgroup.tre GAGA-0002 GAGA-0055 GAGA-0220 GAGA-0234 GAGA-0331 GAGA-0366 GAGA-0386 GAGA-0387 GAGA-0389 GAGA-0392 GAGA-0406 GAGA-0408 GAGA-0495 GAGA-0535 GAGA-0536 GAGA-0537 GAGA-0553 GAGA-0554 GAGA-0577 > GAGA_dated_phylogeny_newick_woutgroup_nostLFR.tre
```

### Running CAFE
We used the script [run_cafe_analyses.sh](run_cafe_analyses.sh) to run the CAFE analyses. The script generates a number of independent runs in order to ensure convergence. It contains the code to run all the following steps that we conducted, just comment/uncomment the respective lines to run each step. 

- I - First, we used all orthogroups to run CAFE with a global lambda model and using a poisson distribution for the root frequency distribution. However, given the number of orthogroups with a high variable number of genes across species, and the large number of species in the phylogeny, the program did not converge. 

- II - Then, we filtered gene families with large gene copy number variance, as these can cause parameter estimates to be non-informative and prevent convergence in CAFE. We tried different values, finally filtering gene families for which the difference between minimum and maximum gene count was higher than 20. This number can be modified and adjusted in the script used for filtering [filter_genefams_forcafe.pl](filter_genefams_forcafe.pl), which also excludes gene families containing TE-related genes (based on the functional annotations). Then, we run CAFE using the ortholog table filtering these orthogroups and ensured that CAFE converged. We checked and evaluated the Maximum-likelihood and estimated lambda(s) of each run using the following command:
```
grep 'L' CAFE*/*/Base_results.txt > Results_lnL_lambda.txt
```

- III - Third, we use the error model in CAFE to estimate a global lambda accounting for non-biological factors. The estimated lambda using the error model is slightly lower, which is in accordance with the fact that errors in genome assembly and gene annotation should artefactually inflate the rate of gene family evolution, and therefore controlling for it should lead to smaller estimates of lambda. In our case, the difference is really small suggesting the good quality of the genomes and annotation (since we filtered short-read based assemblies). 

- IV - Fourth, we tested different lambda models (e.g. we compared a two lambda model, having separate rates for the outgroup species and the ants) using an error model. Note that the models that included more than three lambdas never converged (e.g. one separate rate per each ant subfamily). In this step, we created a tree with the number of lambdas specified within the newick file. The best fitting model was assessed by conducting likelihood ratio tests between the Maximum likelihood scores of the models, and the difference in the parameters is the different number of lambda (function lr.test() in the R library library("extRemes")).

- V - Finally, we used the lambda estimated under the best fitting model, using the error model, to analyze all gene families, including the high variance ones that were filtered to estimate the lambda in the previous steps. 


### Visualizing CAFE results


We plotted the number of expansions and contractions in the phylogeny using the output generated by CAFE with the tool [CafePlotter](https://github.com/moshi4/CafePlotter) following the instructions.



