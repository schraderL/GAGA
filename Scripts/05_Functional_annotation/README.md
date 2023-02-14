# GAGA gene functional annotation

We combined several similarity-based searches to conduct the functional annotation of the protein-coding genes annotated in all ant GAGA genomes. We first used BLASTP to retrieve the 5? best hits against the following databases: ***Zijun check this and add the command/scripts and link to the databases you used, also e-value to filter hits***
- Swissprot
- Trembl
- KEGG
- KOG

In addition, we searched for the specific protein-domain signatures in protein-coding gene sequences using [InterProScan](https://www.ebi.ac.uk/interpro/download/InterProScan/). The [script run_all_interpro.sh](run_all_interpro.sh) creates a job array to run InterPro in the protein sequences for each GAGA genome.

[EggNOG-mapper v2.1.6](https://github.com/eggnogdb/eggnog-mapper) is also used to retrieve the functional annotation using precomputed orthologous groups and phylogenies from the eggNOG database (http://eggnog5.embl.de). The [script run_all_eggnog.sh](run_all_eggnog.sh) creates a job array to run EggNOG-mapper in the protein sequences for each GAGA genome.


Finally, all functional annotations are combined generating a general functional annotation for each ant genome, that includes the best hit with each database, as well as the Gene ontology (GO) and KEGG terms transferred from significant hits. First, we combine homology-based and protein-domain searches, retrieving a summary table with the functional annotations and GO associated terms (Zijun check this line). 

***ZIJUN SCRIPT*** Zijun describe and add here the script you use to generate the annotation.xls and the wego files (i.e.: OUT-0001_gene_merged.wego), as I used them as input in the next script. 

These annotations are then combined with EggNOG to generate a summary table containing the protein identifier, best hit from SwissProt, EggNOG, KEGG, COG, TrEMBL and InterPro, and the number of GOs. In addition, we create a table and a ".annot" file containing the GO and KEGG terms. These steps are run in the [get_all_functional_annotation.pl script](get_all_functional_annotation.pl). All input paths and files need to be specified accordingly in the script (lines 13-23).
```
perl get_all_functional_annotation.pl
```


### Functional annotation of ortholog groups across the 163 ant GAGA genomes

After conducting the [orthology assessment](../06_Analyses/Orthology) across the protein-coding genes of our 163 ant genomes, we retrieved the consensus functional annotation for each orthogroup. Briefly, we used the previously generated functional annotation for each genome, and counted the number of genes with the same annotations in each orthogroup, retrieving the consensus functional annotation. In addition, we transferred for each orthogroup the GO and KEGG terms that were annotated at least in 1/3 of the genes in the orthogroup. 
```
perl get_functional_annotations_per_orthogroup.pl
```  



