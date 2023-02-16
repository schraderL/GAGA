# GAGA orthology assessment

We used [OrthoFinder v2.5.4](https://github.com/davidemms/OrthoFinder) to conduct the orthology assessment across the annotated proteins in our 163 genomes, using one representative isoform per gene. Our comprehensive phylogeny, inferred using non-coding regions retrieved from the whole-genome alignments of the 163 ant genomes using Cactus, was used in OrthoFinder to increase accuracy by retrieving orthogroups at each hierarchical level. The commands to run OrthoFinder can be found in:

***[01_submit_orthofinder.sh](01_submit_orthofinder.sh)***


The resulting file ***N0.tsv***, in which each row contains the genes belonging to a single orthogroup, was used to retrieve the different types or orthologous groups: single copy, multi copies, clade/species specific or others; using the script [orthogroups_v5.py](orthogroups_v5.py).
The script requires as input the N0.tsv output from OrthoFinder, a text file containing the list of species to retrieve the orthogroups (i.e.: it can be all species input in orthofinder, or a subset of those), and the percentage of species needed to classify the orthogroups (i.e.: 90% means that orthogroups with one copy in at least 90% of the genomes will be classified as single copy).
```
pyhton3 orthogroups_v5.py N0.tsv species_list_all.txt 0.75
```
See the following script to submit this step:

***[02_get_orthogroup_classification.sh](02_get_orthogroup_classification.sh)***


Finally, we retrieve all sequences, both protein and cds, for each orthogroup that will be used in the comparative genomics analyses. Here, we use the script [get_orthogroup_sequences](). To retrieve the sequences for the specific list of orthogroups, that it can just be the single copy or all orthogroups. 

***[03_retrieve_orthogroup_sequences.sh](03_retrieve_orthogroup_sequences.sh)***




## Orthology assessment using outgroup non-ant species

We also conducted the orthology assessment using eight outgroup non-ant high-quality genomes from the Apoidea superfamily, that we used in the comparative genomics analyses. 

# Anthophila (bees)
1 - Apis mellifera: https://www.ncbi.nlm.nih.gov/genome/48
2 - Bombus hortorum: https://portal.darwintreeoflife.org/data/root/details/Bombus%20hortorum
3 - Nomada fabriciana: https://portal.darwintreeoflife.org/data/root/details/Nomada%20fabriciana
4 - Osmia bicornis: https://www.ncbi.nlm.nih.gov/genome/76072
5 - Lasioglossum lativentre: https://portal.darwintreeoflife.org/data/root/details/Lasioglossum%20lativentre 

# Crabronidae (digger wasps)
6 - Cerceris rybyensis: https://portal.darwintreeoflife.org/data/root/details/Cerceris%20rybyensis
7 - Nysson spinosus: https://portal.darwintreeoflife.org/data/root/details/Nysson%20spinosus

# Ampulicidae
8 - Ampulex compressa: https://www.ncbi.nlm.nih.gov/genome/?term=Ampulicidae


We re-annotated these genomes using our pipeline described in [Gene_annotation](../../03_Gene_annotation), confirmed the expected phylogeny relationships between these species and the ants based on [Peters et al. 2017](http://dx.doi.org/10.1016/j.cub.2017.01.027), and used the proteomes and phylogeny in OrthoFinder and following steps as detailed above. 






