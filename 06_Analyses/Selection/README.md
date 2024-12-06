# Selective constraint analyses


Scripts used in the selection analyses


1. The protein coding sequence (cds) for all genes in each orthogroup were aligned using PRANK.

Script to create an array of jobs to submit the PRANK alignments (edit the inputs inside the script):
```bash
bash 01_run_all_alignments_prank.sh
```

Then, check if all alignments were generated successfully, or rerun them if they failed, allowing for more memory and time in the cluster. 
```bash
bash 01b_run_check_alignments_scriptrerun_prank.sh
```


2. The quality of the multiple sequence alignments (MSAs) were evaluated using Zorro.

Script to create an array of jobs to submit Zorro (edit the inputs inside the script):
```bash
bash 02_run_all_zorro.sh
```

Then, the orthogroups with average quality below 4 in more than 75% of the unaligned sequence lengths were filtered for the selection analyses using the following scripts (edit the inputs inside the script).
```bash
perl 02b_mask_stop_zorro_codon_alignments.pl
perl 02c_mask_stop_zorro_codon_alignments_part2_filterbadaln.pl
```


3. The Genetic Algorithm for Recombination Detection (GARD) was used in HyPhy to screen the alignments for recombination breakpoints and to split them into separate partitions.

Script to create an array of jobs to submit HyPhy GARD (edit the input and command paths inside the script):
```bash
bash 03_run_all_hyphy_gard.sh
```

Then, the alignments with recombination breakpoint are split into separate partitions, using the following script:
```bash
perl 03b_analyze_gard_v2.pl
```

The following script can be used for a specific orthogroup (HOG), instead of iterating through a whole folder of alignments [03c_analyze_gard_v2_single_HOG.pl](03c_analyze_gard_v2_single_HOG.pl).


4. A maximum likelihood tree was reconstructed for each MSA partition using iq-tree.

Script to create the commands and run iqtree for each MSA (edit the input inside the script):
```bash
bash 04_run_all_genetree_codon_dna_gard.sh
```


5. The alignments for each partition (or entire orthogroup sequences when no recombination was identified) were further evaluated using HmmCleaner, a segment-filtering software that detects putative errors in MSAs. The identified putative misaligned segments were masked as gaps, and blocks with gappy regions in more than 50% of the species were removed, similar to sequences with fewer than 15 unmasked codons.

Script to run HmmCleaner:
```bash
bash 05_run_all_hmmcleaner_cleanaln.sh
```


6. The adaptive branch-site random effects likelihood (aBSREL) model implemented in the HyPhy package was used to detect hallmarks of positive selection (PS) and infer dN/dS across all branches in each partitioned-orthogroup gene tree, using the high-confidence codon alignments.

Script to create the jobs to submit HyPhy aBSREL for all orthogroups (edit the input and command paths inside the script):
```bash
bash 06_run_all_hyphy_absrel.sh
```
In addition, the script [06b_run_all_hyphy_absrel_onlysinglecopy.pl](06b_run_all_hyphy_absrel_onlysinglecopy.pl) allows to run aBSREL for a subset of orthogroups provided as a list in the input file specified in line 21. 

Then, the output for each orthogroup partition is parsed to create summary tables including the positive selection and dN/dS for each branch:
```bash
perl 06c_analyze_absrel.pl
```


7. Inferring general positive selection patterns in the ant phylogeny using the output tables from HyPhy aBSREL.

We used the following script to map the nodes in the gene tree to the species tree conducting a gene-tree versus species-tree reconciliation. For this analysis, we used the 8384 single-copy orthogroups with high-quality codon alignments in ≥ 80% of the sequenced ants. Each node in a gene tree was assigned to the corresponding node in the species tree that included all species present in a gene-tree clade. To avoid incorrectly assigning internal nodes in the gene tree when their evolutionary history deviated from the species tree (e.g. because of incomplete lineage sorting), reconciliations between gene tree nodes and species tree nodes were retained only when at least 60% of the species descending from the species tree node were represented in the gene tree node. The resulting reconciliation allowed us to retrieve the number of genes under positive selection across all nodes of the species tree, and the proportion of positive selection events among them.

```bash
perl 07_analyze_absrel_omegasummary_pernode_relax2withexceptions.pl species_all_80percsp_orthogroups_singlecopy_counts.tsv Hyphy_absrel_withgard_hmmclean_50_rooted_omega_pernode_80pc1to1_fdr001_strict60_relax2_withexcept.txt 0.6 None 0.02 0.001
```
The inputs and further parameters are detailed inside the script (lines 11-50).



8. Identifying signatures of relaxed or intensified selection associated with phenotypic traits.

To assess shifts in the strength of selection on genes associated with phenotypic traits, we used the RELAX model106 implemented in HyPhy. The following script labels tips belonging to species expressing a “test” or “reference” trait in the orthogroup gene trees (using separate partitions when recombination was identified as described above), and create the jobs to submit HyPhy RELAX for all orthogroups (edit the inputs, commands and paths inside the script):
```bash
bash 08_run_all_hyphy_relax.sh ../Species_selection_Queen_worker_Dimorphism_HighVsLow_clades_all.txt

```
The input provided in the command line (e.g. [Species_selection_Queen_worker_Dimorphism_HighVsLow_clades_all.txt](Species_selection_Queen_worker_Dimorphism_HighVsLow_clades_all.txt)) contains the list of species with the "test" or "reference" trait to use for the analysis. 


Then, we use the following script to get gene candidates with intensified or relaxed selection associated with the analyzed trait from the RELAX analyses:
```bash
perl 08b_get_candidates_relax.pl runRELAX/relax/ Relax_results
```
Note that "runRELAX/relax/" is the output folder from hyphy RELAX containing the JSON files, and "Relax_results" is the folder to store the results from this script. 



