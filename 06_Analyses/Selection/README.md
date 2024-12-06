# GAGA Selective constraint analyses


Scripts used in the selection analyses. 


1. The protein coding sequence (cds) for all genes in each orthogroup were aligned using PRANK.

Script to create an array of jobs to submit the PRANK alignments (edit the inputs inside the script):
```bash
    bash 01_run_all_alignments_prank.sh
```

Then, check if all alignments were generated succesfully, or rerun them if they failed, allowing for more memory and time in the cluster. 
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

Script to create an array of jobs to submit HyPhy aBSREL (edit the input and command paths inside the script):
```bash

```

Then, the output is parsed to create summary tables including the positive selection and dN/dS for each branch:


The following script can be used just to retrieve branches under significant positive selection providing the specific p-value or FDR. 



The output was further used to do the Figure 3



7. The using RELAX.

Script to :
```bash

```




----





8. Run Hyphy absrel, fitmg for dn/ds values, fubar... 

/home/projects/ku_00039/people/joeviz/Suz/ortholog_alignments/run_all_hyphy_absrel.sh


8b. Get results table. Check Ignasi scripts (/home/projects/ku_00039/people/joeviz/../igngod/scripts/aBSREL_2_table.py)


pyhton3 /home/projects/ku_00039/people/joeviz/../igngod/scripts/aBSREL_2_table.py


9. Hyphy relax. Rename labels with the groups and run relax

Ignasi script to rename the tree, and use the script in Suz folder to submit relax






### Identifying signatures of selection associated with phenotypic traits

To understand trait evolution and retrieve gene candidates under positive selection, as well as relaxed or intensified selection associated with our analyzed traits, we followed the next steps:




