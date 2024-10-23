# GAGA Selection pipeline

We used the following described pipeline to infer dn/ds, and positive selection using the coding alignment of each orthologous group and their respective gene trees. 


# Pipeline for selection analyses


1. Retrieve orthogroups from species of interest and categorize them

python3 /home/projects/ku_00039/people/joeviz/orthofinder/orthogroups_v5.py N0.tsv species.txt 0.8(Percentage of species included in the orthogroup as single copy)

# Check the orthogroups and select those interesting for analyses: 1to1, gene families, specific in group of species...


2. Retrieve sequences from orthogroups

perl get_orthogroup_sequences.pl ortologtable(Format 1 column OG, rest are gene names; Ex: tables from previous script) outputdir

    2a. Retrieve already aligned sequences from orthogroups. Avoid repeating the codon/protein alignment
perl get_orthogroup_sequences_fromaln.pl ortologtable(Format 1 column OG, rest are genes) dir_with_aln extension outputdir
#Ej: perl get_orthogroup_sequences_fromaln.pl /home/projects/ku_00039/people/joeviz/orthofinder/run_allGAGA_final_annotations/script_orthogroups_subset/species_subset_80percsp_orthogroups_multicopy_genes.tsv all_orthogroups/orthogroups_seqs/ pep.fasta test_subset_fromaln/


3. Align sequences (Maybe change the script to run independent jobs as the job array give problems in computerome2)

# Codon: edit directory with sequences inside the script
bash /home/projects/ku_00039/people/joeviz/orthology_alignments/run_all_alignments_prank.sh 

# Protein 
bash /home/projects/ku_00039/people/joeviz/orthology_alignments/run_all_alignments_prot.sh 


    3b. Check all alignments finished, and re-run if some failed

bash /home/projects/ku_00039/people/joeviz/orthology_alignments/run_check_alignments.sh
bash /home/projects/ku_00039/people/joeviz/orthology_alignments/run_check_alignments_scriptrerun_prank.sh



4. Run zorro in codon alignments

bash /home/projects/ku_00039/people/joeviz/orthology_alignments/all_orthogroups/codon_alignments/run_all_zorro.sh


    4b. Check if all run

bash /home/projects/ku_00039/people/joeviz/orthology_alignments/all_orthogroups/codon_alignments/run_check_zorro_scriptrerun_no2seqs.sh


5. Remove stop codons, and get zorro statistics and (NOT mask alignment)

#module load ngs tools
#module load anaconda3/4.4.0
perl /home/projects/ku_00039/people/joeviz/orthology_alignments/all_orthogroups/codon_alignments/mask_stop_zorro_codon_alignments.pl


    5b. Check the table, and filter all those orthogroups with bad alignments

perl /home/projects/ku_00039/people/joeviz/orthology_alignments/all_orthogroups/codon_alignments/mask_stop_zorro_codon_alignments.pl


6. GARD here, analyze and create partitions

#run GARD
#bash /home/projects/ku_00039/people/joeviz/Suz/ortholog_alignments/run_all_hyphy_gard_nozorromasked.sh

bash /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/run_all_hyphy_gard.sh

#analyze GARD, get summary table and create partition files
/home/projects/ku_00039/people/joeviz/Suz/ortholog_alignments/analyze_gard.pl


7. Run gene trees. USE DNA model with the codon alignments. Decide about using codon or protein trees. (probably protein, are faster, check topologies and support in both to decide). 

bash /home/projects/ku_00039/people/joeviz/orthology_trees/run_all_genetree_codon_dna.sh

    7b. Check if all gene trees finished ok

bash /home/projects/ku_00039/people/joeviz/orthology_trees/run_check_trees.sh
bash /home/projects/ku_00039/people/joeviz/orthology_trees/run_check_trees_scriptrerun.sh



# MOVED TO A PREVIOUS STEP #####7. Run GARD to check for recombination
#####7b. Check results, and re-run gene trees for orthogroups with recombintation using GARD output


8. Run Hyphy absrel, fitmg for dn/ds values, fubar... 

/home/projects/ku_00039/people/joeviz/Suz/ortholog_alignments/run_all_hyphy_absrel.sh


8b. Get results table. Check Ignasi scripts (/home/projects/ku_00039/people/joeviz/../igngod/scripts/aBSREL_2_table.py)


pyhton3 /home/projects/ku_00039/people/joeviz/../igngod/scripts/aBSREL_2_table.py


9. Hyphy relax. Rename labels with the groups and run relax

Ignasi script to rename the tree, and use the script in Suz folder to submit relax

10. Get GOs and do enrichments. 

see scripts in ~/orthology_trees/all_orthologs/









### Trait evolution pipeline

To understand trait evolution and retrieve gene candidates under positive selection, as well as relaxed or intensified selection associated with our analyzed traits, we followed the next steps:




