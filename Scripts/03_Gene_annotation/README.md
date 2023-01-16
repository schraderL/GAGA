### General gene annotation used in GAGA

The folder Gemoma_input_genome_data contains publicly available genomes that we used in GeMoMa for the homology-based annotation. 
Depending on the availability of RNA-seq for each species, we followed two pipelines:

# RNA-seq available

1. We map the RNA-seq and run GeMoMa.
	01_run_gemoma_annotation_andRNAmapping.pl Script generated the submission commands for each species (using our GAGA-IDs).

2. We retrieve a gtf using Stringtie and the aligned RNA-seq reads.
	02_run_stringtie_bam_to_gtf.pl

Pending from here: Zijun

3. Stringtie to generate RNA-seq based annotations
4. Augustus to get de-novo annotations
5. Script to combine all annotations
6. Script to generate a gff3 with one representative isoform (longer) per gene


7. Script to filter all transposable and repeat proteins from the gene annotation
	07_run_gene_annotation_filtering_all.pl Script will search for all genome files and run the main script. See the script to individually run this step in one genome: "07a_run_gene_annotation_filtering.pl" 


# No RNA-seq available

1. We run GeMoMa including the gene annotations from the closest species with RNA-seq available. 
	01_run_gemoma_annotation_noRNAspecies_with_close_relative.pl


Step 2 and 3 are skipped because the absence of RNA-seq; and Step 4, 5, 6 and 7 are the same as above.



