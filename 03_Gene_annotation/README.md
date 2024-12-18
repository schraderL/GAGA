### GAGA General gene annotation

General genome annotation was conducted by combining gene annotation from several sources using a pipeline optimized for the ant genomes generated by the GAGA project. Depending on the availability of RNA-seq for each species, we followed two pipelines:

# RNA-seq available

1.  First, short-read RNA-seq are aligned to the reference repeat soft-masked genome assembly using STAR v2.7.2b default options. In the case of hifi long-reads RNA-seq, we used bbmap 38.90 to align the reads to the genome. Then, we retrieved the publicly available gene annotations from the fruit fly Drosophila melanogaster, the red flour beetle Tribolium castaneum, the parasitoid wasp Nasonia vitripennis, the honeybee Apis mellifera, the clonal raider ant Ooceraea biroi and the Florida carpenter ant Camponotus floridanus, that can be found in the file "GeMoMa_input/Whole_genomes.zip" in [here](https://sid.erda.dk/sharelink/EJZrYWKPrj). The annotations from these insect species are used to conduct homology-based gene predictions using GeMoMa v1.7.1, which also incorporates the RNA-seq evidences for splice site prediction.
	
	**01_run_gemoma_annotation_andRNAmapping.pl** Script generates the submission commands for each species (using our GAGA-IDs [file](GAGA-ID_all.txt)).


2. The independent RNA-seq alignments are merged creating a consensus GTF (Gene transfer format) using Stringtie v2.1.5.
	
	**02_run_stringtie_bam_to_gtf.pl**


3. BestORF (Molquest package, Softberry) is used to identify open reading frame (ORF) in the transcript sequences. Then, the identified ORF, genome assembly and the generated gtf in the Step 2 are input in the following pipeline to retrieve gene models based on RNA-seq evidences.

	**03_runRNAseq_prediction.sh**
  
4. We randomly select \~1,000 high–quality genes from GeMoMa prediction to train Augustus v3.2.2. The de novo gene prediction is then performed using Augustus with the repeat-masked genome, filtering out genes with lower length than 150bp or incomplete ORF.  Transposon-related proteins were identified and filtered using a BLASTP search against Swissprot database with matches proteins containing any one of the following keywords were filtered: transpose, transposon, retro-transposon, retrovirus, retrotransposon, reverse transcriptase, transposase, and retroviral.
	**04_Get_trainingDataSet.sh**
	**04_Run_augustus.sh**

5. Gene annotations from the three evidences are combined generating the final gene annotation for each genome. from which we also generate an annotation with a single representative isoform per gene (i.e.: longest isoform is kept as the representative). 

	**05_combine_gene_predictions.sh**

6. We generate an annotation with a single representative isoform per gene (i.e.: longest isoform is kept as the representative).

	**06_Merge_GeneIsoforms.sh**

7. Transposon-related proteins were identified and filtered using a BLASTP search against Swissprot database and the transposable element protein database from RepeatMasker. Script to filter all transposable and repeat proteins from the gene annotation:
	
	**07_run_gene_annotation_filtering_all.pl** Script will search for all genome files and run the main script. See the script to individually run this step in one genome: "07a_run_gene_annotation_filtering.pl" 


# No RNA-seq available

1. We run GeMoMa including the gene annotations from the closest species with RNA-seq available. 
	
	**01_run_gemoma_annotation_noRNAspecies_with_close_relative.pl**

Step 2 and 3 are skipped because of the absence of RNA-seq; Steps 4 and 7 are the same as above, and Steps 5 and 6 are similar but excluding the lines involving RNA-seq based annotations.

	**05_combine_gene_predictions_noRNAseq.sh**
	**06_Merge_GeneIsoforms_noRNAseq.sh**


