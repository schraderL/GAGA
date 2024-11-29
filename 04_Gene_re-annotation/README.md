# GAGA gene re-annotation
After finishing the general gene annotation described in [03_Gene_annotation](../03_Gene_annotation), we thoroughly evaluated the annotation of specific genes and gene families of interest in ant biology. The full list of assessed genes can be found in Table S1E or in [here](https://docs.google.com/spreadsheets/d/1EI8pShvL_YlbxYEyrlkrod7a-a7W58lR3rWOqlQffgA/edit?usp=sharing).

We set a pipeline to evaluate our initial annotations, annotate previously missed genes, and label putative erroneous annotations to be manually reviewed and fixed if necessary using an Apollo web browser set for all GAGA genomes in the [GAGA webpage](https://db.cngb.org/antbase/project) in "Visualization". 

The inputs required for our [pipeline](https://github.com/IAndreu/Gene-annotation-pipeline) are:
- Genome files
   - Genome file in fasta format.
   - General feature format file (GFF) with the initial general gene annotations.
   - Fasta file containing the initial annotated protein sequences.
   - TSV file with protein domain annotation from InterProScan (see [run_interpro.sh](run_interpro.sh) for the command to retrieve the tsv file).

- Gene/Gene family information to conduct the re-annotation
   - We require a table containing the gene and gene families to re-annotate, in addition to further information and parameters used in the pipeline, such as the InterPro or Pfam domains. The table we used in GAGA can be found [here](https://docs.google.com/spreadsheets/d/1EI8pShvL_YlbxYEyrlkrod7a-a7W58lR3rWOqlQffgA/edit?usp=sharing) (Note: Make sure not to edit it in Excel, or format the cells as text as range numbers are often changed to dates and the pipeline will report an error).
   - If the 3rd column "Blast" is set to Yes, a fasta file containing protein sequences (ideally from closely related organisms) from that gene is expected to be included in the input with format "GENENAME_db.fasta". The protein sequences used in GAGA are available in the file [Gene_family_database](Gene_family_database.zip).

In addition, the pipeline requires the following prerequisites to be installed, requiring for some of them to specify the paths in the main script. See a full detailed manual in the github repository of our pipeline: https://github.com/IAndreu/Gene-annotation-pipeline 
- Perl: Perl is installed by default in most operating systems. See https://learn.perl.org/installing/ for installation instructions.
- Python: Download the newest version available from https://www.python.org/downloads/
- Biopython: Download from https://biopython.org/wiki/Download
- BLAST: Download blast executables from: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
- HMMER: "conda install -c bioconda hmmer" 
- GeMoMa: GeMoMa is implemented in Java using Jstacs and can be downloaded from: http://www.jstacs.de/index.php/GeMoMa. The version of GeMoMa used in our pipeline is v1.7.1.
- MAFFT: Download the newest version available from https://mafft.cbrc.jp/alignment/software/
- InterProScan: Download it from https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html#obtaining-a-copy-of-interproscan
- BITACORA: Download from https://github.com/molevol-ub/bitacora


To run the pipeline, first clone the [GitHub repository](https://github.com/IAndreu/Gene-annotation-pipeline) containing the scripts, and then run the [gene family pipeline script](submit_run_genefamilypipeline.sh). For a more detailed parameters, see the following usage:
```
usage: run_analysis.py [-h]
                       pipeline_dir gene_families_info gene_families_db
                       proteome interpro gff genome bitacora name out_dir
                       threads

Pipeline for gene family re-annotation across hundreds of genomes.

positional arguments:
  pipeline_dir        Main directory of the pipeline that contains the
                      README.md.
  gene_families_info  Excel file containing all the gene families information.
  gene_families_db    Directory containing the query protein databases
                      (GENEFAMILY-NAME_db.fasta), where the “GENEFAMILY-NAME”
                      label is a gene family name from the Excel file. The
                      addition of ”_db” to the database name with its proper
                      extension is mandatory. ej: PATH/TO/DB
  proteome            File with predicted proteins in FASTA format.
  interpro            File with predicted domains from InterPro in TSV format.
  gff                 File with structural annotations in GFF3 format.
  genome              File with genomic sequences in FASTA format.
  bitacora            Path to script runBITACORA_command_line.sh. (edit the script and add the location of the BITACORA scripts folder and GeMoMa path)
  name                Name/ID of the genome i.e.: GAGA-0001.
  out_dir             Output directory. ej: PATH/TO/OUT_DIR
  threads             Number of threads to be used.

optional arguments:
  -h, --help          show this help message and exit
```

Once the pipeline finishes, it will generate a folder containing the curated gene models for each gene or gene family, and a summary table named "table_results.txt" that contains the number of annotated genes through each step, as well as further information such as if the gene is complete, partial, chimeric, or if some gene has been newly annotated (i.e: previously missed in the initial general gene annotation). Importantly, the 13rd column contains the folder containing the best gene models, or alternatively, it will say REVIEW for cases that need manual checking. These cases need to be manually checked, and the curated or checked models added to a folder (in the same directory that contains Step1 and Step2 folders) named "Manual" with the GFF named as "GENENAME_manualmodels.gff3"; and finally edit column 13rd as "Manual". These cases would be, in general terms:
- Missing genes that are expected to be present in the species genome. The action required would be to confirm that the gene is indeed absent, or otherwise the gene model could have been missed by the annotation pipeline, as some factor could be preventing its annotation, such as errors in the start/stop codons, or in the splicing sites. In such cases, the gene can be annotated in Apollo, and the generated GFF file included in the output folder to be integrated into the annotation (Creating a folder named "Manual", and changing column 13 in the "table_results.txt" with "Manual").
- More genes annotated than expected. The action required would be to check that indeed all retrieved annotations are duplications from the same gene, or otherwise the program might be annotating other different genes that, for instance, only share a protein domain but do not present similarity in the other gene regions. Incorrect genes can be deleted from the GFF and added in a newly created "Manual" folder with the file named as "GENENAME_manualmodels.gff3", and change column 13 in the "table_results.txt" with Manual. 
- Putative chimeric genes. Genes that are longer than expected, and have similarity with distinct genes, evidencing that the annotated gene might be an artefactual chimera merging two separate gene models. These cases can be validated and curated using Apollo, and the curated GFF should be added in the Manual folder as explained above. 

Finally, when all cases are reviewed (Note again that it requires to create a folder named "Manual" for each reviewed gene, and the file "table_results.txt" edited accordingly in column 13), the script gen_final_output.py will create an output folder with the final gene models of each gene family. See the [script to generate the final output](submit_run_genefamilypipepine_generatefinaloutput.sh) to run this final step.


## Chemosensory gene families
We used a different pipeline for the chemosensory receptors described in [Chemosensory_gene_families folder](Chemosensory_gene_families), including the olfactory (OR) and gustatory (GR) receptors, given the repetitive nature of these gene families (i.e.: many copies arranged in tandem) and fast evolution in the ants. 
We also use homology-based searches to validate and rename the annotations, and we also classify the chemosensory genes family members as described in the scripts included in this folder. 



## Merging all re-annotations in a single annotation file
Once all gene re-annotations are curated and completed, we combine the GFF files generated with the initial general gene annotation, replacing curated gene models and adding newly annotated genes. Therefore, a final GFF file is generated including all isoforms, and also the GFF containing a representative isoform per gene. We used the following pipeline to retrieve the final GFF files:

It first compares the overlapping positions between curated gene annotations and original annotations. If a new curated gene model overlappes with a model in original annotations (overlapping length >100 bp), the latter is removed. For the representative isoform, curated gene models are prioritized as representantatives.

```
Usage : sh pipeline_representative.sh <GAGA_ID> <prefix>
```

