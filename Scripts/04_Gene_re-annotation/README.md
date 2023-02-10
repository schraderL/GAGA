# GAGA gene re-annotation
After finishing the general gene annotation described in [03_Gene_annotation](../03_Gene_annotation), we thoroughly evaluated the annotation of specific genes and gene families of interest in ant biology. The full list of assessed genes can be found [here](https://docs.google.com/spreadsheets/d/1EI8pShvL_YlbxYEyrlkrod7a-a7W58lR3rWOqlQffgA/edit?usp=sharing).

We set a pipeline to evaluate our initial annotations, annotate previously missed genes, and label putative erroneus annotations to be manually reviewed and fixed if necessary using an Apollo web browser set for all GAGA genomes in the [GAGA webpage](https://db.cngb.org/antbase/project) in Visualization. 

The inputs required for our pipeline are:
- Genome files
   - Genome file in fasta format.
   - General feature format file (GFF) with the intital general gene annotations.
   - Fasta file containing the initial annotated protein sequences.
   - TSV file with protein domain annotation from InterProScan (see [run_interpro.sh](run_interpro.sh) for the command to retrieve the tsv file).

- Gene/Gene family information to conduct the re-annotation
   - We require a table containing the gene and gene families to re-annotate, in addition to further information and parameters used in the pipeline, such as the InterPro or Pfam domains. The table we used in GAGA can be found [here](https://docs.google.com/spreadsheets/d/1EI8pShvL_YlbxYEyrlkrod7a-a7W58lR3rWOqlQffgA/edit?usp=sharing).
   - If the 3rd column "Blast" is set to Yes, a fasta file containing protein sequences (ideally from close organisms) from that gene is expected to be included in the input with format "GENENAME_db.fasta". The protein sequences used in GAGA are available in the file [Gene_family_database](Gene_family_database.zip).

In addition, the pipeline requires the following pre-requisites to be installed, requiring for some of them to specify the paths in the main script. See a full detailed manual in the github repository of our pipeline: https://github.com/IAndreu/Gene-annotation-pipeline 
- Perl: Perl is installed by default in most operating systems. See https://learn.perl.org/installing/ for installation instructions.
- Python: Download the newest version available from https://www.python.org/downloads/
- Biopython: Download from https://biopython.org/wiki/Download
- BLAST: Download blast executables from: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
- HMMER: "conda install -c bioconda hmmer" 
- GeMoMa: GeMoMa is implemented in Java using Jstacs and can be downloaded from: http://www.jstacs.de/index.php/GeMoMa. The version of GeMoMa used in our pipeline is v1.7.1.
- MAFFT: Download the newest version available from https://mafft.cbrc.jp/alignment/software/
- InterProScan: Download it from https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html#obtaining-a-copy-of-interproscan
- BITACORA: Download from https://github.com/molevol-ub/bitacora


To run the pipeline, see the following usage:
```
usage: run_analysis.py [-h]
                       pipeline_dir gene_families_info gene_families_db
                       proteome interpro gff genome bitacora name out_dir
                       threads

Pipeline for gene family re-annotation accross hundreds of genomes.

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
  bitacora            Path to script runBITACORA_command_line.sh.
  name                Name/ID of the genome ej. GAGA-0001.
  out_dir             Output directory. ej: PATH/TO/OUT_DIR
  threads             Number of threads to be used.

optional arguments:
  -h, --help          show this help message and exit
```


## Chemosensory gene families
We used a different pipeline for the chemosensory receptors described in [Chemosensory_gene_families folder](Chemosensory_gene_families), including the olfactory (OR) and gustatory (GR) receptors, given the repetitive nature of these gene families (i.e.: many copies arranged in tandem) and fast evolution in the ants. 



## Merging all re-annotations in a single annotation file
Once all gene re-annotations are curated and completed, we combine the GFF files generated with the initial general gene annotation, replacing curated gene models and adding newly annotated genes. Therefore, a final GFF file is generated including all isoforms, or a representative isoform per gene. 

ZIJUN SCRIPT

