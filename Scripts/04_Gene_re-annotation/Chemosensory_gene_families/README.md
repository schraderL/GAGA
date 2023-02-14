# GAGA gene re-annotation of chemosensory-related gene families

## Odorant receptors (ORs)

The annotation of the odorant receptors is conducted using [HAPpy-ABCENTH pipeline](https://github.com/biorover/HAPpy-ABCENTH). Specifically, we use ABCENTH (Annotation Based on Conserved Exons Noticed Through Homology), a gene finder designed for multigene families with extremely high sequence divergence but highly conserved exon structure, such as ant ORs. It is also designed to avoid gene fusion in tandem arrays, see the [github repository](https://github.com/biorover/HAPpy-ABCENTH) for further details.

We use the [run_OR_annotation_happy_abcenth_finalpipeline.sh script](run_OR_annotation_happy_abcenth_finalpipeline.sh) to first run ABCENTH to annotate the odorant receptor genes using the HMM profiles built from manually curated OR models in five ant genomes (Harpegnathos saltator , Ooceraea biroi, Camponotus floridanus, Solenopsis invicta and Stigmatomma sp. (GAGA-0391)). The HMM files can be found in [RefHMMs_HsSrObCfSi_0.45dist.tar.g](RefHMMs_HsSrObCfSi_0.45dist.tar.gz).
Next, the OR annotations are converted to GFF3 format, and are evaluated using homology-based searches. We use insect OR sequences to identify the OR co-receptor (ORCO), rename the OR genes, and finally we classify them as complete sequences, pseudogenes, partial or fragment annotations. The sequence datasets can be found in [Chemo_db.zip file](Chemo_db.zip).

The pipeline requires the following programs:
- HAPpy-ABCENTH
- [Genome tools](http://genometools.org/)
- InterProScan
- Blast


## Gustatory receptors (GRs)

The annotation of the gustatory receptors is also conducted using [HAPpy-ABCENTH pipeline](https://github.com/biorover/HAPpy-ABCENTH). However, here we conducted both HAPpy and ABCENTH methods in HAPpy-ABCENTH pipeline to annotate the GRs, and combined the annotations from both methods to generate the final gene models. Finally, the annotations are also evaluated, identifying conserved GR members and classifying genes as complete, pseudogenes, partial or fragment models. 

First, we used high-quality GR gene models to create a dataset that was used further to conduct the GR annotations. This dataset consists in GRs from Drosophila melanogaster, Daphnia pulex, Strigamia maritima, Ixodes scapularis, and the following ant genomes: 
GAGA-0063 GAGA-0198 GAGA-0199 GAGA-0245 GAGA-0306 GAGA-0340 GAGA-0365 GAGA-0378 GAGA-0391 GAGA-0392 GAGA-0521 GAGA-0534 GAGA-0552 GAGA-0580 NCBI-0001 NCBI-0002. The protein sequences are in [GR_ant_genewise_db.fasta file](GR_ant_genewise_db.fasta). The species names can be found in the [GAGA species list file](GAGA_species_list.txt). 
In addition, we used the ant annotations to create the exon HMM profiles using HAPpy-ABCENTH, that are in the [RefHMMs_ant_GRs_0.40dist.zip file](RefHMMs_ant_GRs_0.40dist.zip). 


***01_run_GR_annotation_abcenth.sh*** 
This first pipeline will conduct the GR annotation using ABCENTH. It requires the genome assembly and the folder containing the above explained HMM profiles. It also requires to specify the path to the [chemo_db folder](Chemo_db.zip) (lines 76-83) and the OR annotations (these are used to confirm that no OR has been annotated as a GR, given the similarity of these two gene families). All necessary scripts to run the pipeline can be found in the [scripts folder](scripts).

***02_run_GR_annotation_happy.sh***
Here, GRs will be annotated using HAPpy homology-based prediction, which requires GeneWise. The input is the genome assembly, and the (above mentioned file)[GR_ant_genewise_db.fasta file] containing GR proteins. It is also necessary to specify the path to the [chemo_db folder](Chemo_db.zip) (lines 76-83) and the OR annotations. All necessary scripts to run the pipeline can be found in the [scripts folder](scripts).

***03_run_GR_annotation_happy_abcenth_combined.sh*** 
The GR models predicted in the previous two steps are combined generating the final gene models. 


The pipeline requires the following programs (edit the path in the scripts, as well as all input files):
- HAPpy-ABCENTH
- [Genome tools](http://genometools.org/)
- InterProScan
- Blast
- [GeMoMa v1.7.1](http://www.jstacs.de/index.php/GeMoMa) (Note that newer versions will not work in this pipeline, as the module CompareTranscripts has been modified)


# Other chemosensory gene families

For the remaining chemosensory-related gene families, we used mostly the same pipeline, which changes the protein dataset for each gene family. We conducted the annotation of these gene families using the general re-annotation pipeline (they are included in the [excel table](https://docs.google.com/spreadsheets/d/1EI8pShvL_YlbxYEyrlkrod7a-a7W58lR3rWOqlQffgA/edit?usp=sharing)). Then, we curate putative chimeric models, and identify conserved members (and subfamilies such as in the case of IR and iGluRs) and classify them into complete, pseudogene or partial copies. 

- Ionotropic receptors (IRs): [run_IR_annotation.sh](run_IR_annotation.sh)
- CD36-SNMPs: [run_CD36SNMP_annotation.sh](run_CD36SNMP_annotation.sh)
- Pickpocket proteins (PPKs, also known as DEG-ENAG channels): [run_PPK_annotation.sh](run_PPK_annotation.sh)
- Odorant binding proteins (OBPs): [run_OBP_annotation.sh](run_OBP_annotation.sh)
- Chemosensory proteins (CSPs): [run_CSP_annotation.sh](run_CSP_annotation.sh)
- NPC2: [run_NPC2_annotation.sh](run_NPC2_annotation.sh)


We provide a script to create a job array that runs the pipeline in all genomes in [scripts/create_submit_script_allgenomes.sh](scripts/create_submit_script_allgenomes.sh). The code inside can be changed to be adjusted for each gene family, and the necessary paths also have to be modified accordingly. 
