# GAGA Genome assembly
Folder with the scripts used to conduct the genome assembly and quality evaluation.

## Genome assembly
We aimed to generate both long-read PacBio and single-tube long fragment read (stLFR) sequencing data for all of the ant genomes generated in the project. However, stringent DNA-quality requirements and limited biomass availability were challenging, preventing us to generate PacBio reads in 19 species for which we only have stLFR short-reads (and one with 10X genomics). In addition, some species were sequenced by subprojects under the GAGA project, and WGS Illumina reads were generated instead of stLFR. The sequencing data used for each species can be checked in our [summary table](GAGA_genome_stats.xlsx). Then, we classify the genome assembly procedures into three pipelines according to the sequencing data available:  

- **A: Genomes generated with both PacBio long- and stLFR short-reads**
   - A.1. The PacBio SMRT-Analysis package was used for processing polymerase reads, removing sequencing adapters and filtering reads with low quality and short length (parameters: minSubReadLength:500). Note this step is done in the sequencing company (Novogene). 
   - A.2. The clean reads are then assembled using Wtdbg2 v2.5.
   ```
   perl wtdbg2.pl -t 16 -x sq -g 300m -o GAGA_ant_id CLR_reads.fa.gz
   perl wtdbg2.pl -t 16 -x ccs -g 300m -o GAGA_ant_id Hifi_reads.fa.gz
   ```
   - A.3. The assembled contigs from Wtdbg2 were then scaffolded using SSPACE-LongRead. 
   ```
   perl SSPACE-LongRead.pl -c assembled_genome.fa -p reads.fasta -b output
   ```
   - A.4. We performed an additional round of gap filling to eliminate the gaps within scaffolds using [LR_Gapcloser](https://github.com/CAFS-bioinformatics/LR_Gapcloser) with PacBio subreads.
   ```
   sh LR_Gapcloser.sh -i scaffolds.fasta -l input.fasta
   ```
   - A.5. To further improve the accuracy of the genome assembly, two-step polishing was performed on the initial assembly. In the first round, [Arrow](https://github.com/skoren/ArrowGrid) software was used to map the PacBio sequences to the genome assembly. The high coverage PacBio sequencing data could efficiently correct the small indels and substitutions in the initial assembly, obtaining the consensus sequences.
   ```
   /usr/local/smrtlink/smrtcmds/bin/pbcromwell run pb_resequencing  --entry GAGA-id_subSet.xml --entry GAGA-id_RefSet.xml  --config $PWD/cromwell.conf  --backend Local  --output-dir pb_resequencing --nproc 30  >run.log
   ```
   - A.6. Because of the high error rate of PacBio raw reads, the consensus sequences were also subjected to a further polishing step using the short reads with [NextPolish v1.3.0](https://github.com/Nextomics/NextPolish) . 
   - A.7. Finally, the polished scaffolds were further scaffolded using the barcoding information from stLFR reads with [SLR-superscaffolder pipeline](https://github.com/BGI-Qingdao/SLR-superscaffolder).


- **B: Genomes generated with both PacBio long- and WGS Illumina short-reads**
   - The same steps were used as those mentioned in A, except for A.7. (additional scaffolding using stLFR data).

- **C: Genomes generated with only stLFR short-reads**
   - C.1. Raw stLFR reads were first cleaned from adaptors and PCR duplicates, and the barcode IDs were assigned in the read names using [stlfr2supernova_pipeline](https://github.com/BGI-Qingdao/stlfr2supernova_pipeline). 
   - C.2. We used two different genome assembly pipelines and evaluated the genome assembly quality, given the low contiguity metrics retrieved in these short-read genome assemblies:
   - C.2a. MaSuRCa-SLRsuperscaffolder pipeline
      - C.2a.1 The clean reads were then assembled using MaSuRCA v3.3.0 (“JF_SIZE = 2500000000”). 
      - C.2a.2 The resulting assembly was further scaffolded using the barcoding information from stLFR reads to assemble contigs into scaffolds with [SLR-superscaffolder pipeline](https://github.com/BGI-Qingdao/SLR-superscaffolder) (Same script as in previous step A.7). 
   - C.2b. Supernova pipeline
      - C.2b.1 The clean reads were assembled using Supernova pipeline (v2.1.1). See the [github repository](https://github.com/BGI-Qingdao/stlfr2supernova_pipeline) for specific details about the pipeline. 
      - C.2b.2 Then, a scaffolding step was conducted using [Scaf10X](https://github.com/wtsi-hpag/Scaff10X).
   - C.3. Then, both assemblies were compared to retrieve the most contiguous and complete genome assembly (see quality evaluation below).
   - C.4. In some cases, we used [ntJoin](https://github.com/bcgsc/ntJoin) to further scaffold the final genome assembly from MaSuRCa pipeline using the Supernova assembly, generating the final version. 


## Chromosome level assemblies
We generated HiC libraries for 15 ant species, which we further used to curate, scaffold and ultimately resolve the expected number of chromosomes in these genomes. The methods and scripts used are in [HiC_scaffolding](HiC_scaffolding) folder.


## Quality evaluation
Genome assembly quality was evaluated using contiguity metrics, and gene completeness with BUSCO v5.1.2, for each of the above mentioned steps to ensure that the genome assembly quality improved. The scripts used can be found in [Quality_evaluation folder](Quality_evaluation). 

The final genome assemblies were screened and filtered for putative duplicated scaffolds using Funannotate “clean” pipeline v1.8.3. In addition, [purge_haplotigs](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/) was used to purge haplotigs in assemblies that showed high duplication rates (Note that it was only the case in three assemblies). 

Next, the genome assemblies were screened for putative contaminations from other organisms using a pipeline established and optimized for the ant genomes in the GAGA project described here: https://github.com/dinhe878/GAGA-Metagenome-LGT.  In brief, we compiled separate databases containing 1908 complete bacterial genome sequences, 43 complete insect genome sequences, as well as databases containing corresponding bacterial or insect CDS sequences. We divided the genome assembly in 2000 bp sliding windows (overlap 500 bp) and searched each window against the different insect and bacterial databases using mmseqs (release_12-113e3) and identified the single best hit (according to bitscore) for each sliding window against each database. For each scaffold, we calculated the ratio of windows showing higher similarity to bacterial than to eukaryotic databases and used this, along with coverage and GC content information, as lines of evidence to identify contaminant scaffolds.

Finally, the genome assemblies were renamed and validated using again contiguity metrics, gene completeness, and consensus quality (QV) and k-mer completeness using Merqury as shown in the [summary table](GAGA_genome_stats.xlsx). 


### Species validation and contamination checking

Finally, some mitochondrial (CO1 and CytB) and autosomal gene markers (Wingless, LwRh, AbdA, ArgK, CAD) were annotated in the genome assembly after each of the above mentioned steps using [BITACORA](https://github.com/molevol-ub/bitacora) to confirm the species identity of the sequenced data, which allowed to identify a few contaminations and swaps in the sequenced DNA. See [here](Species_barcoding/GAGA_barcoding_species_confirmation.xlsx) the table with the species name validation for the GAGA sequenced ants.



