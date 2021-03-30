#Gene Annotation of *P. californicus* (pleometrotic population) reference genome assembly


<!-- TOC START min:1 max:3 link:true asterisk:false update:true -->
- [Installation of funannotate](#installation-of-funannotate)
  - [install miniconda on jgpogo](#install-miniconda-on-jgpogo)
  - [Start up funannotate and see if everything works](#start-up-funannotate-and-see-if-everything-works)
- [Run Funannotate](#run-funannotate)
  - [Setup environment](#setup-environment)
  - [Preprocess genome assembly](#preprocess-genome-assembly)
  - [Soft Mask genome](#soft-mask-genome)
  - [Train gene annotation programs](#train-gene-annotation-programs)
    - [Prepare GeMoMa input](#prepare-gemoma-input)
    - [Run GeMoMa](#run-gemoma)
    - [Combine GeMoMa predictions](#combine-gemoma-predictions)
    - [Clean gemoma joint gff](#clean-gemoma-joint-gff)
  - [Annotate with funannotate predict](#annotate-with-funannotate-predict)
  - [Processing annotation with funannotate update & fix](#processing-annotation-with-funannotate-update--fix)
    - [funannotate update](#funannotate-update)
    - [funannotate fix](#funannotate-fix)
  - [Functional Annotation](#functional-annotation)
    - [InterProScan](#interproscan)
    - [Orthogroup inference with eggnog](#orthogroup-inference-with-eggnog)
  - [Run funannotate annotate](#run-funannotate-annotate)
  - [Filter Annotations with funannotate](#filter-annotations-with-funannotate)
    - [Funannotate fix](#funannotate-fix-1)
- [Final annotation](#final-annotation)
  - [Overview](#overview)
  - [BUSCOv3 runs](#buscov3-runs)
    - [Output for busco against *unfiltered* annotation](#output-for-busco-against-unfiltered-annotation)
    - [Output for busco against *filtered* annotation](#output-for-busco-against-filtered-annotation)
  - [GAG](#gag)
    - [Output](#output)
- [Supplementary](#supplementary)
  - [calculate RNAseq coverage](#calculate-rnaseq-coverage)
  - [Find homology against uniref90 with mmseqs2](#find-homology-against-uniref90-with-mmseqs2)
  - [Count functionally annotated genes](#count-functionally-annotated-genes)
<!-- TOC END -->



We used funannotate (https://funannotate.readthedocs.io) v1.7.1 to annotate protein coding genes in the minION long-read sequencing-based genome assembly of *Pogonomyrmex californicus*. Details on the assembly can be found here: [../assembly/assembly.md](../assembly/assembly.md)

# Installation of funannotate

## install miniconda on jgpogo
```bash

cd
wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
/bin/bash ~/miniconda.sh -b -p /global/homes/jg/schradel/conda/
# reset PYTHONPATH
PYTHONPATH=""

/global/homes/jg/schradel/conda/bin/conda config --add channels defaults
/global/homes/jg/schradel/conda/bin/conda config --add channels bioconda
/global/homes/jg/schradel/conda/bin/conda config --add channels conda-forge
/global/homes/jg/schradel/conda/bin/conda create -n funannotate python=2.7 funannotate
/global/homes/jg/schradel/conda/bin/conda install -c bioconda eggnog-mapper
```

## Start up funannotate and see if everything works
```bash
#start up conda ENV
/global/homes/jg/schradel/conda/bin/conda init bash
# then restart shell
conda install -c etetoolkit ete3
conda install -c bioconda eggnog-mapper
conda install proteinortho
conda activate funannotate

#check that all modules are installed
funannotate check --show-versions

#download/setup databases to a writable/readable location
funannotate setup -d $HOME/funannotate_db
funannotate setup -b hymenoptera

#set ENV variable for $FUNANNOTATE_DB
echo "export FUNANNOTATE_DB=$HOME/funannotate_db" > /global/homes/jg/schradel/conda/envs/funannotate/etc/conda/activate.d/funannotate.sh
echo "unset FUNANNOTATE_DB" > /global/homes/jg/schradel/conda/envs/funannotate/etc/conda/deactivate.d/funannotate.sh

#run tests -- requires internet connection to download data
export FUNANNOTATE_DB=$HOME/funannotate_db
export GENEMARK_PATH=/global/homes/jg/schradel/software/gmes_linux_64/
export PATH=$PATH:/global/homes/jg/schradel/software/gmes_linux_64/
#https://pypi.org/project/funannotate/
cd /global/homes/jg/schradel/software/gmes_linux_64/
sed -i 's@/usr/bin/perl@/usr/bin/env perl@g' *.pl

ln -s /global/homes/jg/schradel/software/signalp-4.1/signalp /global/homes/jg/schradel/conda/envs/funannotate/bin/
ln -s /global/homes/jg/schradel/software/gmes_linux_64/gmes_petap.pl /global/homes/jg/schradel/conda/envs/funannotate/bin/
funannotate test -t all --cpus 20
funannotate test -t compare --cpus 20
```

# Run Funannotate
-------------------------------------------------------
## Setup environment
Define some basic variables and select the genome assembly fasta file to process. Here, I used the final assembly produced by combining two assemblies with ```ntJoin```.

```bash
base=/global/homes/jg/schradel/data/Pcal/annotation/funannotate

# raw genome fasta
genome=/global/homes/jg/schradel/data/Pcal/minION/assemblyJoined/1/target.fa.k32.w500.n1.all.scaffolds.fa
ln -s $genome ./Pcal.genome.fa
```
-------------------------------------------------------
## Preprocess genome assembly
Here, duplicate scaffolds are removed and the scaffolds are renamed and sorted by size.
```bash
# remove duplicate scaffolds
funannotate clean -i Pcal.genome.fa -o Pcal.genome.clean.fa
# rename and order scaffolds
funannotate sort -i Pcal.genome.clean.fa -o Pcal.genome.clean.sort.fa
```
-------------------------------------------------------
## Soft Mask genome
We used a custom TE annotation pipeline, combining RepeatMasker, RepeatModeler2 and EDTA to soft mask the genome with ```maskFastaFromBed```. Instead, we could also have used funannote mask, but this only superficially masks repeats.
```bash
TEgff=/global/homes/jg/schradel/data/Pcal/annotation/repeat/Pcal/RepeatMasker/Pcal.TEs.gff3
maskFastaFromBed -fi Pcal.genome.clean.sort.fa -bed $TEgff -soft -fo Pcal.genome.clean.sort.masked.fa
readlink -f Pcal.genome.clean.sort.masked.fa
```

-------------------------------------------------------

## Train gene annotation programs
Funannotate was trained with short read, paired-end RNASeq data generated from heads of workers. We also used long-read RNAseq data generated in a MinION.
```bash

# prepare folders
mkdir $base/data
mkdir $base/data/R1
mkdir $base/data/R2
mkdir $base/data/minIonRNAseq

# retrieve RNAseq data
datadir=/global/homes/jg/schradel/data/Pcal/sequencingData/download191029/Schwitte
cd $base/data/R1/
ln -s $datadir/RNA*R1*.fastq.gz .
cd $base/data/R2/
ln -s $datadir/RNA*R2*.fastq.gz .

seqkit stats R1/*
# file                           format  type  num_seqs     sum_len  min_len  avg_len  max_len
# R1/RNA_10_S21_R1_001.fastq.gz  FASTQ   DNA    584,295  43,822,125       75       75       75
# R1/RNA_11_S22_R1_001.fastq.gz  FASTQ   DNA    542,525  40,689,375       75       75       75
# R1/RNA_12_S23_R1_001.fastq.gz  FASTQ   DNA    556,809  41,760,675       75       75       75
# R1/RNA_13_S24_R1_001.fastq.gz  FASTQ   DNA    623,289  46,746,675       75       75       75
# R1/RNA_14_S25_R1_001.fastq.gz  FASTQ   DNA    756,414  56,731,050       75       75       75
# R1/RNA_15_S26_R1_001.fastq.gz  FASTQ   DNA    497,156  37,286,700       75       75       75
# R1/RNA_16_S27_R1_001.fastq.gz  FASTQ   DNA    457,711  34,328,325       75       75       75
# R1/RNA_17_S28_R1_001.fastq.gz  FASTQ   DNA    596,453  44,733,975       75       75       75
# R1/RNA_18_S29_R1_001.fastq.gz  FASTQ   DNA    571,454  42,859,050       75       75       75
# R1/RNA_1_S12_R1_001.fastq.gz   FASTQ   DNA    628,961  47,172,075       75       75       75
# R1/RNA_2_S13_R1_001.fastq.gz   FASTQ   DNA    488,873  36,665,475       75       75       75
# R1/RNA_3_S14_R1_001.fastq.gz   FASTQ   DNA    399,425  29,956,875       75       75       75
# R1/RNA_4_S15_R1_001.fastq.gz   FASTQ   DNA    480,928  36,069,600       75       75       75
# R1/RNA_5_S16_R1_001.fastq.gz   FASTQ   DNA    493,156  36,986,700       75       75       75
# R1/RNA_6_S17_R1_001.fastq.gz   FASTQ   DNA    538,655  40,399,125       75       75       75
# R1/RNA_7_S18_R1_001.fastq.gz   FASTQ   DNA    535,120  40,134,000       75       75       75
# R1/RNA_8_S19_R1_001.fastq.gz   FASTQ   DNA    443,175  33,238,125       75       75       75
# R1/RNA_9_S20_R1_001.fastq.gz   FASTQ   DNA    585,936  43,945,200       75       75       75

# retrieve minIon RNAseq data
cd $base/data/minIonRNAseq
scp -P 50055 juergeng@pyra.uni-muenster.de:/bioinf/seq_data_chestnut/P.californicus/RNAseq_data/Minion/Pcal_RNA_MINION_pass/Pcal_RNAseq_MinION.fastq $base/data/minIonRNAseq

# setup variables
R1=$(echo data/R1/*)
R2=$(echo data/R2/*)
nanopore=$base/data/minIonRNAseq/Pcal_RNAseq_MinION.fastq

# run funannotate train on softmasked genome
funannotate train -i Pcal.genome.clean.sort.masked.fa -o training \
    --left  $R1\
    --right $R2 \
    --nanopore_cdna $nanopore \
    --stranded no --species "Pogonomyrmex californicus" \
    --strain pleo --cpus 20 \
    --memory 25G \
    --max_intronlen 5000
```

-------------------------------------------------------
##Gene Prediction

Homology-based gene prediction was done with GeMoMa against a set of five ant genomes.

### Prepare GeMoMa input
```bash

# download reference genomes and gff

cd /global/homes/jg/schradel/data/genomes/NCBI/
#Pbar
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/187/915/GCF_000187915.1_Pbar_UMD_V03/GCF_000187915.1_Pbar_UMD_V03_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/187/915/GCF_000187915.1_Pbar_UMD_V03/GCF_000187915.1_Pbar_UMD_V03_genomic.gff.gz
#Obir
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/672/135/GCF_003672135.1_Obir_v5.4/GCF_003672135.1_Obir_v5.4_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/672/135/GCF_003672135.1_Obir_v5.4/GCF_003672135.1_Obir_v5.4_genomic.gff.gz
#Cflo
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/227/725/GCF_003227725.1_Cflo_v7.5/GCF_003227725.1_Cflo_v7.5_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/227/725/GCF_003227725.1_Cflo_v7.5/GCF_003227725.1_Cflo_v7.5_genomic.gff.gz
#Sinv
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/075/GCF_000188075.2_Si_gnH/GCF_000188075.2_Si_gnH_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/075/GCF_000188075.2_Si_gnH/GCF_000188075.2_Si_gnH_genomic.gff.gz
#Aech
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/204/515/GCF_000204515.1_Aech_3.9/GCF_000204515.1_Aech_3.9_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/204/515/GCF_000204515.1_Aech_3.9/GCF_000204515.1_Aech_3.9_genomic.fna.gz

#define variables
#Pbar
ref1fa=/global/homes/jg/schradel/data/genomes/NCBI/Pbar/GCF_000187915.1_Pbar_UMD_V03_genomic.fna
ref1gff=/global/homes/jg/schradel/data/genomes/NCBI/Pbar/GCF_000187915.1_Pbar_UMD_V03_genomic.gff
ref1ors=/global/homes/jg/r_stef04/Master/ORgenes/Pbar/final/finalSet/Pbar.OR.gff3

# Aech
ref2fa=/global/homes/jg/schradel/data/genomes/NCBI/Aech/GCF_000204515.1_Aech_3.9_genomic.fna
ref2gff=/global/homes/jg/schradel/data/genomes/NCBI/Aech/GCF_000204515.1_Aech_3.9_genomic.gff

#Cflo
ref3fa=/global/homes/jg/schradel/data/genomes/NCBI/Cflo/GCF_003227725.1_Cflo_v7.5_genomic.fna
ref3gff=/global/homes/jg/schradel/data/genomes/NCBI/Cflo/GCF_003227725.1_Cflo_v7.5_genomic.gff

#Obir
ref4fa=/global/homes/jg/schradel/data/genomes/NCBI/Obir/GCF_003672135.1_Obir_v5.4_genomic.fna
ref4gff=/global/homes/jg/schradel/data/genomes/NCBI/Obir/GCF_003672135.1_Obir_v5.4_genomic.gff

#Sinv
ref5fa=/global/homes/jg/schradel/data/genomes/NCBI/Sinv/GCF_000188075.2_Si_gnH_genomic.fna
ref5gff=/global/homes/jg/schradel/data/genomes/NCBI/Sinv/GCF_000188075.2_Si_gnH_genomic.gff


# RNAseq alignments
rnaseq=/global/homes/jg/schradel/data/Pcal/annotation/funannotate/training/training/transcript.alignments.bam

# which GeMoMa version to use?
gemoma=/global/homes/jg/schradel/software/GeMoMa-1.7/GeMoMa-1.7.jar

# which genome to annotate
genome=/global/homes/jg/schradel/data/Pcal/annotation/funannotate/Pcal.genome.clean.sort.fa

```
### Run GeMoMa
#### Annotate from full reference genome annotation
```bash
#Pbar
ref=Pbar
reffa=$ref1fa
refgff=$ref1gff

base=/global/homes/jg/schradel/data/Pcal/annotation/GeMoMa/$ref/all/
mkdir $base
cd $base

java -Xmx25G -jar $gemoma CLI GeMoMaPipeline threads=25 outdir=$base GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=$genome i=$ref.NCBI a=$refgff g=$reffa GeMoMa.pa=false r=MAPPED ERE.m=$rnaseq

#Aech
ref=Aech
reffa=$ref2fa
refgff=$ref2gff

base=/global/homes/jg/schradel/data/Pcal/annotation/GeMoMa/$ref/all/
mkdir $base
cd $base

java -Xmx25G -jar $gemoma CLI GeMoMaPipeline threads=25 outdir=$base GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=$genome i=$ref.NCBI a=$refgff g=$reffa GeMoMa.pa=false r=MAPPED ERE.m=$rnaseq

#Cflo
ref=Cflo
reffa=$ref3fa
refgff=$ref3gff

base=/global/homes/jg/schradel/data/Pcal/annotation/GeMoMa/$ref/all/
mkdir $base
cd $base

java -Xmx25G -jar $gemoma CLI GeMoMaPipeline threads=25 outdir=$base GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=$genome i=$ref.NCBI a=$refgff g=$reffa GeMoMa.pa=false r=MAPPED ERE.m=$rnaseq

#Obir
ref=Obir
reffa=$ref4fa
refgff=$ref4gff

base=/global/homes/jg/schradel/data/Pcal/annotation/GeMoMa/$ref/all/
mkdir $base
cd $base

java -Xmx25G -jar $gemoma CLI GeMoMaPipeline threads=25 outdir=$base GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=$genome i=$ref.NCBI a=$refgff g=$reffa GeMoMa.pa=false r=MAPPED ERE.m=$rnaseq

#Sinv
ref=Sinv
reffa=$ref5fa
refgff=$ref5gff

base=/global/homes/jg/schradel/data/Pcal/annotation/GeMoMa/$ref/all/
mkdir $base
cd $base

java -Xmx25G -jar $gemoma CLI GeMoMaPipeline threads=25 outdir=$base GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=$genome i=$ref.NCBI a=$refgff g=$reffa GeMoMa.pa=false r=MAPPED ERE.m=$rnaseq

```
### Combine GeMoMa predictions
```bash
base=/global/homes/jg/schradel/data/Pcal/annotation/GeMoMa/GAF/
mkdir $base
cd $base
java -Xmx25G -jar $gemoma CLI GAF \
g=/global/homes/jg/schradel/data/Pcal/annotation/GeMoMa/Aech/all/final_annotation.gff \
g=/global/homes/jg/schradel/data/Pcal/annotation/GeMoMa/Cflo/all/final_annotation.gff \
g=/global/homes/jg/schradel/data/Pcal/annotation/GeMoMa/Obir/all/final_annotation.gff \
g=/global/homes/jg/schradel/data/Pcal/annotation/GeMoMa/Pbar/all/final_annotation.gff \
g=/global/homes/jg/schradel/data/Pcal/annotation/GeMoMa/Sinv/all/final_annotation.gff

readlink -f filtered_predictions.gff
#/global/homes/jg/schradel/data/Pcal/annotation/GeMoMa/GAF/filtered_predictions.gff
```
### Clean gemoma joint gff
```bash
# convert to proper GFF
agat_convert_sp_gxf2gxf.pl -g /global/homes/jg/schradel/data/Pcal/annotation/GeMoMa/GAF/filtered_predictions.gff  -o Pcal.GeMoMa.gff3
a1=$(readlink -f Pcal.GeMoMa.gff3)
```


## Annotate with funannotate predict
`Funannotate predict` uses different programs to predict protein-coding genes.

```bash
base=/global/homes/jg/schradel/data/Pcal/annotation/funannotate
cd $base

funannotate predict -i Pcal.genome.clean.sort.fa --organism other --max_intronlen 10000 \
                --GENEMARK_PATH /global/homes/jg/schradel/software/gmes_linux_64/ \
                -o training -s "Pogonomyrmex californicus" --strain pleo --cpus 20 \
                --other_gff $a1:10 --repeats2evm --busco_db hymenoptera \
                --min_protlen 75 --name PCAL --busco_seed_species hymenoptera


```
--------------------------
## Processing annotation with funannotate update & fix

### funannotate update

```bash
cd $base
funannotate update -i training --cpus 10 --species "Pogonomyrmex californicus"  --strain pleo
```

**Output:**

>Total Gene Models:      26,408
>Total transcripts:      26,940
>New Gene Models:        61
>No Change:              20,219
>Update UTRs:            6,099
>Exons Changed:          29
>Exons/CDS Changed:      0
>Dropped Models:         1
>CDS AED:                0.002
>mRNA AED:               0.032




### funannotate fix

```bash
funannotate fix -i training/update_results/Pogonomyrmex_californicus_pleo.gbk -t training/update_results/Pogonomyrmex_californicus_pleo.tbl
```
-------------------------------------------------------
## Functional Annotation
###  InterProScan
```bash
#@jgpogo
mkdir $base/interpro
pep=$(readlink -f $base/training/update_results/Pogonomyrmex_californicus_pleo.proteins.fa)
cd /global/homes/jg/schradel/dbs/interpro/interproscan-5.44-79.0

./interproscan.sh \
--input $pep \
--disable-precalc \
--output-dir $base/interpro \
--formats TSV,XML,GFF3 \
--goterms

cd /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/interpro

```

### Orthogroup inference with eggnog
```bash
cd $base
download_eggnog_data.py --data_dir ./data/
# create sbt file for NCBI submission
# https://submit.ncbi.nlm.nih.gov/genbank/template/submission/
cd $base
diamondDB=/global/homes/jg/schradel/data/Pcal/annotation/funannotate/data/eggnog_proteins.dmnd
eggnogDataDir=/global/homes/jg/schradel/data/Pcal/annotation/funannotate/data/

emapper.py -m diamond -i training/update_results/Pogonomyrmex_californicus_pleo.proteins.fa --output eggnogOut --cpu 12 --data_dir $eggnogDataDir
```

## Run funannotate annotate
```bash

xml=/global/homes/jg/schradel/data/Pcal/annotation/funannotate/interpro/Pogonomyrmex_californicus_pleo.proteins.fa.xml

funannotate annotate -i training --cpus 20 --sbt Pcal.template.sbt \
--iprscan $xml --busco_db hymenoptera --cpus 20 --eggnog eggnogOut.emapper.annotations

```
## Filter Annotations with funannotate
I filtered all annotations based on RNAseq coverage and functional annotation. Any gene without RNAseq support and without any functional annotation is removed.

### Funannotate fix
```
cd $base
mkdir training/filtering_misc
# get all genes without annotation
cat ./training/annotate_results/Pogonomyrmex_californicus_pleo.annotations.txt|awk -F $'\t' 'BEGIN {OFS = FS} {if ($10=="" && $11=="" && $12=="" && $13=="" && $14=="" && $15=="" && $16=="" && $17=="" &&$18=="" && $19=="" &&$2=="CDS") print $1}'|uniq > training/filtering_misc/no.functionalAnnotation.lst


cat  training/filtering_misc/no.coverage.lst  training/filtering_misc/no.functionalAnnotation.lst|sort|uniq -d > training/filtering_misc/filter.lst

# remove genes that have neither RNAseq coverage nor a functional annotation.
funannotate fix -o training/filtering_results/ -i training/annotate_results/Pogonomyrmex_californicus_pleo.gbk -t training/annotate_results/Pogonomyrmex_californicus_pleo.tbl -d training/filtering_misc/filter.lst
```
--------------------------
# Final annotation
## Overview
The final annotation including peptide, cds, mRNA fastas, a genbank file and a gff is found here:
`/global/homes/jg/schradel/data/Pcal/annotation/funannotate/training/filtering_results`

```bash
cd /global/homes/jg/schradel/data/Pcal/annotation/funannotate/training/filtering_results

Pogonomyrmex_californicus.cds-transcripts.fa
Pogonomyrmex_californicus.gbk
Pogonomyrmex_californicus.gff3
Pogonomyrmex_californicus.mrna-transcripts.fa
Pogonomyrmex_californicus.proteins.fa

```
--------------------------
## BUSCOv3 runs
```bash
base=/global/homes/jg/schradel/data/Pcal/annotation/funannotate
mkdir $base/BUSCO
cd $base/BUSCO
export PATH=$PATH:"/global/projects/programs/source/augustus-3.3.1/bin/"
cp -r /global/projects/programs/source/augustus-3.3.1/config/ .
export AUGUSTUS_CONFIG_PATH=$(readlink -f ./config/)
BUSCO_CONFIG_FILE=/global/projects/programs/source/busco/config/config.ini

# vs unfiltered annotation
python /global/projects/programs/source/busco/scripts/run_BUSCO.py --in $base/training/update_results/Pogonomyrmex_californicus_pleo.proteins.fa --out Pcal.protein.raw --lineage /global/projects/programs/source/busco/hymenoptera_odb9/ --mode prot --cpu 10


# vs filtered annotation
python /global/projects/programs/source/busco/scripts/run_BUSCO.py --in $base/training/filtering_results/Pogonomyrmex_californicus.proteins.fa --out Pcal.protein --lineage /global/projects/programs/source/busco/hymenoptera_odb9/ --mode prot --cpu 10

```
### Output for busco against *unfiltered* annotation

**Target fasta:**
`$base/training/update_results/Pogonomyrmex_californicus_pleo.proteins.fa`

>INFO    Results:
**INFO    C:96.3%[S:91.0%,D:5.3%],F:2.0%,M:1.7%,n:4415**
INFO    4251 Complete BUSCOs (C)
INFO    4016 Complete and single-copy BUSCOs (S)
INFO    235 Complete and duplicated BUSCOs (D)
INFO    89 Fragmented BUSCOs (F)
INFO    75 Missing BUSCOs (M)
INFO    4415 Total BUSCO groups searched
INFO    BUSCO analysis done. Total running time: 270.022828817 seconds
INFO    Results written in /global/homes/jg/schradel/data/Pcal/annotation/funannotate/BUSCO/run_Pcal.protein.raw/

### Output for busco against *filtered* annotation
**Target fasta:**
`$base/training/filtering_results/Pogonomyrmex_californicus.proteins.fa`
>INFO    Results:
**INFO    C:96.1%[S:90.8%,D:5.3%],F:1.9%,M:2.0%,n:4415**
INFO    4244 Complete BUSCOs (C)
INFO    4009 Complete and single-copy BUSCOs (S)
INFO    235 Complete and duplicated BUSCOs (D)
INFO    86 Fragmented BUSCOs (F)
INFO    85 Missing BUSCOs (M)
INFO    4415 Total BUSCO groups searched
INFO    BUSCO analysis done. Total running time: 5720.32516813 seconds
INFO    Results written in /global/homes/jg/schradel/data/Pcal/annotation/funannotate/BUSCO/run_Pcal.protein/

--------------------------

## GAG

```bash
#optimize gff for GAG
agat_convert_sp_gxf2gxf.pl -g Pogonomyrmex_californicus.gff3 -o Pcal.geneAnnotation.v1.0.gff3
gag.py -f $base/Pcal.genome.clean.sort.masked.fa -g Pcal.geneAnnotation.v1.0.gff3 --fix_start_stop
```
### Output
||Reference Genome|Modified Genome|
|-|-|-|
|Total sequence length|252297138|252297138|
|**Number of genes**|**16501**|16501|
|**Number of mRNAs**|**17033**|17033|
|Number of exons|85140|85140|
|Number of introns|68107|68107|
|Number of CDS|15899|15899|
|Overlapping genes|698|698|
|Contained genes|48|48|
|CDS: complete|0|15675|
|CDS: start, no stop|0|110|
|CDS: stop, no start|0|111|
|CDS: no stop, no start|17033|1137|
|Total gene length|68158939|68158939|
|Total mRNA length|71265321|71265321|
|Total exon length|26579488|26579488|
|Total intron length|44822047|44822047|
|Total CDS length|22229451|22229451|
|Shortest gene|63|63|
|Shortest mRNA|63|63|
|Shortest exon|3|3|
|Shortest intron|11|11|
|Shortest CDS|135|135|
|Longest gene|94868|94868|
|Longest mRNA|94868|94868|
|Longest exon|34511|34511|
|Longest intron|56892|56892|
|Longest CDS|40338|40338|
|mean gene length|4131|4131|
|mean mRNA length|4184|4184|
|mean exon length|312|312|
|mean intron length|658|658|
|mean CDS length|1398|1398|
|% of genome covered by genes|27.0|27.0|
|% of genome covered by CDS|8.8|8.8|
|mean mRNAs per gene|1|1|
|mean exons per mRNA|5|5|
|mean introns per mRNA|4|4|

# Supplementary
## calculate RNAseq coverage
```bash

# install tools
# http://lindenb.github.io/jvarkit/

cd ~/software/
git clone "https://github.com/lindenb/jvarkit.git"
cd jvarkit
./gradlew bamstats05

cd /global/homes/jg/schradel/data/Pcal/annotation/coverage



bam=/global/homes/jg/schradel/data/Pcal/annotation/funannotate/training/training/hisat2.coordSorted.bam
gff=/global/homes/jg/schradel/data/Pcal/annotation/funannotate/training/predict_results/Pogonomyrmex_californicus_pleo.gff3
awk '{if ($3=="exon") print $1"\t"$4"\t"$5"\t"$9}' $gff|perl -pe 's/(.*\t).*;Parent=(.*);/$1$2/g' > exons.bed
#bedtools coverage -hist  -a exons.bed -b $bam > out.coverage.tsv

java -jar /global/projects/programs/bin/picard.jar AddOrReplaceReadGroups \
      I=$bam \
      O=output.bam \
      RGPL=Illumina \
      RGLB=lib1 \
      RGPU=xxx \
      RGSM=head

samtools index output.bam
java -jar  ~/software/jvarkit/dist/bamstats05.jar -B exons.bed --merge output.bam > exons.coverage.tsv
cat exons.coverage.tsv|awk '{if ($9==0) print $4}'|cut -f 1 -d "-" > $base/training/filtering_misc/no.coverage.lst

cat exons.coverage.tsv|sed 1d| awk '{if ($9!=0) print $4}'|wc -l
#>9435
cat exons.coverage.tsv|sed 1d| awk '{if ($9!=0) print $4}' > genes.w.coverage.lst
#/global/homes/jg/schradel/data/Pcal/annotation/coverage/genes.w.coverage.lst
```

## Find homology against uniref90 with mmseqs2
mmseqs2 is a brilliant alternative to slow blast. Running it can be a bit mysterious and prone to failure, often due to segmentation faults caused by excessive RAM consumption.

```bash
cd /global/scratch/schradel/BLAST/dbs
mmseqs databases UniRef90 mmseq.UniRef90.db tmp

cd /global/homes/jg/schradel/data/Pcal/annotation/homology
ln -s /global/homes/jg/schradel/data/Pcal/annotation/funannotate/training/predict_results/Pogonomyrmex_californicus_pleo.proteins.fa .
mmseqs createdb Pogonomyrmex_californicus_pleo.proteins.fa mmseq.Pogonomyrmex_californicus_pleo.proteins.db

mmseqs search mmseq.Pogonomyrmex_californicus_pleo.proteins.db /global/scratch/schradel/BLAST/dbs/mmseq.UniRef90.db Pogonomyrmex_californicus_pleo.proteins.UniRef90.results.db tmp -s 1 --split-memory-limit 24G --max-accept 1 --threads 4

mmseqs convertalis mmseq.Pogonomyrmex_californicus_pleo.proteins.db /global/scratch/schradel/BLAST/dbs/mmseq.UniRef90.db Pogonomyrmex_californicus_pleo.proteins.UniRef90.results.db Pogonomyrmex_californicus_pleo.uniref90.m8

# Filter based on e-value
cat Pogonomyrmex_californicus_pleo.proteins.m8|awk '($11 + 0) < 1E-5' |sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > Pogonomyrmex_californicus_pleo.proteins.bh

wc -l Pogonomyrmex_californicus_pleo.proteins.bh
#>15110

# how many genes with RNAseq support had no blast hit?
cut -f 1 Pogonomyrmex_californicus_pleo.proteins.bh|grep -f - /global/homes/jg/schradel/data/Pcal/annotation/coverage/genes.w.coverage.lst -v |wc -l
#> 1279

```


## Count functionally annotated genes
```bash
cd ~/sciebo/Projects/Pcalifornicus/minION/results/annotation/geneAnnotation
awk '{if ($3=="mRNA") print $0}' Pogonomyrmex_californicus.gff3|grep InterPro |perl -pe 's/.*Parent=(.*?)\;.*/$1/g'|sort|uniq |wc -l
#> 10863
awk '{if ($3=="mRNA") print $0}' Pogonomyrmex_californicus.gff3|grep Ontology_term |perl -pe 's/.*Parent=(.*?)\;.*/$1/g'|sort|uniq |wc -l
#> 7500
awk '{if ($3=="mRNA") print $0}' Pogonomyrmex_californicus.gff3|grep PFAM |perl -pe 's/.*Parent=(.*?)\;.*/$1/g'|sort|uniq |wc -l
#> 9532
awk '{if ($3=="mRNA") print $0}' Pogonomyrmex_californicus.gff3|grep EggNog |perl -pe 's/.*Parent=(.*?)\;.*/$1/g'|sort|uniq |wc -l
#> 6218
```
