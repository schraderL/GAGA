# Genome evolution

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Build Status](https://travis-ci.org/username/AwesomeProject.svg?branch=main)](https://travis-ci.org/username/AwesomeProject)

## Outline

- [Chromosome-level synteny](#Chromosome-level)
- [Conserved microsynteny](#Conserved)
- [Rearrangement breakpoint rate](#Rearrangement)
- [Gene co-expression network analysis](#Gene)

## Chromosome-level synteny

[Lasz](https://github.com/lastz/lastz) was utilized to conduct pairwise whole genome alignments, while the [genome alignment tool](https://github.com/hillerlab/GenomeAlignmentTools) was employed to filter and chain these alignments. The Lastz pipeline accommodates various parameters based on phylogenetic distance to achieve optimal alignments.
1. Run Lastz pipeline：

    ```bash
    perl step1_lacnem_v2.pl target.fa query.fa
    cd output/5.net
    ```

2. Filter syntenic net with length < 500kb：

    ```bash
    mkdir netFilter; cd netFilter;
    sh step2_Net_filter.sh ../UCSC.target.filter.net ../../query.sizes
    ```

3. Macro-synteny plot with ideogram.R:

    ```bash
    Rscript ideogram.R # the inputs are from step2.
    ```

## Conserved microsynteny

[Synphoni](https://github.com/nsmro/synphoni) was used to detect conserved  microsyntenic blocks, in which we performed 4 steps:
1. Run Synphoni pipeline：

    ```bash
    sh Synphoni_step1/step2/step2.5/step3/step4.sh
    ```

## Rearrangement breakpoint rate

[AGORA](https://github.com/DyogenIBENS/Agora) was used to re-construct gene order in ancestors, and identify rearrangement breakpoints.

1. Orthology re-construction. Please refer to [our GAGA orthology analysis](https://github.com/schraderL/GAGA/tree/main/06_Analyses/Orthology) for details


2. running AGORA pipeline:

    ```bash
    python3 software/agora/src/agora-generic.py GAGA_phylogeny.nhx Orthogroups/orthologyGroups.%s.list genes/genesSTE.%s.list.bz2 -workingDir=output
    ```

3. detect rearrangement breakpoints and calucate the breakpoint rate:

    ```bash
    sh step3_breakpoints.sh; sh step4_timeNorm.sh
    ```

## Gene co-expression network analysis

[WGCNA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559) was used in gene co-expression network analysis.

1. Running WGCNA on the Monomorium pharaonis developmental transcriptome [Qu et al., 2022](https://www.nature.com/articles/s41559-022-01884-y), we plotted the general expression dynamics of gene co-expression modules across developmental stages using Figure2C.R, the expression dynamics of conserved syntenic blocks using Figure2F.R, and the co-expression networks using Figure2D.R.

    ```bash
    Rscript WGCNA.R 
    ```
