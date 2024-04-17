
# Pathway Analysis of Microbial Function 
## Table of Contents
- [Introduction](#introduction)
- [Dependencies](#dependencies)
- [Usage](#usage)
    - [0. Data Generation](#0-data-generation)
    - [1. Profile Comparison](#1-profile-comparison)
    - [2. Kaiju Annotation](#2-kaiju-annotation)
    - [3. Mifaser Annotation](#3-mifaser-annotation)
    - [4. MCL Clustering](#4-mcl-clustering)
- [Data](#data)
- [Citation](#citation)

## Introduction
This repository is to document all code used in **Pathway Analysis of Microbial Function** or [pamfun](), a research project that generates and analyzes putatative microbial pathways using [fusion](https://doi.org/10.1093/nar/gkad757). Analysis comprises data generation, descriptive statistics, figure generation, and pathway evaluation. Analysis is intended to run in sequential order from 0-2.

## Abstract
Functional metagenomics enables the study of unexplored bacterial diversity, gene families, and pathways essential to microbial communities. Here, we use a co-occurrence-based analysis of predicted microbial protein functions to uncover pathways in genomic and metagenomic biological systems. Our approach, using phylogenetic profiles, improves identification of functional relationships, or participation in the same biochemical pathway, between enzymes over a comparable homology-based approach. We optimized the design of our profiles to identity potential pathways using minimal data, we clustered functionally related enzyme pairs into multi-enzymatic pathways, and evaluated our predictions against reference pathways in KEGG. We then demonstrated a novel extension of this approach to predict inter-bacterial protein interactions amongst members of a marine microbiome. Our work establishes a basis for identifying the potential functional capacities of entire metagenomes, capturing previously unknown functions into putative pathways.


## Dependencies
All scripts are written in R for version 4.0.3., and require the following packages if used as is. Scripts were written for use in an HPC environment, and may require modification to run on other systems.

### R Packages
- [R](https://www.r-project.org/)
- [tidyverse](https://www.tidyverse.org/)
- [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
- [stringr](https://cran.r-project.org/web/packages/stringr/index.html)
- [ggplot2](https://ggplot2.tidyverse.org/)
- [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html)
- [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)
- [RcppParallel](https://cran.r-project.org/web/packages/RcppParallel/index.html)
- [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf)
- [igraph](https://r.igraph.org/)

### External Software
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- [Kaiju](https://bioinformatics-centre.github.io/kaiju/)
- [Mifaser](https://bromberglab.org/project/mifaser/)
- [MCL](https://micans.org/mcl/)
- [seqkit](https://bioinf.shenwei.me/seqkit/)

## Usage
Project scripts are separated into groups, numbered 0-2, representing the order in which they should be run. A description of each group is provided below.

All scripts are designed to be run with R, in the following manner:
```R
Rscript 0_fetch_kegg.R
```


### 0. Data Generation
Scripts in this group are used to generate data for analysis. This includes downloading and processing data from [KEGG](https://www.genome.jp/kegg/), [GTDB](https://gtdb.ecogenomic.org/), and [NCBI](https://ftp.ncbi.nlm.nih.gov/). Scripts in this group are numbered 0-4.

- ``0_fetch_kegg_helper.R`` - Helper functions for ``0_fetch_kegg.R``
- ``0_fetch_kegg.R`` - Fetches KEGG data, including bacterial organism list, enzyme lists, module and pathway lists.
- ``1_assembly.R`` - Fetches and processes NCBI assembly data for bacterial organisms in KEGG.
- ``2_describe.R`` - Generates descriptive statistics for KEGG and fusion data.
- ``3_taxonomy.R`` - Annotates KEGG/NCBI data with taxonomic information from the GTDB.

### 1. Profile Comparison
Scripts in this group are used to build fusion and mmseqs2 profiles, calculate profile saturation and similarity, and evaluate predictive performance against KEGG. Scripts in this group are numbered 0-10.

- ``0_fusion_link_helper.R`` - Helper functions for ``0_fusion_link_prep.R`` and ``0_fusion_link_base.R``.
- ``0_fusion_link_prep.R`` - Prepares data for fusion profile generation.
- ``0_fusion_link_base.R`` - Adds labels to fusion interactions and calculates the saturation of the fusion profiles.
- ``1_mmseqs_prep.sh`` - Generates mmseqs2 databases for KEGG organisms.
- ``1_mmseq_pp.sh`` - Generates mmseqs2 profiles for KEGG organisms.
- ``2_mmseqs_link_pp_merge.R`` - Parses the mmseqs2 search results and converts them into mmseqs2 occurence matrices.
- ``3_fusion_links.R`` - Generated pathway predictions from fusion profiles, including predictions using profiles of varying lengths.
- ``4_blast_compare.R``- Compare best performing fusion and mmseqs2 profiles.
- ``5_variable_lengths.R`` - Calculate the precision and recall of mmseqs/fusion performance at variable lengths.
- ``6_variable_lengths_plot.R`` - Compare and plot the precision recall performance of fusion and mmseqs at variable lengths and saturations.
- ``7_rand_index.R`` - Calculate rand index of fusion and other functional groupings.
- ``8_sequence_identity.R`` - Calculate the sequence identity between proteins in the same phylogenetic profiles.
- ``9_distributions.R`` - Plot the distribution of sequence identity values for phylogenetic profiles.
- ``10_all_comparison.R`` - Compare and plot the distribution of Jaccard Similarity scores between profiles from different methods.
- ``precision_recall.cpp`` - helper function to calculate precision recall with C++ in R.

### 2. Kaiju Annotation
Scripts in this groups are used to annotate metagenome samples with taxonomy information using [kaiju](https://bioinformatics-centre.github.io/kaiju/).

- ``0_kaiju_analysis.R`` - Parse kaiju output files of metagenome samples.
- ``1_kaiju_metadata.R`` - Parse marine metagenome metadata.

### 3. Mifaser Annotation
Script in this group are used to evaluate the performance of [mifaser](https://bromberglab.org/project/mifaser/) in annotation fusion functions and running mifaser on metagenome samples.

- ``0_test_mifaser.R`` - Sets up a preliminary evaluation of mifaser's ability to annotate fusions functions.
- ``1_full_mifaser_setup.R`` - Builds the mifaser reference database from fusion functions.
- ``2_evaluate_mifaser.R`` - evaluates mifaser performance on annotation fusion functions.
### 4. MCL Clustering
Scripts in this group are used to generate and evaluate mcl clusters of fusion functions based on the Jaccard similarity of their profiles.

- ``1_external_fusion.R`` - Creates fusion mcl clustering with parallelization.
- ``1_external_mmseq.R`` - Creates mmseq mcl clustering with parallelization.
- ``2_external_mcl_parse.R`` - Analyze external pathway predictions generated from mcl.
- ``2_external_mcl_prep.R`` - Generate external pathway predictions with mcl.
- ``3_mcl.R`` - Generate internal pathway predictions with mcl.
- ``3_mcl_plot.R`` - Plot mcl results.
- ``4_all_mcl_parse.R`` - Generate putative marine metagenome pathways over all fusions (known and unknown).
- ``jaccard.cpp`` - Calculate Jaccard Similarity with C++ in R.
- ``jaccard_threshold.cpp`` - Calculate Jaccard Similarity with C++ R, skipping similarities below a specified threshold.

## Data

Bacterial genomes, KEGG, and taxonomy data can be downloaded from their respective sources using the scripts provided. Fusion and intermediate data files can be found at the following figshare links;

Henri Chung "Fusion Phylogenetic Profiling of Microbes", Figshare, 2024,

Fusion Pathway Data. doi:10.6084/m9.figshare.25025231

KEGG Reference Data. doi:10.6084/m9.figshare.25066163

Fusion Predictions. doi:10.6084/m9.figshare.25066478

### Citation
Chung H, Bromberg Y, Friedberg I, Assembling bacterial puzzles: piecing together functions into microbial pathways (2024) bioRxiv 2024.03.27.587058; [doi](https://doi.org/10.1101/2024.03.27.587058)

