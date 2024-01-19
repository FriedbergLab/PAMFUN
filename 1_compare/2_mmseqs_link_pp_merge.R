## ---------------------------
## Purpose of script: Parses the mmseqs2 search results and converts them into mmseqs2 dtms.
## Author: Henri Chung
## Date Created: 2022-06-26
## Date Modified: 2024-01-04
## ---------------------------

# load libraries
library(data.table)
library(tidyverse)
# clear working directory
rm(list = ls())
source("compare/code/fusion_link_helper.R")
# Read in fusion data
balanced_organism_accessions <- paste("GCA_", readLines("data/fusion/balanced_organism_accessions"), sep = "")
# fusion data
balanced_data <- data.table::fread("compare/data/balanced_data.csv")
module_proteins <- readLines("data/kegg/module_proteins.txt") 

# combine mmseq2 results
mmseq_folders <- list.files(path = "compare/data/mmseqs_results", pattern = "GCA", full.names = TRUE);
mmseq_files <- lapply(mmseq_folders, list.files, pattern = "mmseq$|mmseq\\.[0-9]$", full.names = TRUE) %>% unlist()
names(temp_wide_list) <-  stringr::str_extract(mmseq_files, pattern = "GCA_[0-9]+")
saveRDS(compact(temp_wide_list), "compare/data/balanced_mmseq_dtms.RDS")
subset_data <- balanced_data[ncbi_accession2 %in% rownames(temp_wide_list)]
# order subset data by rownames(temp_wide_list5)
subset_data <- subset_data[match(rownames(temp_wide_list), ncbi_accession2)]
