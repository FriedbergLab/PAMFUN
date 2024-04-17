## ---------------------------
## Purpose of script: Calculates the sequence identity between proteins in the same phylogenetic profiles
## Author: Henri Chung
## Date Created: 2022-09-21
## Date Modified: 2024-01-04
## ---------------------------

# load libraries
library(data.table)
library(tidyverse)
library(parallel)
# clear working directory
rm(list = ls())

# read in required data
message(Sys.time(), " reading in data")
balanced_fusion_dtms <- readRDS("compare/data/balanced_fusion_dtms.RDS")

# list of balanced accessions with kegg annotations
balanced_organism_accessions <- paste("GCA_", readLines("data/fusion/balanced_organism_accessions"), sep = "")
# fusion data
balanced_data <- data.table::fread("compare/data/balanced_data.csv")
setkeyv(balanced_data, c("ncbi_accession2", "fusion_lvl_1"))

fusions <- lapply(balanced_fusion_dtms, function(.x){rownames(.x)}) %>% unlist() %>% unname()
# group fusions by their fusion_lvl_1 feature in balanced data
fusion_groups <- fusion_data %>%
  filter(ncbi_accession2 %in% fusions) %>% 
  select(fusion_lvl_1, ncbi_accession) %>% 
  group_by(fusion_lvl_1) %>% 
  group_split() %>% lapply(function(.x){.x$ncbi_accession}) %>%
  unlist()

name = 1
system("rm -rf sequence_identity/fusion_tmp/*")
while(length(fusion_groups) > 0){
  message(Sys.time(), " ", length(fusion_groups))
  temp_fusions <- fusion_groups[1:500]
  writeLines(temp_fusions, "identity_test/test.txt")
  temp_data <- balanced_data[fusion_lvl_1 %in% unique(balanced_data[ncbi_accession %in% temp_fusions]$fusion_lvl_1)]
  filename <- paste("sequence_identity/results/", name, "fusion.tsv", sep = "")
  if(!file.exists(filename)){
    p1 <- temp_fusions
    p2 <- temp_data$ncbi_accession
    
    writeLines(unique(p1),"sequence_identity/fusion_tmp/fusionqueryseqs1.txt")
    system("seqkit --quiet grep -f sequence_identity/fusion_tmp/fusionqueryseqs1.txt data/protein_fa/all_seqs.fa > sequence_identity/fusion_tmp/fusionquery1.fasta")
    system("mmseqs createdb sequence_identity/fusion_tmp/fusionquery1.fasta sequence_identity/fusion_tmp/fusionqueryDB1 -v 0")
    #
    writeLines(unique(p2),"sequence_identity/fusion_tmp/fusionqueryseqs2.txt")
    system("seqkit --quiet grep -f sequence_identity/fusion_tmp/fusionqueryseqs2.txt data/protein_fa/all_seqs.fa > sequence_identity/fusion_tmp/fusionquery2.fasta")
    system("mmseqs createdb sequence_identity/fusion_tmp/fusionquery2.fasta sequence_identity/fusion_tmp/fusionqueryDB2 -v 0")
    #
    system("sequence_identity/code/fake_pref.sh sequence_identity/fusion_tmp/fusionqueryDB1 sequence_identity/fusion_tmp/fusionqueryDB2 sequence_identity/fusion_tmp/allvsallpref")
    system("mmseqs align sequence_identity/fusion_tmp/fusionqueryDB1 sequence_identity/fusion_tmp/fusionqueryDB2 sequence_identity/fusion_tmp/allvsallpref sequence_identity/fusion_tmp/allvsallaln -e inf --alignment-mode 3 -v 0")
    system(paste0("mmseqs convertalis sequence_identity/fusion_tmp/fusionqueryDB1 sequence_identity/fusion_tmp/fusionqueryDB2 sequence_identity/fusion_tmp/allvsallaln sequence_identity/fusion_tmp/results.tsv --format-output 'query,target,pident' -v 0"))

    temp_results <- data.table::fread("sequence_identity/fusion_tmp/results.tsv")
    new_results <- list()
    for(i in 1:length(temp_fusions)){
      relevant_fusions <- balanced_data[fusion_lvl_1 %in% unique(balanced_data[ncbi_accession %in% temp_fusions[i]]$fusion_lvl_1)]
      new_results[[i]] <- temp_results[V1 == temp_fusions[i] & V2 %in% relevant_fusions$ncbi_accession]
    }
    new_results2 <- rbindlist(new_results)
    data.table::fwrite(new_results2, filename)
  }
  rm(temp_results)
  fusion_groups <- fusion_groups[!fusion_groups %in% temp_fusions]
  name = name + 1
  system("rm -rf sequence_identity/fusion_tmp/*")
}

