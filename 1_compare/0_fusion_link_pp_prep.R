## ---------------------------
## Purpose of script: This script generates the fusion dtms for the balanced dataset, and random sub samples.
## Author: Henri Chung
## Date Created: 2023-06-21
## Date Modified: 2024-01-04
## ---------------------------

# Load libaries
library(data.table)
library(tidyverse)

# clear working directory
rm(list = ls())

# load helper functions
source("compare/code/fusion_link_helper.R")

# list of balanced accessions with kegg annotations
balanced_organism_accessions <- paste("GCA_", readLines("data/fusion/balanced_organism_accessions"), sep = "")
# fusion data
fusion_data <- data.table::fread("data/fusion/fusion_data.tsv")[, assembly_accession2 := gsub("\\.[0-9]", "", assembly_accession)][, ncbi_accession2 := gsub("\\.[0-9]", "", ncbi_accession)]
all_proteins <- readLines("compare/data/balanced_assemblies_ids")

# subset to fusions that are not singletons, greater than minimum count, and module proteins.
filtered_fusions <- (fusion_data
  [assembly_accession2 %in% balanced_organism_accessions]
  [ncbi_accession2 %in% all_proteins]
  [fusion_lvl_1 != "S" & fusion_lvl_1 != "SU"] 
  [, by = fusion_lvl_1, .(count = .N) ]
  [count >= 5])

# subset to balanced dataset
balanced_data <- (fusion_data
  [fusion_lvl_1 %in% filtered_fusions$fusion_lvl_1]
  [assembly_accession2 %in% balanced_organism_accessions]
  [ncbi_accession2 %in% all_proteins])
  
# write out balanced data
data.table::fwrite(balanced_data, "compare/data/balanced_data.csv", sep = ",")
setkey(balanced_data, "assembly_accession2")
writeLines(balanced_data$ncbi_accession2, "compare/data/balanced_fusion_ncbi_accession.txt")

# list of proteins involved in kegg modules by ncbi accession
module_proteins <- readLines("data/kegg/module_proteins.txt")
module_data <- balanced_data[ncbi_accession2 %in% module_proteins]

# read in fasta files for organisms in balanced dataset.
kegg_balanced_accessions <- list.files("data/protein_fa", pattern = "*.gz$") %>% gsub("\\..*", "", .)

# loop through each organisms by assembly and calculate the distance between proteins.
results_list <- list()
for(i in 1:length(kegg_balanced_accessions)){

    # subset to module data for current assembly
    query_data <- module_data[assembly_accession2 == kegg_balanced_accessions[i]]
    if(nrow(query_data) == 0){results_list[[i]] <- NULL; next}

    # create hash table of ncbi_accession2 to fusion_lvl_1
    ncbi_to_fusion <- new.env(hash=TRUE)
    hash_key <- lapply(1:nrow(query_data), function(.x){ ncbi_to_fusion[[query_data$ncbi_accession2[.x]]] <- query_data$fusion_lvl_1[.x];})
    rm(hash_key)

    # create dtm for current assembly
    current_dtm <- return_dtm(fusion_rows = query_data$fusion_lvl_1, f_data = balanced_data)
    current_dtm <- current_dtm[, colnames(current_dtm) %in% kegg_balanced_accessions]
    current_bindtm <- current_dtm
    current_bindtm[current_dtm > 0] <- 1
    
    # expand out fusion dtm to ncbi protein rows
    fusion_dtm <- lapply(query_data$ncbi_accession2, function(.y){
        row <- ncbi_to_fusion[[.y]]
        if(!(row %in% rownames(current_bindtm))){
            res <-  as.data.frame(matrix(0, ncol = ncol(current_bindtm)))
            colnames(res) <- colnames(current_bindtm)
        }else{
            res <- current_bindtm[row,]
        }
        return(res)}) %>% 
    bind_rows() %>% as.data.frame()
    rownames(fusion_dtm) <- query_data$ncbi_accession2

    message(Sys.time(), " ", i, " fusion")  
    results_list[[i]] <- fusion_dtm
}
names(results_list) <- kegg_balanced_accessions
results_list <- plyr::compact(results_list)
saveRDS(results_list, "compare/data/balanced_fusion_dtms.RDS")

# write out list of assemblies in final dtm
dtm_assemblies <- colnames(results_list[[1]])
writeLines(dtm_assemblies, "compare/data/balanced_fusion_dtm_assemblies.txt")

protein_names <- unlist(lapply(results_list, rownames))
writeLines(protein_names, "compare/data/balanced_fusion_dtm_protein_names.txt")
# create random sub samples of dtms
dtm_assemblies <- readLines("compare/data/balanced_fusion_dtm_assemblies.txt")
profile_sizes <- c(50, 100, 300, 500, 700, 900, 1100, length(dtm_assemblies))
set.seed(123)
sample_assemblies <- sapply(profile_sizes, function(.x){replicate(sample(dtm_assemblies, .x, replace = FALSE), n = 100, simplify = FALSE)}, simplify = FALSE)
names(sample_assemblies) <- as.character(profile_sizes)
saveRDS(sample_assemblies, "compare/data/balanced_fusion_sample_assemblies2.RDS")

