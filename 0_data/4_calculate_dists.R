
## ---------------------------
## Purpose of script: Calculate distance between adjacenet proteins in a fusion link
## Author: Henri Chung
## Date Created: 2022-06-26
## Date Modified: 2024-01-04
## ---------------------------

# load libraries
library(tidyverse)
library(data.table)
library(parallel)
rm(list = ls())

# load data
fusion_long <- readRDS("compare/outputs/fusion_links_long.RDS")
kegg_assemblies <- data.table::fread("data/kegg/kegg_assemblies.csv")
# kegg_balanced_accessions <- readLines("compare/data/balanced_kegg_fusion_accessions")
balanced_organism_accessions <- paste("GCA_", readLines("data/fusion/balanced_organism_accessions"), sep = "")
kegg_balanced_assemblies <- kegg_assemblies[assembly_id2 %in% balanced_organism_accessions]


shared_ind <- unique(fusion_long$assembly_accession2) %in% unique(fusion_long$assembly_accession2)
shared_assemblies <- unique(fusion_long$assembly_accession2)[shared_ind]

all_data <- (fusion_long
  [assembly_accession2 %in% shared_assemblies]
  [!is.na(sim)]
  [, c("Var1", "Var2") := tstrsplit(link, "-", fixed=TRUE)])
setkey(all_data, "link")


# load input data
ft_files <- list.files("data/feature_tables", full.names = TRUE)
ft_ind <- (str_extract(ft_files, "GCA_[0-9]+") %in% kegg_balanced_assemblies$assembly_id2) 
assembly_stats <- data.table::fread("data/assembly_stats/assembly_lengths.tsv")
colnames(assembly_stats) <- c("assembly_accession2", "length")
kegg_assemblies <- data.table::fread("data/kegg/kegg_assemblies.csv")

# get chr distance between any two proteins
chr_dist <- function(.x){
    split <- strsplit(.x, split = "-")[[1]]
    if(split[1] == "NA"){return(NA)}
    if(protein_to_strand[[split[1]]] != protein_to_strand[[split[[2]]]]){return(NA)}
    dist <- abs(protein_to_pos[[split[1]]] - protein_to_pos[[split[2]]])
    length <- protein_to_chr_length[[split[1]]]
    if(dist > length/2){dist <- length-dist}
    return(dist)
}

ft_ind <- (str_extract(ft_files, "GCA_[0-9]+") %in% kegg_balanced_assemblies$assembly_id2) 
ft_list <- list()
for(i in 1:length(ft_files[ft_ind])){ft_list[[i]] <- data.table::fread(ft_files[ft_ind][i], fill = TRUE); message(i)}
to_fix <- ft_files[ft_ind][which(sapply(ft_list, nrow) == 0)]

feature_table <- rbindlist(ft_list)
feature_table <- (feature_table
  [,product_accession2 := gsub("\\.[0-9]", "", product_accession)]
  [product_accession2 != ""]
  [,assembly_accession2 := str_extract(assembly, "GCA_[0-9]+")]
  [assembly_stats, length := i.length, on = "assembly_accession2"]
)
setkey(feature_table, product_accession2)
unique(feature_table[is.na(length)]$assembly_accession2)

# create a hash table for the protein positions
protein_to_pos <- new.env(hash=TRUE)
hash_key <- lapply(1:nrow(feature_table), function(.y){ protein_to_pos[[feature_table$product_accession2[.y]]] <- feature_table$start[.y];})
rm(hash_key)

# create a hash table for the protein strands (+/-)
protein_to_strand <- new.env(hash=TRUE)
hash_key <- lapply(1:nrow(feature_table), function(.y){ protein_to_strand[[feature_table$product_accession2[.y]]] <- feature_table$strand[.y];})
rm(hash_key)

protein_to_chr_length <- new.env(hash=TRUE)
hash_key <- lapply(1:nrow(feature_table), function(.y){ protein_to_chr_length[[feature_table$product_accession2[.y]]] <- feature_table$length[.y];})
rm(hash_key)

Var1 <- as.character(all_data$Var1); head(Var1) 
Var2 <- as.character(all_data$Var2); head(Var2)
pos1 <- sapply(Var1, function(.x){protein_to_pos[[.x]]}); head(pos1)
pos2 <- sapply(Var2, function(.x){protein_to_pos[[.x]]}); head(pos2)
strand1 <- sapply(Var1, function(.x){protein_to_strand[[.x]]}); head(strand1)
strand2 <- sapply(Var2, function(.x){protein_to_strand[[.x]]}); head(strand2)
l <- sapply(Var1, function(.x){protein_to_chr_length[[.x]]}); head(l)
dist <- abs(pos1 - pos2)
dist[dist > l/2] <- abs(l[dist > l/2] - dist[dist > l/2])
ind = strand1 == strand2
all_data$dist <- dist
all_data$valid = ind
saveRDS(all_data, "compare/data/protein_links_chr_dist.RDS")
# quit()



