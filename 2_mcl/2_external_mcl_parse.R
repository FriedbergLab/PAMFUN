## ---------------------------
## Purpose of script: analyze external pathway predictions generated from mcl
## Author: Henri Chung
## Date Created: 2023-05-26
## Date Modified: 2024-01-04
## ---------------------------
library(furrr)
library(tidyverse)
library(data.table)
library(fastmap)
rm(list = ls())

# read in reference data
balanced_data <- data.table::fread("compare/data/balanced_data.csv")
# read in list and remove NA elements
protein_interaction_list <- readRDS("data/kegg/protein_interaction_list.rds") 
protein_interaction_list <- protein_interaction_list[!unlist(lapply(protein_interaction_list, is.logical))]
# reshape table to proteins for each module and org
proteins <- lapply(protein_interaction_list, function(.x){
    proteins <- .x %>% group_by(module_id, org) %>% group_split() %>% lapply(FUN = function(.y){
        as.character(unique(c(.y$Var1, .y$Var2)))
    })
    return(proteins)
}) %>% unlist(recursive = FALSE)

# create new df with module_id, org, and component proteins
protein_interaction_df <- protein_interaction_list %>%
    rbindlist() %>%
    select(module_id, org) %>%
    unique()
protein_interaction_df$proteins <- proteins
protein_interaction_df <- unnest(protein_interaction_df, proteins) %>%
    filter(!is.na(proteins))

# generate protein to fusion fastmap
protein_to_fusion_fastmap <- fastmap::fastmap()
values <- balanced_data$fusion_lvl_1
names(values) <- balanced_data$ncbi_accession2
protein_to_fusion_fastmap$mset(.list = as.list(values))

# generate protein to module fastmap
protein_to_module_fastmap <- fastmap::fastmap()
values <- protein_interaction_df %>% group_by(proteins) %>% group_split() %>% lapply(FUN = function(.x){
    as.character(unique(.x$module_id))
})
names(values) <- unique(protein_interaction_df$proteins)
values2 <- lapply(values, function(.x).x[[1]])
names(values2) <- names(values)
protein_to_module_fastmap$mset(.list = values2)

# read in clusters
system("tr -cd '[:print:]\n' < /work/idoerg/hchung/pamfun2/mcl/outputs/results_dt.hipmcl > /work/idoerg/hchung/pamfun2/mcl/outputs/results_dt2.hipmcl")
system("tr -cd '[:print:]\n' < /work/idoerg/hchung/pamfun2/mcl/outputs/mmseq_dt.hipmcl > /work/idoerg/hchung/pamfun2/mcl/outputs/mmseq_dt2.hipmcl")

fusion_external_clusters_list <- readLines("/work/idoerg/hchung/pamfun2/mcl/outputs/results_dt2.hipmcl") %>% sapply(strsplit, split = "\\t")
fusion_external_clusters <- lapply(fusion_external_clusters_list, function(.x){
    res <- strsplit(.x, split = " ")[[1]]
    res <- res[grepl("^[a-zA-Z]", res)]
    return(res)
    })

mmseq_external_clusters_list <- readLines("/work/idoerg/hchung/pamfun2/mcl/outputs/mmseq_dt2.hipmcl") %>% sapply(strsplit, split = "\\t")
mmseq_external_clusters <- lapply(mmseq_external_clusters_list, function(.x){
    if(length(.x) == 0){
        return(NULL)
    }
    res <- strsplit(.x, split = " ")[[1]]
    res <- res[grepl("^[a-zA-Z]", res)]
    return(res)
    })
# remove clusters with only one protein
fusion_external_clusters <- fusion_external_clusters[sapply(fusion_external_clusters, length) != 1] %>% unname()
mmseq_external_clusters <- mmseq_external_clusters[sapply(mmseq_external_clusters, length) != 1] %>% unname()
external_modules <- readRDS("mcl/data/external_mcl_targets.RDS") 

# Function to calculate overlap
calculate_jaccard <- function(ref_cluster, input_cluster) { 
  if(length(input_cluster) == 1){return(0)}
  intersect_len <- length(intersect(ref_cluster, input_cluster)); 
  union_len <- length(union(ref_cluster, input_cluster)); 
  return(intersect_len / union_len) 
}

fusion_result <- lapply(external_modules, function(ref_cluster){
    overlaps <- sapply(fusion_external_clusters, calculate_jaccard, ref_cluster = ref_cluster)
    max_overlap_index <- which.max(overlaps)
    max_overlap_value <- max(overlaps)
    res <- tibble::tibble("cluster_id" = max_overlap_index, "jaccard" = max_overlap_value)
    return(res)}) %>% bind_rows(.id = "module_id")

mmseq_result <- lapply(external_modules, function(ref_cluster){
    overlaps <- sapply(mmseq_external_clusters, calculate_jaccard, ref_cluster = ref_cluster)
    max_overlap_index <- which.max(overlaps)
    max_overlap_value <- max(overlaps)
    res <- tibble::tibble("cluster_id" = max_overlap_index, "jaccard" = max_overlap_value)
    return(res)}) %>% bind_rows(.id = "module_id")
summary(fusion_result$jaccard); sd(fusion_result$jaccard)
summary(mmseq_result$jaccard); sd(mmseq_result$jaccard)

compare_df <- fusion_result %>% left_join(mmseq_result, by = c("module_id")) %>% rename("fusion" = jaccard.x, "mmseq" = jaccard.y, "fusion_cluster" = cluster_id.x, "mmseq_cluster" = cluster_id.y)
compare_df %>% filter(fusion > mmseq) %>% arrange(desc(fusion))
