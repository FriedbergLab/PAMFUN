## ---------------------------
## Purpose of script: Generate internal mcl clusters
## Author: Henri Chung
## Date Created: 2023-05-26
## Date Modified: 2024-01-04
## ---------------------------

library(igraph)
library(tidyverse)
library(data.table)
rm(list = ls())
# Define a list of reference clusters
module_clusters <- readRDS("mcl/data/mcl_targets.RDS")
module_size <- lapply(module_clusters, length) %>% stack() %>% separate(ind, into = c("module_id", "assembly_accession"), sep = "\\.") %>% rename("size" = values)

create_mcl_clusters <- function(.x, i = 2){
  temp_infilename <- paste0("mcl/temp_", .x[1]$Var1[1], ".csv")
  temp_outfilename <- paste0("mcl/temp_", .x[1]$Var1[1], ".mcl")
  data.table::fwrite(.x[,-1], temp_infilename, sep = "\t", col.names = FALSE)
  system(paste0("mcl ",temp_infilename," -q x -V all --abc -I ", i, " -o ", temp_outfilename))
  mcl_clusters <- readLines(temp_outfilename) %>% sapply(strsplit, split = "\\t")
  names(mcl_clusters) <- as.character(1:length(mcl_clusters))
  system(paste0("rm ", temp_infilename))
  system(paste0("rm ", temp_outfilename))
  return(mcl_clusters)
}

# Function to calculate overlap
calculate_jaccard <- function(ref_cluster, input_cluster) { 
  if(length(input_cluster) == 1){return(0)}
  intersect_len <- length(intersect(ref_cluster, input_cluster)); 
  union_len <- length(union(ref_cluster, input_cluster)); 
  return(intersect_len / union_len) 
}

calculate_result <- function(edge_list, n_cores = 4, res = 1){
  result <- parallel::mclapply(edge_list, function(.x){
    input_clusters <- create_mcl_clusters(.x, i = res)
    assembly <- .x$assembly_accession2[1]
    reference_clusters <- module_clusters[grepl(assembly, names(module_clusters))]
    temp_result <- lapply(reference_clusters, function(ref_cluster){
      overlaps <- sapply(input_clusters, calculate_jaccard, ref_cluster = ref_cluster)
      max_overlap_index <- which.max(overlaps)
      max_overlap_value <- max(overlaps)
      return(tibble::tibble("cluster_id" = max_overlap_index, "jaccard" = max_overlap_value))
    }) %>% bind_rows(.id = "reference") %>%
      separate(reference, into = c("module_id", "assembly_accession"), sep = "\\.")
  }, mc.cores = n_cores)
  return(result)
}



seed <- 123
sim_floors <- seq(0, 0.9, 0.1)
ncores = 4
fusion_data <- data.table::fread("mcl/data/fusion_mcl_data.csv")
mmseqs_data <- data.table::fread("mcl/data/mmseqs_mcl_data.csv")
for(i in 1:length(sim_floors)){
  sim_floor <- sim_floors[i]
  fusion_edges <- fusion_data[fusion_sim > sim_floor]
  mmseqs_edges <- mmseqs_data[mmseqs_sim > sim_floor]
  # make random comparisons
  set.seed(seed)
  random_fusion_edges <- fusion_data
  random_fusion_values <- sample(fusion_data$fusion_sim, nrow(fusion_data), replace = FALSE)
  random_fusion_edges$fusion_sim <- random_fusion_values
  set.seed(seed)
  random_mmseqs_edges <- mmseqs_data
  random_mmseqs_values <- sample(mmseqs_data$mmseqs_sim, nrow(mmseqs_data), replace = FALSE)
  random_mmseqs_edges$mmseqs_sim <- random_mmseqs_values
  
  # split into lists
  fusion_edge_list <- split(fusion_edges, by = "assembly_accession2")
  mmseqs_edge_list <- split(mmseqs_edges,  by = "assembly_accession2")
  random_fusion_edge_list <- split(random_fusion_edges, by = "assembly_accession2")
  random_mmseqs_edge_list <- split(random_mmseqs_edges, by = "assembly_accession2")

  louvain_values <- seq(1.2, 4, 0.4)
  edge_lists <- list(fusion_edge_list, mmseqs_edge_list, random_fusion_edge_list, random_mmseqs_edge_list)
  names(edge_lists) <- c("fusion", "mmseqs", "random_fusion", "random_mmseqs")

  message(Sys.time(), " start ", sim_floor, " ", seed)
  results_list <- lapply(edge_lists, function(edge_list){
    lapply(louvain_values, function(louvain_value){
      calculate_result(edge_list, res = louvain_value, n_cores = ncores) %>% bind_rows() %>% mutate(louvain_value = louvain_value)
    }) %>% bind_rows()
  }) %>% bind_rows(.id = "method") %>% mutate(sim_floor = sim_floor)
  message(Sys.time(), " end ", sim_floor, " ", seed)
  saveRDS(results_list, paste0("mcl/outputs/results_list_",seed,"_",sim_floor, ".RDS"))
}

quit()

