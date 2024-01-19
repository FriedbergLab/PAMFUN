## ---------------------------
## Purpose of script: Command line arguments for mmseq mcl clustering with parallelization
## Author: Henri Chung
## Date Created: 2023-05-26
## Date Modified: 2024-01-04
## ---------------------------
# load libraries
message(Sys.time(), " - Loading libraries")
library(data.table)
library(tidyverse)
library(parallel)
library(tictoc)
library(Rcpp)
library(RcppArmadillo)
# clear working directory
rm(list = ls())
sourceCpp("mcl/code/jaccard.cpp")
if(!file.exists("mcl/data/mmseqs_data.RDS")){
  # read in fusion profiles
  balanced_mmseq_dtms <- readRDS("compare/data/balanced_mmseq_dtms.RDS")
  # create vector of unique proteins that have profiles
  unique_proteins <- unlist(lapply(balanced_mmseq_dtms, rownames))
  # bind fusion dtms to one data table
  balanced_dt <- rbindlist(balanced_mmseq_dtms, fill = TRUE)
  balanced_dt[balanced_dt < 40] <- 0
  balanced_dt[balanced_dt >= 40] <- 1
  balanced_dt[is.na(balanced_dt)] <- 0
  protein_saturation <- rowSums(balanced_dt)/ncol(balanced_dt)
  names(protein_saturation) <- unique_proteins
  low_sat_proteins <- names(protein_saturation)[protein_saturation < 0.9]
  all_mat <- as.matrix(balanced_dt)
  rownames(all_mat) <- unique_proteins
  subset_mat <- all_mat[low_sat_proteins,]
  mmseqs_data <- list(subset_mat, low_sat_proteins)
  names(mmseqs_data) <- c("subset_mat", "unique_proteins")
  saveRDS(mmseqs_data, "mcl/data/mmseqs_data.RDS")
  }else{
  message(Sys.time(), " - Loading mat data")
  mmseqs_data <- readRDS("mcl/data/mmseqs_data.RDS")
  all_mat <- mmseqs_data$subset_mat
  unique_proteins <- mmseqs_data$unique_proteins
}

k = 10
l = 1000
floors <- seq(1, length(unique_proteins), by = ceiling(length(unique_proteins)/k))
ceilings <- c(floors[-1], length(unique_proteins))
# command line argument
args <- commandArgs(trailingOnly = TRUE)
ind <- as.numeric(args[[1]])
subset_proteins <- unique_proteins[floors[ind]:ceilings[ind]]
subset_mat <- all_mat[subset_proteins,]

res_list <- list()
input_list <- seq(1, length(unique_proteins), by = l) 
res_list <- parallel::mclapply(input_list, function(.x){
  lower <- .x
  upper <- .x + (l - 1)
  temp_proteins <- unique_proteins[lower:upper]
  temp_proteins <- temp_proteins[!is.na(temp_proteins)]
  temp_mat <- all_mat[temp_proteins,]
  temp_dist <- Jaccard_cpp(temp_mat, subset_mat, rownames(temp_mat), rownames(subset_mat))
  temp_long <- mefa4::Melt(temp_dist) %>% filter(value > 0.9)
  return(temp_long)
}, mc.cores = 36)
# find errors
missing_ind <- sapply(res_list, function(.x){if(is.character(.x)|is.null(.x)){return(TRUE)}else{return(FALSE)}})
missing_list <- input_list[missing_ind]
while(length(missing_list) > 0){
  message(Sys.time(), " - ", length(missing_list))
  res_list[missing_ind] <- parallel::mclapply(missing_list, function(.x){
    message(.x)
    lower <- .x
    upper <- .x + (l - 1)
    temp_proteins <- unique_proteins[lower:upper]
    temp_proteins <- temp_proteins[!is.na(temp_proteins)]
    temp_mat <- all_mat[temp_proteins,]
    temp_dist <- Jaccard_cpp(temp_mat, subset_mat, rownames(temp_mat), rownames(subset_mat))
    temp_long <- mefa4::Melt(temp_dist) %>% filter(value > 0.9)
    return(temp_long)
  }, mc.cores = min(36, length(missing_list)))
  missing_ind <- sapply(res_list, is.character)
  missing_list <- input_list[missing_ind]
}
saveRDS(res_list, paste0("mcl/data/mmseqs_list_", ind, ".RDS"))


