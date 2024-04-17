## ---------------------------
## Purpose of script: Calculates the precision and recall of mmseqs/fusion performance at variable lengths.
## Author: Henri Chung
## Date Created: 2022-06-26
## Date Modified: 2024-01-04
## ---------------------------

# load libraries
library(data.table)
library(tidyverse)
library(Rcpp)
# clear working directory
rm(list = ls())
source("compare/code/fusion_link_helper.R")
rm(precision_recall)
sourceCpp("compare/code/precision_recall.cpp")
args <- commandArgs(trailingOnly = TRUE)
ind <- as.numeric(args[1])

# read in balanced dtms
balanced_fusion_dtms <- readRDS("compare/data/balanced_fusion_dtms.RDS")
balanced_fusion_dt <- rbindlist(balanced_fusion_dtms, use.names = TRUE)
rownames(balanced_fusion_dt) <- unlist(lapply(balanced_fusion_dtms, rownames))

# read in samples assemblies for varying lengths
sample_assemblies0 <- readRDS("compare/data/balanced_fusion_sample_assemblies.RDS")

name <- names(sample_assemblies0)[[ind]]
sample_assemblies <- list(sample_assemblies0[[ind]])
message(name)
pred_name <- paste0("compare/data/fusion_variable_length_preds_", name, ".RDS")
sat_name <- paste0("compare/data/fusion_variable_length_sat_", name, ".RDS")
# calculate fusion preds at varying lengths
if (!file.exists(pred_name) | !file.exists(sat_name)) {
  # Your code here
  preds_list <- list()
  sat_list <- list()
  for(i in 1:length(sample_assemblies)){
      temp_preds_list <- list()
      temp_sat_list <- list()
      for(j in 1:length(sample_assemblies[[i]])){
          temp_assemblies <- sample_assemblies[[i]][[j]]
          if(length(temp_assemblies) == 1393 & j != 1){
              temp_preds_list[[j]] <- NA
              temp_sat_list[[j]] <- NA
          }else{
              temp_preds_list[[j]] <- generate_preds(sample_assemblies[[i]][[j]], dtm = balanced_fusion_dtms, preds_only = TRUE)
              message(Sys.time(), " size ", i, " sample ", j)
              #
              temp_mat <- balanced_fusion_dt[,..temp_assemblies]
              temp_res <- rowSums(temp_mat)/ncol(temp_mat)
              names(temp_res) <- rownames(balanced_fusion_dt)
              temp_sat_list[[j]] <- temp_res
          }
      }
      names(temp_preds_list) <- paste("sample", 1:length(temp_preds_list), sep = "_")
      preds_list[[i]] <- temp_preds_list
      names(temp_sat_list) <- paste("sample", 1:length(temp_sat_list), sep = "_")
      sat_list[[i]] <- temp_sat_list
  }
  # save the predictions using fusion variable lengths
  names(preds_list) <- lapply(sample_assemblies, function(.x){length(.x[[1]])})
  saveRDS(preds_list, pred_name)
  # save the protein saturation using fusion variable lengths
  names(sat_list) <- lapply(sample_assemblies, function(.x){length(.x[[1]])})
  saveRDS(sat_list, sat_name)
}



compare_variable_lengths <- function(preds, saturation, long){
  long <- long[order(link)]
  results_list <- list(); l = 1
  names_list <- list()
  lengths <- names(preds)
  for(i in 1:length(lengths)){
    for(j in 1:length(preds[[i]])){
      long$sim <- preds[[i]][[j]][1:nrow(long)]
      profile_saturation <- saturation[[i]][[j]]
      for (k in c(0.5, 0.6, 0.7, 0.8, 0.9, 1)) {
        low_sat_proteins <- profile_saturation[profile_saturation <= k] %>% names()
        temp <- copy(long)
        temp[, sim2 := ifelse(!(Var1 %in% low_sat_proteins) | !(Var2 %in% low_sat_proteins), 0, sim)]
        results_list[[l]] <- precision_recall(as.logical(temp[["truth"]]), temp[["sim2"]]) %>% as.data.frame()
        names_list[[l]] <- paste(i, j, k)
        message(Sys.time(), " ", l)
        l = l+1
      }
    }
  }
  names(results_list) <- names_list[1:length(results_list)]
  return(results_list)
}

# Calculate performance of the fusion and mmseqs profiles at different lengths.
message(Sys.time(), " loading data.")
fusion_long <- readRDS("compare/outputs/fusion_links_long.RDS")
fusion_variable_length_preds <- readRDS(pred_name)
fusion_variable_length_sat <- readRDS(sat_name)
fusion_vl_results <- compare_variable_lengths(preds = fusion_variable_length_preds, saturation = fusion_variable_length_sat, long = fusion_long)
rm(fusion_variable_length_preds, fusion_variable_length_sat, fusion_long); gc()
saveRDS(fusion_vl_results, paste0("compare/outputs/fusion_vl_results_", name,".RDS"))
quit()
