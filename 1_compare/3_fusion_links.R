## ---------------------------
## Purpose of script: Generated Fusion predictions from phylogenetic profiles of varying lengths
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
source("/work/idoerg/hchung/pamfun2/compare/code/fusion_link_helper.R")
args <- commandArgs(trailingOnly = TRUE)
ind <- as.numeric(args[1])

# read in balanced dtms
balanced_fusion_dtms <- readRDS("/work/idoerg/hchung/pamfun2/compare/data/balanced_fusion_dtms.RDS")
balanced_fusion_dt <- rbindlist(balanced_fusion_dtms, use.names = TRUE)
rownames(balanced_fusion_dt) <- unlist(lapply(balanced_fusion_dtms, rownames))

# read in samples assemblies for varying lengths
sample_assemblies0 <- readRDS("/work/idoerg/hchung/pamfun2/compare/data/balanced_fusion_sample_assemblies2.RDS")

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
