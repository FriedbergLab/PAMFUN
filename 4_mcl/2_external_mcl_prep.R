## ---------------------------
## Purpose of script: generate external pathway predictions generated from mcl
## Author: Henri Chung
## Date Created: 2023-05-26
## Date Modified: 2024-01-04
## ---------------------------
library(furrr)
library(tidyverse)
library(data.table)
rm(list = ls())

# bind dist calculations into single dt
if(!file.exists("mcl/data/results_dt.csv")){ # fusion
    results_list_files <- list.files("mcl/data", pattern = "fusion_list*", full.names = TRUE)
    results_list <- lapply(results_list_files, function(.x){
        rbindlist(readRDS(.x))
    })
    results_dt <- rbindlist(results_list)
    data.table::fwrite(results_dt, "mcl/data/results_dt.csv", sep = "\t", col.names = FALSE)
    rm(results_list, results_list_files, results_dt); gc()
}

if(!file.exists("mcl/data/mmseq_dt.csv")){ # mmseq
    mmseq_list_files <- list.files("mcl/data", pattern = "mmseqs_list*", full.names = TRUE)
    mmseq_list <- lapply(mmseq_list_files, function(.x){
        data.table::rbindlist(readRDS(.x))
    })
    mmseq_dt <- data.table::rbindlist(mmseq_list)
    data.table::fwrite(mmseq_dt, "mcl/data/mmseq_dt.csv", sep = "\t", col.names = FALSE)
    rm(mmseq_list, mmseq_list_files, mmseq_dt); gc()
}

module_clusters <- readRDS("mcl/data/mcl_targets.RDS")
module_list <- module_clusters %>%
    stack() %>%
    separate(ind, into = c("module_id", "assembly_accession"), sep = "\\.") %>%
    select(-assembly_accession) %>%
    group_by(module_id) %>%
    group_split()
external_modules <- module_list %>% lapply(function(.x){.x$values})
names(external_modules) <- unlist(lapply(module_list, function(.x){.x$module_id[1]}))
saveRDS(external_modules, "mcl/data/external_mcl_targets.RDS")