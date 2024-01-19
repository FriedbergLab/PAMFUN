## ---------------------------
## Purpose of script: Adds labels to fusion interactions and calculates the saturation of the fusion profiles.
## Author: Henri Chung
## Date Created: 2022-06-26
## Date Modified: 2024-01-04
## ---------------------------

# Load libaries
library(data.table)
library(tidyverse)
# clear working directory
rm(list = ls())
source("compare/code/fusion_link_helper.R")

# read in dtm data for each organism
balanced_fusion_dtms <- readRDS("compare/data/balanced_fusion_dtms.RDS")

# read in kegg links
module_links <- data.table::fread("data/kegg/protein_links_df.csv")
setkey(module_links, "string")

# convert to dtm to long format
internal_long <- lapply(balanced_fusion_dtms, dist_to_long, cols = colnames(balanced_fusion_dtms[[1]]))
fusion_long_dt <- rbindlist(internal_long, idcol = "assembly_accession2", use.names = TRUE)

# annotate long data with correct pathway labels
long_data <- fusion_long_dt[,c("assembly_accession2", "link", "sim")]
setkey(long_data, "link")
long_data_truth <- long_data[link %in% module_links$string]
long_data_false <- long_data[!(link %in% module_links$string)]
long_data_truth$truth <- TRUE
long_data_false$truth <- FALSE
long_data <- rbind(long_data_truth, long_data_false)
links <- strsplit(long_data$link, split = "-")
Var1 <- sapply(links, function(x) x[1])
Var2 <- sapply(links, function(x) x[2])
long_data[, `:=`(Var1 = Var1, Var2 = Var2)]
long_data <- long_data[order(link)]
long_data$truth <- factor(long_data$truth, levels = c(TRUE, FALSE))
saveRDS(long_data, "compare/outputs/fusion_links_long.RDS")

# calcaulate the saturation of the fusion profiles for the full dataset.
calculate_saturation <- function(cols, dt){
    temp <- as.matrix(dt)[,cols]
    rownames(temp) <- rownames(dt)
    sat <- rowSums(temp)/ncol(temp)
    names(sat) <- rownames(dt)
    return(sat)
}

balanced_fusion_dt <- rbindlist(balanced_fusion_dtms, use.names = TRUE)
rownames(balanced_fusion_dt) <- unlist(lapply(balanced_fusion_dtms, rownames))
rownames <- unlist(lapply(balanced_fusion_dtms, rownames))
profile_saturation <- calculate_saturation(colnames(balanced_fusion_dt), dt = balanced_fusion_dt)
saveRDS(profile_saturation, "compare/outputs/fusion_profile_saturation.RDS")