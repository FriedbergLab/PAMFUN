## ---------------------------
## Purpose of script: General marine metagenome clusters over all fusions.
## Author: Henri Chung
## Date Created: 2023-05-26
## Date Modified: 2024-01-04
## ---------------------------
library(tidyverse)
library(data.table)
rm(list = ls())
# fusion to ec numbers
fusion_data <- data.table::fread("data/fusion/fusion_data.tsv")[,fusion_lvl_1 := as.character(fusion_lvl_1)]

# read in mifaser files
mifaser_files <- list.files("/work/idoerg/hchung/pamfun2/metagenome/mifaser_output", pattern = "analysis.tsv", full.names = TRUE, recursive = TRUE)
mifaser_list <- lapply(mifaser_files, function(.x){
    temp <- (data.table::fread(.x, header = TRUE, col.names = c("fusion_lvl_1", "count"))
        [,fusion_lvl_1 := gsub("^0*", "", gsub("\\.", "",fusion_lvl_1))]
        [,.(count = sum(count)), by = fusion_lvl_1]
        [,percent := count/sum(count) * 100]
        [percent >= 0.01]
        [rev(order(percent))]) 
    return(temp)
})
names(mifaser_list) <- sapply(mifaser_files, function(.x){strsplit(.x, split = "\\/")[[1]][8]})
mifaser_list <- compact(mifaser_list)


#####################################################
# cluster the EC's detected by mifaser
mifaser_dt <- rbindlist(mifaser_list, idcol = "ncbi_sra_accession")
mifaser_dtm <- dcast(mifaser_dt, fusion_lvl_1 ~ ncbi_sra_accession, value.var = "count", fill = 0)
mifaser_mat <- as.matrix(mifaser_dtm[,-c("fusion_lvl_1")])
mifaser_mat[mifaser_mat > 0] <- 1
rownames(mifaser_mat) <-mifaser_dtm$fusion_lvl_1

mifaser_freqs <- rownames(mifaser_mat)[(rowSums(mifaser_mat) > 1)]
mifaser_freqs <- rownames(mifaser_mat)[(rowSums(mifaser_mat) < 0.9*ncol(mifaser_mat))]
mifaser_mat <- mifaser_mat[mifaser_freqs,]
dim(mifaser_mat)

# cluster mifaser dist using mcl
mifaser_dist <- dist(mifaser_mat, method = "binary")
mifaser_pairs <- as.data.frame(as.matrix(mifaser_dist)) %>% 
    mutate(fusion = rownames(mifaser_mat)) %>%
    gather(key = "fusion2", value = "dist", -fusion) %>%
    filter(fusion != fusion2) %>%
    mutate(sim = 1 - dist) %>%
    mutate(combo = paste(fusion, fusion)) %>%
    unique()

data.table::fwrite(mifaser_pairs[,c("fusion", "fusion2", "sim")], "kaiju/outputs/all_mifaser_pairs.csv", sep = "\t")
system("module load mcl/14-137-c2ddjpo")
system("mcl kaiju/outputs/all_mifaser_pairs.csv -I 4 --abc -o kaiju/outputs/all_mcl_mifaser")

all_fusion_clusters <- readLines("kaiju/outputs/all_mcl_mifaser") %>%
    sapply(function(.x){strsplit(.x, split = "\t")[[1]]}) %>%
    unname()

# remove list element that contains "fusion"
all_fusion_clusters <- all_fusion_clusters[!grepl("fusion", all_fusion_clusters)]
lapply(all_fusion_clusters, length) %>% unlist()
# bind list into dataframe with cluster id
all_fusion_clusters_list <- lapply(seq_along(all_fusion_clusters), function(i) {
  data.frame(value = all_fusion_clusters[[i]], cluster = i)
})

# Combine all dataframes into one dataframe
all_fusion_clusters_df <- do.call(rbind, all_fusion_clusters_list) %>%
    rename(fusion_lvl_1 = value) %>%
    as.data.table()

clusters_df <- fusion_data[,c("fusion_lvl_1", "ncbi_accession")][all_fusion_clusters_df, on = "fusion_lvl_1"][!is.na(ncbi_accession)]
data.table::fwrite(clusters_df, "kaiju/outputs/all_mifaser_clusters.tsv", sep = "\t")
table(clusters_df$cluster) %>% length()