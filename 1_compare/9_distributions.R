## ---------------------------
## Purpose of script: Plots the distribution of sequence identity values for phylogenetic profiles
## Author: Henri Chung
## Date Created: 2022-09-21
## Date Modified: 2024-01-04
## ---------------------------

# load libraries
library(data.table)
library(tidyverse)
library(parallel)
# clear working directory
rm(list = ls())


# read in required data
message(Sys.time(), " reading in data")
balanced_fusion_dtms <- readRDS("compare/data/balanced_fusion_dtms.RDS")
balanced_mmseqs_dtms <- readRDS("compare/data/balanced_mmseq_dtms.RDS") 


# are fusion profiles more or less filled than mmseqs dtms?
fusion_rowsums <- lapply(balanced_fusion_dtms, function(.x){rowSums(.x)}) %>% unlist() 
mmseqs40_rowsums <- lapply(balanced_mmseqs_dtms, function(.x){if(is.null(.x)){return(NULL)}else{rowSums(.x > 40)}}) %>% unlist() 
mmseqs60_rowsums <- lapply(balanced_mmseqs_dtms, function(.x){if(is.null(.x)){return(NULL)}else{rowSums(.x > 60)}}) %>% unlist()
mmseqs80_rowsums <- lapply(balanced_mmseqs_dtms, function(.x){if(is.null(.x)){return(NULL)}else{rowSums(.x > 80)}}) %>% unlist()
mmseqs100_rowsums <- lapply(balanced_mmseqs_dtms, function(.x){if(is.null(.x)){return(NULL)}else{rowSums(.x > 100)}}) %>% unlist()

missing <- names(fusion_rowsums)[!(names(fusion_rowsums) %in% names(mmseqs40_rowsums))]
strsplit(missing, split = "\\.") %>% lapply(function(.x){.x[[1]]}) %>% unlist() %>% unique()
summary(fusion_rowsums)
summary(mmseqs40_rowsums)
summary(mmseqs60_rowsums)
summary(mmseqs80_rowsums)
summary(mmseqs100_rowsums)

# plot as histograms
p1_data <- data.frame(fusion = fusion_rowsums, mmseqs40 = mmseqs40_rowsums) %>% pivot_longer(everything(), names_to = "method", values_to = "rowsums")
pdf("compare/outputs/rowsums_histogram.pdf")
p1 <- ggplot(p1_data, aes(x = rowsums, fill = method)) + geom_histogram(position = "identity", alpha = 0.5, bins = 50) + facet_wrap(~method) + theme_bw(); p1
dev.off()
#
fusion_df <- fusion_rowsums %>%
    stack() %>%
    rename(fusion_rowsums = values)

mmseqs40_df <- mmseqs40_rowsums %>%
    stack() %>%
    rename(mmseqs40_rowsums = values) 

compare_df <- left_join(fusion_df, mmseqs40_df, by = "ind")

p2 <- ggplot(compare_df, aes(x = mmseqs40_rowsums, y = fusion_rowsums)) + geom_point() + theme_bw() + xlab("MMseqs2 - 40") + ylab("Fusion") + labs(title = "Filled elements in profile")
ggsave("compare/outputs/rowsums_scatter.pdf", p2, width = 10, height = 10, units = "in")