# Script Name: 4_blast_compare.R
# Description: Calculates pair vs group recall evaluation metrics for fusion and mmseqs
# Author: HC
# Date: Feb 14, 2024
####################################################################################################

# load libraries
library(data.table)
library(tidyverse)
# clear working directory
rm(list = ls())

if(!file.exists("compare/outputs/compare_list.RDS")){
  fusion_sims <- readRDS("compare/data/fusion_variable_length_preds_900.RDS")[[1]][[1]]
  fusion_sat <- readRDS("compare/data/fusion_variable_length_sat_900.RDS")[[1]][[1]]
  mmseqs_sims <- readRDS("compare/data/mmseqs_variable_length_preds_40_900.RDS")[[1]][[1]]
  mmseqs_sat <- readRDS("compare/data/mmseqs_variable_length_sat_40_900.RDS")[[1]][[1]]
  compare_list <- list(fusion_sims, fusion_sat, mmseqs_sims, mmseqs_sat)
  names(compare_list) <- c("fusion_sims", "fusion_sat", "mmseqs_sims", "mmseqs_sat")
  saveRDS(compare_list, "compare/outputs/compare_list.RDS")
}else{
  compare_list <- readRDS("compare/outputs/compare_list.RDS")
  fusion_sims <- compare_list$fusion_sims
  fusion_sat <- compare_list$fusion_sat
  mmseqs_sims <- compare_list$mmseqs_sims
  mmseqs_sat <- compare_list$mmseqs_sat
  
}
fusion_long <- readRDS("compare/data/protein_links_chr_dist.RDS")

# read in KEGG data
only_single_step_modules <- readLines("data/kegg/only_single_step_modules.txt")
module_links <- data.table::fread("data/kegg/protein_links_df.csv")[!(module_id %in% only_single_step_modules)]
setkey(module_links, "string")
module_proteins <- module_links %>% group_split(module_id) %>% lapply(function(.x){unique(c(.x$Var1, .x$Var2))})
names(module_proteins) <- module_links %>% group_split(module_id) %>% lapply(function(.x){unique(.x$module_id)})
fusion_data <- data.table::fread("data/fusion/fusion_data.tsv")[, assembly_accession2 := gsub("\\.[0-9]", "", assembly_accession)][, n := .N, by = fusion_lvl_1][,ncbi_accession2 := gsub("\\.[0-9]", "", ncbi_accession)]
# filter proteins by saturation
mmseqs_proteins <- names(mmseqs_sat)[mmseqs_sat < 0.9]
fusion_proteins <- names(fusion_sat)[fusion_sat < 0.9]
intersection <- intersect(mmseqs_proteins, fusion_proteins)

if(!("fusion_sim" %in% colnames(fusion_long))){
compare_df <- (fusion_long
  [,fusion_sim := fusion_sims]
  [,mmseqs_sim := mmseqs_sims]
  [,sim := NULL])
}

# determine precision and recall at f1 max
compare_df <- compare_df[Var1 %in% intersection & Var2 %in% intersection]
fusion_preds <- compare_df[fusion_sim >= 0.93][Var1 %in% fusion_proteins & Var2 %in% fusion_proteins][module_links, org := i.org, on = .(link = string)][module_links, module_id := i.module_id, on = .(link = string)]
mmseqs_preds <- compare_df[mmseqs_sim >= 0.93][Var1 %in% mmseqs_proteins & Var2 %in% mmseqs_proteins][module_links, org := i.org, on = .(link = string)][module_links, module_id := i.module_id, on = .(link = string)]

table(fusion_preds$truth)/nrow(fusion_preds)
table(mmseqs_preds$truth)/nrow(mmseqs_preds)

# custom function to split into groups
module_grouper <- function(.x){
    temp_list <- tibble(
        proteins = list(unique(c(.x$Var1, .x$Var2))), 
        module = .x$module_id[1], 
        org = .x$org[1])
    return(temp_list)
}
modules_groups <- module_links %>% filter(Var1 %in% intersection & Var2 %in% intersection) %>% group_by(module_id, org) %>% group_split() %>% lapply(module_grouper) %>% bind_rows()
fusion_groups <- fusion_preds %>% group_by(module_id, org) %>% group_split() %>% lapply(module_grouper) %>% bind_rows() %>% rename(fusion_proteins = proteins)
mmseqs_groups <- mmseqs_preds %>% group_by(module_id, org) %>% group_split() %>% lapply(module_grouper) %>% bind_rows() %>% rename(mmseqs_proteins = proteins)

compare_groups <- modules_groups %>% 
    left_join(fusion_groups, by = c("module", "org")) %>% 
    left_join(mmseqs_groups, by = c("module", "org"))

for(i in 1:nrow(compare_groups)){
    message(i/nrow(compare_groups))
    temp_fusion <- compare_groups$fusion_proteins[[i]]
    temp_mmseq <- compare_groups$mmseqs_proteins[[i]]
    temp_ref <- compare_groups$proteins[[i]]
    
    fusion_recall <- length(intersect(temp_ref, temp_fusion)) / length(temp_ref)
    mmseqs_recall <- length(intersect(temp_ref, temp_mmseq)) / length(temp_ref)
    compare_groups$fusion_recall[[i]] <- fusion_recall
    compare_groups$mmseqs_recall[[i]] <- mmseqs_recall
}
compare_groups <- compare_groups %>% unnest(fusion_recall) %>% unnest(mmseqs_recall)
compare_groups$l <- lapply(compare_groups$proteins, length) %>% unlist()
compare_groups <- filter(compare_groups, l > 1) 
summary(compare_groups$fusion_recall)
summary(compare_groups$mmseqs_recall)


p1_data <- compare_groups %>%
    pivot_longer(cols = c(fusion_recall, mmseqs_recall), names_to = "method", values_to = "recall") 
p1 <- ggplot(aes(x = l, y = recall, color = method), data = p1_data) + geom_point() + facet_wrap(~method) + labs(title = "Group Recall") + xlab("Module Length") + theme_bw()
ggsave("compare/outputs/groups_recall.png", p1, width = 6, height = 6, units = "in")

p1b_data <- p1_data %>% group_by(l, method) %>% summarise(mean_recall = mean(recall), n = n())
p1b<- ggplot(p1b_data, aes(x = l, y = mean_recall, color = method, size = n)) + geom_point() + facet_wrap(~method) + labs(title = "Group Recall") + xlab("Module Length") + theme_bw()
ggsave("compare/outputs/groups_recall_mean.png", p1b, width = 6, height = 6, units = "in")

p1c_data <- p1b$data %>% filter(method == "fusion_recall") %>% ungroup()
p1c<- ggplot(p1c_data, aes(x = l, y = mean_recall, color = method, size = n)) + ylim(0, 1) + geom_point() + facet_wrap(~method) + labs(title = "Group Recall") + xlab("Module Length") + theme_bw()
ggsave("compare/outputs/groups_recall_mean_fusion_only.png", p1c, width = 6, height = 6, units = "in")

module_links_ref <- module_links %>% filter(Var1 %in% intersection & Var2 %in% intersection & Var1 != Var2) %>% group_by(module_id, org) %>% summarize(module_length = length(unique(c(Var1, Var2))), num_pairs = n()/2)
fusion_preds2 <- fusion_preds %>% filter(truth == TRUE) %>% group_by(org, module_id) %>% summarize(n_preds = n()) %>% right_join(module_links_ref, by = c("module_id", "org")) %>% mutate(recall = n_preds / num_pairs)
mmseqs_preds2 <- mmseqs_preds %>% filter(truth == TRUE) %>% group_by(org, module_id) %>% summarize(n_preds = n()) %>% right_join(module_links_ref, by = c("module_id", "org")) %>% mutate(recall = n_preds / num_pairs)

p2_data <- rbind(fusion_preds2 %>% mutate(method = "fusion"), mmseqs_preds2 %>% mutate(method = "mmseqs")) %>% replace(is.na(.),0) 
p2 <- ggplot(aes(x = num_pairs, y = recall, color = method), data = p2_data) + geom_point() + facet_wrap(~method) + labs(title = "Pair Recall") + xlab("Number of Pairs") + theme_bw()
ggsave("compare/outputs/pairs_recall.png", p2, width = 6, height = 6, units = "in")

p2b_data <- p2_data %>% group_by(module_length, method) %>% summarise(mean_recall = mean(recall), num_modules = n())
p2b <- ggplot(p2b_data, aes(x = module_length, y = mean_recall, color = method, size = num_modules)) + geom_point() + facet_wrap(~method) + labs(title = "Pair Recall") + xlab("Number of Pairs") + theme_bw()
ggsave("compare/outputs/pairs_recall_mean.png", p2b, width = 6, height = 6, units = "in")

p2c_data <- p2b_data %>% filter(method == "fusion") 
p2c <- ggplot(p2c_data, aes(x = module_length, y = mean_recall, color = method, size = num_modules)) + geom_point() + facet_wrap(~method) + labs(title = "Pair Recall") + xlab("Number of Pairs") + ylim(0, 1) + theme_bw()
ggsave("compare/outputs/pairs_recall_mean_fusion_only.png", p2c, width = 6, height = 6, units = "in")

p3_data_a <- p1c_data %>% 
  rename(module_length = l, recall = mean_recall, num_modules = n) %>%
  mutate(method = "group") %>%
  ungroup() %>% select(-.group)

p3_data_b <- p2c_data %>%
  mutate(method = "pair") %>%
  rename(recall = mean_recall) %>%
  ungroup()

colnames(p3_data_a); colnames(p3_data_b)
p3_data <- rbind(p3_data_a, p3_data_b)

p3 <- ggplot(aes(x = module_length, y = recall, shape = method, color = method), data = p3_data) + 
  geom_point() + 
  labs(title = "Recall Comparison of Modules") + 
  xlab("Module Proteins") + ylab("Mean Recall") +
  theme_bw() 
ggsave("compare/outputs/pair_vs_group_recall.png", p3, width = 6, height = 6, units = "in")

p3a <- ggplot(aes(x = module_length, y = recall, size = num_modules, color = method), data = p3_data) + 
  geom_point() + 
  labs(title = "Recall Comparison of Modules") + 
  xlab("Module Proteins") + ylab("Mean Recall") +
  theme_bw() 
ggsave("compare/outputs/pair_vs_group_recall_size.png", p3a, width = 6, height = 6, units = "in")

p3b <- ggplot(aes(x = module_length, y = recall, color = num_modules, shape = method), data = p3_data) + 
  geom_point() + 
  labs(title = "Recall Comparison of Modules") + 
  xlab("Module Proteins") + ylab("Mean Recall") +
  theme_bw() 
ggsave("compare/outputs/pair_vs_group_recall_color.png", p3b, width = 6, height = 6, units = "in")

p3c <- p3_data %>%
  mutate(Method = ifelse(method == "pair", "Pair", "Group")) %>%
  mutate(Method = factor(Method, levels = c("Pair", "Group"))) %>%
  ggplot(aes(x = module_length, y = recall, shape = Method, color = Method)) + 
  geom_point() + 
  labs(title = "", color = "Recall Method", shape = "Recall Method") + 
  xlab("Number of Module Proteins") + ylab("Mean Recall") +
  theme_bw() 
ggsave("compare/outputs/pair_vs_group_recall2.png", p3c, width = 5, height = 5, units = "in")

p3d <- p3_data %>%
  mutate(Method = ifelse(method == "pair", "Pair", "Group")) %>%
  mutate(Method = factor(Method, levels = c("Pair", "Group"))) %>%
  ggplot(aes(y = module_length, x = recall, shape = Method, color = Method)) + 
  geom_point() + 
  labs(title = "", color = "Recall Method", shape = "Recall Method") + 
  ylab("Number of Module Proteins") + xlab("Mean Recall") +
  theme_bw() 
ggsave("compare/outputs/pair_vs_group_recall4.png", p3d, width = 5, height = 3, units = "in")

p3e <- p3_data %>%
  mutate(Method = ifelse(method == "pair", "Pair", "Group")) %>%
  mutate(Method = factor(Method, levels = c("Pair", "Group"))) %>%
  ggplot(aes(y = module_length, x = recall, shape = Method, color = Method)) + 
  geom_point() + 
  labs(title = "", color = "Recall Method", shape = "Recall Method") + 
  ylab("Number of Module Proteins") + xlab("Mean Recall") +
  theme_bw() +
  theme(legend.position = c(0.9, 0.1), 
        legend.justification = c(0.9, 0.1), 
        legend.box.just = "right",
        legend.margin = margin(t = 5, r = 5, b = 5, l = 5), 
        legend.box.margin = margin(), 
        legend.background = element_rect(linetype = "solid", color = "black", fill = "white"))
ggsave("compare/outputs/pair_vs_group_recall5.png", p3e, width = 5, height = 3, units = "in")
