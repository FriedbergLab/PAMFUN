## ---------------------------
## Purpose of script: Generate plots for mcl results
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

fusion_edges <- data.table::fread("mcl/data/fusion_mcl_data.csv")[fusion_sim > 0]
mmseqs_edges <- data.table::fread("mcl/data/mmseqs_mcl_data.csv")[mmseqs_sim > 0]
random_fusion_edges <- fusion_edges; random_fusion_edges$fusion_sim <- sample(fusion_edges$fusion_sim, nrow(fusion_edges), replace = FALSE)
random_mmseqs_edges <- mmseqs_edges; random_mmseqs_edges$mmseqs_sim <- sample(mmseqs_edges$mmseqs_sim, nrow(mmseqs_edges), replace = FALSE)

fusion_edge_list <- split(fusion_edges, by = "assembly_accession2")
mmseqs_edge_list <- split(mmseqs_edges,  by = "assembly_accession2")
random_fusion_edge_list <- split(random_fusion_edges,  by = "assembly_accession2")
random_mmseqs_edge_list <- split(random_mmseqs_edges, by = "assembly_accession2")

# Create an igraph graph object from the edge list
results_files_list <- list.files("mcl/outputs", pattern = "results.*RDS", full.names = TRUE); results_files_list
results_list_raw <- lapply(results_files_list[-length(results_files_list)], readRDS)
results_list <- lapply(results_list_raw, function(.x){
  res <- .x %>% 
    left_join(module_size, by = c("module_id", "assembly_accession")) %>%
    mutate(size_group = case_when(
      size < 5 ~ "2-5",
      size >= 5 & size < 20 ~ "6-20",
      size >= 20 ~ "20+")) %>%
    mutate(louvain_value = factor(louvain_value), 
         size_group = factor(size_group, levels = c("2-5", "6-20", "20+"))) %>%
    mutate(method_label = case_when(
      method == "fusion" ~ "Fusion",
      method == "mmseqs" ~ "MMseqs2",
      method == "random_fusion" ~ "Random\nFusion",
      method == "random_mmseqs" ~ "Random\nMMseqs2"
    )) %>%
    mutate(random_label = ifelse(grepl("Random", method_label), "Random", "Non-random"))
  return(res)}) %>% bind_rows()

####################
# Plotting

p0 <- results_list  %>%
  group_by(method, louvain_value, sim_floor, size_group) %>% 
  summarize(median = median(jaccard), n = n(), .groups = "drop") %>%
  ggplot(aes(x = louvain_value, y = median, fill = method)) + 
  geom_bar(stat = "identity", position = "dodge")  +
  facet_grid(sim_floor~size_group) +
  theme_bw()
png("mcl/outputs/p0.png", width = 10, height = 10, units = "in", res = 300)
print(p0)
dev.off()

p1 <- results_list %>%
  group_by(method, sim_floor, louvain_value, size_group) %>%
  arrange(desc(jaccard), .by_group = TRUE) %>%
  filter(louvain_value == 4) %>%
  slice_head(n = 10000) %>%
  ggplot(aes(x = method_label, y = jaccard, fill = method_label)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(sim_floor~size_group) +
  theme_bw() +
  labs(fill = "Method") + ylab("Jaccard Index") + xlab("")
png("mcl/outputs/p1.png", width = 10, height = 10, units = "in", res = 300)
print(p1)
dev.off()


p2 <- results_list %>%
  group_by(method, sim_floor, louvain_value, size_group) %>%
  filter(louvain_value == 4) %>%
  ggplot(aes(x = method_label, y = jaccard, fill = method_label)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(sim_floor~size_group) +
  theme_bw() +
  labs(fill = "Method") + ylab("Jaccard Index") + xlab("")
png("mcl/outputs/p2.png", width = 10, height = 10, units = "in", res = 300)
print(p2)
dev.off()

p3 <- p2$data %>%
  filter(sim_floor == 0.9) %>%
  filter(louvain_value == 4) %>%
  ggplot(aes(x = method_label, y = jaccard, fill = method_label)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~size_group) +
  theme_bw() +
  labs(fill = "Method") + ylab("Jaccard Index") + xlab("")
png("mcl/outputs/p3.png", width = 12, height = 4, units = "in", res = 300)
print(p3)
dev.off()

p4 <- p1$data %>%
  filter(sim_floor == 0.9) %>% 
  filter(louvain_value == 4) %>%
  ggplot(aes(x = method_label, y = jaccard, fill = method_label)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~size_group) +
  theme_bw() +
  labs(fill = "Method") + ylab("Jaccard Index") + xlab("") 
png("mcl/outputs/p4.png", width = 12, height = 4, units = "in", res = 300)
print(p4)
dev.off()

fusion_jaccard <- filter(p4$data, method == "fusion")$jaccard
mmseq_jaccard <- filter(p4$data, method == "mmseqs")$jaccard
t.test(fusion_jaccard, mmseq_jaccard, paired = TRUE)
summary(fusion_jaccard); sd(fusion_jaccard)
summary(mmseq_jaccard); sd(mmseq_jaccard)
p4b <- p4$data %>%
  ggplot(aes(x = size_group, y = jaccard, fill = method_label)) + 
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  labs(fill = "Method") + ylab("Jaccard Index") + xlab("Module Size")
png("mcl/outputs/p4b.png", width = 6, height = 4, units = "in", res = 300)
print(p4b)
dev.off()