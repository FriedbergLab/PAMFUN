## ---------------------------
## Purpose of script: compare saturation levels of best performing fusion and mmseqs2 profiles.
## Author: Henri Chung
## Date Created: 2022-06-26
## Date Modified: 2024-01-04
## ---------------------------

# load libraries
library(data.table)
library(tidyverse)
# clear working directory
rm(list = ls())

# read in required data
message(Sys.time(), " loading data.")
fusion_data <- data.table::fread("data/fusion/fusion_data.tsv")[, assembly_accession2 := gsub("\\.[0-9]", "", assembly_accession)][, n := .N, by = fusion_lvl_1][,ncbi_accession2 := gsub("\\.[0-9]", "", ncbi_accession)]
fusion_long <- readRDS("compare/data/protein_links_chr_dist.RDS")# readRDS("compare/outputs/fusion_long.RDS")
profile_saturation <- readRDS("compare/outputs/fusion_profile_saturation.RDS")
low_sat_proteins <- names(profile_saturation[profile_saturation <= 0.9])
# read in kegg links
only_single_step_modules <- readLines("data/kegg/only_single_step_modules.txt")
module_links <- data.table::fread("data/kegg/protein_links_df.csv")[!(module_id %in% only_single_step_modules)]
setkey(module_links, "string")
module_proteins <- module_links %>% group_split(module_id) %>% lapply(function(.x){unique(c(.x$Var1, .x$Var2))})
module_lengths <- data.table::fread("data/kegg/module_lengths.csv")
names(module_proteins) <- module_links %>% group_split(module_id) %>% lapply(function(.x){unique(.x$module_id)})
module_df <- stack(module_proteins) %>% rename(ncbi_accession2 = values, module_id = ind) %>% left_join(fusion_data[,c("ncbi_accession2", "n")], by = "ncbi_accession2") %>% filter(!is.na(n))
long_data_list <- readRDS("compare/data/filtered_long_data_list.RDS")

# Compare saturation of profiles
balanced_fusion_dtms <- readRDS("compare/data/balanced_fusion_dtms.RDS")
balanced_dt <- rbindlist(balanced_fusion_dtms, idcol = "id", use.names = TRUE)
rownames(balanced_dt) <- unlist(lapply(balanced_fusion_dtms, function(.x){rownames(.x)}))
fusion_profile_saturation <- unlist(lapply(balanced_fusion_dtms, function(.x){rowSums(.x > 0)}))/ncol(balanced_fusion_dtms[[1]])
names(fusion_profile_saturation) <- gsub("GCA.*\\.", "", names(fusion_profile_saturation))
rm(balanced_fusion_dtms)

results_list <- readRDS("compare/data/balanced_mmseq_dtms.RDS")
results_dt <- rbindlist(results_list, , idcol = "id", use.names = TRUE)
rownames(results_dt) <- unlist(lapply(results_list, function(.x){rownames(.x)}))
mmseqs_profile_saturation <- unlist(lapply(results_list, function(.x){rowSums(.x >= 40)}))/ncol(results_list[[1]])
names(mmseqs_profile_saturation) <- gsub("compare.*mmseq\\.", "", names(mmseqs_profile_saturation))
rm(results_list)

fusion_profile_saturation_df <- stack(fusion_profile_saturation) %>% rename(f_sat = "values")
mmseqs_profile_saturation_df <- stack(mmseqs_profile_saturation) %>% mutate(ind = gsub("GCA_[0-9]+\\.", "", ind)) %>% rename(m_sat = "values")

fm_sat <- left_join(fusion_profile_saturation_df, mmseqs_profile_saturation_df , by = "ind") %>% rename(Var1 = "ind") %>% as.data.table()
wilcox.test(fm_sat$f_sat, fm_sat$m_sat)
summary(fm_sat$f_sat); sd(fm_sat$f_sat)
summary(fm_sat$m_sat); sd(fm_sat$m_sat, na.rm = TRUE)

# compare similarities
fusion_stats <- long_data_list[["fusion"]] %>% 
  select(link, sim, truth, Var1, Var2) %>% 
  left_join(module_links, by = c("Var1" = "Var1", "Var2" = "Var2", "link" = "string")) %>%
  mutate(module_id = ifelse(is.na(module_id), "No Module", module_id)) %>%
  left_join(fm_sat, by = "Var1") %>%
  rename(f_sat1 = f_sat, m_sat1 = m_sat) %>%
  left_join(fm_sat, by = c("Var2" = "Var1")) %>%
  rename(f_sat2 = f_sat, m_sat2 = m_sat)  %>%
  mutate(improv = ifelse(module_id %in% fusion_sig_modules, "fusion_improv", "base")) %>%
  mutate(improv = ifelse(improv == "base" & module_id != "No Module", "mmseqs_improv", improv))

fusion_improv_stats <- fusion_stats %>%  filter(module_id %in% fusion_sig_modules) 
fusion_worse_stats <- fusion_stats %>% filter(!(module_id %in% fusion_sig_modules)) 
summary(fusion_improv_stats$sim); summary(fusion_worse_stats$sim)

mmseq_stats <- long_data_list[["mmseqs_40"]] %>% 
  select(link, sim, truth, Var1, Var2) %>% 
  left_join(module_links, by = c("Var1" = "Var1", "Var2" = "Var2", "link" = "string")) %>%
  mutate(module_id = ifelse(is.na(module_id), "No Module", module_id)) %>%
  left_join(fm_sat, by = "Var1") %>%
  rename(f_sat1 = f_sat, m_sat1 = m_sat) %>%
  left_join(fm_sat, by = c("Var2" = "Var1")) %>%
  rename(f_sat2 = f_sat, m_sat2 = m_sat) %>%
  mutate(improv = ifelse(module_id %in% fusion_sig_modules, "fusion_improv", "base")) %>%
  mutate(improv = ifelse(improv == "base" & module_id != "No Module", "mmseqs_improv", improv))

fusion_stats %>% filter(improv == "base") %>% pull(sim) %>% summary()
fusion_stats %>% filter(improv == "fusion_improv") %>% pull(sim) %>% summary()
fusion_stats %>% filter(improv == "mmseqs_improv") %>% pull(sim) %>% summary()

mmseq_stats %>% filter(improv == "base") %>% pull(sim) %>% summary()
mmseq_stats %>% filter(improv == "fusion_improv") %>% pull(sim) %>% summary()
mmseq_stats %>% filter(improv == "mmseqs_improv") %>% pull(sim) %>% summary()

#
fusion_all <- fusion_stats$sim
fusion_modules <- fusion_stats %>% filter(module_id != "No Module") %>% pull(sim)
fusion_fusion_improv <- filter(fusion_stats, improv == "fusion_improv")$sim
fusion_mmseqs_improv <- filter(fusion_stats, improv == "mmseqs_improv")$sim

mmseq_all <- mmseq_stats$sim
mmseq_modules <- mmseq_stats %>% filter(module_id != "No Module") %>% pull(sim)
mmseq_fusion_improv <- filter(mmseq_stats, improv == "fusion_improv")$sim
mmseq_mmseqs_improv <- filter(mmseq_stats, improv == "mmseqs_improv")$sim

p1_values <- list(fusion_all, fusion_modules, fusion_fusion_improv, fusion_mmseqs_improv, mmseq_all, mmseq_modules, mmseq_fusion_improv, mmseq_mmseqs_improv)
names(p1_values) <- c("fusion_all", "fusion_modules", "fusion_fusion_improv",  "fusion_mmseqs_improv", "mmseq_all", "mmseq_modules", "mmseq_fusion_improv", "mmseq_mmseqs_improv")
p1_hist <- lapply(p1_values, function(.x){tabulate(findInterval(.x*100, vec=seq(1,100,length.out=100)), nbins=100)})
p1_prop <- lapply(p1_hist, function(.x){
    .y <- .x / sum(.x)
    names(.y) <- seq(1,100,length.out=100)
    res <- stack(.y) %>%
      rename(Jaccard_Similarity = ind, Frequency = values) %>%
      mutate(method = names(.x))
    return(res)})


p1_data <- bind_rows(p1_prop, .id = "method")
p1 <- p1_data %>% 
  mutate(method =  factor(method, levels =  c("fusion_all", "fusion_modules", "fusion_fusion_improv",  "fusion_mmseqs_improv", "mmseq_all", "mmseq_modules", "mmseq_fusion_improv", "mmseq_mmseqs_improv"))) %>%
  mutate(method_label = case_when(method == "fusion_all" ~ "Fusion",
                                  method == "fusion_modules" ~ "Fusion",
                                  method == "fusion_fusion_improv" ~ "Fusion" ,
                                  method == "fusion_mmseqs_improv" ~ "Fusion",
                                  method == "mmseq_all" ~ "MMseqs2 - 40",
                                  method == "mmseq_modules" ~ "MMseqs2 - 40",
                                  method == "mmseq_fusion_improv" ~ "MMseqs2 - 40",
                                  method == "mmseq_mmseqs_improv" ~ "MMseqs2 - 40")) %>%
  mutate(class_label = case_when(method == "fusion_all" ~ "All",
                                 method == "fusion_modules" ~ "Modules",
                                 method == "fusion_fusion_improv" ~ "Fusion Improved" ,
                                 method == "fusion_mmseqs_improv" ~ "MMseqs2 Improved",
                                 method == "mmseq_all" ~ "All",
                                  method == "mmseq_modules" ~ "Modules",
                                 method == "mmseq_fusion_improv" ~ "Fusion Improved",
                                 method == "mmseq_mmseqs_improv" ~ "MMseqs2 Improved")) %>%
  mutate(class_label = factor(class_label, levels = c("All", "Modules", "Fusion Improved", "MMseqs2 Improved"))) %>%
  mutate(Jaccard_Similarity = as.numeric(Jaccard_Similarity)) %>% 
  ggplot(aes(x = Jaccard_Similarity/100, y = Frequency, fill = method_label)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.5) + 
  facet_wrap(~class_label, scales = "free") +
  theme_bw() + 
  labs(x = "Jaccard Similarity", y = "Frequency", title = "Jaccard Similarity of Identified Proteins", fill = "Method") 
ggsave("compare/outputs/overlapping_histogram_all.png", p1, width = 10, height = 5)

mmseq_improv_stats <- mmseq_stats %>% filter(module_id %in% fusion_sig_modules) 
mmseq_worse_stats <- mmseq_stats %>% filter(!(module_id %in% fusion_sig_modules) & module_id != "No Module")
summary(mmseq_improv_stats$sim); summary(mmseq_worse_stats$sim)

fusion_sim <- fusion_stats %>% select(sim, improv) %>% mutate(method = "fusion")
mmseq_sim <- mmseq_stats %>% select(sim, improv) %>% mutate(method = "mmseqs_40")

pdf("compare/outputs/sim_boxplot.pdf", width = 6, height = 6)
compare_sim <- rbind(fusion_sim, mmseq_sim) %>%
  mutate(improv = ifelse(improv == "fusion_improv", "Fusion Improved", "MMseqs Improved"), method = ifelse(method == "fusion", "Fusion", "MMseqs (40)")) %>%
  ggplot(aes(x = method, y = sim, fill = improv)) + geom_boxplot() +
  theme_bw() + ylim(0, 1) +
  xlab("") + ylab("Jaccard Similarity") + guides(fill = guide_legend(title = "Module Type")) +
  labs(title = "Difference in Jaccard Similarity")
compare_sim
dev.off()
wilcox.test(filter(fusion_sim, improv == "base")$sim, filter(fusion_sim, improv != "base")$sim)
wilcox.test(filter(mmseq_sim, improv == "base")$sim, filter(mmseq_sim, improv != "base")$sim)
tiff("compare/outputs/sim_boxplot.tiff", width = 6, height = 6, units = "in", res = 300)
compare_sim
dev.off()

pdf("compare/outputs/sat_boxplot.pdf", width = 4, height = 6)
compare_sat <- fm_sat %>% 
  pivot_longer(cols = c("f_sat", "m_sat"), names_to = "method", values_to = "sat") %>%
  mutate(method = ifelse(method == "f_sat", "Fusion", "MMseqs (40)")) %>%
  ggplot(aes(x = method, y = sat, fill = method)) + geom_boxplot() +
  ylab("Saturation") + labs(title = "Profile Saturation Comparison") + 
  theme_bw() + guides(fill = FALSE) 
compare_sat
dev.off()
tiff("compare/outputs/sat_boxplot.tiff", width = 6, height = 6, units = "in", res = 300)
compare_sat
dev.off()
wilcox.test(fm_sat$f_sat, fm_sat$m_sat)
