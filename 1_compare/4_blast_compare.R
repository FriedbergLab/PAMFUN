## ---------------------------
## Purpose of script: Detailed comparison of best performing fusion and mmseqs2 profiles.
## Author: Henri Chung
## Date Created: 2022-06-26
## Date Modified: 2024-01-04
## ---------------------------


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
module_df <- stack(module_proteins) %>% rename(ncbi_accession2 = values, module_id = ind) %>% left_join(fusion_data[,c("ncbi_accession2", "n")], by = "ncbi_accession2") %>% filter(!is.na(n))
length(unique(filter(module_df, ncbi_accession2 %in% names(mmseqs_sat))$module_id))

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

table(fusion_preds$truth) / nrow(fusion_preds)
table(mmseqs_preds$truth) / nrow(mmseqs_preds)

# Group predictions based on module 
module_lengths <- data.table::fread("data/kegg/module_lengths.csv")
module_links_count <- unique(module_links[Var1 != "" & Var2 != ""][, count := .N, by = module_id][, .(module_id, count)])
module_classes <- read.csv("data/kegg/module_classes.csv") %>% 
  left_join(module_links_count, by = "module_id") %>% 
  group_by(subclass) %>% 
  mutate(subclass_n = length(unique(module_id))) %>%
  left_join(module_lengths, by = "module_id") %>%
  mutate(class = ifelse(is.na(class), "Signature modules", class)) %>%
  mutate(subclass = ifelse(is.na(subclass), "Signature modules", subclass)) 

# are there any modules that are not predicted by either method?
fusion_modules <- module_links[string %in% fusion_preds[truth == TRUE]$link]$module_id %>% table() %>% stack() %>% rename(module_id = "ind", fusion_values = "values")
mmseqs_modules <- module_links[string %in% mmseqs_preds[truth == TRUE]$link]$module_id %>% table() %>% stack() %>% rename(module_id = "ind", mmseqs_values = "values") 

module_dists <- (fusion_long
  [module_links, module_id := i.module_id, on = .(link = string)]
  [truth == TRUE & valid == TRUE]
  [as.data.table(module_classes), on = "module_id"]
  [, .(mean_dist = mean(dist, na.rm = TRUE), median_dist = median(dist, na.rm = TRUE)), by = module_id]
  [!is.na(median_dist)]
)

# how many modules are predicted by both methods?
unique_modules <- unique(module_links$module_id) 
# how many modules are missing from both methods?
temp_fusion_missing_modules <- unique_modules[!unique_modules %in% unique(fusion_modules$module_id)]
temp_mmseqs_missing_modules <- unique_modules[!unique_modules %in% unique(mmseqs_modules$module_id)]
both_missing <- table(c(temp_fusion_missing_modules, temp_mmseqs_missing_modules)) %>% stack() %>% filter(values == 2) %>% pull(ind) %>% as.character()
fusion_missing <- temp_fusion_missing_modules[!temp_fusion_missing_modules %in% both_missing]
mmseqs_missing <- temp_mmseqs_missing_modules[!temp_mmseqs_missing_modules %in% both_missing]

# sort by number of proteins in module
filter(module_classes, module_id %in% both_missing)$subclass %>% table() %>% stack() %>% arrange(values) 
filter(module_classes, module_id %in% fusion_missing)$subclass %>% table() %>% stack() %>% arrange(values) 
filter(module_classes, module_id %in% mmseqs_missing)$subclass %>% table() %>% stack() %>% arrange(values) 

# any modules found exclusively by fusion? (equivalent to mmseqs_missing)
fusion_exclusive <- unique_modules[unique_modules %in% unique(fusion_modules$module_id) & !unique_modules %in% unique(mmseqs_modules$module_id)]
# how many predicted by fusion exclusively?
fusion_preds[module_id %in% fusion_exclusive] %>% nrow()
# how many are not already in the same fusion?
fusion_preds[module_id %in% fusion_exclusive][fusion_sim != 1]
# list of proteins in fusion exclusive modules
fusion_exclusive_proteins <- unique(c(fusion_preds[module_id %in% fusion_exclusive]$Var1, fusion_preds[module_id %in% fusion_exclusive]$Var2))
# protein name data of fusion exclusive proteins
fusion_data[ncbi_accession2 %in% fusion_exclusive_proteins] 
fusion_exclusive_data <- fusion_data[ncbi_accession2 %in% fusion_exclusive_proteins]

# names of proteins involved in fusion exclusive but arent in the same fusions
apply(fusion_preds[module_id %in% fusion_exclusive][fusion_sim != 1], 1, function(.x){
  res <- fusion_exclusive_data[ncbi_accession2 %in% c(.x[["Var1"]], .x[["Var2"]]),]$protein_name
  sort(res)
}, simplify = FALSE)

#
mmseq_files <- list.files("compare/data/mmseqs_results", full.names = TRUE)
mmseq_assembly_files <- mmseq_files[grepl(paste(fusion_exclusive_not_1_assemblies, collapse = "|"), mmseq_files)]
mmseqs_matches <- lapply(mmseq_assembly_files, function(.x){
  temp_file <- data.table::fread(list.files(.x, full.names = TRUE))[,V1 := gsub("\\.[0-9]", "", V1)][,V2 := gsub("\\.[0-9]", "", V2)][V3 > 40]
  temp_file2 <- temp_file[V1 %in% fusion_exclusive_proteins_not_1] 
  res <- split(temp_file2, temp_file2$V1) %>% lapply(function(.y){.y$V2})
  return(res)
}) %>% unlist(recursive = FALSE)

balanced_data <- data.table::fread("compare/data/balanced_data.csv")
fusion_matches <- lapply(names(mmseqs_matches), function(.x){
  balanced_data[fusion_lvl_1 == balanced_data[ncbi_accession2 == .x]$fusion_lvl_1]$ncbi_accession2
}); names(fusion_matches) <- names(mmseqs_matches)

mmseq_exclusive_proteins <- lapply(names(fusion_matches), function(.x){ mmseqs_matches[.x][[1]][!(mmseqs_matches[.x][[1]] %in% fusion_matches[.x][[1]])]}) 
names(mmseq_exclusive_proteins) <- names(fusion_matches)

mmseq_fusions <- balanced_data[ncbi_accession2 %in% mmseq_exclusive_proteins[[9]]]$fusion_lvl_1 %>% unique()
fusion_fusions <- balanced_data[ncbi_accession2 == names(mmseq_exclusive_proteins)[[9]]]$fusion_lvl_1

balanced_data[ncbi_accession2 %in% names(mmseq_exclusive_proteins)][protein_name == "tocopherol cyclase"]
balanced_data[ncbi_accession2 == names(mmseq_exclusive_proteins)[[9]]]$protein_name
a <- balanced_data[fusion_lvl_1 %in% mmseq_fusions]$assembly_accession2 %>% trimws() %>% unique()
b <- balanced_data[fusion_lvl_1 %in% fusion_fusions]$assembly_accession2 %>% trimws() %>% unique()

# compare taxonomy of fusion exclusive proteins
taxonomy <- (data.table::fread("data/taxonomy/gtdb_reference.csv")
  [, assembly_accession2 := gsub("\\..*", "", assembly_accession)]
  [assembly_accession2 %in% unique(fusion_preds$assembly_accession)]
)
taxonomy[assembly_accession2 %in% a]$phylum
taxonomy[assembly_accession2 %in% b]$phylum

balanced_data[fusion_lvl_1 %in% mmseq_fusions]$protein_name %>% trimws() %>% unique()
balanced_data[fusion_lvl_1 %in% fusion_fusions]$protein_name %>% trimws() %>% unique()
balanced_data[fusion_lvl_1 %in% mmseq_fusions][protein_name == "Tocopherol cyclase"]
# mmseqs exclusive 
mmseqs_exclusive <- unique_modules[unique_modules %in% unique(mmseqs_modules$module_id) & !unique_modules %in% unique(fusion_modules$module_id)]
mmseqs_preds[module_id %in% mmseqs_exclusive] %>% nrow()
mmseqs_preds[module_id %in% mmseqs_exclusive]


# which modules are predicted as "most complete"
# pathway analysis
module_links_org_count <- unique(module_links[, count := .N, by = .(module_id, org)][, .(module_id, org, count)][!is.na(module_id)])
fusion_links_org_count <- unique(fusion_preds[, count := .N, by = .(module_id, org)][, .(module_id, org, count)][!is.na(module_id)])
mmseqs_links_org_count <- unique(mmseqs_preds[, count := .N, by = .(module_id, org)][, .(module_id, org, count)][!is.na(module_id)])

compare_links_org_count <- (module_links_org_count
  [fusion_links_org_count, fusion_count := i.count, on = .(module_id, org)]
  [mmseqs_links_org_count, mmseqs_count := i.count, on = .(module_id, org)]
  [,fusion_perc := fusion_count / count]
  [,mmseqs_perc := mmseqs_count / count]
  [!is.na(fusion_perc) & !is.na(mmseqs_perc)]
  [count > 1]
)

fusion_completeness_improv <- compare_links_org_count[fusion_perc > mmseqs_perc]
mmseqs_completeness_improv <- compare_links_org_count[mmseqs_perc > fusion_perc]
nrow(fusion_completeness_improv)
nrow(mmseqs_completeness_improv)
summary(compare_links_org_count$fusion_perc)
summary(compare_links_org_count$mmseqs_perc)
#
joined_dt <- right_join(mmseqs_modules, fusion_modules, by = "module_id") %>% 
  left_join(module_classes, by = "module_id") %>%
  left_join(module_dists, by = "module_id") %>% 
  mutate(mmseqs_values = ifelse(is.na(mmseqs_values), 0, mmseqs_values),
         fusion_values = ifelse(is.na(fusion_values), 0, fusion_values)) %>%
  mutate(perc = fusion_values / mmseqs_values, fusion_complete = fusion_values / count, mmseqs_complete = mmseqs_values / count) %>%
  filter(!is.na(perc))

# group data table by module id class and count rows
fusion_mmseqs_comparison_class <- joined_dt %>%
  group_by(class) %>%
  summarize(mmseqs_values = sum(mmseqs_values), fusion_values = sum(fusion_values)) %>%
  mutate(perc = fusion_values / mmseqs_values) %>%
  arrange(desc(perc))

# plot increase in class predictions as barplot
comparison_plot <- fusion_mmseqs_comparison_class %>%
  filter(!is.na(class) & !is.infinite(perc)) %>%
  mutate(color = ifelse(perc > 1, "red", "black")) %>%
  mutate(class = reorder(class, perc)) %>%
  ggplot(aes(x = class, y = perc, fill = color)) + 
  geom_text(aes(label = signif(perc, 2)), hjust = -0.5) +
  geom_bar(stat = "identity") +
  xlab("KEGG Module Class") + ylab("Fold Change") + labs(title = "Fold Change in Class Link Prediction") +
  theme_bw() +
  coord_flip() +
  scale_fill_manual(values = c("black", "red1")) + 
  guides(fill = FALSE) 

pdf("compare/outputs/fusion_mmseqs_comparison.pdf", height = 5, width = 7)
comparison_plot
dev.off()
png("compare/outputs/fusion_mmseqs_comparison.png", height = 5, width = 7, units = "in", res = 300)
comparison_plot
dev.off()
tiff("compare/outputs/fusion_mmseqs_comparison.tiff", height = 5, width = 7, units = "in", res = 300)
comparison_plot
dev.off()

comparison_class_table <- comparison_plot$data %>% select(-color) %>% rename("Module Class"  = "class", "MMseqs2 - 40" = "mmseqs_values", "Fusion" = "fusion_values", "Fold Change" = "perc") %>% mutate(`Fold Change` = signif(`Fold Change`, 3)) 
comparison_class_table %>% write.csv("compare/outputs/fusion_mmseqs_comparison_class.csv")

#subclass analysis
fusion_mmseqs_comparison_subclass <- joined_dt %>%
  group_by(subclass) %>%
  summarize(mmseqs_values = sum(mmseqs_values), fusion_values = sum(fusion_values)) %>%
  mutate(perc = fusion_values / mmseqs_values) %>%
  arrange(desc(perc))

comparison_plot_subclass <- fusion_mmseqs_comparison_subclass %>%
  filter(perc >= 1 & !is.na(subclass) & !is.infinite(perc)) %>%
  mutate(color = ifelse(perc > 1, "red", "black")) %>%
  mutate(subclass = reorder(subclass, perc)) %>%
  ggplot(aes(x = subclass, y = perc, fill = color)) + 
  geom_text(aes(label = signif(perc, 3)), hjust = -0.5) +
  geom_bar(stat = "identity") +
  xlab("KEGG Module Subclass") + ylab("Fold Change") + labs(title = "Fold Change in Subclass Link Prediction") +
  theme_bw() +
  coord_flip() +
  scale_fill_manual(values = c("black", "red1")) + 
  guides(fill = FALSE) 

pdf("compare/outputs/fusion_mmseqs_comparison_subclass.pdf", height = 5, width = 7)
comparison_plot_subclass
dev.off()
png("compare/outputs/fusion_mmseqs_comparison_subclass.png", height = 5, width = 7, units = "in", res = 300)
comparison_plot_subclass
dev.off()
tiff("compare/outputs/fusion_mmseqs_comparison_subclass.tiff", height = 5, width = 7, units = "in", res = 300)
comparison_plot_subclass
dev.off()

comparison_subclass_table <- comparison_plot_subclass$data %>% 
  select(-color) %>% 
  rename("Module Subclass"  = "subclass", "MMseqs2 - 40" = "mmseqs_values", "Fusion" = "fusion_values", "Fold Change" = "perc")
comparison_subclass_table %>% write.csv("compare/outputs/fusion_mmseqs_comparison_subclass.csv")

# size of fusions in module
module_fusion_size <- module_df %>% 
  left_join(module_classes, by = "module_id") %>% 
  group_by(module_id) %>% 
  summarize(median_fusion_size = median(n), mean_fusion_size = mean(n)) 

# is it related to the number of proteins in the module?
module_parts <- data.table::fread("data/kegg/module_parts.csv")

# shannon diversity and count
shannon <- function(counts) {
  p <- counts / sum(counts)
  -sum(p * log(p))
}

module_ec_count <- module_parts %>%
  select(module_id, ec_number) %>%
  left_join(module_classes, by = "module_id") %>%
  unique() %>%
  mutate(ec_number3 = gsub("(\\d+)$", "-", ec_number)) %>%
  filter(ec_number != "") %>%
  group_by(module_id) %>%
  nest() %>%
  mutate(n = map_dbl(data, ~length(unique(.x$ec_number))),
         n3 = map_dbl(data, ~length(unique(.x$ec_number3))),
         shannon = map_dbl(data, ~shannon(table(.x$ec_number))),
         shannon3 = map_dbl(data, ~shannon(table(.x$ec_number3))))

# number of unique orthologs, or steps in in each module
module_pages <- readRDS("data/kegg/module_pages.RDS")

module_names <- lapply(module_pages, function(.x){.x$name}) %>% 
  stack() %>%
  rename(module_id = ind, module_name = values) 

module_orthologs <- lapply(module_pages, function(.x){
  unique_rxns <- length(.x$orthology)
  n_kos <- sapply(.x$orthology, function(.x){str_extract_all(.x, pattern = "K[0-9]{5}")}) %>% unname() %>% unlist() %>% unique() %>% length()
  res <- tibble(unique_rxns = unique_rxns, n_kos = n_kos)
  return(res)
  }) %>% bind_rows(.id = "module_id")

# separate based on organisms
organism_list <- readRDS("data/kegg/organism_list.RDS")

taxonomy <- (data.table::fread("data/taxonomy/gtdb_reference.csv")
  [, assembly_accession2 := gsub("\\..*", "", assembly_accession)]
  [assembly_accession2 %in% unique(fusion_preds$assembly_accession)]
)

org_to_assembly <- organism_list %>% 
  lapply(function(.x){head(.x$protein, 1)}) %>% 
  rbindlist() %>% 
  mutate(org := lapply(.$kegg_protein, function(.x){strsplit(.x, split = ":")[[1]][1]}), ncbi_accession2 := gsub("ncbi-proteinid:", "", ncbi_protein)) %>% 
  unnest(org) %>%
  left_join(fusion_data, by = "ncbi_accession2") %>%
  select(assembly_accession2, org)

module_org_count <- module_links %>% 
  select(module_id,org) %>% 
  left_join(org_to_assembly, by = "org") %>%
  left_join(select(taxonomy, c("assembly_accession2", "phylum", "species")), by = "assembly_accession2") %>%
  left_join(module_classes, by = "module_id") %>% 
  unique() %>% 
  group_by(module_id) %>% 
  summarize(n_org = length(unique(org)), n_phyla = length(unique(phylum))) %>%
  mutate(org_per_phyla = n_org / n_phyla) %>% 
  left_join(module_orthologs)

# Does fusion identify pathways with more 
module_subunits <- lapply(module_pages, function(.x){length(.x$orthology[grepl("subunit", .x$orthology)])}) %>%
  stack() %>%
  rename(module_id = ind, n_subunits = values)

# bind module metadata to fusion and mmseqs2 aggregate
joined_dt2 <- joined_dt %>%
  left_join(module_fusion_size, "module_id") %>%
  left_join(module_ec_count, "module_id") %>%
  left_join(module_org_count, "module_id") %>%
  left_join(module_subunits, "module_id") %>%
  mutate(perc = ifelse(is.infinite(perc), fusion_values, perc))  %>%
  mutate(n = ifelse(is.na(n), 0, n),
         n3 = ifelse(is.na(n3), 0, n3),
         shannon = ifelse(is.na(shannon), 0, shannon),
         shannon3 = ifelse(is.na(shannon3), 0, shannon3)) %>%
  mutate(perc = ifelse(is.infinite(perc), fusion_values, perc))


# hypergeometric test for modules
hypergeo_link_test <- joined_dt2 %>% select(-data) %>% filter(perc > 1) %>% select(module_id, fusion_values, mmseqs_values, perc) 
hypergeo_link_test$pval.adj <- apply(hypergeo_link_test, 1, function(.x){
  m <- filter(module_links_count, module_id == .x[[1]])$count
  # calculate hypergeometric test
  phyper(as.numeric(.x[[2]]), m = m, n = sum(module_links_count$count)-m, k = nrow(fusion_preds), lower.tail = FALSE, log.p = FALSE)
}) %>% p.adjust(method = "BH") 
fusion_sig_modules <- hypergeo_link_test %>% left_join(module_classes, by = "module_id") %>% filter(pval.adj <= 0.05) %>% pull(module_id)
joined_dt2 <- joined_dt2 %>% mutate(method = ifelse(module_id %in% fusion_sig_modules, "fusion", "mmseqs"))


comparison_table <- joined_dt %>%
  filter(module_id %in% fusion_sig_modules) %>%
  left_join(module_names) %>%
  select(module_id, module_name, mmseqs_values, fusion_values, count,  perc) %>%
  mutate(perc = signif(perc, 3)) %>%
  arrange(desc(perc)) %>%
  rename("Module ID" = module_id, "Name" = module_name, "Fold Change" = perc, "MMseqs2" = mmseqs_values, "Fusion" = fusion_values, "Total" = count) 
comparison_table
write.csv(comparison_table, "compare/outputs/fusion_mmseqs_comparison.csv", row.names = FALSE)


comparison_table2 <- joined_dt %>%
  filter(module_id %in% fusion_sig_modules) %>%
  group_by(class) %>%
  summarize(mmseqs_values = sum(mmseqs_values), fusion_values = sum(fusion_values), total = sum(count)) %>%
  mutate(perc = signif(fusion_values/mmseqs_values, 3)) %>%
  arrange(desc(perc)) %>%
  rename("Module Class" = class, "Fold Change" = perc, "MMseqs2" = mmseqs_values, "Fusion" = fusion_values, "Total" = total) 
comparison_table2
write.csv(comparison_table2, "compare/outputs/fusion_mmseqs_comparison2.csv", row.names = FALSE)

# wilcox test on each column
wilcox_variables <- c("count", "n", "n3", "shannon", "shannon3", "n_org", "median_dist", "median_length", "median_fusion_size", "n_phyla", "org_per_phyla", "unique_rxns", "n_kos")#, "n_subunits")
wilcox_results <- lapply(joined_dt2[,wilcox_variables], function(.x){
  wilcox.test(.x ~ method, data = joined_dt2)
}) %>%  lapply(function(.x){.x$p.value}) %>% 
  unlist() %>% 
  sort() %>%
  p.adjust(method = "BH") %>% 
  stack() %>% 
  as_tibble() %>%
  rename("adj.pval" = values, "variable" = ind) %>% 
  arrange(adj.pval)
wilcox_results

# plot values against each other 
# calculate median value for each variable by method
joined_dt2 %>% select(c("method", wilcox_variables)) %>% group_by(method) %>% summarize_all(median, na.rm = TRUE) %>% t() %>% as.data.frame() 

p2_data <- joined_dt2 %>%
  select(method, wilcox_variables) %>%
  pivot_longer(cols = -method, names_to = "variable", values_to = "value") %>%
  mutate(
    variable_label = case_when(
      variable == "count" ~ "Protein Pairs",
      variable == "n" ~ "Unique EC4",
      variable == "n3" ~ "Unique EC3",
      variable == "shannon" ~ "Shannon Diversity (EC4)",
      variable == "shannon3" ~ "Shannon Diversity (EC3)",
      variable == "n_org" ~ "Unique Organisms",
      variable == "median_dist" ~ "Median Chromosome Distance",
      variable == "median_fusion_size" ~ "Median Fusion Size",
      variable == "n_phyla" ~ "Unique Phyla",
      variable == "org_per_phyla" ~ "Organisms per Phyla",
      variable == "unique_rxns" ~ "Unique Reactions",
      variable == "median_length" ~ "Module Steps",
      variable == "n_kos" ~ "Unique KOs",
      variable == "n_ortho" ~ "Unique Orthologs",
      TRUE ~ as.character(variable)  # Default case, retains the original value if no condition is met
    )
  ) %>%
  mutate(variable = factor(variable, levels = wilcox_variables)) %>%
  mutate(method = ifelse(method == "fusion", "Fusion", "MMSeqs2")) 

pdf("compare/outputs/fusion_mmseqs_variable_comparison.pdf", height = 8, width = 8)
p2 <- p2_data %>%
  filter(!(variable %in% c("n_ortho", "org_per_phyla"))) %>%
  mutate(value = ifelse(variable == "count", log10(value), value), variable_label = ifelse(variable == "count", "Protein Pairs (log10)", variable_label)) %>%
  ggplot(aes(x = variable_label, y = value, fill = method)) +
  geom_boxplot() + 
  theme_bw() +
  xlab("") + 
  labs(title = "Fusion vs MMseqs2 - Variable Comparison") +
  scale_fill_manual(values = c("red1", "black")) +
  facet_wrap(~variable_label, scales = "free") +
  theme(strip.text = element_blank()) +
  labs(fill = "Method", title = "") +
  xlab("") + ylab("")
p2
dev.off()
png("compare/outputs/fusion_mmseqs_variable_comparison.png", height = 8, width = 8, units = "in", res = 300)
p2
dev.off()
tiff("compare/outputs/fusion_mmseqs_variable_comparison.tiff", height = 8, width = 8, units = "in", res = 300)
p2
dev.off()

# compare taxonomy of differences

# group fusion preds by assembly_accession2 and count rows
fusion_preds2 <- fusion_preds[truth == TRUE][, .(fusion_n = .N), by = assembly_accession2][,c("assembly_accession2", "fusion_n")]
mmseqs_preds2 <- mmseqs_preds[truth == TRUE][, .(mmseqs_n = .N), by = assembly_accession2][,c("assembly_accession2", "mmseqs_n")]

# join with fusion_preds2 and mmseqs_preds2
joined_taxa <- (left_join(mmseqs_preds2, fusion_preds2, by = "assembly_accession2")
  [taxonomy, on = "assembly_accession2"]
  [,c("assembly_accession2", "phylum", "fusion_n", "mmseqs_n")]
  [,perc_diff := fusion_n/mmseqs_n]
  [perc_diff > 1])

taxa_count <- taxonomy %>% group_by(phylum) %>% summarize(n = n())
# calculate hypergeometric test to determine significant of phylum counts with perc > 1
draws <- nrow(joined_taxa)
hypergeo_tests_phylum <- joined_taxa %>% group_by(phylum) %>% summarize(n = n()) %>% arrange(desc(n))
hypergeo_tests_phylum$pval_adj <- apply(hypergeo_tests_phylum, 1, function(.x){
  m <- filter(taxa_count, phylum == .x[[1]])$n
  # calculate hypergeometric test
  phyper(as.numeric(.x[[2]]), m = m, n = sum(taxa_count$n)-m, k = draws, lower.tail = FALSE, log.p = FALSE)
}) %>% p.adjust(method = "BH")
hypergeo_tests_phylum %>% filter(pval_adj <= 0.05)

# compare similarities
fusion_stats <- compare_df %>% 
  select(link, fusion_sim, truth, Var1, Var2) %>% 
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
summary(fusion_improv_stats$f_sat1); summary(fusion_worse_stats$f_sat2)
summary(fusion_improv_stats$m_sat1); summary(fusion_worse_stats$m_sat2)

mmseq_stats <- compare_df %>% 
  select(link, mmseqs_sim, truth, Var1, Var2) %>% 
  left_join(module_links, by = c("Var1" = "Var1", "Var2" = "Var2", "link" = "string")) %>%
  mutate(module_id = ifelse(is.na(module_id), "No Module", module_id)) %>%
  left_join(fm_sat, by = "Var1") %>%
  rename(f_sat1 = f_sat, m_sat1 = m_sat) %>%
  left_join(fm_sat, by = c("Var2" = "Var1")) %>%
  rename(f_sat2 = f_sat, m_sat2 = m_sat) %>%
  mutate(improv = ifelse(module_id %in% fusion_sig_modules, "fusion_improv", "base")) %>%
  mutate(improv = ifelse(improv == "base" & module_id != "No Module", "mmseqs_improv", improv))

mmseq_improv_stats <- mmseq_stats %>% filter(module_id %in% fusion_sig_modules) 
mmseq_worse_stats <- mmseq_stats %>% filter(!(module_id %in% fusion_sig_modules) & module_id != "No Module")
summary(mmseq_improv_stats$f_sat1); summary(mmseq_worse_stats$f_sat2)
summary(mmseq_improv_stats$m_sat1); summary(mmseq_worse_stats$m_sat2)

fusion_all <- fusion_stats$fusion_sim
fusion_modules <- fusion_stats %>% filter(module_id != "No Module") %>% pull(fusion_sim)
fusion_fusion_improv <- filter(fusion_stats, improv == "fusion_improv")$fusion_sim
fusion_mmseqs_improv <- filter(fusion_stats, improv == "mmseqs_improv")$fusion_sim

mmseq_all <- mmseq_stats$mmseqs_sim
mmseq_modules <- mmseq_stats %>% filter(module_id != "No Module") %>% pull(mmseqs_sim)
mmseq_fusion_improv <- filter(mmseq_stats, improv == "fusion_improv")$mmseqs_sim
mmseq_mmseqs_improv <- filter(mmseq_stats, improv == "mmseqs_improv")$mmseqs_sim

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
