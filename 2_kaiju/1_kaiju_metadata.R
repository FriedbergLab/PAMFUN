## ---------------------------
## Purpose of script: Parse marine metagenome metadata 
## Author: Henri Chung
## Date Created: 2022-09-21
## Date Modified: 2024-01-04
## ---------------------------
library(tidyverse)
library(rvest)
rm(list = ls())


# Download study metadata
url <- "https://www.nature.com/articles/sdata2018176/tables/4"
webpage <- read_html(url)
tables <- html_nodes(webpage, "table")
dataframes <- lapply(tables, html_table)

PRJNA385854_metadata <- dataframes[[1]] %>%
    janitor::clean_names()
write.csv(PRJNA385854_metadata, "kaiju/outputs/PRJNA385854_metadata.csv", row.names = FALSE)

# export dataframe connecting taxid to modules
pathways_list <- readRDS("data/kegg/full_pathways_list.RDS")
complete_modules <- stack(unlist(lapply(pathways_list, function(.x){.x$complete}), recursive = FALSE)) %>%
   rename("module_id" = values, "org_code" = ind) %>%
   group_by(org_code) %>% 
   nest()

# export dataframe connecting modules to ec numbers
kegg_assemblies <- data.table::fread("data/kegg/kegg_assemblies.csv") %>%
    left_join(complete_modules, by = "org_code")
saveRDS(kegg_assemblies, "kaiju/outputs/org_module.RDS")

# export dataframe connecting fusions to ec numbers
system("grep '>' mifaser/data/fusion_db_small.fa > mifaser/data/fusion_db_small_headers.txt")
fusion_headers <- readLines("mifaser/data/fusion_db_small_headers.txt") %>% substr(., start = 2, stop = 9)
fusion_data <- data.table::fread("data/fusion/fusion_data.tsv")[, assembly_accession2 := gsub("\\.[0-9]", "", assembly_accession)][, ncbi_accession2 := gsub("\\.[0-9]", "", ncbi_accession)]
setkey(fusion_data, ncbi_accession2)
ec_data <- data.table::fread("data/fusion/uniprot-pe1-exp-ec.extended_by_hfsp.mapping")[hfsp >= 20][,ec_number := gsub("n", "", ec_number)]

fusion_ec_data <- fusion_data[ncbi_accession2 %in% fusion_headers][ec_data, ec_number := i.ec_number, on = "seguid"][,c("fusion_lvl_1", "ec_number", "ncbi_accession2", "assembly_accession2")]
data.table::fwrite(fusion_ec_data, "kaiju/data/fusion_ec_data.tsv", sep = "\t", row.names = FALSE)


# compare overlaps between metagenome preds and modules
metagenome_preds <- data.table::fread("kaiju/outputs/metagenome_preds.csv")
modules_ec <- data.table::fread("data/kegg/module_parts.csv")

modules_ec_list <- split(modules_ec$ec_number, modules_ec$module) 
overlaps <- apply(metagenome_preds,1, function(.x){
    temp <- c(.x[["EC1"]], .x[["EC2"]])
    matches <- lapply(modules_ec_list, function(.y){
        sum(temp %in% .y)
    }) %>% unlist()
    max_match <- matches[which(matches == max(matches))]
    res <- tibble("module_id"= names(max_match), "overlap" = max_match, pair = paste(.x[["EC1"]], .x[["EC2"]], sep = "_"))
    return(res)
}) %>% bind_rows()
arrange(overlaps, desc(overlap))


quit()

