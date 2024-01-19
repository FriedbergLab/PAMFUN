## ---------------------------
## Purpose of script: Calculate rand index of fusion and other functional groupings.
## Author: Henri Chung
## Date Created: 2022-06-26
## Date Modified: 2024-01-04
## ---------------------------

# load libraries
library(RColorBrewer)
library(data.table)
library(tidyverse)
library(parallel)
# clear working directory
rm(list = ls())

# read in required data
message(Sys.time(), " reading in data")
balanced_fusion_dtms <- readRDS("compare/data/balanced_fusion_dtms.RDS") 
balanced_mmseqs_dtms <- readRDS("compare/data/balanced_mmseq_dtms.RDS") 
# add names to mmseqs dtms
mmseq_folders <- list.files(path = "compare/data/mmseqs_results", pattern = "GCA", full.names = TRUE);
mmseq_files <- lapply(mmseq_folders, list.files, pattern = "mmseq$|mmseq\\.[0-9]$", full.names = TRUE) 

# list of balanced accessions with kegg annotations
balanced_organism_accessions <- paste("GCA_", readLines("data/fusion/balanced_organism_accessions"), sep = "")
# fusion data
balanced_data <- data.table::fread("compare/data/balanced_data.csv")
setkey(balanced_data, "assembly_accession2")

# create reference table for protein id to ec annotation
ec_data <- data.table::fread("data/fusion/uniprot-pe-exp_seguid_to_ec_mapping.tsv")
balanced_ec_data <- balanced_data[ec_data, ec_number := i.ec_number, on = "seguid"][,c("fusion_lvl_1", "ec_number", "ncbi_accession2", "assembly_accession2")]
setkeyv(balanced_ec_data, c("assembly_accession2", "fusion_lvl_1"))

module_proteins <- lapply(balanced_fusion_dtms, rownames) %>% unlist() %>% unname()

# calculate Rand Index for proteins in the same profile to see if they have the same EC number
rand_index <- function(.x){
  .x <- .x[!is.na(.x)]
  counts <- table(.x)
  unique_terms <- unique(.x)
  pairs <- data.table::CJ(unique_terms, unique_terms)
  pair_counts <- apply(pairs, 1, function(.y){
    counts[.y[1]] * counts[.y[2]]
  })
  same_pairs <- apply(pairs, 1, function(.y){
    return(.y[1] == .y[2])
  })
  same_counts <- sum(pair_counts[same_pairs])
  num <- (sum(same_counts) - length(.x))/2
  den <- (sum(pair_counts) - length(.x))/2
  #res <- paste0(sum(same_counts), "/", sum(pair_counts))
  res <- paste0(num, "/",den)
  return(res)
}

is_same_ec_mmseqs <- function(.x, threshold){
  if(is.null(.x) | .x == "NULL"){return(NULL)}
  message(Sys.time(), " ", which(names(balanced_mmseqs_dtms) == .x))
  temp_file <- mmseq_files[grepl(.x, mmseq_files)][[1]]
  temp_dt <- (data.table::fread(temp_file)[V3 > threshold]
    [,query2 := gsub("\\.[0-9]", "", V1)]
    [query2 %in% module_proteins]
    [,ncbi_accession2 := gsub("\\.[0-9]", "", V2)]
    [fusion_ec_data, ec_number := i.ec_number, on = "ncbi_accession2"]
    [!is.na(ec_number)]
    [ , .SD[.N > 1], by = V1]
  )

  result <- split(temp_dt, temp_dt$query2) %>% lapply(function(.x){
    res <- rand_index(.x$ec_number)
    return(res)
  })

  if(length(result) == 0){return(NULL)}else{return(stack(result))}
}

summarize_fractions <- function(fractions) {
  total_numerator <- 0
  total_denominator <- 0
  for(frac in fractions) {
    parts <- as.numeric(unlist(strsplit(frac, "/")))
    total_numerator <- total_numerator + parts[1]
    total_denominator <- total_denominator + parts[2]
  }
  combined_decimal <- total_numerator / total_denominator
  percentage <- 100 * combined_decimal
  return(percentage)
}

# calculate the Rand Index for proteins in the fusion
fusion_ec_data <- balanced_ec_data[!is.na(ec_number)][assembly_accession2 %in% names(balanced_mmseqs_dtms)][ , .SD[.N > 1], by = fusion_lvl_1]
fusion_ec_fractions <- split(fusion_ec_data, fusion_ec_data$fusion_lvl_1) %>% lapply(function(.x){
  res <- rand_index(.x$ec_number)
  res <- as.list(rep(res, sum(.x$ncbi_accession2 %in% module_proteins)))
  return(res)
})
summarize_fractions(unlist(fusion_ec_fractions)) 

# calculate the Rand Index for proteins in the mmseqs 
mmseqs_assemblies <- unique(fusion_ec_data[ncbi_accession2 %in% module_proteins]$assembly_accession2)
mmseqs_ec_fractions <- lapply(mmseqs_assemblies, is_same_ec_mmseqs, threshold = 40)
mmseqs_ec_fractions2 <- bind_rows(mmseqs_ec_fractions, .id = "assembly_accession2")
summarize_fractions(mmseqs_ec_fractions2$values)

result <- tibble(method = c("fusion", "mmseqs"), rand_index = c(summarize_fractions(unlist(fusion_ec_fractions)), summarize_fractions(mmseqs_ec_fractions2$values)), definition = "ec")
saveRDS(result, "compare/outputs/rand_index_ec.RDS")

#################################################################


# create reference teable for protein id to module annotation
protein_links0 <- data.table::fread("data/kegg/protein_links_df.csv")
protein_links_a <- protein_links0[, .(module_id, Var1)]
protein_links_b <- protein_links0[, .(module_id, Var2)][, Var1 := Var2][, Var2 := NULL]
protein_links <- (unique(rbind(protein_links_a, protein_links_b))
  [, ncbi_accession2 := Var1]
  [,Var1 := NULL]
  [balanced_ec_data, assembly_accession2 := i.assembly_accession2, on = "ncbi_accession2"]
  [balanced_ec_data, fusion_lvl_1 := i.fusion_lvl_1, on = "ncbi_accession2"]
  [!is.na(fusion_lvl_1)]
  [!is.na(module_id)]
)
setkeyv(protein_links, c("assembly_accession2", "fusion_lvl_1"))

balanced_module_data <- balanced_data[protein_links, module_id := i.module_id, on = "ncbi_accession2"][!is.na(module_id)]
module_proteins <- lapply(balanced_fusion_dtms, rownames) %>% unlist() %>% unname()
# calculate Rand Index for proteins in the same profile to see if they have the same EC number
rand_index <- function(.x){
  .x <- .x[!is.na(.x)]
  counts <- table(.x)
  unique_terms <- unique(.x)
  pairs <- data.table::CJ(unique_terms, unique_terms)
  pair_counts <- apply(pairs, 1, function(.y){
    counts[.y[1]] * counts[.y[2]]
  })
  same_pairs <- apply(pairs, 1, function(.y){
    return(.y[1] == .y[2])
  })
  same_counts <- sum(pair_counts[same_pairs])
  num <- (sum(same_counts) - length(.x))/2
  den <- (sum(pair_counts) - length(.x))/2
  #res <- paste0(sum(same_counts), "/", sum(pair_counts))
  res <- paste0(num, "/",den)
  return(res)
}

summarize_fractions <- function(fractions) {
  total_numerator <- 0
  total_denominator <- 0
  for(frac in fractions) {
    parts <- as.numeric(unlist(strsplit(frac, "/")))
    total_numerator <- total_numerator + parts[1]
    total_denominator <- total_denominator + parts[2]
  }
  combined_decimal <- total_numerator / total_denominator
  percentage <- 100 * combined_decimal
  return(percentage)
}

# calculate the Rand Index for proteins in the fusion
fusion_module_data <- balanced_module_data[!is.na(module_id)][assembly_accession2 %in% names(balanced_mmseqs_dtms)][ , .SD[.N > 1], by = fusion_lvl_1]
fusion_module_fractions <- split(fusion_module_data, fusion_module_data$fusion_lvl_1) %>% lapply(function(.x){
  res <- rand_index(.x$module_id)
  res <- as.list(rep(res, sum(.x$ncbi_accession2 %in% module_proteins)))
  return(res)
})
summarize_fractions(unlist(fusion_module_fractions)) 


is_same_module_mmseqs <- function(.x, threshold){
  if(is.null(.x) | .x == "NULL"){return(NULL)}
  message(Sys.time(), " ", which(names(balanced_mmseqs_dtms) == .x))
  temp_dtm <- balanced_mmseqs_dtms[[.x]]
  temp_file <- mmseq_files[grepl(.x, mmseq_files)][[1]]
  temp_dt <- (data.table::fread(temp_file)[V3 > threshold]
    [,query2 := gsub("\\.[0-9]", "", V1)]
    [query2 %in% module_proteins]
    [,ncbi_accession2 := gsub("\\.[0-9]", "", V2)]
    [balanced_module_data, module_id := i.module_id, on = "ncbi_accession2"]
    [!is.na(module_id)]
    [ , .SD[.N > 1], by = V1]
  )

  result <- split(temp_dt, temp_dt$query2) %>% lapply(function(.x){
    res <- rand_index(.x$module_id)
    return(res)
  })

  if(length(result) == 0){return(NULL)}else{return(stack(result))}
}


# calculate the Rand Index for proteins in the mmseqs 
mmseqs_assemblies <- unique(fusion_module_data$assembly_accession2)
mmseqs_module_fractions <- lapply(mmseqs_assemblies, is_same_module_mmseqs, threshold = 40)
mmseqs_module_fractions2 <- bind_rows(mmseqs_module_fractions, .id = "assembly_accession2")
summarize_fractions(mmseqs_module_fractions2$values)

result <- tibble(method = c("fusion", "mmseqs"), rand_index = c(summarize_fractions(unlist(fusion_module_fractions)), summarize_fractions(mmseqs_module_fractions2$values)), definition = "module")
saveRDS(result, "compare/outputs/rand_index_module.RDS")

#################################################################



#  calculate Rand Index for proteins in the same profile to see if they have the same KO
orthology_to_protein <- (data.table::fread("data/kegg/orthology_to_protein.csv")
  [balanced_ec_data, fusion_lvl_1 := i.fusion_lvl_1, on = "ncbi_accession2"]
  [balanced_ec_data, assembly_accession2 := i.assembly_accession2, on = "ncbi_accession2"])
setkeyv(orthology_to_protein, c("fusion_lvl_1", "assembly_accession2"))
setkey(orthology_to_protein, "ncbi_accession2")  


balanced_ortholog_data <- balanced_data[orthology_to_protein, ortholog := i.ortholog, on = "ncbi_accession2"][!is.na(ortholog)]
module_proteins <- lapply(balanced_fusion_dtms, rownames) %>% unlist() %>% unname()
# calculate Rand Index for proteins in the same profile to see if they have the same EC number
rand_index <- function(.x){
  .x <- .x[!is.na(.x) & .x != ""]
  counts <- table(.x)
  unique_terms <- unique(.x)
  pairs <- data.table::CJ(unique_terms, unique_terms)
  pair_counts <- apply(pairs, 1, function(.y){
    counts[.y[1]] * counts[.y[2]]
  })
  same_pairs <- apply(pairs, 1, function(.y){
    return(.y[1] == .y[2])
  })
  same_counts <- sum(pair_counts[same_pairs])
  num <- (sum(same_counts) - length(.x))/2
  den <- (sum(pair_counts) - length(.x))/2
  res <- paste0(num, "/",den)
  return(res)
}

summarize_fractions <- function(fractions) {
  total_numerator <- 0
  total_denominator <- 0
  for(frac in fractions) {
    parts <- as.numeric(unlist(strsplit(frac, "/")))
    total_numerator <- total_numerator + parts[1]
    total_denominator <- total_denominator + parts[2]
  }
  combined_decimal <- total_numerator / total_denominator
  percentage <- 100 * combined_decimal
  return(percentage)
}

# calculate the Rand Index for proteins in the fusion
fusion_ortholog_data <- balanced_ortholog_data[!is.na(ortholog)][ncbi_accession2 %in% module_proteins][assembly_accession2 %in% names(balanced_mmseqs_dtms)][ , .SD[.N > 1], by = fusion_lvl_1]

fusion_ortholog_fractions <- split(fusion_ortholog_data, fusion_ortholog_data$fusion_lvl_1) %>% lapply(function(.x){
  res <- rand_index(.x$ortholog)
  res <- as.list(rep(res, sum(.x$ncbi_accession2 %in% module_proteins)))
  return(res)
})
summarize_fractions(unlist(fusion_ortholog_fractions)) 


is_same_ortholog_mmseqs <- function(.x, threshold){
  if(is.null(.x) | .x == "NULL"){return(NULL)}
  message(Sys.time(), " ", which(names(balanced_mmseqs_dtms) == .x))
  temp_dtm <- balanced_mmseqs_dtms[[.x]]
  temp_file <- mmseq_files[grepl(.x, mmseq_files)][[1]]
  temp_dt <- (data.table::fread(temp_file)[V3 > threshold]
    [,query2 := gsub("\\.[0-9]", "", V1)]
    [query2 %in% module_proteins]
    [,ncbi_accession2 := gsub("\\.[0-9]", "", V2)]
    [balanced_ortholog_data, ortholog := i.ortholog, on = "ncbi_accession2"]
    [!is.na(ortholog)]
    [ , .SD[.N > 1], by = V1]
  )

  result <- split(temp_dt, temp_dt$query2) %>% lapply(function(.x){
    res <- rand_index(.x$ortholog)
    return(res)
  })

  if(length(result) == 0){return(NULL)}else{return(stack(result))}
}

# calculate the Rand Index for proteins in the mmseqs 
mmseqs_assemblies <- unique(fusion_ortholog_data$assembly_accession2)
mmseqs_ortholog_fractions <- lapply(mmseqs_assemblies, is_same_ortholog_mmseqs, threshold = 40)
mmseqs_ortholog_fractions2 <- bind_rows(mmseqs_ortholog_fractions, .id = "assembly_accession2")
result <- tibble(method = c("fusion", "mmseqs"), rand_index = c(summarize_fractions(unlist(fusion_ortholog_fractions)), summarize_fractions(mmseqs_ortholog_fractions2$values)), definition = "ortholog")
saveRDS(result, "compare/outputs/rand_index_ortholog.RDS")


#################################################################



balanced_ortholog_data <- balanced_data[orthology_to_protein, ortholog := i.ortholog, on = "ncbi_accession2"][!is.na(ortholog)]
module_proteins <- lapply(balanced_fusion_dtms, rownames) %>% unlist() %>% unname()
module_steps <- data.table::fread("data/kegg/non_orthologous_replacements.csv") %>% mutate(combos = paste(a, b, sep = "_")) %>% filter(a != b) %>% unique()
# calculate Rand Index for proteins in the same profile to see if they have the same non-orthologous KOs
rand_index_nors <- function(.x){
  .x <- .x[!is.na(.x) & .x != ""]
  counts <- table(.x)
  unique_terms <- unique(.x)
  pairs <- data.table::CJ(unique_terms, unique_terms)
  pair_counts <- apply(pairs, 1, function(.y){
    counts[.y[1]] * counts[.y[2]]
  })
  same_pairs <- paste(pairs[[1]], pairs[[2]], sep = "_")
  ind <- same_pairs %in% module_steps$combos
  same_counts <- sum(pair_counts[ind])
  num <- (sum(same_counts) - length(.x))/2
  den <- (sum(pair_counts) - length(.x))/2
  res <- paste0(num, "/",den)
  return(res)
}


summarize_fractions <- function(fractions) {
  total_numerator <- 0
  total_denominator <- 0
  for(frac in fractions) {
    parts <- as.numeric(unlist(strsplit(frac, "/")))
    total_numerator <- total_numerator + parts[1]
    total_denominator <- total_denominator + parts[2]
  }
  combined_decimal <- total_numerator / total_denominator
  percentage <- 100 * combined_decimal
  return(percentage)
}

# calculate the Rand Index for proteins in the fusion
fusion_ortholog_data <- balanced_ortholog_data[!is.na(ortholog)][ncbi_accession2 %in% module_proteins][assembly_accession2 %in% names(balanced_mmseqs_dtms)][ , .SD[.N > 1], by = fusion_lvl_1]
fusion_ortholog_fractions <- split(fusion_ortholog_data, fusion_ortholog_data$fusion_lvl_1) %>% lapply(function(.x){
  res <- rand_index_nors(.x$ortholog)
  return(res)
})
summarize_fractions(unlist(fusion_ortholog_fractions)) 


is_same_ortholog_mmseqs <- function(.x, threshold){
  if(is.null(.x) | .x == "NULL"){return(NULL)}
  message(Sys.time(), " ", which(names(balanced_mmseqs_dtms) == .x))
  temp_dtm <- balanced_mmseqs_dtms[[.x]]
  temp_file <- mmseq_files[grepl(.x, mmseq_files)][[1]]
  temp_dt <- (data.table::fread(temp_file)[V3 > threshold]
    [,query2 := gsub("\\.[0-9]", "", V1)]
    [query2 %in% module_proteins]
    [,ncbi_accession2 := gsub("\\.[0-9]", "", V2)]
    [balanced_ortholog_data, ortholog := i.ortholog, on = "ncbi_accession2"]
    [!is.na(ortholog)]
    [ , .SD[.N > 1], by = V1]
  )

  result <- split(temp_dt, temp_dt$query2) %>% lapply(function(.x){
    res <- rand_index_nors(.x$ortholog)
    res <- as.list(rep(res, sum(.x$ncbi_accession2 %in% module_proteins)))
    return(res)
  })

  if(length(result) == 0){return(NULL)}else{return(stack(result))}
}

# calculate the Rand Index for proteins in the mmseqs 
mmseqs_assemblies <- unique(fusion_ortholog_data$assembly_accession2)
mmseqs_ortholog_fractions <- lapply(mmseqs_assemblies, is_same_ortholog_mmseqs, threshold = 40)
mmseqs_ortholog_fractions2 <- bind_rows(mmseqs_ortholog_fractions, .id = "assembly_accession2")

result <- tibble(method = c("fusion", "mmseqs"), rand_index = c(summarize_fractions(unlist(fusion_ortholog_fractions)), summarize_fractions(mmseqs_ortholog_fractions2$values)), definition = "nors")
saveRDS(result, "compare/outputs/rand_index_nors.RDS")

#################################################################


# calculate Rand Index for proteins in the same profile to see if they have the same EC number
rand_index <- function(.x){
  .x <- .x[!is.na(.x)]
  counts <- table(.x)
  unique_terms <- unique(.x)
  pairs <- data.table::CJ(unique_terms, unique_terms)
  pair_counts <- apply(pairs, 1, function(.y){
    counts[.y[1]] * counts[.y[2]]
  })
  same_pairs <- apply(pairs, 1, function(.y){
    return(.y[1] == .y[2])
  })
  same_counts <- sum(pair_counts[same_pairs])
  num <- (sum(same_counts) - length(.x))/2
  den <- (sum(pair_counts) - length(.x))/2
  #res <- paste0(sum(same_counts), "/", sum(pair_counts))
  res <- paste0(num, "/",den)
  return(res)
}

is_same_ec_mmseqs <- function(.x, threshold){
  if(is.null(.x) | .x == "NULL"){return(NULL)}
  message(Sys.time(), " ", which(names(balanced_mmseqs_dtms) == .x))
  temp_dtm <- balanced_mmseqs_dtms[[.x]]
  temp_file <- mmseq_files[grepl(.x, mmseq_files)][[1]]
  temp_dt <- (data.table::fread(temp_file)[V3 > threshold]
    [,query2 := gsub("\\.[0-9]", "", V1)]
    [query2 %in% module_proteins]
    [,ncbi_accession2 := gsub("\\.[0-9]", "", V2)]
    [balanced_ec_data, ec_number := i.ec_number, on = "ncbi_accession2"]
    [!is.na(ec_number)]
    [ , .SD[.N > 1], by = V1]
  )

  result <- split(temp_dt, temp_dt$query2) %>% lapply(function(.x){
    res <- rand_index(.x$ec_number)
    res <- as.list(rep(res, sum(.x$ncbi_accession2 %in% module_proteins)))
    return(res)
  })

  if(length(result) == 0){return(NULL)}else{return(stack(result))}
}

summarize_fractions <- function(fractions) {
  total_numerator <- 0
  total_denominator <- 0
  for(frac in fractions) {
    parts <- as.numeric(unlist(strsplit(frac, "/")))
    total_numerator <- total_numerator + parts[1]
    total_denominator <- total_denominator + parts[2]
  }
  combined_decimal <- total_numerator / total_denominator
  percentage <- 100 * combined_decimal
  return(percentage)
}

# calculate the Rand Index for proteins in the fusion
fusion_ec_data <- balanced_ec_data[!is.na(ec_number)][ncbi_accession2 %in% module_proteins][assembly_accession2 %in% names(balanced_mmseqs_dtms)][ , .SD[.N > 1], by = fusion_lvl_1]
fusion_ec_fractions <- split(fusion_ec_data, fusion_ec_data$fusion_lvl_1) %>% lapply(function(.x){
  res <- rand_index(.x$ec_number)
  res <- as.list(rep(res, sum(.x$ncbi_accession2 %in% module_proteins)))
  return(res)
})
summarize_fractions(unlist(fusion_ec_fractions)) 

# calculate the Rand Index for proteins in the mmseqs 
mmseqs_assemblies <- unique(fusion_ec_data$assembly_accession2)
mmseqs_ec_fractions <- lapply(mmseqs_assemblies, is_same_ec_mmseqs, threshold = 40)
mmseqs_ec_fractions2 <- bind_rows(mmseqs_ec_fractions, .id = "assembly_accession2")
summarize_fractions(mmseqs_ec_fractions2$values)

result <- tibble(method = c("fusion", "mmseqs"), rand_index = c(summarize_fractions(unlist(fusion_ec_fractions)), summarize_fractions(mmseqs_ec_fractions2$values)), definition = "ec3")
saveRDS(result, "compare/outputs/rand_index_ec3.RDS")

#################################################################

# read in data and compare
compare_ec <- readRDS("compare/outputs/rand_index_ec.RDS")
compare_ec3 <- readRDS("compare/outputs/rand_index_ec3.RDS") 
compare_modules <- readRDS("compare/outputs/rand_index_module.RDS") 
compare_kos <- readRDS("compare/outputs/rand_index_ortholog.RDS") 
compare_nors <- readRDS("compare/outputs/rand_index_nors.RDS")

compare_list <- list(compare_ec, compare_ec3, compare_modules, compare_kos, compare_nors)

p1_data <- bind_rows(compare_list) %>%
  rename(variable = "definition") %>%
  mutate(rand_index = signif(rand_index, 3)) %>%
  mutate(method = factor(method, levels = c("fusion", "mmseqs"))) %>%
  mutate(variable2  = case_when(
    variable == "ec" ~ "EC4",
    variable == "ec3" ~ "EC3",
    variable == "module" ~ "Module ID",
    variable == "ortholog" ~ "KO",
    variable == "nors" ~ "NOR"
  )) 

table1 <- p1_data %>% pivot_wider(names_from = method, values_from = rand_index) %>%
  select(variable2, fusion, mmseqs) %>%
  rename("Fusion" = fusion, "MMseqs" = mmseqs, "Definition" = variable2)
table1
write_csv(table1, "compare/outputs/compare_rand_index1.csv")
table2 <- p1_data %>% select(method, rand_index, variable2) %>% 
  pivot_wider(names_from = variable2, values_from = rand_index) %>%
  mutate(method = case_when(
    method == "fusion" ~ "Fusion",
    method == "mmseqs" ~ "MMseqs"
  ))
write_csv(table2, "compare/outputs/compare_rand_index2.csv")
table2