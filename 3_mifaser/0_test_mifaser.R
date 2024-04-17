## ---------------------------
## Purpose of script: Set up mifaser testing
## Author: Henri Chung
## Date Created: 2022-09-21
## Date Modified: 2024-01-04
## ---------------------------

library(data.table)
library(tidyverse)
# clear working directory
rm(list = ls())

balanced_data <- data.table::fread("compare/data/balanced_data.csv")
balanced_organism_accessions <- unique(balanced_data$assembly_accession2)
# load input data for feature tables
ft_files <- list.files("data/feature_tables", full.names = TRUE)
ft_ind <- (str_extract(ft_files, "GCA_[0-9]+") %in% balanced_organism_accessions) 
assembly_stats <- data.table::fread("data/assembly_stats/assembly_lengths.tsv")
colnames(assembly_stats) <- c("assembly_accession2", "length")
kegg_assemblies <- data.table::fread("data/kegg/kegg_assemblies.csv")

# read in feature tables (chromosome positions)
ft_ind <- (str_extract(ft_files, "GCA_[0-9]+") %in% balanced_organism_accessions) 
ft_list <- list()
for(i in 1:length(ft_files[ft_ind])){ft_list[[i]] <- data.table::fread(ft_files[ft_ind][i], fill = TRUE)}
to_fix <- ft_files[ft_ind][which(sapply(ft_list, nrow) == 0)]

# format feature table as chr position of proteins
feature_table <- rbindlist(ft_list)
feature_table <- (feature_table
  [,product_accession2 := gsub("\\.[0-9]", "", product_accession)]
  [product_accession2 != ""]
  [,assembly_accession2 := str_extract(assembly, "GCA_[0-9]+")]
  [assembly_stats, length := i.length, on = "assembly_accession2"]
)

# load list of balanced genomes
gn_files <- list.files("data/genomes", full.names = TRUE)
gn_file_assemblies <- str_extract(gn_files, "GCA_[0-9]+")
ind <- gn_file_assemblies %in% balanced_organism_accessions
balanced_gn_files <- gn_files[ind]
balanced_gn_assemblies <- gn_file_assemblies[ind]

# split balanced genomes into training and test set
set.seed(123) # for reproducibility
split_index <- sample(c(TRUE, FALSE), size = length(balanced_gn_assemblies), replace = TRUE, prob = c(0.7, 0.3))
train_gn <- balanced_gn_assemblies[split_index]
test_gn <- balanced_gn_assemblies[!split_index]
train_gn_files <- balanced_gn_files[split_index]
test_gn_files <- balanced_gn_files[!split_index]

if(!file.exists("mifaser/data/train_gn.txt")){
  writeLines(train_gn, "mifaser/data/train_gn.txt")
}

# extract proteins from balanced protein fastas into fusion train db
count <- system("wc -l mifaser/data/fusion_train_db.fa", intern = TRUE) %>%
  strsplit(" ")
if(count[[1]][1] == "0") {
  system("touch mifaser/data/fusion_train_db.fa")
  system("while read line; do find data/protein_fa -type f -name \"*$line*\" -exec zcat {} + >> mifaser/data/fusion_train_db.fa ; done < mifaser/data/train_gn.txt")
}

# format training proteins to have fusion annotation as pseudo-EC number
insert_decimals <- function(x) {
  x_str <- sprintf("%04d", as.integer(x))
  split_at <- sample(1:(nchar(x_str)-1), 3)
  split_at <- sort(split_at, decreasing = TRUE)
  for (i in split_at) {
    x_str <- paste0(substr(x_str, 1, i), ".", substr(x_str, i+1, nchar(x_str)))
  }
  return(x_str)
}

# write fusion training fasta file
if(!file.exists("mifaser/data/fusion_train.fa")){
  # read in fasta file of proteins
  train_proteins <- seqinr::read.fasta("mifaser/data/fusion_train_db.fa")
  train_proteins_names <- gsub("\\.[0-9]", "", names(train_proteins))

  # create table of proteins with mifaser names
  mifaser_names <- (balanced_data
    [ncbi_accession2 %in% train_proteins_names] # filter to proteins in fusion_train_db
    [,fake_ec := sapply(fusion_lvl_1, insert_decimals)] # replace fusion_lvl_1 with fake ec number
    [,protein_name := gsub("\\|", ",", protein_name)] # replace | in protein_name 
    [,mifaser_name := paste(ncbi_accession2, "|", assembly_accession2, protein_name, "|", fake_ec, "|", fusion_lvl_1)] # create mifaser names
  )

  # match order of train_proteins to mifaser_names
  mifaser_names <- mifaser_names[match(train_proteins_names, mifaser_names$ncbi_accession2),]
  # assign names
  names(train_proteins) <- mifaser_names$mifaser_name
  # remove NA
  train_proteins <- train_proteins[!is.na(names(train_proteins))]
  # change sequence to uppercase
  train_proteins2 <- lapply(train_proteins, toupper)
  #write to file
  seqinr::write.fasta(train_proteins2, names = names(train_proteins2), "mifaser/data/fusion_train.fa")
}

# generate test sequences
get_genome_ind <- function(.x){
    size <- unlist(lapply(.x, length))
    ind = which(size == max(size))

    return(.x[ind])
}

generate_sequences <- function(vector, N, seed = 123) {
  set.seed(seed)
  sequences <- tibble(index = 1:N) %>%
    mutate(seq_length = map_int(index, ~ sample(50:250, 1)),
           start_index = map_int(seq_length, ~ sample(1:(length(vector) - .x + 1), 1)),
           end_index = start_index + seq_length - 1,
           sequence = map2(start_index, end_index, ~ vector[.x:.y])) %>%
    select(-index)
  
  return(sequences)
}


# create simulated short reads
if(!file.exists("mifaser/data/test_reads_metadata.RDS")){
  fragment_number <- 10000
  sim_list <- list()
  for(i in 1:length(test_gn_files)){
    temp_fasta <- seqinr::read.fasta(test_gn_files[i]) %>% get_genome_ind()
    if(length(temp_fasta) == 0){next}
    sim_reads <- generate_sequences(temp_fasta[[1]], N = fragment_number, seed = 123) %>%
      mutate(name = paste0(names(temp_fasta), " | start = ", start_index, " | end = ", end_index))
    message("simulating reads ", i)
    sim_list[[i]] <- sim_reads
  }
  names(sim_list) <- gsub("\\.[0-9]", "", str_extract( test_gn_files, "GCA_.*\\.[0-9]"))
  sim_list2 <- sim_list[!is.na(sim_list)]
  all_sim_reads <- bind_rows(sim_list2, .id = "assembly_accession2")
  saveRDS(all_sim_reads, "mifaser/data/test_reads_metadata.RDS")
}else{
  all_sim_reads <- readRDS("mifaser/data/test_reads_metadata.RDS") 
}
# assigned short reads to fusions by calculating max percent overlap with known proteins
overlap_percent <- function(df, start, end, N) {
  a_start <- df$start
  a_end <- df$end
  names <- df$product_accession2
  overlap <- pmin(a_end, end-1) - pmax(a_start, start-1)
  overlap[overlap < 0] <- 0
  percent_overlap <- overlap / (a_end - a_start)
  top_n <- order(percent_overlap, decreasing = TRUE)[1:N]
  filtered_top_n <- top_n[overlap[top_n] >= 10]
  
  if (length(top_n) == 0) {
    cat("No significant overlaps could be found.\n")
    return(NULL)
  }

  result <- list(proteins = names[filtered_top_n], percent_overlap = percent_overlap[filtered_top_n])
  return(result)
}

# assign fusion annotations to simulated reads
all_sim_reads$fusions <- NA
n_cores = 4
N <- nrow(all_sim_reads)
all_sim_reads$fusion <- mclapply(1:N, function(i){
  temp_ft <- feature_table[assembly_accession2 == all_sim_reads$assembly_accession2[i]]
  protein <- overlap_percent(start = all_sim_reads$start_index[i], end = all_sim_reads$end_index[i], df = temp_ft, N = 10)
  if(length(protein$proteins) == 0){
    res <- "0"
  }else{
    res <- as.character(balanced_data[ncbi_accession2 %in% protein$proteins]$fusion_lvl_1)
    if(length(res) == 0){res <- "0"}
  }
  message(signif(i/N, 2))
  return(res)
}, m.cores = n_cores)

#
values <- as.list(balanced_data$fusion_lvl_1)
names(values) <- balanced_data$ncbi_accession2
protein_to_fusion <- fastmap::fastmap()
protein_to_fusion$mset(.list = values)
#
feature_list <- split(feature_table, feature_table$assembly_accession2)
assembly_to_feature <- fastmap::fastmap()
assembly_to_feature$mset(.list = feature_list)
all_sim_reads$fusions <- NA
n_cores = 4
N <- nrow(all_sim_reads)
#
all_sim_reads$fusion <- parallel::mclapply(1:N, function(i){
  temp_ft <- assembly_to_feature$get(all_sim_reads$assembly_accession2[i])
  protein <- overlap_percent(start = all_sim_reads$start_index[i], end = all_sim_reads$end_index[i], df = temp_ft, N = 10)
  if(length(protein$proteins) == 0){
    res <- "0"
  }else{
    res <- as.character(protein_to_fusion$mget(protein$proteins))
    if(length(res) == 0){res <- "0"}
  }
  message(signif(i/N, 2))
  return(res)
}, mc.cores = n_cores)
#
missing <- lapply(all_sim_reads$fusions, is.null) %>% unlist()
ind <- which(missing == TRUE)
while(length(ind) > 0){
  all_sim_reads$fusions[ind] <- lapply(1:length(ind), function(i){
    temp_ft <- assembly_to_feature$get(all_sim_reads$assembly_accession2[ind[i]])
    protein <- overlap_percent(start = all_sim_reads$start_index[ind[i]], end = all_sim_reads$end_index[ind[i]], df = temp_ft, N = 10)
    if(length(protein$proteins) == 0){
      res <- "0"
    }else{
      res <- as.character(protein_to_fusion$mget(protein$proteins))
      if(length(res) == 0){res <- "0"}
    }
    message(signif(i/length(ind), 2))
    return(res)
  })
  missing <- lapply(all_sim_reads$fusion, is.null) %>% unlist()
  ind <- which(missing == TRUE)
}

all_sim_reads <- all_sim_reads %>% unnest(fusions)
all_sim_reads$sequence <- lapply(all_sim_reads$sequence, toupper)
all_sim_reads$name <- trimws(all_sim_reads$name)

# write the test short reads to file
seqinr::write.fasta(all_sim_reads$sequence, all_sim_reads$name, "mifaser/data/test_reads.fa")
saveRDS(all_sim_reads, "mifaser/data/test_reads_metadata.RDS")
quit()
# run misafter with custom database
system("mifaser -D mifaser/database/fusion_train mifaser/data/fusion_train.fa")
system("mifaser -n -d mifaser/database/fusion_train -f mifaser/data/test_reads.fa -o mifaser/outputs/fusion_database.out")
