## ---------------------------
## Purpose of script: Download assemblt information from NCBI
## Author: Henri Chung
## Date Created: 2021-02-09
## Date Modified: 2024-01-04
## ---------------------------

# Package names
packages <- c("data.table", "tidyverse", "mefa4", "Matrix", "parallel", "proxyC")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
lapply(packages, function(x){suppressPackageStartupMessages(library(x, character.only = TRUE))})

# clear working directory
rm(list = ls())

# limit to fusion data with kegg annotations
kf_assemblies <- readLines("data/kegg/kegg_fusion_assemblies") %>% gsub("\\.[0-9]", "", .)

# Load Genbank assembly accession 
# Assembly features
dataFolder <- "data/"
assembly_summary_filename <- paste(dataFolder, "fusion/assembly_summary_genbank.txt", sep = "")
assembly_summary <- data.table::fread(assembly_summary_filename) %>%
  janitor::clean_names() %>%
  mutate(aa = gsub("\\.[0-9]+", "", number_assembly_accession)) %>%
  as.data.table()

# remove version number from genbank assembly accessions
subset_assembly_summary <- filter(assembly_summary, aa %in% kf_assemblies) %>%
  rowwise() %>%
  mutate(asm_name = gsub("M2/40", "M2_40", asm_name)) %>%
  mutate(asm_name = gsub(" ", "_", asm_name)) %>%
  filter(ftp_path != "na") %>%
  mutate(download_path = paste(ftp_path, "/", number_assembly_accession, "_", asm_name, sep = ""))
#GCA_001029795
# separate out ftp links
subset_ftps <- gsub("\\(T\\)", "_T", unique(subset_assembly_summary$download_path))
asm_names <- subset_assembly_summary$asm_name
feature_table_suffix <- "_feature_table.txt.gz"
assembly_stats_suffix <- "_assembly_stats.txt"
protein_fa_suffix <- "_protein.faa.gz"
genomic_suffix <- "_genomic.fna.gz"
# loop through ftp links and download assembly stats and feature tables for all assemblies.
message("Starting downloads.")
n = length(subset_ftps)
for(i in 1:length(subset_ftps)){
  message(Sys.time(), " ", i, "/", n) 
  asm <- str_split(subset_ftps[i], pattern = "/")[[1]]
  asm <- asm[length(asm)]

  ft_url <- paste(subset_ftps[i], feature_table_suffix, sep = "")
  ft_des <- paste("data/feature_tables/", asm, feature_table_suffix, sep = "")
  tryCatch({
    download.file(ft_url, destfile = ft_des, method = "wget", extra = "-nc", quiet = TRUE)
  }, error = function(e) {message("")
  })

  as_url <- paste(subset_ftps[i], assembly_stats_suffix, sep = "")
  as_des <- paste("data/assembly_stats/", asm, assembly_stats_suffix, sep = "")
  tryCatch({
    download.file(as_url, destfile = as_des, method = "wget", extra = "-nc", quiet = TRUE)
  }, error = function(e) {message("")
  })

  pf_url <- paste(subset_ftps[i], protein_fa_suffix, sep = "")
  pf_des <- paste("data/protein_fa/", asm, protein_fa_suffix, sep = "")
  tryCatch({
    download.file(pf_url, destfile = pf_des, method = "wget", extra = "-nc", quiet = TRUE)
  }, error = function(e) {message("")
  })

  gn_url <- paste(subset_ftps[i], genomic_suffix, sep = "")
  gn_des <- paste("data/genomes/", asm, genomic_suffix, sep = "")
  tryCatch({
    download.file(gn_url, destfile = gn_des, method = "wget", extra = "-nc", quiet = TRUE)
  }, error = function(e) {message("")
  })
}
quit()

#grep "Primary Assembly[[:space:]]\+[I|na|1|circular|large|eg_1|sequence1|gsn.131|DSM 122|M2/40_rep1|cPNK|MARIT|TJEJU|DPRO|Kuenenia_stuttgartiensis_MBR1]\+[[:space:]]\+Chromosome[[:space:]]\+assembled-molecule[[:space:]]\+total-length"  data/assembly_stats/*.txt > data/assembly_stats/assembly_lengths.tsv
assembly_stats <- data.table::fread("data/assembly_stats/assembly_lengths.tsv")[,c("V1", "V6")]
colnames(assembly_stats) <- c("filename", "length")
assembly_stats[,filename := stringr::str_extract(filename, "GCA_\\d+")]
data.table::fwrite(assembly_stats, "data/assembly_stats/assembly_lengths.tsv")
quit()
# Concatenate feature tables into single file.
feature_table_files <- list.files("data/feature_tables", pattern = "*.txt", full.names = TRUE) # files to search
subset_search <- paste(unique(subset_data$assembly_accession), collapse = "|") # accession numbers to search for in downloaded files
subset_files <- grepl(subset_search, feature_table_files) # list of subset files
feature_table_files2 <- feature_table_files[subset_files] # list of subset feature tables
feature_table_list <- lapply(feature_table_files, data.table::fread) # read in feature tables
names(feature_table_list) <- feature_table_files2 # name
feature_table_list <- feature_table_list[lengths(feature_table_list) != 0] # remove null elements

# join feature tables together, rename columns
feature_table_list <- lapply(feature_table_list, function(.x){
  colnames(.x) <- c("feature", "class", "assembly", "assembly_unit", "seq_type", 
  "chromosome", "genomic_accession", "start", "end", "strand", "product_accession", 
  "non-redundant_refseq", "related_accession", "name", "symbol", "GeneID", "locus_tag",
  "feature_interval_length", "product_length", "attributes")
  return(.x)
  })

# save feature list as RDS objects
saveRDS(feature_table_list, "data/feature_table_list.RDS")

# bind assembly lengths to feature table
#assembly_lengths <- read_csv("data/fusion/assembly_lengths.csv")
assembly_lengths_dt <- as.data.table(assembly_lengths)
feature_table_dt <- mclapply(feature_table_list, function(.x){
  res <- .x[feature == "CDS" & seq_type == "chromosome", ]
  setnames(res, "product_accession", "ncbi_accession")
  setnames(res, "assembly", "assembly_accession")
  res[subset_data, on = "ncbi_accession", fusion_lvl_1 := i.fusion_lvl_1]
  res[assembly_lengths_dt, on = "assembly_accession", length := i.length]
  mycols <- c("assembly_accession", "start", "end", "ncbi_accession", "fusion_lvl_1", "length")
  res <- res[, ..mycols]
  return(res)
  }, mc.cores = 8) %>%
  bind_rows() %>%
  filter(!is.na(fusion_lvl_1)) %>%
  as.data.table()
  
saveRDS(feature_table_dt, "data/feature_table_dt.RDS")

# section to fetch assembly lengths for suppressed assemblies
suppressed <- filter(subset_data, !assembly_accession %in% assembly_lengths_df$assembly_accession) %>%
  pull(assembly_accession) %>%
  unique()
ncbi_prefix <- "https://www.ncbi.nlm.nih.gov/assembly/"

nums <- NULL
for(i in 1:length(suppressed)){
  url <- paste(ncbi_prefix, suppressed[i], sep = "")
  temp <- readLines(con = url)
  line <- temp[grepl(temp, pattern = "Total sequence")] %>% 
   str_extract(pattern = "align_r\\\">.*</td></tr><tr><td>Total ungapped") 
  line <- gsub("align_r\\\">",  "", line)
  line <- gsub("</td></tr><tr><td>Total ungapped", "", line)
  nums[i] <- as.numeric(gsub(",", "", line))
}

suppressed_df <- tibble(
  length = nums,
  assembly_accession = suppressed
  )
