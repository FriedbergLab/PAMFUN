## ---------------------------
## Purpose of script: Generate test set of proteins for mifaser evaluation
## Author: Henri Chung
## Date Created: 2022-09-21
## Date Modified: 2024-01-04
## ---------------------------
library(data.table)
library(tidyverse)
# clear working directory
rm(list = ls())

# read in background data
balanced_data <- fread("compare/data/balanced_data.csv")

# read in mifaser protein fasta
count <- system("wc -l mifaser/data/protein_fa.fasta", intern = TRUE) %>%
  strsplit(" ") 
if(count[[1]][1] == "0") {
  system("touch mifaser/data/full_protein_fa.fasta")
  system("zcat data/protein_fa/* >> mifaser/data/protein_fa.fasta")
}

# read in fasta file of proteins
all_proteins <- seqinr::read.fasta("mifaser/data/protein_fa.fasta")
all_proteins_names <- gsub("\\.[0-9]", "", names(all_proteins))

# find precision/recall of mifaser for fusions on training set
# run on metagenome.

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

# create table of proteins with mifaser names
mifaser_names <- (balanced_data
  [ncbi_accession2 %in% all_proteins_names] # filter to proteins in fu
  [,fake_ec := sapply(fusion_lvl_1, insert_decimals)] # replace fusion_lvl_1 with fake ec number
  [,protein_name := gsub("\\|", ",", protein_name)] # replace | in protein_name 
  [,mifaser_name := paste(ncbi_accession2, "|", assembly_accession2, protein_name, "|", fake_ec, "|", fusion_lvl_1)] # create mifaser names
)
# match order of train_proteins to mifaser_names
mifaser_names <- mifaser_names[match(all_proteins_names, mifaser_names$ncbi_accession2),]
# assign names
names(all_proteins) <- mifaser_names$mifaser_name
# remove NA
all_proteins <- all_proteins[!is.na(names(all_proteins))]
# change sequence to uppercase
all_proteins2 <- lapply(all_proteins, toupper)
#write to file
seqinr::write.fasta(all_proteins2, names = names(all_proteins2), "mifaser/data/fusion_db.fa")

# create mifaser custom database
system("mifaser -D /mifaser/database/fusion_database mifaser/data/fusion_db.fa")
# run misafter with custom database
system("mifaser -d mifaser/fusion_database -f mifaser/data/test_reads.fa -o mifaser/test_out")
quit()

