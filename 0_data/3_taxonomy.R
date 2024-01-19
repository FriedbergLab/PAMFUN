## ---------------------------
## Purpose of script: Maps fusion assemblies to GTDB Taxonomy.
## Author: Henri Chung
## Date Created: 2021-01-12
## Date Modified: 2024-01-04
## ---------------------------

# Package names
packages <- c("data.table", "readr")

# Packages loading
lapply(packages, function(x){suppressPackageStartupMessages(library(x, character.only = TRUE))})

# clear working directory
rm(list = ls())

# Read in fusion data
message(Sys.time(), ": Reading in data")
fusion_data <- data.table::fread("data/fusion/fusion_data.tsv")

# wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/qc_failed_r207.tsv
fusion_assemblies <- (
	unique(fusion_data[,("assembly_accession"), with = FALSE])
	[,accession := gsub("GCA_|\\.[0-9]", "", assembly_accession)]
	)
fusion_accessions <- fusion_assemblies[,accession]


# download the GTDB taxonomy file
# wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/bac120_taxonomy_r207.tsv.gz
# gzip -d bac120_taxonomy_r207.tsv.gz
# List of genome assemblies failing the internal GTDB QC criteria.

# read in GTDB data
taxa_reference <- (data.table::fread("data/taxonomy/bac120_taxonomy_r207.tsv", header = FALSE)
	[, c("gtdb_accession", "domain") := data.table::tstrsplit(V1, "\\t")]
	[, 1:9 := lapply(.SD, function(x) gsub("^.{1,2}_{1,2}", "", x))]
	[, assembly_accession := gsub("GCF_", "GCA_", gtdb_accession)]
	)
colnames(taxa_reference) <- c("x", "phylum", "class", "order", "family", "genus", "species", "gtdb_accession", "domain", "assembly_accession")
taxa_cols <- c("gtdb_accession", "assembly_accession", "domain", "phylum", "class", "order", "family", "genus", "species")
taxa_reference <- taxa_reference[,..taxa_cols]

# load quality control failed taxonomy metadata for assemblies
qc_fail <- (data.table::fread("data/taxonomy/qc_failed_r207.tsv")
	[,c("Accession", "GTDB taxonomy")]
	[, c("domain", "phylum", "class", "order", "family", "genus", "species") := data.table::tstrsplit(`GTDB taxonomy`, "; ")]
	[, 1:9 := lapply(.SD, function(x) gsub("^.{1,2}_{1,2}", "", x))]
	[, assembly_accession := gsub("GCF_", "GCA_", Accession)]
	)
colnames(qc_fail) <- c("gtdb_accession", "x", "domain", "phylum", "class", "order", "family", "genus", "species", "assembly_accession")
qc_cols <- c("gtdb_accession", "assembly_accession", "domain", "phylum", "class", "order", "family", "genus", "species")
qc_reference <- qc_fail[,..qc_cols]

# bind taxonomy data together.
gtdb_reference <- (as.data.table(rbind(taxa_reference, qc_reference))
	[, accession := gsub("GCA_|\\.[0-9]", "", assembly_accession)]
	[accession %in% fusion_accessions]
	)
fwrite(gtdb_reference, "data/taxonomy/gtdb_reference.csv")
