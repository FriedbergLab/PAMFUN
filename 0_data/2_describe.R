## ---------------------------
## Purpose of script: Generates desciptive statistics for fusion dataset.
## Author: Henri Chung
## Date Created: 2021-01-12
## Date Modified: 2024-01-04
## ---------------------------

# Package names
packages <- c("data.table", "tidyverse")

# Packages loading
lapply(packages, function(x){suppressPackageStartupMessages(library(x, character.only = TRUE))})

# clear working directory
rm(list = ls())

# load config file.
source("code/config.R")
# Read in fusion data
message(Sys.time(), ": Reading in data")
my_data <- data.table::fread(subset_filename); output_file <- "outputs/subset_stats.txt"

out <- file(output_file, open = "wt")
sink(out, type = "message")

# Count of objects
n_assemblies <- length(unique(my_data[,assembly_accession]))
message("Number of unique organisms/assemblies: ", n_assemblies)

n_proteins <- length(unique(my_data[,ncbi_accession]))
message("Number of unique proteins: ", n_proteins)

n_fusions <- length(unique(my_data[,fusion_lvl_1]))
message("Number of unique fusions: ", n_fusions)

single_a <- my_data[, if(.N==1) .SD, by=fusion_lvl_1]
message("Number of singletons by count: ", nrow(single_a))

single_b <-  my_data[fusion_lvl_1 %like% "S|SU"][, if(.N>1) .SD, by=fusion_lvl_1]
message("Number of singletons by designation (S|SU): ", nrow(single_b))

single_fusions <- rbind(single_a, single_b)
n_singles <- nrow(single_fusions)
message("Total number of singletons: ", n_singles)

doubles_fusions <- my_data[, if(.N==2) .SD, by=fusion_lvl_1]
n_doubles <- nrow(doubles_fusions)
message("Number of doubletons: ", n_doubles)

multi_fusions <- my_data[, if(.N > 2) .SD, by=fusion_lvl_1][!(fusion_lvl_1 %like% "S|SU")]
n_multi <- nrow(multi_fusions)
message("Number of fusions above singletons/doubletons: ", n_multi)
closeAllConnections() 

# Assembly stats
 
# number of proteins
protein_assemblies <- my_data[, by = assembly_accession, .(n_proteins = .N)][order(-n_proteins)]

# number of fusions
fusion_assemblies_cols <- c("fusion_lvl_1", "assembly_accession")
fusion_assemblies_a <- unique(my_data[, ..fusion_assemblies_cols])
fusion_assemblies <- fusion_assemblies_a[, by = assembly_accession, .(n_fusions = .N) ][order(-n_fusions)]

# number of single fusions
single_fusions_vec <- single_fusions[,fusion_lvl_1]
single_assemblies_a <- fusion_assemblies_a[ fusion_lvl_1 %in% single_fusions_vec] 
single_assemblies <- single_assemblies_a[, by = assembly_accession, .(n_single = .N) ][order(-n_single)]

# number of multi fusions
multi_assemblies_a <- fusion_assemblies_a[ !(fusion_lvl_1 %in% single_fusions_vec)] 
multi_assemblies <- multi_assemblies_a[, by = assembly_accession, .(n_multi = .N) ][order(-n_multi)]

# types of proteins
seqtype_cols <- c("assembly_accession", "ncbi_accession", "seq_type")
seqtype_assemblies <- unique(my_data[, ..seqtype_cols])
seqtype_assemblies[, count := 1]
seqtype_assemblies <- dcast(seqtype_assemblies, assembly_accession ~ seq_type, value.var = "count", fun.agg = sum)

# assembly information
assembly_metadata  <- (single_assemblies
	[fusion_assemblies, on = "assembly_accession"]
	[multi_assemblies, nomatch = 0, on = "assembly_accession"]
	[seqtype_assemblies, nomatch = 0, on = "assembly_accession"]
	)
fwrite(assembly_metadata, gsub("_stats.txt", "_assembly_metadata.tsv", output_file))

# fusion stats

# size of each fusion
protein_fusions_cols <- c("fusion_lvl_1", "ncbi_accession")
protein_fusions_a <- unique(my_data[, ..protein_fusions_cols])
protein_fusions <- protein_fusions_a[, by = fusion_lvl_1, .(n_proteins = .N) ][order(-n_proteins)]


# number of participating assemblies
assemblies_fusions <- fusion_assemblies_a[, by = fusion_lvl_1, .(n_assemblies = .N) ][order(-n_assemblies)]

# EC - Fusion information + analysis

# read in fusion to hfsp mapping
seguid_to_ec <- data.table::fread(seguid_to_ec_filename)[hfsp >= hfsp_threshold][,ec_number := gsub("n", "", ec_number)]

# Filter out EC to fusion annotations 
# Filter out EC annotations to fusion clusters where the EC annotations,
# is under a threshold of the total proportion of EC annotations in a cluster.
fusion_ec_cols <- c("fusion_lvl_1", "ec_number", "seguid")
fusion_to_ec2 <- unique(my_data[seguid_to_ec, on = "seguid", ec_number := i.ec_number]
  [,..fusion_ec_cols])
  [!is.na(fusion_lvl_1)]
  [, by = c("fusion_lvl_1", "ec_number"), .(n_proteins = .N) ]
  [, by = "fusion_lvl_1", .(prop = n_proteins / sum(n_proteins), ec_number = ec_number)]
  [prop > ec_threshold]
  [,prop := NULL]
  )
fwrite(fusion_to_ec, gsub("_stats.txt", "_fusion_to_ec.tsv", output_file))

# written description of ecs
ec_descriptions <- data.table::fread(ecfunction_raw_filename)

# Calculate how many EC numbers are assigned to multiple fusion clusters.
ec_duplicates <- fusion_to_ec[, if(.N > 1) .SD, by=ec_number]

# Calculate how many fusion clusters have more than one EC numbers.
fusion_duplicates <- fusion_to_ec[, if(.N > 1) .SD, by=fusion_lvl_1]

# seqtype breakdown of fusions;
seqtype_cols <- c("fusion_lvl_1", "ncbi_accession", "seq_type")
seqtype_fusions <- unique(my_data[, ..seqtype_cols])
seqtype_fusions[, count := 1]
seqtype_fusions <- dcast(seqtype_fusions, fusion_lvl_1 ~ seq_type, value.var = "count", fun.agg = sum)


# fusion information
fusion_metadata <- (protein_fusions
	[assemblies_fusions, on = "fusion_lvl_1"]
	[fusion_to_ec, on = 'fusion_lvl_1', ec_number := i.ec_number]
	[ec_descriptions, on = "ec_number", ec_desc := i.desc]
	[seqtype_fusions, nomatch = 0, on = "fusion_lvl_1"]
	)

fwrite(fusion_metadata, gsub("_stats.txt", "_fusion_metadata.tsv", output_file))

# fusion descriptions
message("Number of fusions with EC annotations: ", length(unique(fusion_metadata$fusion_lvl_1)))
message("Number of unique ECs: ", length(unique(fusion_metadata$ec_number)))

ec_freq <- fusion_metadata[, by = ec_number, .(n_fusions = .N) ][order(-n_fusions)]

message("Most annotated ECs:")
head(ec_freq)


# fusions with the most and least copies
fusion_ratio <- fusion_metadata[, ratio := n_proteins / n_assemblies]


