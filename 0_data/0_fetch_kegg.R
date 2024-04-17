## ---------------------------
## Purpose of script: Download pathway information from KEGG
## Author: Henri Chung
## Date Created: 2021-01-12
## Date Modified: 2024-01-04
## ---------------------------

# Package names
packages <- c("rvest", "rjson", "tidyverse", "data.table", "stringr", "parallel")

# Packages loading
lapply(packages, function(x){suppressPackageStartupMessages(library(x, character.only = TRUE))})

# clear working directory
rm(list = ls())
source("code/0_fetch_kegg_helper.R")
outputFolder <- "data/kegg/"

#####

# Identify Modules specific to bacteria
# --------------------------- 

# read in assembly file 
assemblies_df <- data.table::fread("./data/balanced_organism_set.treemmer-RTL-0.9.gtdb.incl_taxonomy.tsv", keepLeadingZeros = TRUE)

# First, fetch all KEGG organism information 
message(Sys.time(), " COLLECTING MODULE DATA")
kegg_taxa <- read_html("https://rest.kegg.jp/list/organism") %>%
	html_text() %>%
	strsplit(split = "\n") %>%
	as.data.frame(col.names = "text") %>%
	separate(text, into = c("kegg_id", "kegg_abbr", "species", "taxonomy"), sep = "\t") %>%
	as.data.table()

# filter to bacteria
kegg_bacteria <- kegg_taxa[taxonomy %like% "Bacteria"][,kegg_id]

# split the bacterial species in querable chunks in groups of 10
bac_list <- split(kegg_bacteria, ceiling(seq_along(kegg_bacteria)/10))
bac_queries <- lapply(bac_list, function(.x){
	combo <- paste(paste("gn:", .x, sep = ""), collapse = "+")
	paste("https://rest.kegg.jp/get/", combo, sep = "")
	}) 

# fetch the bacteria organism entry pages 
bac_pages <- list()
for(i in 1:length(bac_queries)){
	bac_pages[[i]] <- multi_kegg(bac_queries[[i]], pattern = "T[0-9]{5}")
	Sys.sleep(1)
	message(i)
}
bac_pages2 <- unlist(bac_pages, recursive = FALSE)
saveRDS(bac_pages2, "./data/kegg/bacterial_kegg_pages.RDS")

# filter the bacterial entries in KEGG to our dataset using assembly accession
kegg_assemblies <- lapply(bac_pages2, function(.x){
	string <- trimws(.x$data_source[1])
 	assemblies <- regmatches(string, regexpr("GC[A|F]_[0-9]{9}\\.\\d", string))
	if(length(assemblies) == 0){assemblies <- "NO REF"}
 	org_code <- trimws(.x$org_code)
	taxonomy <- gsub("TAX:", "", trimws(.x$taxonomy)[1]) 
 	return(c(org_code, assemblies, taxonomy))
	}) %>%
	bind_rows() %>%
	t() %>% 
	as.data.frame() %>%
	rownames_to_column("kegg_id") %>%
	as.data.table()

# assemblies in kegg with the org code annotation
colnames(kegg_assemblies) <- c("kegg_id", "org_code", "assembly_id", "tax_id")
kegg_assemblies[,assembly_id2 := gsub("\\.[0-9]", "", assembly_id)] 
data.table::fwrite(kegg_assemblies, "./data/kegg/kegg_assemblies.csv")


fusion_data <- data.table::fread("./data/fusion/fusion_data.tsv")
balanced_organism_accessions <- paste("GCA_", readLines("data/fusion/balanced_organism_accessions"), sep = "")
kegg_fusion_data <- (fusion_data
	[,assembly_accession2 := gsub("\\.[0-9]", "", assembly_accession)]
	[,ncbi_accession2 := gsub("\\.[0-9]", "", ncbi_accession)]
	[assembly_accession2 %in% kegg_assemblies$assembly_id2]
	[assembly_accession2 %in% balanced_organism_accessions]
)
# fusion data subset to only include kegg entries
data.table::fwrite(kegg_fusion_data, "data/fusion/kegg_fusion_data.tsv")
kegg_fusion_assemblies <- unique(kegg_fusion_data[,assembly_accession])
# the assemblies in both the fusion and kegg dataset
write.table(kegg_fusion_assemblies, "data/kegg/kegg_fusion_assemblies", col.names = FALSE, row.names = FALSE, quote = FALSE)

# loop through the organisms from KEGG and retrieve revalent information
org_codes <- as.character(kegg_assemblies[assembly_id2 %in% kegg_fusion_data$assembly_accession2][,org_code])

pathways_list <- lapply(org_codes, retry_kegg)
names(pathways_list) <- org_codes

# manually check for errors
pathway_errors <- pathways_list[unlist(lapply(pathways_list, function(x) "simpleError" %in% class(x)))] %>% names()
while(length(pathway_errors) > 0){
	pathways_list[pathway_errors] <- lapply(pathway_errors, retry_kegg)
	pathway_errors <- pathways_list[unlist(lapply(pathways_list, function(x) "simpleError" %in% class(x)))] %>% names()
}
saveRDS(pathways_list, "data/kegg/organism_list.RDS")

full_org_codes <- as.character(kegg_assemblies[,org_code])
full_pathways_list <- parallel::mclapply(full_org_codes, retry_kegg, mc.cores = 4)
names(full_pathways_list) <- full_org_codes
pathways_list <- full_pathways_list[names(pathways_list)]
# manually check for errors
full_pathway_errors <- full_pathways_list[unlist(lapply(full_pathways_list, function(x) "simpleError" %in% class(x)))] %>% names()
while(length(full_pathway_errors) > 0){
	full_pathways_list[full_pathway_errors] <- lapply(full_pathway_errors, retry_kegg)
	full_pathway_errors <- full_pathways_list[unlist(lapply(full_pathways_list, function(x) "simpleError" %in% class(x)))] %>% names()
}
saveRDS(full_pathways_list, paste0(outputFolder, "full_pathways_list.RDS"))

# list out all KEGG protein ids involved.
protein_list <- lapply(pathways_list, function(.x){.x$protein})
protein_df <- data.table::rbindlist(protein_list)
protein_dt <- (protein_df
	[,org_code := gsub(":.*", "", kegg_protein)]
	[kegg_assemblies, tax_id := i.tax_id, on = "org_code"]
	[, kegg_protein2 := gsub(".*:", "", kegg_protein)]
	[,string_id := paste(tax_id, kegg_protein2, sep = ".")] # this doesnt work need to find aliases ecj:JW0001
	[,ncbi_protein := gsub("ncbi-proteinid:", "", ncbi_protein)]
	[(ncbi_protein %in% gsub("\\.[0-9]", "", kegg_fusion_data$ncbi_accession2))]
	)
data.table::fwrite(protein_dt, "data/kegg/kegg_proteins.csv")

unique(protein_df$org_code)[!(unique(protein_df$org_code) %in% unique(protein_dt$org_code))]
# Fetch list of all modules
module_list <- read_html("https://rest.kegg.jp/list/module") %>%
		html_text() %>%
		read.table(text = ., sep = "\t") %>%
		pull(V1) %>%
		gsub("md:", "", .) %>%
		unique() %>%
		as.list()

# subset modules into querable chunks of size 10
module_queries <- split(module_list, ceiling(seq_along(module_list)/10)) %>%
	lapply(., function(.x){
		combo <- paste(.x, collapse = "+")
		paste("https://rest.kegg.jp/get/", combo, sep = "")
		})

# fetch modules
module_pages <- list()
for(i in 1:length(module_queries)){
	module_pages[[i]] <- multi_kegg(module_queries[[i]], pattern = "M[0-9]{5}")
	message(i)
}
module_pages2 <- unlist(module_pages, recursive = FALSE)
names(module_pages2) <- module_list
saveRDS(module_pages2, "data/kegg/module_pages.RDS")

module_classes <- lapply(module_pages2, function(.x){
	res <- trimws(gsub("CLASS", "", .x$class )) %>% paste(collapse = " ")
	if(res == ""){res <- names(.x)[grepl("class", names(.x))] %>% strsplit(split = "_") %>% unlist() %>% toupper(); return(res[2])}else{return(res)}
	}) %>% 
	stack() %>% 
	filter(values != "") %>% 
	separate(values, into = c("type", "class", "subclass"), sep = ";") %>% 
	select(-type) %>%
	rename("module_id" = ind)

fill_module_classes <- tibble(
	module_id = names(module_pages2)[!(names(module_pages2) %in% module_classes$module_id)],
	class = "Signature modules",
	subclass = "Signature modules"
)
data.table::fwrite(rbind(module_classes,fill_module_classes),  "data/kegg/module_classes.csv")

# extract the definition string of modules.
module_definitions <- lapply(module_pages2, function(.x){
	res <- trimws(gsub("DEFINITION", "", .x$definition )) %>% paste(collapse = " ")
	if(res == ""){res <- names(.x)[grepl("definition", names(.x))] %>% strsplit(split = "_") %>% unlist() %>% toupper(); return(res[2])}else{return(res)}
	}) %>% stack() %>% filter(values != "")
module_definitions_list <- as.list(module_definitions$values)
names(module_definitions_list) <- module_definitions$ind

# calculate the possible paths for each module
module_permutations <- lapply(module_definitions_list, decipher)
names(module_permutations) <- names(module_definitions_list)

# convert the paths to long format
module_permutations_long <- module_permutations %>%
	stack() %>%
	as_tibble() %>%
	rowwise() %>%
	mutate(steps = purrr::map(values, function(.x){
		steps <- strsplit(.x, split = " ")[[1]]
		return(steps)})) %>% as.data.table()

# separate long paths into list
module_permutations_list <- split(module_permutations_long, by = "ind")

# reformat pathways_list data
pathways_dt <- lapply(pathways_list, function(.x){
	if((is.character(.x$message)) | (length(.x$complete[[1]])== 0)){
		return(NULL)
	}else{
	temp <- (.x$module
		[.x$ortholog, ortholog := i.ortholog, on = "kegg_protein"]
		[.x$protein, ncbi_protein := i.ncbi_protein, on = "kegg_protein"]
		[,module_id := gsub("md:[a-z]+_", "", module_id)]
		[module_id %in% .x$complete[[1]]]
		[,ncbi_protein := gsub("ncbi-proteinid:", "", ncbi_protein)]
		[,ortholog := gsub("ko:", "", ortholog)]
		[, c("org", "kegg_protein") := tstrsplit(kegg_protein, ":", fixed=TRUE)]
	) 
	return(temp)
	}
	}) %>% purrr::compact()

# using pathways_dt and module_permutations_list
# convert the possible module permutations (in KOs) into 
# gene permutations.
permutations_list <- parallel::mclapply(pathways_dt, function(temp){
	message(Sys.time()," ", temp$org[1])
	temp_kos <- split(temp, by = "module_id") %>% lapply(function(.y){
		module <- .y$module_id[[1]]
		paths <- lapply(module_permutations_list[[module]]$steps, function(.z){
			if(sum(.z %in% .y$ortholog) == length(.z)){
					res <- .y[ortholog %in% .z] %>% 
						split(by = "ortholog") %>% 
						lapply(function(.a){pull(.a, ncbi_protein)}) %>% 
							expand.grid() %>% 
							apply(MARGIN = 1,FUN = paste, simplify = FALSE) %>%
							lapply(paste, collapse = " ")
					return(res)
				}else{
					return(NULL)}
			}) %>% unlist(recursive = FALSE)
		return(paths)
		}) %>% unlist(recursive = FALSE) %>% stack()
	colnames(temp_kos) <- c("protein_path", "module_id")
	return(temp_kos)
	}, mc.cores = 4) 
permutations_dt <- rbindlist(permutations_list, idcol = "org_code")[,module_id := substr(module_id, 1, 6)]
permutations_dt$length <- lapply(permutations_dt$protein_path, function(.x){length(unlist(strsplit(.x, " ")))})
data.table::fwrite(permutations_dt, "data/kegg/permutations_dt.csv")

# extract lengths by module_id group
module_lengths_list <- split(permutations_dt, by = "module_id")
module_lengths <- lapply(module_lengths_list, function(.x){median(unlist(.x$length))}) %>% stack() %>% rename("module_id" = ind, "median_length" = values)
write_csv(module_lengths, "data/kegg/module_lengths.csv")

# Parse module pages 
names(module_pages2) <- lapply(module_pages2, function(.x){name <- str_extract(.x$entry, "[M]\\d{5}")})
module_parts <- lapply(module_pages2, function(.x){
	temp <- trimws(.x$orthology)
	name <- str_extract(.x$entry, "[M]\\d{5}")
	res <- bind_rows(lapply(temp, separate_string))
	return(res)
})
module_parts <- bind_rows(module_parts, .id = "module_id")
# dataframe of modules, associated kegg orthologs, and their equivalent ec number
fwrite(module_parts, "data/kegg/module_parts.csv")
module_parts <- data.table::fread("data/kegg/module_parts.csv")

# Map modules to pathways (which modules belong to what pathway)
module_pathways <- read_html("https://rest.kegg.jp/link/pathway/module") %>%
	html_text() %>%
	read.table(text = ., sep = "\t") %>% 
	as_tibble() %>%
	mutate(V1 = gsub("md:", "", V1), V2 = gsub("path:", "", V2)) %>%
	rename("module_id" = V1, pathway = "V2") %>%
	as.data.table()
data.table::fwrite(module_pathways, "data/kegg/module_pathways.csv")


# Identify EC Information from KEGG
# --------------------------- 

# acquire list of EC numbers used in kegg db
kegg_ec <- read_html("https://rest.kegg.jp/list/ec") %>%
	html_text() %>%
	read.table(text = ., sep = "\t") %>%
	as.data.table()

# split the KEGG EC numbers into query-able splits
ec_numbers <- unique(c(
	unlist(str_extract_all(kegg_ec[,V1], "[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+")), 
	unlist(str_extract_all(kegg_ec[,V2], "[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+"))
))
kb_list <- split(ec_numbers, ceiling(seq_along(ec_numbers)/10))
# use the KEGG rest API to query the entries for these EC numbers.
ec_queries <- lapply(kb_list, function(.x){
	combo <- paste(.x, collapse = "+")
	paste("https://rest.kegg.jp/get/", combo, sep = "")
	})

# fetch ec numbers
ec_pages <- list(); sink("ec_page_progress")
for(i in 1:length(ec_queries)){ec_pages[[i]] <- multi_kegg(ec_queries[[i]], pattern = "[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+"); print(i)}

ec_pages <- parallel::mclapply(ec_queries, function(.x){
	res <- multi_kegg(.x, pattern = "[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+")
	print(which(ec_queries == .x))
	return(res)
	}, mc.cores = 4)

ec_pages2 <- unlist(ec_pages, recursive = FALSE)
names(ec_pages2) <- lapply(ec_pages, names) %>% unlist()
# reshape ec number with functional descriptions
ec_descs <- lapply(ec_pages2, function(.x){trimws(.x$name)})
names(ec_descs) <- names(ec_pages2)
ec_df <- as.data.table(stack(ec_descs))
colnames(ec_df) <- c("ec_desc", "ec_number")
data.table::fwrite(ec_df, "data/kegg/ec_descriptions.csv")

# Write dataframe of all protein-protein interactions
protein_interaction_list <- list()
for(i in 1:length(pathways_list)){
    temp <- pathways_list[[i]]  
    if(length(temp) == 1){message(i); protein_interaction_list[[i]] <- NA; next}
	complete_modules <- temp$complete[[1]]
    if(length(complete_modules) < 1){message(i); protein_interaction_list[[i]] <- NA; next}
	if(is.null(temp$module)){missing <- c(missing, names(pathways_list)[[i]]); protein_interaction_list[[i]] <- NA; next}
    modules_to_proteins <- temp$module[, module_id := gsub(".*_", "", module_id)][module_id %in% complete_modules]
    temp_df <- (temp$protein
        [,ncbi_protein := gsub(".*:", "", ncbi_protein)]
        [modules_to_proteins, on = .(kegg_protein)]
        [, c("org", "kegg_protein") := tstrsplit(kegg_protein, ":", fixed=TRUE)])
    protein_interaction_list[[i]] <- temp_df %>% 
        group_by(module_id, org) %>% 
        nest() %>% 
        mutate(links = map(data, function(.x){
            res <- expand.grid(.x$ncbi_protein, .x$ncbi_protein) %>%
                mutate(string = paste(Var1, Var2, sep = "-"))
            })) %>%
		select(-data) %>%
        unnest(links) #%>%
        #unnest(data)
}
names(protein_interaction_list) <- names(pathways_list)
saveRDS(protein_interaction_list, "data/kegg/protein_interaction_list.rds")
protein_interaction_list2 <- protein_interaction_list[!is.na(protein_interaction_list)]
protein_links_df <- bind_rows(protein_interaction_list2)
data.table::fwrite(protein_links_df, "data/kegg/protein_links_df.csv")

# write proteins in a module to file
module_proteins <- as.character(unique(c(protein_links_df$Var1, protein_links_df$Var2)))
writeLines(module_proteins, "data/kegg/module_proteins.txt")
# write module id and descriptions
module_data <- read_html("https://rest.kegg.jp/list/module") %>%
	html_text() %>%
	read.table(text = ., sep = "\t")
colnames(module_data) <- c("module_id", "module_desc")
write.csv(module_data, "data/kegg/module_data.csv", row.names = FALSE)

pathway_data <- read_html("https://rest.kegg.jp/list/pathway") %>%
	html_text() %>%
	read.table(text = ., sep = "\t")
colnames(pathway_data) <- c("pathway_id", "pathway_desc")
write.csv(pathway_data, "data/kegg/pathway_data.csv", row.names = FALSE)

permutations_dt <- data.table::fread("data/kegg/permutations_dt.csv")
permutations_dt[, n := lapply(protein_path, function(.x){stringr::str_count(.x,  " ")})]
setkey(permutations_dt, "module_id")
# filter to only module_ids with only 1 unique n value
single_step_modules <- unique(permutations_dt[n == 0]$module_id)
unique_steps_modules <- permutations_dt[, .(n2 = length(unique(n))), by = module_id]
only_single_step_modules <- unique(unique_steps_modules[n2 == 1][module_id %in% single_step_modules]$module_id)
writeLines(only_single_step_modules, "data/kegg/only_single_step_modules.txt")

# identify non-orthologous replacements
pathways_list <- readRDS("data/kegg/organism_list.RDS")
orthology_to_protein <- list()
for(i in 1:length(pathways_list)){
	temp_ortholog <- pathways_list[[i]]$ortholog
	temp_protein <- pathways_list[[i]]$protein
	temp <- right_join(temp_ortholog, temp_protein, by = "kegg_protein")
	orthology_to_protein[[i]] <- temp[, ortholog := gsub("ko:", "", ortholog)][, ncbi_accession2 := gsub(".*:", "", ncbi_protein)][, .(ortholog, ncbi_accession2)]
}
orthology_to_protein <- rbindlist(orthology_to_protein)
data.table::fwrite(orthology_to_protein, "data/kegg/orthology_to_protein.csv")

setkey(orthology_to_protein, "ortholog")
module_steps <- lapply(module_definitions_list, function(.x){
	step_list <- decipher2(.x)[[1]]
	})  %>% unlist(recursive = FALSE) %>% unique() 
# remove list elements of length <= 1
module_steps2 <- module_steps[lengths(module_steps) > 1] %>% lapply(function(.x){
	# remove elements that contain " "
	.y <- .x[!grepl(" ", .x)]
	res <- data.table::CJ(.y, .y)
	if(!is.data.table(res)){return(NULL)}
	colnames(res) <- c("a", "b")
	return(res)
	}) %>% rbindlist() 
data.table::fwrite(module_steps2, "data/kegg/unique_non_orthologous_replacements.csv")

module_steps3 <- parallel::mclapply(module_definitions_list, function(.x){
	step_list <- decipher2(.x)[[1]]
	nor <- lapply(step_list, function(.y){
		# remove elements that contain " "
		.y <- .y[!grepl(" ", .y)]
		res <- data.table::CJ(.y, .y)
		if(!is.data.table(res)){return(NULL)}
		colnames(res) <- c("a", "b")
		return(res)
		}) %>% rbindlist() 
	nor_sum <- nor[, combo := paste(a, b, sep = "_")]
	return(nor_sum)
	}, mc.cores = 8) %>% rbindlist(id = "module_id")
data.table::fwrite(module_steps3, "data/kegg/non_orthologous_replacements.csv")