## ---------------------------
## Purpose of script: Parse kaiju output files of metagenome samples.
## Author: Henri Chung
## Date Created: 2022-09-21
## Date Modified: 2024-01-04
## ---------------------------

library(tidyverse)
library(data.table)
rm(list = ls())

# read in files
balanced_organism_accessions <- paste("GCA_", readLines("data/fusion/balanced_organism_accessions"), sep = "")

# taxon_id to modules
org_modules <- readRDS("kaiju/outputs/org_module.RDS") %>%
    rename("taxon_id" = tax_id) %>%
    filter(assembly_id2 %in% balanced_organism_accessions)

# fusion to ec numbers
fusion_ec_data <- data.table::fread("kaiju/data/fusion_ec_data.tsv")[,fusion_lvl_1 := as.character(fusion_lvl_1)]

# modules to ec numbers
module_parts <- data.table::fread("data/kegg/module_parts.csv") %>%
    select(module_id, ec_number) %>%
    unique() %>%
    filter(ec_number %in% unique(fusion_ec_data$ec_number))

# sample 
PRJNA385854_metadata <- read.csv("kaiju/outputs/PRJNA385854_metadata.csv")

#####################################################

# read in kaiju taxonomy reference file
options(warn = 1)
taxa_reference <- data.table::fread("data/taxonomy/gtdb_reference.csv")
# read in kaiju taxonomy summary files
kaiju_files <- list.files("kaiju/outputs", pattern = "*.summary", full.names = TRUE) 
kaiju_taxa <- lapply(kaiju_files, function(.x){
    taxa <- data.table::fread(.x, col.names = c("taxon_id", "count", "percent"))
    if(nrow(taxa) == 0){
        return(c())
    }
    temp <- org_modules[taxa, on =   "taxon_id"][!is.na(kegg_id)]
})

# sort kaiju by phyla
kaiju_taxa_data <- lapply(kaiju_taxa, function(.x){
    temp <- (.x
        [, assembly_accession := assembly_id]
        [taxa_reference, phylum := i.phylum, on = "assembly_accession"]
        [,c("phylum", "percent")]
        [,.(percent = sum(percent)), by = phylum]
        [percent < 1, phylum := "Other"]
        [!is.na(phylum)]
        [,.(percent = sum(percent)), by = phylum])
    # add new row phylum is "unknown" and percent is the 1 - the sum of all other phyla
    temp <- rbind(temp, data.table(phylum = "Unknown", percent = 100 - sum(temp$percent)))
    return(temp)
})
names(kaiju_taxa_data) <- gsub("\\.filtered\\.out\\.summary", "", basename(kaiju_files))
kaiju_taxa_df <- bind_rows(kaiju_taxa_data, .id = "ncbi_sra_accession") %>%
    left_join(PRJNA385854_metadata, by = "ncbi_sra_accession") %>%
    mutate(collection_date = as.POSIXct(collection_date, format="%Y-%m-%dT%H:%M:%S"))

# Plot phyla distribution of kaiju results
# create a vector of unique kaiju phylums with Other and Unknown last
kaiju_phylums <- kaiju_taxa_df %>%
    filter(phylum != "Other" & phylum != "Unknown") %>%
    pull(phylum) %>%
    unique() %>%
    c("Other", "Unknown")
kaiju_pallete <-  c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#999999", "#CCCCCC")


# Analyze EC and module functionality of kaiju results
# list EC numbers for each sample
kaiju_list <- lapply(kaiju_taxa, function(temp){
    kaiju_ec <- lapply(compact(temp[percent >= 1]$data), function(.y){
        modules_ec <- unique(left_join(.y, module_parts, by = "module_id"))
        ec_number <- unique(modules_ec$ec_number) 
        return(ec_number)
        }) %>%
        unlist() %>%
        unique()
    return(kaiju_ec)
    })
names(kaiju_list) <- gsub("\\.filtered\\.out\\.summary", "", basename(kaiju_files))
kaiju_list <- compact(kaiju_list)
# list modules for each sample
kaiju_modules <- lapply(kaiju_taxa, function(temp){
    res <- lapply(compact(temp[percent >= 1]$data), function(.y){return(.y$module_id)})%>%
        unlist() %>%
        unique()
})
names(kaiju_modules) <- gsub("\\.filtered\\.out\\.summary", "", basename(kaiju_files))

kaiju_list %>% lapply(length) %>% unlist() %>% summary()
kaiju_modules %>% lapply(length)  %>% unlist() %>% summary()

kaiju_list %>% unlist() %>% unique() %>% length()
kaiju_modules %>% unlist() %>% unique() %>% length()

module_classes <- data.table::fread("data/kegg/module_classes.csv")
kaiju_modules_df <- kaiju_modules %>% 
    unlist() %>% 
    table() %>% 
    stack() %>% 
    rename("module_id" = ind, "count" = values) %>% 
    left_join(module_classes, by = "module_id") %>%
    filter(!is.na(class)) %>%
    mutate(percent = count / 472 * 100) %>%
    mutate(quantile = case_when(
        percent <= 10 ~ "0-10%",
        percent > 10 & percent <= 50 ~ "10-50%",
        percent > 50 & percent <= 90 ~ "50-90%",
        percent > 90 ~ "90-100%"
    )) %>%
    group_by(quantile, class) %>%
    summarize(n = n()) %>%
    group_by(quantile) %>%
    mutate(perc = n / sum(n) * 100) 

# list of modules in each sample based on species
kaiju_modules
# list of EC numbers in each sample based on species   
kaiju_list
# module metadata
module_classes

#####################################################

# read in mifaser files
mifaser_files <- list.files("metagenome/mifaser_output", pattern = "analysis.tsv", full.names = TRUE, recursive = TRUE)
mifaser_list <- lapply(mifaser_files, function(.x){
    temp <- (data.table::fread(.x, header = TRUE, col.names = c("fusion_lvl_1", "count"))
        [,fusion_lvl_1 := gsub("^0*", "", gsub("\\.", "",fusion_lvl_1))]
        [fusion_ec_data, ec_number := i.ec_number, on = "fusion_lvl_1"]
        [,.(count = sum(count)), by = ec_number]
        [,percent := count/sum(count) * 100]
        [percent >= 0.1]
        [rev(order(percent))]) 
    return(temp)
})
names(mifaser_list) <- sapply(mifaser_files, function(.x){strsplit(.x, split = "\\/")[[1]][8]})
mifaser_list <- compact(mifaser_list)

# synce mifaser and kaiju samples 
ind <- names(kaiju_list)[(names(kaiju_list) %in% names(mifaser_list))]
kaiju_list <- kaiju_list[ind]
kaiju_modules <- kaiju_modules[ind]
mifaser_list <- mifaser_list[ind]

# test how many modules are detected in an unfiltered mifaser analysis that should be there according to a kaiju analysis.
mifaser_modules <- lapply(names(mifaser_list), function(.x){
    mifaser_temp <- mifaser_list[[.x]]
    kaiju_temp <- kaiju_modules[[.x]]
    module_temp <- module_parts[module_id %in% kaiju_temp]
    res <- module_temp %>%
        group_by(module_id) %>%
        group_split() %>%
        lapply(function(.y){
            temp <- sum(.y$ec_number %in% mifaser_temp$ec_number)/nrow(.y)
            names(temp) <- .y$module_id[1]
            return(temp)
        }) %>% 
        unlist() %>%
        stack() %>% 
        arrange(desc(values))
})
names(mifaser_modules) <- names(mifaser_list)

# bind detected mifaser modules with the number of steps in each module
module_count <- module_parts %>% group_by(module_id) %>% summarize(n = n())
mifaser_modules_df <- bind_rows(mifaser_modules, .id = "sample_id") %>% 
    rename("module_id" = ind, "percent" = values) %>%
    left_join(module_count, by = "module_id")

# filter list of all detectable modules to only kaiju predicted modules
filtered_kaiju_modules <- kaiju_modules %>%
    unlist() %>%
    unname()
filtered_kaiju_modules <- filtered_kaiju_modules[filtered_kaiju_modules %in% mifaser_modules_df$module_id]

#####################################################


# randomly sample EC numbers from the fusion data and calculate the percent of modules that are detected
set.seed(123)
random_mifaser_modules <- lapply(names(mifaser_list), function(.x){
    mifaser_temp <- mifaser_list[[.x]]
    kaiju_temp <- kaiju_modules[[.x]]
    module_temp <- module_parts[module_id %in% kaiju_temp]
    random_mifaser <- sample(unique(fusion_ec_data$ec_number), size = nrow(mifaser_temp), replace = TRUE)
    res <- module_temp %>%
        group_by(module_id) %>%
        group_split() %>%
        lapply(function(.y){
            temp <- sum(.y$ec_number %in% random_mifaser )/nrow(.y)
            names(temp) <- .y$module_id[1]
            return(temp)
        }) %>% 
        unlist() %>%
        stack() %>%
        arrange(desc(values))
})
random_mifaser_modules_df <- bind_rows(random_mifaser_modules, .id = "sample_id") %>% 
    rename("module_id" = ind, "percent" = values) %>%
    left_join(module_count, by = "module_id")


# compare mifaser with random
mifaser_detected_modules <- filter(mifaser_modules_df, percent == 1 & n > 1) %>% pull(module_id) %>% table() 
random_mifaser_detected_modules <- filter(random_mifaser_modules_df, percent == 1 & n > 1) %>% pull(module_id) %>% table()   

mifaser_detected_modules
random_mifaser_detected_modules

sum(mifaser_detected_modules) # 2509
sum(random_mifaser_detected_modules) # 84

################################
# test how many modules are detected in an unfiltered mifaser analysis that should be there according to a kaiju analysis.
unexpected_mifaser_modules <- lapply(names(mifaser_list), function(.x){
    mifaser_temp <- mifaser_list[[.x]]
    kaiju_temp <- kaiju_modules[[.x]]
    module_temp <- module_parts[!(module_id %in% kaiju_temp)]
    res <- module_temp %>%
        group_by(module_id) %>%
        group_split() %>%
        lapply(function(.y){
            temp <- sum(.y$ec_number %in% mifaser_temp$ec_number)/nrow(.y)
            names(temp) <- .y$module_id[1]
            return(temp)
        }) %>% 
        unlist() %>%
        stack() %>% 
        arrange(desc(values))
})
names(unexpected_mifaser_modules) <- names(mifaser_list)

# bind detected mifaser modules with the number of steps in each module
unexpected_mifaser_modules_df <- bind_rows(unexpected_mifaser_modules, .id = "sample_id") %>% 
    rename("module_id" = ind, "percent" = values) %>%
    left_join(module_count, by = "module_id")

unexpected_mifaser_detected_modules <- filter(unexpected_mifaser_modules_df, percent == 1 & n > 1) %>% pull(module_id) %>% table() 
sum(unexpected_mifaser_detected_modules) # 2226

set.seed(123)
random_unexpected_mifaser_modules <- lapply(names(mifaser_list), function(.x){
    mifaser_temp <- mifaser_list[[.x]]
    kaiju_temp <- kaiju_modules[[.x]]
    module_temp <- module_parts[!(module_id %in% kaiju_temp)]
    random_mifaser <- sample(unique(fusion_ec_data$ec_number), size = nrow(mifaser_temp), replace = TRUE)
    res <- module_temp %>%
        group_by(module_id) %>%
        group_split() %>%
        lapply(function(.y){
            temp <- sum(.y$ec_number %in% random_mifaser )/nrow(.y)
            names(temp) <- .y$module_id[1]
            return(temp)
        }) %>% 
        unlist() %>%
        stack() %>%
        arrange(desc(values))
})
names(random_unexpected_mifaser_modules) <- names(mifaser_list)

# bind detected mifaser modules with the number of steps in each module
random_unexpected_mifaser_modules_df <- bind_rows(random_unexpected_mifaser_modules, .id = "sample_id") %>% 
    rename("module_id" = ind, "percent" = values) %>%
    left_join(module_count, by = "module_id")

random_unexpected_mifaser_detected_modules <- filter(random_unexpected_mifaser_modules_df, percent == 1 & n > 1) %>% pull(module_id) %>% table()
sum(random_unexpected_mifaser_detected_modules) # 590


# count the unique number of ec modules detected in each sample
ec_list <- list()
for(i in 1:length(mifaser_list)){
    temp_sample <- names(mifaser_list)[i]
    temp_kaiju <- kaiju_list[[temp_sample]]
    temp_mifaser <- mifaser_list[[temp_sample]]$ec_number

    ec_list[[i]] <- list(mifaser_only = temp_mifaser[!(temp_mifaser %in% temp_kaiju)],
                kaiju_only = temp_kaiju[!(temp_kaiju %in% temp_mifaser)],
                both = temp_kaiju[temp_kaiju %in% temp_mifaser])
}
lapply(ec_list, function(.x){length(.x$mifaser_only)}) %>% unlist() %>% summary()
lapply(ec_list, function(.x){length(.x$kaiju_only)}) %>% unlist() %>% summary()
lapply(ec_list, function(.x){length(.x$both)}) %>% unlist() %>% summary()

#####################################################
# cluster the EC's detected by mifaser
mifaser_dt <- rbindlist(mifaser_list, idcol = "ncbi_sra_accession")
mifaser_dtm <- dcast(mifaser_dt, ec_number ~ ncbi_sra_accession, value.var = "count", fill = 0)
mifaser_mat <- as.matrix(mifaser_dtm[,-c("ec_number")])
mifaser_mat[mifaser_mat > 0] <- 1
rownames(mifaser_mat) <-mifaser_dtm$ec_number

mifaser_freqs <- rownames(mifaser_mat)[(rowSums(mifaser_mat) > 1)]

# cluster mifaser dist using mcl
mifaser_dist <- dist(mifaser_mat, method = "binary")
mifaser_pairs <- as.data.frame(as.matrix(mifaser_dist)) %>% 
    mutate(ec_number = rownames(mifaser_mat)) %>%
    gather(key = "ec_number2", value = "dist", -ec_number) %>%
    filter(ec_number != ec_number2) %>%
    mutate(sim = 1 - dist) %>%
    filter(sim >= 0.9) %>%
    mutate(combo = paste(ec_number, ec_number2)) %>%
    unique()
data.table::fwrite(mifaser_pairs[,c("ec_number", "ec_number2", "sim")], "kaiju/outputs/mifaser_pairs.csv", sep = "\t")
system("mcl kaiju/outputs/mifaser_pairs.csv -I 4 --abc -o kaiju/outputs/mcl_mifaser")

# calculate list of modules by EC numbers in mifaser mat
module_list <- module_parts %>% 
    filter(ec_number %in% rownames(mifaser_mat)) %>% 
    group_by(module_id) %>% filter(n() > 1) %>% 
    group_split() %>% 
    lapply(function(.x){res <- list(.x$ec_number); names(res) <- .x$module_id[1]; return(res)}) %>% 
    unlist(recursive = FALSE) 
# read in mcl clustering
mcl_mifaser <- readLines("kaiju/outputs/mcl_mifaser") %>% sapply(function(.x){strsplit(.x, split = "\t")[[1]]}) %>% unname()
mcl_mifaser <- mcl_mifaser[!grepl("ec_number", mcl_mifaser)]
metagenome_mcl <- lapply(mcl_mifaser, function(.x){
    res <- lapply(module_list, function(.y){
        temp <- sum(.y %in% .x)/length(.y)
        names(temp) <- names(.y)
        return(temp)
    }) %>% stack() %>%
    arrange(desc(values))
}) 
metagenome_results <- bind_rows(metagenome_mcl, .id = "cluster") %>% mutate(cluster = as.numeric(as.factor(cluster))) %>% filter(values == 1)
metagenome_kaiju_pred <- metagenome_results  %>% filter(ind %in% unique(unlist(kaiju_modules))); metagenome_kaiju_pred %>% rename("module_id" = ind) %>% left_join(module_classes, by = "module_id") %>% arrange(class)
metagenome_kaiju_unpred <- metagenome_results  %>% filter(!(ind %in% unique(unlist(kaiju_modules)))); metagenome_kaiju_unpred %>% rename("module_id" = ind) %>% left_join(module_classes, by = "module_id") %>% arrange(class)

#####################################################
# set up random comparison
set.seed(123)
random_mifaser_pairs <- as.data.frame(as.matrix(mifaser_dist)) %>% 
    mutate(ec_number = rownames(mifaser_mat)) %>%
    gather(key = "ec_number2", value = "dist", -ec_number) %>%
    filter(ec_number != ec_number2) %>%
    mutate(sim = 1 - dist) %>%
    mutate(sim = sample(sim)) %>%
    filter(sim >= 0.9) %>%
    mutate(combo = paste(ec_number, ec_number2)) %>%
    unique()
data.table::fwrite(random_mifaser_pairs[,c("ec_number", "ec_number2", "sim")], "kaiju/outputs/random_mifaser_pairs.csv", sep = "\t")
system("mcl kaiju/outputs/random_mifaser_pairs.csv -I 4 --abc -o kaiju/outputs/random_mcl_mifaser")

random_mcl_mifaser <- readLines("kaiju/outputs/random_mcl_mifaser") %>% sapply(function(.x){strsplit(.x, split = "\t")[[1]]}) %>% unname()
random_mcl_mifaser <- random_mcl_mifaser[!grepl("ec_number", random_mcl_mifaser)]

random_metagenome_mcl <- lapply(random_mcl_mifaser, function(.x){
    res <- lapply(module_list, function(.y){
        temp <- sum(.y %in% .x)/length(.y)
        names(temp) <- names(.y)
        return(temp)
    }) %>% stack() %>%
    arrange(desc(values))
})
random_metagenome_results <- bind_rows(random_metagenome_mcl, .id = "cluster") %>% mutate(cluster = as.numeric(as.factor(cluster))) %>% filter(values == 1)
random_metagenome_kaiju_pred <- random_metagenome_results  %>% filter(ind %in% unique(unlist(kaiju_modules))); random_metagenome_kaiju_pred
random_metagenome_kaiju_unpred <- random_metagenome_results  %>% filter(!(ind %in% unique(unlist(kaiju_modules)))); random_metagenome_kaiju_unpred


#####################################################

kaiju_pairs <- lapply(kaiju_list, function(.x){
    if(is.null(.x)){return(NULL)}
    temp <- data.table::CJ(.x, .x)
    colnames(temp) <- c("V1", "V2")
    return(temp[!is.na(V1) & !is.na(V2) & V1 != V2])
}) %>% rbindlist() %>% unique() %>% mutate(combo = paste(V1, V2))

mifaser_only <- mifaser_pairs[!(mifaser_pairs$combo %in% kaiju_pairs$combo),]; mifaser_only
# identify unique pairs
unique_pairs <- apply(mifaser_only, 1, function(.x){
    temp <- sort(unlist(strsplit(c(.x[1], .x[2]), split = " ")))
    return(paste(temp, collapse = " "))
}) %>% unique()

preds <- mifaser_pairs %>% 
    filter(combo %in% unique_pairs) %>% 
    rename("EC1" = ec_number, "EC2" = ec_number2) %>% 
    mutate(Jaccard = 1 - dist) %>% 
    select(-c("dist", "combo", "sim"))
data.table::fwrite(preds, "kaiju/outputs/metagenome_preds.csv", sep = "\t")

# mcl the mifaser only preds
system("mcl kaiju/outputs/metagenome_preds.csv --abc -I 4 -o kaiju/outputs/mcl_metagenome_preds")
mifaser_only_mcl <- readLines("kaiju/outputs/mcl_metagenome_preds") %>% sapply(function(.x){strsplit(.x, split = "\t")[[1]]}) %>% unname()
mifaser_only_mcl <- mifaser_only_mcl[!grepl("EC1", mifaser_only_mcl)]

mifaser_only_metagenome_mcl <- lapply(mifaser_only_mcl, function(.x){
    res <- lapply(module_list, function(.y){
        temp <- sum(.y %in% .x)/length(.y)
        names(temp) <- names(.y)
        return(temp)
    }) %>% stack() %>%
    arrange(desc(values))
}) 
mifaser_only_metagenome_results <- bind_rows(mifaser_only_metagenome_mcl, .id = "cluster") %>% mutate(cluster = as.numeric(as.factor(cluster))) %>% filter(values == 1)
mifaser_only_metagenome_results$ind[mifaser_only_metagenome_results$ind %in% unique(unlist(kaiju_modules))] %>% as.character()

mifaser_only_mcl_overlap <- lapply(mifaser_only_mcl , function(.x){
    res <- lapply(module_list, function(.y){
        temp <- sum(.y %in% .x)/length(.y)
        names(temp) <- names(.y)
        return(temp)
    }) %>% stack() %>%
    arrange(desc(values))
}) 

module_metadata <- data.table::fread("data/kegg/module_classes.csv")

mifaser_only_modules <- mifaser_only_mcl_overlap %>% 
    bind_rows() %>% 
    filter(values == 1) %>% 
    rename("module_id" = ind) %>%
    left_join(module_metadata, by = "module_id") %>%
    arrange(class)
mifaser_only_modules$module_id %in% unlist(kaiju_list)