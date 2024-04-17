## ---------------------------
## Purpose of script: Helper functions to build fusion profiles.
## Author: Henri Chung
## Date Created: 2021-01-12
## Date Modified: 2024-01-04
## ---------------------------

# function to return cooccurance matrix
return_dtm <- function(fusion_rows, f_data){
    subset_data <- (f_data[fusion_lvl_1 %in% fusion_rows])
    cast_data <- subset_data[,.(val=.N),.(assembly_accession2, fusion_lvl_1)]
    DTM <- dcast(cast_data, fusion_lvl_1 ~ assembly_accession2, value.var = "val", fill = 0)[,fusion_lvl_1 := as.character(fusion_lvl_1)] %>%
        column_to_rownames("fusion_lvl_1")
    all_assemblies <-  unique(f_data$assembly_accession2)
    empty_assemblies <- all_assemblies[!(all_assemblies %in% unique(subset_data$assembly_accession2))]
    empty_df <- as.data.frame(matrix(0, nrow = nrow(DTM), ncol = length(empty_assemblies)))
    colnames(empty_df) <- empty_assemblies
    DTM <- cbind(empty_df, DTM)
    DTM_mat <- as.matrix(DTM)
    return(DTM_mat)
}

# convert cooccurance matrix to long format
dist_to_long <- function(.x, cols){
  temp <- .x[sort(rownames(.x)),]
  temp <- temp[,cols]
  temp_dist <- parallelDist::parDist(as.matrix(temp), method = "binary")
  temp_mat <- as.matrix(temp_dist)
  temp_mat[upper.tri(temp_mat, diag = TRUE)] <- NA
  temp_long <- reshape2::melt(temp_mat) %>% 
    filter(!is.na(value)) %>% 
    mutate(sim = 1-value) %>%
    mutate(link = paste(Var1, Var2, sep = "-")) %>%
    as.data.table()
  return(temp_long)
}

# generate predictions
generate_preds <- function(cols, dtms, preds_only = FALSE){
    fusion_long <- lapply(dtms, dist_to_long, cols = cols)
    fusion_long_dt <- rbindlist(fusion_long, idcol = "assembly_accession2", use.names = TRUE)
    res <- fusion_long_dt[order(link)]
    if(preds_only != TRUE){
        return(res)
    }else{
        return(res$sim)
    }
}

# Compare precision
precision_recall <- function(df, truth_col, estimate_col) {
  thresholds <- seq(0, 1, by = 0.05)
  df[[truth_col]] <- as.logical(df[[truth_col]])
  pr <- sapply(thresholds, function(threshold) {
    predicted_positive <- df[[estimate_col]] >= threshold
    true_positive <- sum(predicted_positive & df[[truth_col]])
    false_positive <- sum(predicted_positive & !df[[truth_col]])
    false_negative <- sum(!predicted_positive & df[[truth_col]])
    precision <- true_positive / (true_positive + false_positive)
    recall <- true_positive / (true_positive + false_negative)
    c(precision = precision, recall = recall)
  })  
  data.frame(threshold = thresholds, precision = pr["precision",], recall = pr["recall",])
}


# compare profiles of different lengths
compare_variable_lengths <- function(preds, saturation, long, k_vals = c(0.5, 0.6, 0.7, 0.8, 0.9, 1)){
  long <- long[order(link)]
  results_list <- list(); l = 1
  names_list <- list()
  lengths <- names(preds)
  for(i in 1:length(lengths)){
    for(j in 1:length(preds[[i]])){
      long$sim <- preds[[i]][[j]][1:nrow(long)]
      profile_saturation <- saturation[[i]][[j]]
      for (k in k_vals) {
        low_sat_proteins <- profile_saturation[profile_saturation <= k] %>% names()
        temp <- copy(long)
        temp[, sim2 := ifelse(!(Var1 %in% low_sat_proteins) | !(Var2 %in% low_sat_proteins), 0, sim)]
        results_list[[l]] <- precision_recall(temp, "truth", "sim2")
        names_list[[l]] <- paste(i, j, k)
        message(Sys.time(), " ", l)
        l = l+1
        rm(temp); gc()
      }
    }
  }
  names(results_list) <- names_list[1:length(results_list)]
  return(results_list)
}