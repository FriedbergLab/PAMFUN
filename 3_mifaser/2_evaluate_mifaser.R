## ---------------------------
## Purpose of script: Evaluate mifaser results
## Author: Henri Chung
## Date Created: 2022-09-21
## Date Modified: 2024-01-04
## ---------------------------
# evaluate results
library(tidyverse)
library(data.table)
rm(list = ls())

all_sim_reads <- readRDS("mifaser/data/test_reads_metadata.RDS")
mifaser_Res <- read_tsv("mifaser/outputs/fusion_database_cdhit.out/ec_count.tsv", col_names = FALSE)
collapse_columns <- function(df) {
  if (ncol(df) > 2) {
    first_col <- df[, 1, drop = FALSE]
    other_cols <- df[, -1]
    collapsed_col <- apply(other_cols, 1, paste, collapse = "\t")
    new_df <- cbind(first_col, collapsed_col)
    return(new_df)
  } else {
    return(df)
  }
}

evaluate_preds <- function(.x, .y){
    preds <- .x %>% left_join(.y, by = c("name")) %>% 
        select(seq_length, name, fusions, data) %>%
        mutate(data = purrr::map(data, ~ifelse(is.null(.x), "0", unique(.x$ec)))) %>%
        mutate(truth = purrr::map2(fusions, data, ~ sum(.x %in% .y) > 0)) %>%
        unnest(c(data, truth)) %>% rename(ec = "data") %>%
        select(name, fusions, ec, truth, seq_length)
    # need to create negative, or proteins not connected to any fusion
    P = filter(preds, ec != "0")
    TP = sum(P$truth)
    FP = nrow(P) - TP
    N = filter(preds, ec == "0")
    TN = sum(N$truth)
    FN = nrow(N) - TN
    Precision = TP / (TP + FP); Precision # 0.945
    Recall = TP / (TP + FN); Recall #0.571
    Sensitivity = TP / (TP + FN); Sensitivity #0.571
    Specificity = TN / (TN + FP); Specificity # 0.829
    F1 = 2 * Precision * Recall / (Precision + Recall); F1 # 0.712
    Accuracy = (TP + TN) / (TP + TN + FP + FN); Accuracy # 0.613
    res_df <- data.frame(Precision, Recall, Sensitivity, Specificity, F1, Accuracy)
    return(res_df)
}

res_parse <- function(.x){
    temp <- read_delim(.x, col_names = FALSE, delim = "\t") %>% collapse_columns()
    #temp <- data.table::fread(.x, header = FALSE, sep = "\t", fill = TRUE) %>% collapse_columns()
    colnames(temp) <- c("name", "ec")
    res <- temp %>%
        #mutate(name =  gsub("\\|.*", "", name)) %>%
        separate_rows(ec, sep = "\t") %>%
        filter(!is.na(ec)) %>%
        mutate(ec = gsub("^0+", "",gsub("\\.", "", ec)))  %>%
        group_by(name) %>% nest(data = "ec")
    return(res)
}

#mifaser_output_files <- list.files(pattern = "*ec_count.tsv", path = "mifaser/outputs/db_samples_out", full.names = TRUE, recursive = TRUE)[1]
#mifaser_output_files <- list.files(pattern = "*ec_count.tsv", path = "mifaser/outputs/fusion_database.out", full.names = TRUE, recursive = TRUE)
mifaser_output_files <- list.files(pattern = "*ec_count.tsv", path = "mifaser/outputs/fusion_database_cdhit.out", full.names = TRUE, recursive = TRUE)


res_list <- list()
for(i in 1:length(mifaser_output_files)){
    res <- res_parse(mifaser_output_files[i])
    res_eval <- evaluate_preds(all_sim_reads, res)
    res_list[[i]] <- res_eval
}
names(res_list) <- mifaser_output_files


res_eval
#   Precision    Recall Sensitivity Specificity        F1  Accuracy
# 1 0.9446623 0.5254617   0.5254617   0.7076516 0.6752953 0.5428171