## ---------------------------
## Purpose of script: Compare the plot the distribution of Jaccard Similarity scores between profiles from different methods.
## Author: Henri Chung
## Date Created: 2022-09-21
## Date Modified: 2024-01-04
## ---------------------------

# load libraries
library(data.table)
library(tidyverse)
library(parallel)
# clear working directory
rm(list = ls())

# read in required data
message(Sys.time(), " reading in data")
balanced_fusion_dtms <- readRDS("compare/data/balanced_fusion_dtms.RDS")
balanced_mmseqs_dtms <- readRDS("compare/data/balanced_mmseq_dtms.RDS") 
profile_saturation <- readRDS("compare/outputs/fusion_profile_saturation.RDS")
low_sat_proteins <- names(profile_saturation[profile_saturation <= 0.9])

# read in long form data
fusion_long <- readRDS("compare/outputs/fusion_links_long.RDS")
mmseqs_long_40 <- readRDS("compare/data/mmseqs_links_long_40.RDS")

# generate a random comparison of mmseqs profile saturations
dist_to_long_random <- function(.x, cols){
  temp <- .x[sort(rownames(.x)),]
  temp <- temp[,cols]
  temp <- t(apply(temp, 1, sample))
  temp_dist <- parallelDist::parDist(as.matrix(temp), method = "binary")
  temp_mat <- as.matrix(temp_dist)
  temp_mat[upper.tri(temp_mat, diag = TRUE)] <- NA
  temp_long <- reshape2::melt(temp_mat) %>% 
    filter(!is.na(value)) %>% 
    mutate(sim = 1-value) 
  return(temp_long$sim)
}

# convert to dtm to long 
message(Sys.time(), " ", "converting to long")
if(!file.exists("compare/outputs/sum_fusion_long_random.RDS")){
  fusion_long_random_list <- list()
  for(i in 1:100){
    message(Sys.time(), " ", i)
    fusion_long_random <- lapply(balanced_fusion_dtms, dist_to_long_random, cols = colnames(balanced_fusion_dtms[[1]])) %>% unlist()
    bins <- tabulate(findInterval(fusion_long_random*100, vec=seq(1,100,length.out=100)), nbins=100)
    names(bins) <- as.character(seq(1,100,length.out=100)/100)
    fusion_long_random_list[[i]] <- bins
  }
  sum_fusion_long_random <-  colSums(bind_rows(fusion_long_random_list))
  saveRDS(sum_fusion_long_random, "compare/outputs/sum_fusion_long_random.RDS")
  message(Sys.time(), " finished")
  }else{
  sum_fusion_long_random <- readRDS("compare/outputs/sum_fusion_long_random.RDS")
}

message(Sys.time(), " ", "mmseqs_bin_dtms created")
if(!file.exists("compare/outputs/sum_mmseqs_long_random.RDS")){
  mmseqs_bin_dtms <- lapply(balanced_mmseqs_dtms, FUN = function(.y){
    .y[.y <= 40] <- 0
    .y[.y > 40] <- 1
    return(.y)
  })

  mmseqs_long_random_list <- list()
  for(i in 1:100){
    message(Sys.time(), " ", i)
    mmseqs_long_random <- lapply(mmseqs_bin_dtms, dist_to_long_random, cols = colnames(balanced_fusion_dtms[[1]])) %>% unlist()
    bins <- tabulate(findInterval(mmseqs_long_random*100, vec=seq(1,100,length.out=100)), nbins=100)
    names(bins) <- as.character(seq(1,100,length.out=100)/100)
    mmseqs_long_random_list[[i]] <- bins
  }
  sum_mmseqs_long_random <-  colSums(bind_rows(mmseqs_long_random_list))
  saveRDS(sum_mmseqs_long_random, "compare/outputs/sum_mmseqs_long_random.RDS")
  message(Sys.time(), " finished")
}else{
  sum_mmseqs_long_random <- readRDS("compare/outputs/sum_mmseqs_long_random.RDS")
}

# similarity histogram data
sum_intervals <- function(.x){
  bins <- tabulate(findInterval(.x*100, vec=seq(1,100,length.out=100)), nbins=100)
  names(bins) <- as.character(seq(1,100,length.out=100)/100)
  return(bins)
}

sim_data <- list(sum_intervals(fusion_long$sim), sum_intervals(mmseqs_long_40[[1]]$sim), sum_fusion_long_random, sum_mmseqs_long_random)
names(sim_data) <- c("fusion", "mmseqs", "fusion_random", "mmseqs_random")
sim_data_all <- lapply(sim_data, function(bins){
  prop <- bins / sum(bins)
  names(prop) <- seq(1,100,length.out=100)
  res <- stack(prop) %>%
    rename(Jaccard_Similarity = ind, Frequency = values) 
})

# plot
p1_data <- bind_rows(sim_data_all, .id = "method") 
p1 <- p1_data %>% 
  mutate(method =  factor(method, levels =  c("fusion", "fusion_random", "mmseqs",  "mmseqs_random"))) %>%
  mutate(method_label = case_when(method == "fusion" ~ "Fusion",
                                  method == "fusion_random" ~ "Fusion (Random)",
                                  method == "mmseqs" ~ "MMseqs2 - 40", 
                                  method == "mmseqs_random" ~ "MMseqs2 (Random)")) %>%
  mutate(class_label = ifelse(grepl("Fusion", method_label), "Fusion", "MMseqs2")) %>%
  mutate(class_label = factor(class_label, levels = c("Fusion", "MMseqs2"))) %>%
  mutate(Jaccard_Similarity = as.numeric(Jaccard_Similarity)) %>% 
  ggplot(aes(x = Jaccard_Similarity/100, y = Frequency, fill = method_label)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.5) + 
  facet_wrap(~class_label, scales = "free") +
  theme_bw() + 
  labs(x = "Jaccard Similarity", y = "Frequency", title = "Jaccard Similarity of Identified Proteins", fill = "Method") 
ggsave("compare/outputs/jaccard_similarity_histogram.png", p1, width = 10, height = 5)


summary_stats_from_freq_table <- function(data_frame, value = "ind", frequency = "values") {
  data_frame[[value]] <- as.numeric(as.character(data_frame[[value]]))
  data_frame[[frequency]] <- as.numeric(data_frame[[frequency]])
  # Compute mean
  mean_val <- sum(data_frame[[value]] * data_frame[[frequency]]) / sum(data_frame[[frequency]])
  
  # Helper function to find the value at a given percentile
  find_percentile <- function(p) {
    cumulative_freq <- cumsum(data_frame[[frequency]])
    total_freq <- tail(cumulative_freq, n=1)
    target_freq <- p * total_freq
    idx <- which(cumulative_freq >= target_freq)[1]
    return(data_frame[[value]][idx])
  }
  
  # Compute summary statistics
  stats <- list(
    min = min(data_frame[[value]]),
    first_quartile = find_percentile(0.25),
    median = find_percentile(0.5),
    mean = mean_val,
    third_quartile = find_percentile(0.75),
    max = max(data_frame[[value]])
  )
  
  return(stats)
}

sim_boxplot_data <- lapply(sim_data, function(.x){stack(summary_stats_from_freq_table(stack(.x)))})
names(sim_boxplot_data) <- c("fusion", "mmseqs", "fusion_random", "mmseqs_random")
p2_data <- bind_rows(sim_boxplot_data, .id = "method") %>% pivot_wider(names_from = "ind", values_from = "values")

p2 <- p2_data %>%
  mutate(method_label = case_when(method == "fusion" ~ "Fusion",
                                method == "fusion_random" ~ "Fusion (Random)",
                                method == "mmseqs" ~ "MMseqs2 - 40", 
                                method == "mmseqs_random" ~ "MMseqs2 (Random)")) %>%
  ggplot(aes(x = method_label, ymin=min, lower=first_quartile, middle = median, upper = third_quartile, ymax = max, fill = method_label)) + 
  geom_boxplot(stat="identity") +
  theme_bw() +
  labs(x = "Method", y = "Jaccard Similarity", fill = "Method") 
ggsave("compare/outputs/jaccard_similarity_boxplot.png", p2)
