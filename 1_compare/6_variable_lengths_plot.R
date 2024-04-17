## ---------------------------
## Purpose of script: Compare and plot the precision recall performance of fusion and mmseqs at variable lengths and saturations
## Author: Henri Chung
## Date Created: 2022-06-26
## Date Modified: 2024-01-04
## ---------------------------


library(data.table)
library(tidyverse)
library(zoo)
# clear working directory
rm(list = ls())

# define profile sizes
profile_sizes <- c(50, 100, 300, 500, 700, 900, 1100, 1393) %>% as.character()

# read in results
fusion_vl_files <- list.files("compare/outputs", pattern = "fusion_vl_results_[0-9]*", full.names = TRUE)
fusion_vl_results <- fusion_vl_files %>%
  lapply(readRDS) %>%
  setNames(., fusion_vl_files) %>%
  lapply(rbindlist, idcol = "file") %>%
  rbindlist(idcol = "profile_size") %>%
  separate(file, into = c("empty", "replicate", "filter"), sep = " ") %>%
  mutate(profile_size = str_extract(profile_size, "[0-9]+"))

# write function to parse mmseq results of variable lengths
mmseqs_parse <- function(.x){
  temp_files  <- list.files("compare/outputs", pattern = .x, full.names = TRUE) 
  temp <- temp_files %>%
    lapply(readRDS) %>%
    setNames(temp_files) %>%
    lapply(rbindlist, idcol = "file") %>%
    rbindlist(idcol = "profile_size") %>%
    separate(file, into = c("empty", "replicate", "filter"), sep = " ") %>%
    mutate(profile_size = gsub("_vl", "", str_extract(profile_size, "[0-9]+_vl")))
  return(temp)
}

# apply parsing function to mmseqs results of different bit score thresholds
mmseqs_40_vl_results <- mmseqs_parse("mmseqs_40_[0-9]+_vl_results")
mmseqs_60_vl_results <- mmseqs_parse("mmseqs_60_[0-9]+_vl_results")
mmseqs_80_vl_results <- mmseqs_parse("mmseqs_80_[0-9]+_vl_results")
mmseqs_100_vl_results <- mmseqs_parse("mmseqs_100_[0-9]+_vl_results")

# combine results
vl_results <- list(fusion_vl_results, mmseqs_40_vl_results, mmseqs_60_vl_results, mmseqs_80_vl_results, mmseqs_100_vl_results)
names(vl_results) <- c("fusion", "mmseqs_40", "mmseqs_60", "mmseqs_80", "mmseqs_100")

# identify threshold for max f1
threshold_test <- bind_rows(vl_results, .id = "method")  %>%
  mutate(profile_size = as.character(profile_size)) %>%
  filter((profile_size == "1393" & replicate == "1") | (profile_size != "1393")) %>%
  mutate(method = factor(method, levels = c("fusion", "mmseqs_40", "mmseqs_60", "mmseqs_80", "mmseqs_100"))) %>%
  filter(!is.nan(precision) & thresholds > 0) %>%
  mutate(f1 = 2 * precision * recall / (precision + recall)) %>%
  group_by(method, profile_size,filter, replicate) %>%
  filter(method %in% c("fusion", "mmseqs_40") & replicate == 1 & filter == 0.9 & profile_size == 900)

vl_df <- bind_rows(vl_results, .id = "method")  %>%
  mutate(profile_size = as.character(profile_size)) %>%
  mutate(method = factor(method, levels = c("fusion", "mmseqs_40", "mmseqs_60", "mmseqs_80", "mmseqs_100"))) %>%
  filter(!is.nan(precision) & thresholds > 0) %>%
  mutate(f1 = 2 * precision * recall / (precision + recall)) %>%
  group_by(method, profile_size,filter, replicate) %>%
  filter(f1 == max(f1)) %>%
  group_by(method, profile_size,filter) %>%
  summarize(mean_precision = mean(precision), mean_recall = mean(recall), mean_f1 = mean(f1),
            sd_precision = sd(precision), sd_recall = sd(recall), sd_f1 = sd(f1)) %>%
  mutate(profile_size = factor(profile_size, levels = profile_sizes)) %>%
  mutate(filter = factor(filter, levels = c("0.5", "0.6", "0.7", "0.8", "0.9", "1"))) %>%
  mutate(method =  factor(method, levels = c("fusion", "fusion2", "mmseqs_40", "mmseqs_60", "mmseqs_80", "mmseqs_100"))) 
head(vl_df)
write.csv(vl_df, "compare/outputs/variable_lengths.csv")


# generate plots
dodge_distance = 0.3
p2a <- vl_df %>%
  mutate(filter_label = paste("Max Profile Saturation:", filter)) %>% 
  mutate(profile_label = paste("Profile Length:", profile_size)) %>%
  mutate(profile_label = factor(profile_label, levels = unique(paste("Profile Length:", levels(profile_size))))) %>% 
  ggplot(aes(x = profile_size, y = mean_precision, group = method, color = method, shape = method)) +
  geom_point(position = position_dodge(width = dodge_distance)) + geom_line(position = position_dodge(width = dodge_distance))  + 
  geom_errorbar(aes(ymin=mean_precision-sd_precision, ymax=mean_precision+sd_precision), width=.2, linetype = 1, position = position_dodge(width = dodge_distance)) +
  facet_wrap(~filter_label) +
  theme_bw() + 
  scale_color_manual(
    name = "Method", 
    values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
    labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  labs(x = "Profile Length", y = "Mean Precision") +
  guides(color = guide_legend(override.aes = list(shape = unique(vl_df$method))), shape = "none")

p2b <- vl_df %>%
  mutate(filter_label = paste("Max Profile Saturation:", filter)) %>% 
  mutate(profile_label = paste("Profile Length:", profile_size)) %>%
  mutate(profile_label = factor(profile_label, levels = unique(paste("Profile Length:", levels(profile_size))))) %>% 
  ggplot(aes(x = filter, y = mean_precision, group = method, color = method, shape = method)) +
  geom_point(position = position_dodge(width = dodge_distance)) + geom_line(position = position_dodge(width = dodge_distance))  + 
  geom_errorbar(aes(ymin=mean_precision-sd_precision, ymax=mean_precision+sd_precision), width=.1, linetype = 1, position = position_dodge(width = dodge_distance)) +
  facet_wrap(~profile_label, nrow = 2) +
  theme_bw() + 
  scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  scale_shape_manual(values = c("fusion" = 1, "mmseqs_40" = 3, "mmseqs_60" = 4, "mmseqs_80" = 5, "mmseqs_100" = 6)) +
  labs(x = "Max Profile Saturation", y = "Mean Precision")  +
  guides(color = guide_legend(override.aes = list(shape = unique(vl_df$method))), shape = "none")
pdf("compare/outputs/fusion_mmseqs_variable_lengths_precision.pdf", height = 8, width = 10)
print(p2a)
print(p2b)
dev.off()

png("compare/outputs/fusion_mmseqs_variable_lengths_precision1.png", width = 10, height = 6, units = "in", res = 300)
print(p2a)
dev.off()
png("compare/outputs/fusion_mmseqs_variable_lengths_precision2.png", width = 10, height = 6, units = "in", res = 300)
print(p2b)
dev.off()


p3a <- vl_df %>%
  mutate(filter_label = paste("Max Profile Saturation:", filter)) %>% 
  mutate(profile_label = paste("Profile Length:", profile_size)) %>%
  mutate(profile_label = factor(profile_label, levels = unique(paste("Profile Length:", levels(profile_size))))) %>% 
  ggplot(aes(x = profile_size, y = mean_recall, group = method, color = method, shape = method)) +
  geom_point(position = position_dodge(width = dodge_distance)) + geom_line(position = position_dodge(width = dodge_distance))  + 
  geom_errorbar(aes(ymin=mean_recall-sd_recall, ymax=mean_recall+sd_recall), width=.2, linetype = 1, position = position_dodge(width = dodge_distance)) +
  facet_wrap(~filter_label) +
  theme_bw() +
  scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  scale_shape_manual(values = c("fusion" = 1, "mmseqs_40" = 3, "mmseqs_60" = 4, "mmseqs_80" = 5, "mmseqs_100" = 6)) +
  labs(x = "Profile Length", y = "Mean Recall") +
  guides(color = guide_legend(override.aes = list(shape = unique(vl_df$method))), shape = "none")
p3b <- vl_df %>%
  mutate(filter_label = paste("Max Profile Saturation:", filter)) %>% 
  mutate(profile_label = paste("Profile Length:", profile_size)) %>%
  mutate(profile_label = factor(profile_label, levels = unique(paste("Profile Length:", levels(profile_size))))) %>% 
  ggplot(aes(x = filter, y = mean_recall, group = method, color = method, shape = method)) +
  geom_point(position = position_dodge(width = dodge_distance)) + geom_line(position = position_dodge(width = dodge_distance))  + 
  geom_errorbar(aes(ymin=mean_recall-sd_recall, ymax=mean_recall+sd_recall), width=.1, linetype = 1, position = position_dodge(width = dodge_distance)) +
  facet_wrap(~profile_label, nrow = 2) +
  theme_bw() + 
  scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  scale_shape_manual(values = c("fusion" = 1, "mmseqs_40" = 3, "mmseqs_60" = 4, "mmseqs_80" = 5, "mmseqs_100" = 6)) +
  labs(x = "Max Profile Saturation", y = "Mean Recall") +
  guides(color = guide_legend(override.aes = list(shape = unique(vl_df$method))), shape = "none")

pdf("compare/outputs/fusion_mmseqs_variable_lengths_recall.pdf", height = 8, width = 10)
print(p3a)
print(p3b)
dev.off()
png("compare/outputs/fusion_mmseqs_variable_lengths_recall1.png", width = 10, height = 6, units = "in", res = 300)
print(p3a)
dev.off()
png("compare/outputs/fusion_mmseqs_variable_lengths_recall2.png", width = 10, height = 6, units = "in", res = 300)
print(p3b)
dev.off()


p4a <- vl_df %>%
  mutate(filter_label = paste("Max Profile Saturation:", filter)) %>% 
  mutate(profile_label = paste("Profile Length:", profile_size)) %>%
  mutate(profile_label = factor(profile_label, levels = unique(paste("Profile Length:", levels(profile_size))))) %>% 
  ggplot(aes(x = profile_size, y = mean_f1, group = method, color = method, shape = method)) +
  geom_point(position = position_dodge(width = dodge_distance)) + geom_line(position = position_dodge(width = dodge_distance))  + 
  geom_errorbar(aes(ymin=mean_f1-sd_f1, ymax=mean_f1+sd_f1), width=.2, linetype = 1, position = position_dodge(width = dodge_distance)) +
  facet_wrap(~filter_label) +
  theme_bw() +
  scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  scale_shape_manual(values = c("fusion" = 1, "mmseqs_40" = 3, "mmseqs_60" = 4, "mmseqs_80" = 5, "mmseqs_100" = 6)) +
  labs(x = "Profile Length", y = "Mean F1") +
  guides(color = guide_legend(override.aes = list(shape = unique(vl_df$method))), shape = FALSE)
p4b <- vl_df %>%
  mutate(filter_label = paste("Max Profile Saturation:", filter)) %>% 
  mutate(profile_label = paste("Profile Length:", profile_size)) %>%
  mutate(profile_label = factor(profile_label, levels = unique(paste("Profile Length:", levels(profile_size))))) %>% 
  ggplot(aes(x = filter, y = mean_f1, group = method,  color = method, shape = method)) +
  geom_point(position = position_dodge(width = dodge_distance)) + geom_line(position = position_dodge(width = dodge_distance))  + 
  geom_errorbar(aes(ymin=mean_f1-sd_f1, ymax=mean_f1+sd_f1), width=.1, linetype = 1, position = position_dodge(width = dodge_distance)) +
  facet_wrap(~profile_label, nrow = 2) +
  theme_bw() + 
  scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  scale_shape_manual(values = c("fusion" = 1, "mmseqs_40" = 3, "mmseqs_60" = 4, "mmseqs_80" = 5, "mmseqs_100" = 6)) +
  labs(x = "Max Profile Saturation", y = "Mean F1") +
  guides(color = guide_legend(override.aes = list(shape = unique(vl_df$method))), shape = FALSE)

pdf("compare/outputs/fusion_mmseqs_variable_lengths_f1.pdf", height = 8, width = 10)
print(p4a)
print(p4b)
dev.off()
png("compare/outputs/fusion_mmseqs_variable_lengths_f11.png", width = 10, height = 6, units = "in", res = 300)
print(p4a)
dev.off()
png("compare/outputs/fusion_mmseqs_variable_lengths_f12.png", width = 10, height = 6, units = "in", res = 300)
print(p4b)
dev.off()

p5a <- vl_df %>%
  filter(profile_size == "900") %>%
  mutate(filter_label = paste("Max Saturation:", filter)) %>% 
  mutate(profile_label = paste("Profile Length:", profile_size)) %>%
  mutate(profile_label = factor(profile_label, levels = unique(paste("Profile Length:", levels(profile_size))))) %>% 
  ggplot(aes(x = filter, y = mean_f1, group = method,  color = method, shape = method)) +
  geom_point(position = position_dodge(width = dodge_distance)) + geom_line(position = position_dodge(width = dodge_distance))  + 
  geom_errorbar(aes(ymin=mean_f1-sd_f1, ymax=mean_f1+sd_f1), width=.1, linetype = 1, position = position_dodge(width = dodge_distance)) +
  facet_wrap(~profile_label, nrow = 2) +
  theme_bw() + 
  ylim(0, 0.22) +
  scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  scale_shape_manual(values = c("fusion" = 1, "mmseqs_40" = 3, "mmseqs_60" = 4, "mmseqs_80" = 5, "mmseqs_100" = 6)) +
  labs(x = "Max Profile Saturation", y = "Mean F1") +
  guides(color = guide_legend(override.aes = list(shape = unique(vl_df$method))), shape = FALSE) +
  theme(axis.title.x = element_text(size = 14, vjust = -0.5), 
        axis.title.y = element_text(size = 14, vjust = 3),
        axis.text.x = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.4,0.2), "cm"))

p5b <- vl_df %>%
  filter(filter == "0.9") %>%
  mutate(filter_label = paste("Max Saturation:", filter)) %>% 
  mutate(profile_label = paste("Profile Length:", profile_size)) %>%
  mutate(profile_label = factor(profile_label, levels = unique(paste("Profile Length:", levels(profile_size))))) %>% 
  ggplot(aes(x = profile_size, y = mean_f1, group = method, color = method, shape = method)) +
  geom_point(position = position_dodge(width = dodge_distance)) + geom_line(position = position_dodge(width = dodge_distance))  + 
  geom_errorbar(aes(ymin=mean_f1-sd_f1, ymax=mean_f1+sd_f1), width=.2, linetype = 1, position = position_dodge(width = dodge_distance)) +
  facet_wrap(~filter_label) +
  theme_bw() +
  ylim(0, 0.22) +
  scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  scale_shape_manual(values = c("fusion" = 1, "mmseqs_40" = 3, "mmseqs_60" = 4, "mmseqs_80" = 5, "mmseqs_100" = 6)) +
  labs(x = "Profile Length", y = "Mean F1") +
  guides(color = guide_legend(override.aes = list(shape = unique(vl_df$method))), shape = FALSE) +
  theme(axis.title.x = element_text(size = 14, vjust = -0.5), 
        axis.title.y = element_text(size = 14, vjust = 3),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

p5c <- vl_df %>%
  filter(profile_size == "900") %>%
  mutate(filter_label = paste("Max Saturation:", filter)) %>% 
  mutate(profile_label = paste("Profile Length:", profile_size)) %>%
  mutate(profile_label = factor(profile_label, levels = unique(paste("Profile Length:", levels(profile_size))))) %>% 
  ggplot(aes(x = filter, y = mean_precision, group = method,  color = method, shape = method)) +
  geom_point(position = position_dodge(width = dodge_distance)) + geom_line(position = position_dodge(width = dodge_distance))  + 
  geom_errorbar(aes(ymin=mean_precision-sd_precision, ymax=mean_precision+sd_precision), width=.1, linetype = 1, position = position_dodge(width = dodge_distance)) +
  facet_wrap(~profile_label, nrow = 2) +
  theme_bw() + 
  ylim(0, 0.8) +
  scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  scale_shape_manual(values = c("fusion" = 1, "mmseqs_40" = 3, "mmseqs_60" = 4, "mmseqs_80" = 5, "mmseqs_100" = 6)) +
  labs(x = "Max Profile Saturation", y = "Mean Precision") +
  guides(color = guide_legend(override.aes = list(shape = unique(vl_df$method))), shape = FALSE) +
  theme(axis.title.x = element_text(size = 14, vjust = -0.5), 
        axis.title.y = element_text(size = 14, vjust = 3),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.4,0.2), "cm"))

p5d <- vl_df %>%
  filter(filter == "0.9") %>%
  mutate(filter_label = paste("Max Saturation:", filter)) %>% 
  mutate(profile_label = paste("Profile Length:", profile_size)) %>%
  mutate(profile_label = factor(profile_label, levels = unique(paste("Profile Length:", levels(profile_size))))) %>% 
  ggplot(aes(x = profile_size, y = mean_precision, group = method,  color = method, shape = method)) +
  geom_point(position = position_dodge(width = dodge_distance)) + geom_line(position = position_dodge(width = dodge_distance))  + 
  geom_errorbar(aes(ymin=mean_precision-sd_precision, ymax=mean_precision+sd_precision), width=.1, linetype = 1, position = position_dodge(width = dodge_distance)) +
  facet_wrap(~filter_label) +
  theme_bw() +
  ylim(0, 0.8) +
  scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  scale_shape_manual(values = c("fusion" = 1, "mmseqs_40" = 3, "mmseqs_60" = 4, "mmseqs_80" = 5, "mmseqs_100" = 6)) +
  labs(x = "Profile Length", y = "Mean Precision") +
  guides(color = guide_legend(override.aes = list(shape = unique(vl_df$method))), shape = FALSE) +
  theme(axis.title.x = element_text(size = 14, vjust = -0.5), 
        axis.title.y = element_text(size = 14, vjust = 3),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

p5e <- vl_df %>%
  filter(profile_size == "900") %>%
  mutate(filter_label = paste("Max Saturation:", filter)) %>% 
  mutate(profile_label = paste("Profile Length:", profile_size)) %>%
  mutate(profile_label = factor(profile_label, levels = unique(paste("Profile Length:", levels(profile_size))))) %>% 
  ggplot(aes(x = filter, y = mean_recall, group = method,  color = method, shape = method)) +
  geom_point(position = position_dodge(width = dodge_distance)) + geom_line(position = position_dodge(width = dodge_distance))  + 
  geom_errorbar(aes(ymin=mean_recall-sd_recall, ymax=mean_recall+sd_recall), width=.1, linetype = 1, position = position_dodge(width = dodge_distance)) +
  facet_wrap(~profile_label, nrow = 2) +
  theme_bw() + 
  ylim(0, 0.17) +
  scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  scale_shape_manual(values = c("fusion" = 1, "mmseqs_40" = 3, "mmseqs_60" = 4, "mmseqs_80" = 5, "mmseqs_100" = 6)) +
  labs(x = "Max Profile Saturation", y = "Mean Recall") +
  guides(color = guide_legend(override.aes = list(shape = unique(vl_df$method))), shape = FALSE) +
  theme(axis.title.x = element_text(size = 14, vjust = -0.5), 
        axis.title.y = element_text(size = 14, vjust = 3),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.4,0.2), "cm"))

p5f <- vl_df %>%
  filter(filter == "0.9") %>%
  mutate(filter_label = paste("Max Saturation:", filter)) %>% 
  mutate(profile_label = paste("Profile Length:", profile_size)) %>%
  mutate(profile_label = factor(profile_label, levels = unique(paste("Profile Length:", levels(profile_size))))) %>% 
  ggplot(aes(x = profile_size, y = mean_recall, group = method,  color = method, shape = method)) +
  geom_point(position = position_dodge(width = dodge_distance)) + geom_line(position = position_dodge(width = dodge_distance))  + 
  geom_errorbar(aes(ymin=mean_recall-sd_recall, ymax=mean_recall+sd_recall), width=.1, linetype = 1, position = position_dodge(width = dodge_distance)) +
  facet_wrap(~filter_label) +
  theme_bw() +
  ylim(0, 0.17) +
  scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  scale_shape_manual(values = c("fusion" = 1, "mmseqs_40" = 3, "mmseqs_60" = 4, "mmseqs_80" = 5, "mmseqs_100" = 6)) +
  labs(x = "Profile Length", y = "Mean Recall") +
  guides(color = guide_legend(override.aes = list(shape = unique(vl_df$method))), shape = FALSE)+
  theme(axis.title.x = element_text(size = 14, vjust = -0.5), 
        axis.title.y = element_text(size = 14, vjust = 3),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        #legend.key.size = unit(2, 'cm'), 
        legend.title = element_text(size=14),
        legend.text = element_text(size=14))


library(ggpubr)
shared_legend <- get_legend(p5f)
p5 <- ggarrange(p5a, p5c, p5e, p5b, p5d, p5f, ncol=3, nrow=2, common.legend = TRUE, legend="right", labels = c("A", "B", "C", "D", "E", "F"), legend.grob = shared_legend)
pdf("compare/outputs/fusion_mmseqs_variable_lengths.pdf", height = 5, width = 10)
print(p5)
dev.off()


png("compare/outputs/fusion_mmseqs_variable_lengths.png", width = 12, height = 6, units = "in", res = 300)
print(p5)
dev.off()

p5_vertical <- ggarrange(p5a, p5b, p5c, p5d, p5e, p5f, ncol = 2, nrow = 3, common.legend = TRUE, legend="bottom", labels = c("A", "C", "E", "B", "D", "F"))
png("compare/outputs/fusion_mmseqs_variable_lengths_vertical.png", width = 6, height = 12, units = "in", res = 300)
print(p5_vertical)
dev.off()
# prc curves
prc_plot <- vl_results %>%
  rbindlist(idcol = "method") %>%
  mutate(profile_size = factor(profile_size, levels = profile_sizes)) %>%
  mutate(filter = factor(filter, levels = c("0.5", "0.6", "0.7", "0.8", "0.9", "1"))) %>%
  mutate(method =  factor(method, levels = c("fusion", "fusion2", "mmseqs_40", "mmseqs_60", "mmseqs_80", "mmseqs_100"))) %>%
  mutate(filter_label = paste("Max Profile Saturation:", filter)) %>% 
  mutate(profile_label = paste("Profile Length:", profile_size)) %>%
  mutate(profile_label = factor(profile_label, levels = unique(paste("Profile Length:", levels(profile_size))))) %>% 
  mutate(replicate = paste(replicate, method)) %>%
  ggplot(aes(x = recall, y = precision, color = method, group = replicate)) + geom_line() + facet_grid(profile_size~filter) + theme_bw() +
  scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  xlim(0, 0.3)
png("compare/outputs/prc_vl_plot.png", width = 10, height = 10)
print(prc_plot)
dev.off()

roc_plot <- vl_results %>%
  rbindlist(idcol = "method") %>%
  mutate(profile_size = factor(profile_size, levels = profile_sizes)) %>%
  mutate(filter = factor(filter, levels = c("0.5", "0.6", "0.7", "0.8", "0.9", "1"))) %>%
  mutate(method =  factor(method, levels = c("fusion", "fusion2", "mmseqs_40", "mmseqs_60", "mmseqs_80", "mmseqs_100"))) %>%
  mutate(filter_label = paste("Max Profile Saturation:", filter)) %>% 
  mutate(profile_label = paste("Profile Length:", profile_size)) %>%
  mutate(profile_label = factor(profile_label, levels = unique(paste("Profile Length:", levels(profile_size))))) %>% 
  mutate(replicate = paste(replicate, method)) %>%
  ggplot(aes(x = 1 - specificity, y = recall, color = method, group = replicate)) + geom_line() + facet_grid(profile_size~filter) + theme_bw() +
    scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  xlab("1 - Specificity") + ylab("Sensitivity")
png("compare/outputs/roc_vl_plot.png", width = 10, height = 10)
print(roc_plot)
dev.off()

best_data <- prc_plot$data %>% filter(profile_size == 900 & filter == 0.9)

best_prc_plot <- best_data %>%
  ggplot(aes(x = recall, y = precision, color = method, group = replicate)) + geom_line() + facet_grid(profile_size~filter) + 
  theme_bw() +
  scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  xlim(0, 0.3)
ggsave("compare/outputs/best_prc_plot.png", best_prc_plot, width = 12, height = 8, units = "in", dpi = 300)

best_prc_plot2 <- best_data %>%
  ggplot(aes(x = recall, y = precision, color = method, group = replicate)) + geom_line() + facet_grid(profile_size~filter) + 
  theme_bw() +
  scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) 
ggsave("compare/outputs/best_prc_plot2.png", best_prc_plot2, width = 12, height = 8, units = "in", dpi = 300)

best_prc_plot3 <- best_data %>%
  ggplot(aes(x = recall, y = precision, color = method, group = replicate)) + geom_line() + facet_grid(profile_size~filter) + 
  theme_bw() +
  scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) 
ggsave("compare/outputs/best_prc_plot3.png", best_prc_plot3, width = 12, height = 8, units = "in", dpi = 300)

error_margin <- function(mean, sd, n) {
  z_value <- 1.645
  error_margin <- z_value * (sd / sqrt(n))
  return(error_margin)
}

test_data <- vl_results %>%
  rbindlist(idcol = "method") %>%
  rowwise() %>%
  mutate(recall_ind = signif(recall, 2)) %>%
  group_by(method, profile_size, filter, recall_ind) %>% 
  summarize(mean_precision = mean(precision), sd_precision = sd(precision), recall_ind = mean(recall_ind)) %>%
  rowwise() %>%
  mutate(error_margin = error_margin(mean_precision, sd_precision, n = 10)) 


best_prc_plot <- test_data %>%
  filter(profile_size == 900 & filter == 0.9) %>%
    mutate(method =  factor(method, levels = c("fusion", "fusion2", "mmseqs_40", "mmseqs_60", "mmseqs_80", "mmseqs_100"))) %>%
  ggplot(aes(x = recall_ind, y = mean_precision, fill = method, color = method)) + 
  geom_ribbon(aes(ymin = mean_precision - error_margin, ymax = mean_precision + error_margin), alpha = 0.2, color = NA) + 
  geom_line() + 
  theme_bw() +
    scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  scale_fill_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  xlab("Recall") + ylab("Precision") +
  coord_cartesian(ylim = c(0.5, 1), xlim = c(0, 0.15)) +
png("compare/outputs/best_prc_plot.png", width = 5, height = 3, units = "in", res = 300)
print(best_prc_plot)
dev.off()

best_prc_plot2 <- test_data %>%
  filter(profile_size == 900 & filter == 0.9) %>%
    mutate(method =  factor(method, levels = c("fusion", "fusion2", "mmseqs_40", "mmseqs_60", "mmseqs_80", "mmseqs_100"))) %>%
  ggplot(aes(x = recall_ind, y = mean_precision, fill = method, color = method)) + 
  geom_ribbon(aes(ymin = mean_precision - error_margin, ymax = mean_precision + error_margin), alpha = 0.2, color = NA) + 
  geom_line() + 
  theme_bw() +
    scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  scale_fill_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  xlab("Recall") + ylab("Precision") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 0.25)) +
png("compare/outputs/best_prc_plot2.png", width = 5, height = 3, units = "in", res = 300)
print(best_prc_plot2)
dev.off()
ggsave("compare/outputs/best_prc_plot2.png", best_prc_plot2, width = 8, height = 6, units = "in", dpi = 300)

best_prc_plot3 <- test_data %>%
  filter(profile_size == 900 & filter == 0.9) %>%
    mutate(method =  factor(method, levels = c("fusion", "fusion2", "mmseqs_40", "mmseqs_60", "mmseqs_80", "mmseqs_100"))) %>%
  ggplot(aes(x = recall_ind, y = mean_precision, fill = method, color = method)) + 
  geom_ribbon(aes(ymin = mean_precision - error_margin, ymax = mean_precision + error_margin), alpha = 0.2, color = NA) + 
  geom_line() + 
  theme_bw() +
    scale_color_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  scale_fill_manual(
  name = "Method", 
  values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
  labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  xlab("Recall") + ylab("Precision") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
png("compare/outputs/best_prc_plot3.png", width = 5, height = 3, units = "in", res = 300)
print(best_prc_plot3)
dev.off()
ggsave("compare/outputs/best_prc_plot3.png", best_prc_plot3, width = 8, height = 6, units = "in", dpi = 300)

tile_plot <- function(.x, measure, methodx, start_color = "white", end_color = "red", fill_max = 0.5, fill_min = 0){
  res <- .x %>%
  filter(method == methodx) %>%
  rename(mean = !!paste0("mean_", measure), sd = !!paste0("sd_", measure)) %>%
  ggplot(aes(x = filter, y = profile_size, fill = mean)) +
  geom_tile() +
  scale_fill_gradient(low = start_color, high = end_color) +
  theme_bw() + labs(fill = paste("Mean", stringr::str_to_title(measure))) + 
  ylab("Profile Length") + xlab("Max Profile Saturation") + 
  geom_text(aes(label = paste0(round(mean, 2)," (", round(sd, 3), ")")), size = 3) +
  labs(title = paste0(gsub("mmseqs_", "MMseqs2 ", stringr::str_to_title(methodx), ignore.case = TRUE), " - Mean ", stringr::str_to_title(measure)))
  return(res)
}

pdf("compare/outputs/fusion_mmseqs_variable_lengths.pdf")
p0 <- tile_plot(vl_df, methodx = "fusion", measure = "precision",); print(p0)
p0b <- tile_plot(vl_df, methodx = "fusion", measure = "recall"); print(p0b)
p0c <- tile_plot(vl_df, methodx = "fusion", measure = "f1"); print(p0c)

p1 <- tile_plot(vl_df, methodx = "mmseqs_40", measure = "precision"); print(p1)
p1b <- tile_plot(vl_df, methodx = "mmseqs_40", measure = "recall"); print(p1b)
p1c <- tile_plot(vl_df, methodx = "mmseqs_40", measure = "f1"); print(p1c)
dev.off()

png("compare/outputs/fusion_mmseqs_variable_lengths_precision.png", width = 7, height = 7, units = "in", res = 300)
print(p0)
dev.off()
png("compare/outputs/fusion_mmseqs_variable_lengths_recall.png", width = 7, height = 7, units = "in", res = 300)
print(p0b)
dev.off()
png("compare/outputs/fusion_mmseqs_variable_lengths_f1.png", width = 7, height = 7, units = "in", res = 300)
print(p0c)
dev.off()

png("compare/outputs/mmseqs_40_variable_lengths_precision.png", width = 7, height = 7, units = "in", res = 300)
print(p1)
dev.off()
png("compare/outputs/mmseqs_40_variable_lengths_recall.png", width = 7, height = 7, units = "in", res = 300)
print(p1b)
dev.off()
png("compare/outputs/mmseqs_40_variable_lengths_f1.png", width = 7, height = 7, units = "in", res = 300)
print(p1c)
dev.off()