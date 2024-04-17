library(tidyverse)
library(ggplot2)
library(ggpubr)
rm(list = ls())

Fig2_data <- readRDS("Fig2_data.rds")
dodge_distance = 0.3
Fig2a <- Fig2_data %>%
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
  guides(color = guide_legend(override.aes = list(shape = unique(Fig2_data$method))), shape = FALSE) +
  theme(axis.title.x = element_text(size = 14, vjust = -0.5), 
        axis.title.y = element_text(size = 14, vjust = 3),
        axis.text.x = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.4,0.2), "cm"))

Fig2b <- Fig2_data %>%
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
  guides(color = guide_legend(override.aes = list(shape = unique(Fig2_data$method))), shape = FALSE) +
  theme(axis.title.x = element_text(size = 14, vjust = -0.5), 
        axis.title.y = element_text(size = 14, vjust = 3),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

Fig2c <- Fig2_data %>%
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
  guides(color = guide_legend(override.aes = list(shape = unique(Fig2_data$method))), shape = FALSE) +
  theme(axis.title.x = element_text(size = 14, vjust = -0.5), 
        axis.title.y = element_text(size = 14, vjust = 3),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.4,0.2), "cm"))

Fig2d <- Fig2_data %>%
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
  guides(color = guide_legend(override.aes = list(shape = unique(Fig2_data$method))), shape = FALSE) +
  theme(axis.title.x = element_text(size = 14, vjust = -0.5), 
        axis.title.y = element_text(size = 14, vjust = 3),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

Fig2e <- Fig2_data %>%
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
  guides(color = guide_legend(override.aes = list(shape = unique(Fig2_data$method))), shape = FALSE) +
  theme(axis.title.x = element_text(size = 14, vjust = -0.5), 
        axis.title.y = element_text(size = 14, vjust = 3),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.4,0.2), "cm"))

Fig2f <- Fig2_data %>%
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
  guides(color = guide_legend(override.aes = list(shape = unique(Fig2_data$method))), shape = FALSE) +
  theme(axis.title.x = element_text(size = 14, vjust = -0.5), 
        axis.title.y = element_text(size = 14, vjust = 3),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14))


shared_legend <- get_legend(Fig2f)
Fig2 <- ggarrange(Fig2a, Fig2c, Fig2e, Fig2b, Fig2d, Fig2f, ncol=3, nrow=2, common.legend = TRUE, legend="right", labels = c("A", "B", "C", "D", "E", "F"), legend.grob = shared_legend)
png("Fig2.png", width = 10, height = 5, units = "in", res = 300)
print(Fig2)
dev.off()