library(tidyverse)
library(ggplot2)
rm(list = ls())

system("cp ../large_prc_plot.rds FigS1_data.rds")

FigS1_data <- readRDS("FigS1_data.rds")

FigS1 <- FigS1_data  %>%
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
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) 
png("FigS1.png", width = 5, height = 3, units = "in", res = 300)
print(FigS1)
dev.off()