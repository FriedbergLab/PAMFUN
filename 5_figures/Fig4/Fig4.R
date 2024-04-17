library(tidyverse)
library(ggplot2)
rm(list = ls())

Fig1_data <- readRDS("Fig1_data.rds")
Fig1 <- Fig1_data  %>%
  # plot recall on x axis and precision on y axis
  ggplot(aes(x = recall_ind, y = mean_precision, fill = method, color = method)) + 
  # add error bars
  geom_ribbon(aes(ymin = mean_precision - error_margin, ymax = mean_precision + error_margin), alpha = 0.2, color = NA) + 
  # add lines
  geom_line() + 
  # black and white color scheme
  theme_bw() +
  # define color legend
  scale_color_manual(
    name = "Method", 
    values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
    labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  scale_fill_manual(
    name = "Method", 
    values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
    labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) +
  # label x axis and y axis
  xlab("Recall") + ylab("Precision") +
  # set x and y axis limits
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 0.25)) 

# save to file
png("Fig1.png", width = 5, height = 3, units = "in", res = 300)
print(Fig1)
dev.off()