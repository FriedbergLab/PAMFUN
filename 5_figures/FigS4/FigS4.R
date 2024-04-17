library(tidyverse)
library(ggplot2)
rm(list = ls())

FigS4_data <- readRDS("FigS4_data.rds")

FigS4 <- FigS4_data %>%
  filter(sim_floor == 0.9) %>%
  filter(louvain_value == 4) %>%
  ggplot(aes(x = method_label, y = jaccard, fill = method_label)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~size_group) +
  theme_bw() +
  labs(fill = "Method") + ylab("Jaccard Index") + xlab("")

# save to file
ggsave("FigS4.png", FigS4, width = 10, height = 5)