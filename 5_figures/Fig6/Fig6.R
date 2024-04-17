library(tidyverse)
library(ggplot2)
rm(list = ls())

Fig6_data <- readRDS("Fig6_data.rds")

Fig6 <- Fig6_data  %>%
  ggplot(aes(x = size_group, y = jaccard, fill = method_label)) + 
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  labs(fill = "Method") + ylab("Jaccard Index") + xlab("Module Size")

png("Fig6.png", width = 5, height = 3, units = "in", res = 300)
print(Fig6)
dev.off()

