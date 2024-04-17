library(tidyverse)
library(ggplot2)
rm(list = ls())

FigS2_data <- readRDS("FigS2_data.RDS")

FigS2 <- FigS2_data %>%
  filter(!(variable %in% c("n_ortho", "org_per_phyla"))) %>%
  mutate(value = ifelse(variable == "count", log10(value), value), variable_label = ifelse(variable == "count", "Protein Pairs (log10)", variable_label)) %>%
  ggplot(aes(x = variable_label, y = value, fill = method)) +
  geom_boxplot() + 
  theme_bw() +
  xlab("") + 
  labs(title = "Fusion vs MMseqs2 - Variable Comparison") +
  scale_fill_manual(values = c("red1", "black")) +
  facet_wrap(~variable_label, scales = "free") +
  theme(strip.text = element_blank()) +
  labs(fill = "Method", title = "") +
  xlab("") + ylab("")

pdf("FigS2.pdf", height = 8, width = 8)
FigS2
dev.off()