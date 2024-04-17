library(tidyverse)
library(ggplot2)
rm(list = ls())

FigS3_data <- readRDS("FigS3_data.rds")

FigS3 <- FigS3_data %>% 
  mutate(method =  factor(method, levels = c("fusion","mmseqs_40", "mmseqs_60", "mmseqs_80", "mmseqs_100"))) %>%
  mutate(method_label = case_when(method == "fusion" ~ "Fusion",
                                  method == "mmseqs_40" ~ "MMseqs2 - 40",
                                  method == "mmseqs_60" ~ "MMseqs2 - 60",
                                  method == "mmseqs_80" ~ "MMseqs2 - 80",
                                  method == "mmseqs_100" ~ "MMseqs2 - 100")) %>%
  mutate(Sequence_Identity = as.numeric(Sequence_Identity)) %>% 
  ggplot(aes(x = Sequence_Identity, y = Frequency, fill = method)) +
  scale_fill_manual(
    name = "Method", 
    values = c("fusion" = "#000000", "mmseqs_40" = "#E69F00", "mmseqs_60" = "#56B4E9", "mmseqs_80" = "#009E73", "mmseqs_100" = "#F0E442"),
    labels = c("Fusion", "MMseqs2 - 40", "MMseqs2 - 60", "MMseqs2 - 80", "MMseqs2 - 100")) + 
  geom_bar(stat = "identity", position = "identity", alpha = 0.5) + 
  theme_bw() + ylim(0, 0.06) + 
  labs(fill = "Method",x = "Sequence Identity", y = "Frequency", title = "Sequence Identity of Identified Proteins") +
  facet_wrap(~method_label, scales = "free_x")

# save to file
ggsave("FigS3.png", FigS3, width = 10, height = 5)
