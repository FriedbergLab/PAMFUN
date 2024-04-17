library(tidyverse)
library(ggplot2)
rm(list = ls())

Fig5_data <- readRDS("Fig5_data.rds")

Fig5 <- Fig5_data %>%
  mutate(Method = ifelse(method == "pair", "Pair", "Group")) %>%
  mutate(Method = factor(Method, levels = c("Pair", "Group"))) %>%
  ggplot(aes(y = module_length, x = recall, shape = Method, color = Method)) + 
  geom_point() + 
  labs(title = "", color = "Recall Method", shape = "Recall Method") + 
  ylab("Number of Module Proteins") + xlab("Mean Recall") +
  theme_bw() +
  theme(legend.position = c(0.9, 0.1), # move legend to bottom right corner
        legend.justification = c(0.9, 0.1), # align the legend to the bottom right
        legend.box.just = "right", # right justify legend box
        legend.margin = margin(t = 5, r = 5, b = 5, l = 5), # adjusts legend margin 
        # legend.box.background = element_blank(), # adjusts legend background 
        legend.box.margin = margin(), # adjusts spacing around the legend box 
        legend.background = element_rect(linetype = "solid", color = "black", fill = "white")) # black border

# save to file
ggsave("Fig5.png", Fig5, width = 5, height = 3, units = "in", dpi = 300)