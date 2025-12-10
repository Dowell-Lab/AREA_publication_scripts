# Set working directory
setwd("/Users/bella/Downloads/")

# Set Environment 
library(tidyverse)

library(clusterProfiler)

library(introdataviz)

library(ggplot2)

# Load data in 
# boolean is binary and rank is value 
# Run is same vibe as participant ID 
import.boolean.NB <- read_csv("NB_boolean.csv")

import.rank.NB <- read_csv("NB_rankfile.csv")

import.results.NB <- read_csv("outputarea_scores_20250919-083438.adjpval.csv")

# Merge and prep df for graphng 

NB_merged <- merge (import.boolean.NB, import.rank.NB, by = "Participant")

merged_hypoxia_NB <- NB_merged  %>%
  mutate(group.hypoxia = case_when(
    hypoxia_unfavorable == 1 &  hypoxia_favorable == 0 ~ "Unfavorable Hypoxia Prognosis", 
    hypoxia_unfavorable == 0 & hypoxia_favorable == 1 ~ "Favorable Hypoxia Prognosis"))

merged_hypoxia_NB$GAL <- as.numeric(merged_hypoxia_NB$GAL)

ggplot(merged_hypoxia_NB, aes( x = "NB Hypoxia Prognosis", y = GAL, fill = group.hypoxia, color = group.hypoxia, split = group.hypoxia))  +
  labs (x = "NB Hypoxia Prognosis", y = "GAL", title = "Enrichment of GAL") +
  theme_minimal() + 
  geom_jitter(
    aes(color = group.hypoxia),
    size = 1,
    alpha = 0.5,
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = 0.6))  +
  geom_boxplot(width = .1) +
  geom_split_violin(alpha = 0.7, trim = FALSE) +
  scale_fill_manual(values = c("Unfavorable Hypoxia Prognosis" = "#E0F2F1FF", "Favorable Hypoxia Prognosis" = "#E0F2F1FF")) +
  scale_color_manual(values = c("Unfavorable Hypoxia Prognosis" = "#4DB6ACFF", "Favorable Hypoxia Prognosis" = "#FF6F00FF")) +
  theme(legend.position = "none") + 
  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5))
