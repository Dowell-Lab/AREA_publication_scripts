# Set working directory
setwd("~/DowellLab/raw_data/area_metabolites_data/")

# Set Environment 
library(tidyverse)
library(dplyr)

# Load data in 
import.metabolites <- read_csv("areaoutputarea_scores_20250721-144032.adjpval.csv")

import.binary <- read_csv("metabolomics_binary.csv")

import.values <- read_csv("metabolomics_values.csv")

import.whicht21 <- read_csv("whichT21.csv")
# this file shows the various ways T21 is referred to in the data. Currently only using the complete T21 column for DS classification 

# Filtering based on morbidity (keeping T21, morbidity, and ParticipantID)
binary.filtered.obesity <- import.binary [ ,-c(1:13, 15:30, 32:49)]

binary.filtered.sleepapnea <- import.binary [ ,-c(1:13, 15:45, 47:49)]

binary.filtered.skeletal <- import.binary [ ,-c(1:13,15:43, 45:49)]

binary.filteredcongenitalheartdisease <- import.binary [ ,-c(1:13,15,17:49)]

binary.folliculitis <- import.binary [ ,-c(1:13, 15:21, 23:49)]

binary.just.T21 <- import.binary [ ,-c(1:13, 15:49)]

binary.respiratory_system_disorder <- import.binary [ ,-c(1:13, 15:41, 43:49)]

binary.test <- import.binary [ ,-c(1:13, 15:30, 32:45, 47:49)]
  
# Merge 

test_merge <- merge()
merged_obesity <- merge(binary.filtered.obesity,import.values, by = "ParticipantID")

merged_sleepapnea <- merge(binary.filtered.sleepapnea, import.values, by = "ParticipantID")

merged_skeletal <- merge(binary.filtered.skeletal, import.values, by = "ParticipantID")

merged_cogenitalheartdisease <- merge(binary.filteredcongenitalheartdisease, import.values, by = "ParticipantID")

merged_folliculitis <- merge(binary.folliculitis, import.values, by = "ParticipantID")

merged_just_T21 <- merge (binary.just.T21, import.values, by = "ParticipantID")

merged_respiratory <- merge (binary.respiratory_system_disorder, import.values, by = "ParticipantID")

merged_test <- merge (binary.test, import.values, by = "ParticipantID")

# Set up some groups to be able to do violin plot and also make sure they're in the order you want 
merged_obesity <- merged_obesity %>%
  mutate(group.obesity = case_when(
    obesity_disorder == 0 & complete_trisomy_21 == 0 ~ "D21, Obesity-", 
    obesity_disorder == 1 & complete_trisomy_21 == 0 ~ "D21, Obesity+", 
    obesity_disorder == 0 & complete_trisomy_21 == 1 ~ "T21, Obesity-", 
    obesity_disorder == 1 & complete_trisomy_21 == 1 ~ "T21, Obesity+"))

merged_sleepapnea <- merged_sleepapnea %>%
  mutate(Condition = case_when(
    sleep_apnea_syndrome == 0 & complete_trisomy_21 == 0 ~ "Disomy 21 without Sleep Apnea Syndrome", 
    sleep_apnea_syndrome == 1 & complete_trisomy_21 == 0 ~ "Disomy 21 with Sleep Apnea Syndrome", 
    sleep_apnea_syndrome == 0 & complete_trisomy_21 == 1 ~ "Trisomy 21 without Sleep Apea Syndrome", 
    sleep_apnea_syndrome == 1 & complete_trisomy_21 == 1 ~ "Trisomy 21 with Sleep Apnea Syndrome"))

merged_skeletal <- merged_skeletal %>%
  mutate(group.skeletal = case_when(
    skeletal_system_disorder == 0 & complete_trisomy_21 == 0 ~ "D21, Skeletal-", 
    skeletal_system_disorder == 1 & complete_trisomy_21 == 0 ~ "D21, Skeletal+", 
    skeletal_system_disorder == 0 & complete_trisomy_21 == 1 ~ "T21, Skeletal-", 
    skeletal_system_disorder == 1 & complete_trisomy_21 == 1 ~ "T21, Skeletal+"))

merged_cogenitalheartdisease <- merged_cogenitalheartdisease %>%
  mutate(group.heart = case_when(
    congenital_heart_disease == 0 & complete_trisomy_21 == 0 ~ "D21, Congenital heart-", 
    congenital_heart_disease == 1 & complete_trisomy_21 == 0 ~ "D21, Congenital heart+", 
    congenital_heart_disease == 0 & complete_trisomy_21 == 1 ~ "T21, Congenital heart-", 
    congenital_heart_disease == 1 & complete_trisomy_21 == 1 ~ "T21, Congenital heart+"))

merged_folliculitis <- merged_folliculitis %>%
  mutate(group.folliculitis = case_when(
    folliculitis == 0 & complete_trisomy_21 == 0 ~ "D21, Folliculitis-", 
    folliculitis == 1 & complete_trisomy_21 == 0 ~ "D21, Folliculitis+", 
    folliculitis == 0 & complete_trisomy_21 == 1 ~ "T21, Folliculitis-", 
    folliculitis == 1 & complete_trisomy_21 == 1 ~ "T21, Folliculitis+"))

merged_just_T21 <- merged_just_T21 %>%
  mutate(Condition = case_when(
    complete_trisomy_21 == 0 ~ "Disomy 21", 
    complete_trisomy_21 == 1 ~ "Trisomy 21"))

merged_respiratory <- merged_respiratory %>%
  mutate(group.respiratory = case_when(
    respiratory_system_disorder == 0 & complete_trisomy_21 == 0 ~ "D21, Respiratory_disorder-", 
    respiratory_system_disorder == 1 & complete_trisomy_21 == 0 ~ "D21, Respiratory_disorder+", 
    respiratory_system_disorder == 0 & complete_trisomy_21 == 1 ~ "T21, Respiratory_disorder-", 
    respiratory_system_disorder == 1 & complete_trisomy_21 == 1 ~ "T21, Respiratory_disorder+"))


#-----------------------------------------------------------------------------------
# ----------------------------------Violin plots------------------------------------
#-----------------------------------------------------------------------------------

# violin plot of obesity (remember the metabolite is negatively enriched for obesity, so we expect lower)
ggplot(merged_obesity, aes( x = group.obesity, y = `(5-L-Glutamyl)-L-glutamine`, fill = group.obesity, color = group.obesity))  +
  geom_violin(trim = FALSE) +
  labs (x = "Obesity", y = "(5-L-Glutamyl)-L-glutamine", title = "Enrichment of 5-L-Glutamyl-L-glutamine") +
  theme_minimal() + 
  geom_jitter(alpha = .5, size = 1) +
  geom_boxplot(width = .1) +
  scale_fill_manual(values = c("D21, Obesity-" = "#E0F2F1FF", "D21, Obesity+" = "#E0F2F1FF", "T21, Obesity-" = "#00695CFF", "T21, Obesity+" = "#00695CFF")) +
  scale_color_manual(values = c("D21, Obesity-" = "black", "D21, Obesity+" = "#FF6F00FF", "T21, Obesity-" = "black", "T21, Obesity+" = "#FF6F00FF")) +
  theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5))
 
# stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white", color = "black")
  #scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20))

ggsave('../../figures/area_metabolites/obesity_glutamine_t21_d21_violin.png', 
       height = 6, width = 7)


# violin plot of sleep apnea (positvely enriched)
ggplot(merged_sleepapnea, aes( x = Condition, y = `Taurolithocholic acid`, fill = Condition, color = Condition))  +
  geom_violin(trim = FALSE, color = "black") +
  labs (x = "Sleep Apnea Syndrome", y = " Enrichment of Taurolithocholic acid", title = "Enrichment of Taurolithocholic Acid") +
  theme_minimal() + 
  geom_jitter(alpha = .5, size = 1) +
  #geom_rug() +
  geom_boxplot(width = .1) +
  scale_fill_manual(values = c("Disomy 21 without Sleep Apnea Syndrome" = "tomato1", "Disomy 21 with Sleep Apnea Syndrome" = "tomato1", "Trisomy 21 without Sleep Apea Syndrome" = "cornflowerblue", "Trisomy 21 with Sleep Apnea Syndrome" = "cornflowerblue")) +
  scale_color_manual(values = c("Disomy 21 without Sleep Apnea Syndrome" = "pink", "Disomy 21 with Sleep Apnea Syndrome" = "lightgreen", "Trisomy 21 without Sleep Apea Syndrome" = "pink", "Trisomy 21 with Sleep Apnea Syndrome" = "lightgreen")) +
  theme(legend.position = "top") + 
  theme(plot.title = element_text(hjust = 0.5))
  

ggsave('../../figures/area_metabolites/sleepapnea_taurolithocholicacid_t21_d21_violin.png', 
       height = 6, width = 7)

# violin plot of skeletal system disorder (negatively enriched)
ggplot(merged_skeletal, aes( x = group.skeletal, y = Malate, fill = group.skeletal, color = group.skeletal))  +
  geom_violin(trim = FALSE, color = "black") +
  labs (x = "Skeletal Sleep Disorder", y = "Malate", title = "Enrichment of Malate") +
  theme_minimal() + 
  geom_jitter(alpha = .5, size = 1) +
  geom_boxplot(width = .1) +
  scale_fill_manual(values = c("D21, Skeletal-" = "#E0F2F1FF", "D21, Skeletal+" = "#E0F2F1FF", "T21, Skeletal-" = "#00695CFF", "T21, Skeletal+" = "#00695CFF")) +
  scale_color_manual(values = c("D21, Skeletal-" = "black", "D21, Skeletal+" = "#FF6F00FF", "T21, Skeletal-" = "black", "T21, Skeletal+" = "#FF6F00FF")) +
  theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave('../../figures/area_metabolites/skeletalsystemdisorder_malate_t21_d21_violin.png', 
       height = 6, width = 7)

 # stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white", color = "black")
#scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20))

# violin plot of congenital heart disease (positvely enriched)
ggplot(merged_cogenitalheartdisease, aes( x = group.heart, y = `Taurolithocholic acid`, fill = group.heart))  +
  geom_violin(trim = FALSE, color = "black") +
  labs (x = "Congenital Heart Disease", y = "Taurolithocholic acid", title = "Enrichment of Taurolithocholic Acid") +
  theme_minimal() + 
  geom_jitter(alpha = .5, size = 1) +
  #geom_rug() +
  geom_boxplot() + 
  theme(axis.text.x = element_text(size = 7))

ggsave('../../figures/area_metabolites/congenital_taurolithocholicacid_t21_d21_violin.png', 
       height = 6, width = 7)

# violin plot of folliculitis (skin condition) (negatively enriched)
ggplot(merged_folliculitis, aes( x = group.folliculitis, y = Serine, fill = group.folliculitis))  +
  geom_violin(trim = FALSE, color = "black") +
  labs (x = "Folliculitis", y = "Serine") +
  theme_minimal() + 
  geom_jitter(alpha = .5, size = 1) +
  geom_rug() +
  geom_boxplot() + 
  theme(axis.text.x = element_text(size = 7))

# violin plot of just T21 (negatively enriched)
ggplot(merged_just_T21, aes( x = Condition, y = Hypoxanthine, fill = Condition, , color = Condition))  +
  geom_violin(trim = FALSE, color = "black") +
  labs (y = "Enrichment Level of Hypoxanthine", title = "Enrichment of Hypoxanthine") +
  theme_minimal() + 
  geom_jitter(alpha = .5, size = 1) +
  geom_boxplot(width = .1) +
scale_fill_manual(values = c("Disomy 21" = "tomato1", "Trisomy 21" = "cornflowerblue")) +
 scale_color_manual(values = c("Disomy 21" = "pink", "Trisomy 21" = "lightblue")) +
  theme(legend.position = "top") + 
 theme(plot.title = element_text(hjust = 0.5))

ggsave('../../figures/area_metabolites/justT21_violin.png', 
       height = 6, width = 7)

# violin plot of respiratory system disorder (negatively enriched)
ggplot(merged_respiratory, aes( x = group.respiratory, y = Citrate, fill = group.respiratory))  +
  geom_violin(trim = FALSE, color = "black") +
  labs (x = "Respiratory System Disorder", y = "Citrate") +
  theme_minimal() + 
  geom_jitter(alpha = .5, size = 1) +
  geom_rug() +
  geom_boxplot() + 
  theme(axis.text.x = element_text(size = 7))


# sanity check 

merged_obesity <- merged_obesity %>%
  mutate(group.justobesity = case_when(
    obesity_disorder == 0 ~ "Obesity-", 
    obesity_disorder == 1 ~ "Obesity+"))

  ggplot(merged_obesity, aes( x = group.justobesity, y = `(5-L-Glutamyl)-L-glutamine`, fill = group.justobesity))  +
  geom_violin(trim = FALSE, color = "black") +
  labs (x = "Obesity", y = "(5-L-Glutamyl)-L-glutamine") + 
  #scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20))
  theme_minimal() +
  geom_rug() +
  geom_boxplot()
  
  ggsave('../../figures/area_metabolites/obesity_glutamine_violin.png', 
         height = 6, width = 7)
  
  #-----------------------------------------------------------------------
  #---------------------Bubble Plots---------------------------------------
  #------------------------------------------------------------------------
  
 comorbid_list <- unique(import.metabolites$binary_attribute)
  
  # bubble plot of all metabolites vs complete T21
filtered_metabolites_t21 <- import.metabolites %>% 
  filter(binary_attribute == "complete_trisomy_21") %>% 
  filter(p_value_bonf < .01) %>% 
  filter(abs(NES) > 3) %>% 
  slice_head(n=30)
  #arrange(NES >0, desc(p_value_bonf), NES <0)  

# some of our p bonf are zero --> change them to a super low value so that they show up on the plot (do this after looking at the plot)

filtered_metabolites_t21$pval_fixed <- ifelse(filtered_metabolites_t21$p_value_bonf == 0,
                                             .Machine$double.xmin,
                                             filtered_metabolites_t21$p_value_bonf)


filtered_metabolites_t21$value <- factor(filtered_metabolites_t21$value, levels = unique(filtered_metabolites_t21$value))

ggplot(filtered_metabolites_t21, aes(x = NES, y = reorder(value, NES), size = -log10(pval_fixed), color = NES)) +
  geom_point(alpha = .7) + 
  labs (x = "Normalized Enrichment Score", y = "Metabolite", title = "Top 30 Enriched Metabolites for Complete Trisomy 21") + 
  scale_color_gradient(low = "tomato1", high = "cornflowerblue") +
  theme(axis.text.y = element_text(size = 6)) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

ggsave('../../figures/area_metabolites/t21_metabolites_bubble.png', 
       height = 6, width = 7)

# bubble plot of obesity and its metabolites 
filtered_metabolites_obesity <- import.metabolites %>% 
  filter(binary_attribute == "obesity_disorder") %>% 
  filter(p_value_bonf < .01) %>% 
  filter(abs(NES) > 3) 
 # slice_head(n = 25)

filtered_metabolites_obesity$pval_fixed <- ifelse(filtered_metabolites_obesity$p_value_bonf == 0,
                                              .Machine$double.xmin,
                                              filtered_metabolites_t21$p_value_bonf)

ggplot(filtered_metabolites_obesity, aes(x = NES, y = reorder(value, NES), size = -log10(pval_fixed), color = NES)) +
  geom_point(alpha = .7) + 
  labs (x = "Normalized Enrichment Score", y = "Metabolite", title = "Obesity Disorder") + 
  scale_color_gradient(low = "tomato1", high = "cornflowerblue") +
  theme(axis.text.y = element_text(size = 6)) + 
  theme_minimal()

ggsave('../../figures/area_metabolites/obesity_metabolites_bubble.png', 
       height = 6, width = 7)

# bubble plot of sleep apnea and its metabolites 
filtered_metabolites_sleepapnea <- import.metabolites %>% 
  filter(binary_attribute == "obstructive_sleep_apnea_syndrome") %>% 
  filter(p_value_bonf < .01) %>% 
  filter(abs(NES) > 3) %>% 
  slice_head(n=30)

ggplot(filtered_metabolites_sleepapnea, aes(x = NES, y = reorder(value, NES), size = -log10(p_value_bonf), color = NES)) +
  geom_point(alpha = .7) + 
  labs (x = "Normalized Enrichment Score", y = "Metabolite", title = "Top 30 Metabolites Enriched for Obstructive Sleep Apnea Syndrome") + 
  scale_color_gradient(low = "tomato1", high = "cornflowerblue") +
  theme(axis.text.y = element_text(size = 6)) + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

ggsave('../../figures/area_metabolites/sleepapnea_metabolites_bubble.png', 
       height = 6, width = 7)

# bubble plot for skeletal 
filtered_metabolites_skeletal <- import.metabolites %>% 
  filter(binary_attribute == "skeletal_system_disorder") %>% 
  filter(p_value_bonf < .01) %>% 
  filter(abs(NES) > 3) 

ggplot(filtered_metabolites_skeletal, aes(x = NES, y = reorder(value, NES), size = -log10(p_value_bonf), color = NES)) +
  geom_point(alpha = .7) + 
  labs (x = "Normalized Enrichment Score", y = "Metabolite", title = "Skeletal System Disorder") + 
  scale_color_gradient(low = "tomato1", high = "cornflowerblue") +
  theme(axis.text.y = element_text(size = 6)) + 
  theme_minimal()

ggsave('../../figures/area_metabolites/skeletal_metabolites_bubble.png', 
       height = 6, width = 7)
# bubble plot of congenital heart disease and its metabolites 

filtered_metabolites_congenital <- import.metabolites %>% 
  filter(binary_attribute == "congenital_heart_disease") %>% 
  filter(p_value_bonf < .01) %>% 
  filter(abs(NES) > 3) %>% 
  slice_head(n = 25)

ggplot(filtered_metabolites_congenital, aes(x = NES, y = value, size = -log10(p_value_bonf), color = NES)) +
  geom_point(alpha = .7) + 
  labs (x = "Normalized Enrichment Score", y = "Metabolite", title = "Congenital Heart Disease") + 
  scale_color_gradient(low = "tomato1", high = "cornflowerblue") +
  theme(axis.text.y = element_text(size = 6)) + 
  theme_minimal()

ggsave('../../figures/area_metabolites/congenital_heart_metabolites_bubble.png', 
       height = 6, width = 7)















