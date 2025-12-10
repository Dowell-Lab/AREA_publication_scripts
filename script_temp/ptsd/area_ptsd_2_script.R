# Set working directory
setwd("/Users/ispe1418/area_ptsd_R")

# Set Environment 
library(tidyverse)

library(clusterProfiler)

library(Cairo)

# Load data in 
import.boolean <- read_csv("ptsd.boolean.csv")

import.rank <- read_csv("geneexptpm.rank.csv")

import.ptsd.results <- read_csv("outputarea_scores_20250924-152326.adjpval.csv")

import.gene.list <- read_csv("ptsd.genedf.csv")

filtered.gene.list <- import.gene.list[ ,-c(1:5, 7:11, 13:26)]

# change column name from "Run" to "participant_id" for clarity 
colnames(import.boolean)[2] <- "participant_id"

colnames(import.rank)[1] <- "participant_id"

# Merge and change format of dfs for graphing 

violinplot_merged <- merge (import.boolean, import.rank, by = "participant_id")

violinplot_merged_long <-violinplot_merged %>%
  pivot_longer(cols=c(6:54047), # pick whichever cols are samples with counts 
               names_to='gene_id',
               values_to='gene_expression')

merged_filtered <- merge (filtered.gene.list, violinplot_merged_long, by = "gene_id")

violinplot_merged_filtered_wider <- merged_filtered [ ,c("participant_id", "gene_expression", "gene_name")] %>%
  pivot_wider(
    names_from = gene_name,
    values_from = gene_expression)


ready_df <- merge(import.boolean, violinplot_merged_filtered_wider, by = "participant_id")

#----------------------------------------------------------------------------------------
#--------------------------------------Violin Plots--------------------------------------
#----------------------------------------------------------------------------------------

# NDUFA1, CCDC85B, SNORD54, FKBP5, and SNORD46 were found to be differentially expressed in ppl. with current compared to never PTSD
# all but FKBP5 are down regulated in current 
ready_df <- ready_df %>%
  mutate(group.ptsd = case_when(
    Never == 0 & Past == 0 & Current == 1 ~ " Current ", 
    Never == 0 & Past == 1 & Current == 0 ~ " Past ",
    Never == 1 & Past == 0 & Current == 0 ~ " Never "))

ready_df$FKBP5 <- as.numeric(ready_df$FKBP5)
ready_df$NDUFA1 <- as.numeric(ready_df$NDUFA1)
ready_df$CCDC85B <- as.numeric(ready_df$CCDC85B)
ready_df$SNORD54<- as.numeric(ready_df$SNORD54)
ready_df$SNORD46 <- as.numeric(ready_df$SNORD46)
ready_df$NFIA <- as.numeric(ready_df$NFIA)
ready_df$LOC100128365 <- as.numeric(ready_df$LOC100128365)
ready_df$FOXO3 <- as.numeric(ready_df$FOXO3)


# Graph violin plots
FKBP5_ptsd_violin <- ggplot(ready_df, aes( x = group.ptsd, y = FKBP5, fill = group.ptsd, color = group.ptsd))  +
  geom_violin(trim = FALSE, color = "black") +
  labs (x = "PTSD State", y = "TPM", title = "Expression of FKBP5 in PTSD Study Patients") +
  theme_minimal() + 
  geom_jitter(alpha = .5, size = 1) +
  geom_boxplot(width = .1) + 
  scale_fill_manual(values = c(" Current " = "#cd34b5", " Past " = "#fa8775", " Never " = "#ffd700")) +
  scale_color_manual(values = c(" Current " = "#8B1A7A", " Past " = "#C8432F", " Never " = "#B8860B")) +
  theme(legend.position = "none") + 
  theme(plot.title = element_text(face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.2, color = "grey") +
  theme(plot.title = element_text(hjust = 0.5)) 

print(FKBP5_ptsd_violin)

ggsave('FKBP5_ptsd_violin.png', 
       height = 6, width = 7)

CairoSVG("FKBP5_ptsd_violin.svg")
dev.off()

NDUFA1_ptsd_violin <- ggplot(ready_df, aes( x = group.ptsd, y = NDUFA1, fill = group.ptsd, color = group.ptsd))  +
  geom_violin(trim = FALSE, color = "black") +
  labs (x = "PTSD State", y = "TPM", title = "Expression of NDUFA1 in PTSD Study Patients") +
  theme_minimal() + 
  geom_jitter(alpha = .5, size = 1) +
  geom_boxplot(width = .1) + 
  scale_fill_manual(values = c(" Current " = "#cd34b5", " Past " = "#fa8775", " Never " = "#ffd700")) +
  scale_color_manual(values = c(" Current " = "#8B1A7A", " Past " = "#C8432F", " Never " = "#B8860B")) + 
  theme(legend.position = "none") + 
  theme(plot.title = element_text(face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.2, color = "grey") +
  theme(plot.title = element_text(hjust = 0.5))

print(NDUFA1_ptsd_violin)

ggsave('NDUFA1_ptsd_violin.png', 
       height = 6, width = 7)


CairoSVG("NDUFA1_ptsd_violin.svg")
dev.off()

CCDC85B_ptsd_violin <- ggplot(ready_df, aes( x = group.ptsd, y = CCDC85B, fill = group.ptsd, color = group.ptsd))  +
  geom_violin(trim = FALSE, color = "black") +
  labs (x = "PTSD State", y = "TPM", title = "Expression of CCDC85B in PTSD Study Patients") +
  theme_minimal() + 
  geom_jitter(alpha = .5, size = 1) +
  geom_boxplot(width = .1) + 
  scale_fill_manual(values = c(" Current " = "#cd34b5", " Past " = "#fa8775", " Never " = "#ffd700")) +
  scale_color_manual(values = c(" Current " = "#8B1A7A", " Past " = "#C8432F", " Never " = "#B8860B")) + 
  theme(legend.position = "none") + 
  theme(plot.title = element_text(face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.2, color = "grey") +
  theme(plot.title = element_text(hjust = 0.5))

print(CCDC85B_ptsd_violin)

ggsave('CCDC85B_ptsd_violin.png', 
       height = 6, width = 7)

CairoSVG("CCDC85B_ptsd_violin.svg")
dev.off()

SNORD54_ptsd_violin <- ggplot(ready_df, aes( x = group.ptsd, y = SNORD54, fill = group.ptsd, color = group.ptsd))  +
  geom_violin(trim = FALSE, color = "black") +
  labs (x = "PTSD State", y = "TPM", title = "Expression of SNORD54 in PTSD Study Patients") +
  theme_minimal() + 
  geom_jitter(alpha = .5, size = 1) +
  geom_boxplot(width = .1) + 
  scale_fill_manual(values = c(" Current " = "#cd34b5", " Past " = "#fa8775", " Never " = "#ffd700")) +
  scale_color_manual(values = c(" Current " = "#8B1A7A", " Past " = "#C8432F", " Never " = "#B8860B")) + 
  theme(legend.position = "none") + 
  theme(plot.title = element_text(face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.2, color = "grey") +
  theme(plot.title = element_text(hjust = 0.5))

print(SNORD54_ptsd_violin)

ggsave('SNORD54_ptsd_violin.png', 
       height = 6, width = 7)

CairoSVG("=SNORD54_ptsd_violin.svg")
dev.off()

SNORD46_ptsd_violin <- ggplot(ready_df, aes( x = group.ptsd, y = SNORD46, fill = group.ptsd, color = group.ptsd))  +
  geom_violin(trim = FALSE, color = "black") +
  labs (x = "PTSD State", y = "TPM", title = "Expression of SNORD46 in PTSD Study Patients") +
  theme_minimal() + 
  geom_jitter(alpha = .5, size = 1) +
  geom_boxplot(width = .1) + 
  scale_fill_manual(values = c(" Current " = "#cd34b5", " Past " = "#fa8775", " Never " = "#ffd700")) +
  scale_color_manual(values = c(" Current " = "#8B1A7A", " Past " = "#C8432F", " Never " = "#B8860B")) +
  theme(legend.position = "none") + 
  theme(plot.title = element_text(face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.2, color = "grey") +
  theme(plot.title = element_text(hjust = 0.5))

print(SNORD46_ptsd_violin)

ggsave('SNORD46_ptsd_violin.png', 
       height = 6, width = 7)

CairoSVG("=SNORD46_ptsd_violin.svg")
dev.off()

NFIA_ptsd_violin <- ggplot(ready_df, aes( x = group.ptsd, y = NFIA, fill = group.ptsd, color = group.ptsd))  +
  geom_violin(trim = FALSE, color = "black") +
  labs (x = "PTSD State", y = "TPM", title = "Expression of NFIA in PTSD Study Patients") +
  theme_minimal() + 
  geom_jitter(alpha = .5, size = 1) +
  geom_boxplot(width = .1) +
  scale_fill_manual(values = c(" Current " = "#cd34b5", " Past " = "#fa8775", " Never " = "#ffd700")) +
  scale_color_manual(values = c(" Current " = "#8B1A7A", " Past " = "#C8432F", " Never " = "#B8860B")) +
  theme(legend.position = "none") + 
  theme(plot.title = element_text(face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.2, color = "grey") +
  theme(plot.title = element_text(hjust = 0.5))

print(NFIA_ptsd_violin)

ggsave('NFIA_ptsd_violin.png', 
       height = 6, width = 7)

CairoSVG("=NFIA_ptsd_violin.svg")
dev.off()

LOC100128365_ptsd_violin <- ggplot(ready_df, aes( x = group.ptsd, y = LOC100128365, fill = group.ptsd, color = group.ptsd))  +
  geom_violin(trim = FALSE, color = "black") +
  labs (x = "PTSD State", y = "TPM", title = "Expression of LOC100128365 in PTSD Study Patients") +
  theme_minimal() + 
  geom_jitter(alpha = .5, size = 1) +
  geom_boxplot(width = .1) +
  scale_fill_manual(values = c(" Current " = "#cd34b5", " Past " = "#fa8775", " Never " = "#ffd700")) +
  scale_color_manual(values = c(" Current " = "#8B1A7A", " Past " = "#C8432F", " Never " = "#B8860B")) +
  theme(legend.position = "none") + 
  theme(plot.title = element_text(face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.2, color = "grey") +
  theme(plot.title = element_text(hjust = 0.5))

print(LOC100128365_ptsd_violin)

ggsave('LOC100128365_ptsd_violin.png', 
       height = 6, width = 7)

CairoSVG("=LOC100128365_ptsd_violin.svg")
dev.off()

FOXO3_ptsd_violin <- ggplot(ready_df, aes( x = group.ptsd, y = FOXO3, fill = group.ptsd, color = group.ptsd))  +
  geom_violin(trim = FALSE, color = "black") +
  labs (x = "PTSD State", y = "TPM", title = "Expression of FOXO3 in PTSD Study Patients") +
  theme_minimal() + 
  geom_jitter(alpha = .5, size = 1) +
  geom_boxplot(width = .1) +
  scale_fill_manual(values = c(" Current " = "#cd34b5", " Past " = "#fa8775", " Never " = "#ffd700")) +
  scale_color_manual(values = c(" Current " = "#8B1A7A", " Past " = "#C8432F", " Never " = "#B8860B")) +
  theme(legend.position = "none") + 
  theme(plot.title = element_text(face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.2, color = "grey") +
  theme(plot.title = element_text(hjust = 0.5))

print(FOXO3_ptsd_violin)

ggsave('FOXO3_ptsd_violin.png', 
       height = 6, width = 7)

CairoSVG("=FOXO3_ptsd_violin.svg")
dev.off()

# Bubble Plots 

filtered_current_ptsd_results <- import.ptsd.results %>% 
  filter(boolean_attribute == "Current") %>% 
  filter(p_value_bonf < .01) %>% 
  filter(abs(NES) > 3) %>% 
  slice_head(n = 50)

filtered_never_ptsd_results <- import.ptsd.results %>% 
  filter(boolean_attribute == "Never") %>% 
  filter(p_value_bonf < .01) %>% 
  filter(abs(NES) > 3) #%>% 
 # slice_head(n = 50)

filtered_past_ptsd_results <- import.ptsd.results %>% 
  filter(boolean_attribute == "Past") %>% 
  filter(p_value_bonf < .01) %>% 
  filter(abs(NES) > 3) %>% 
  slice_head(n = 50)

# clean import.gene.list to merge easily and change gene_id to value to match results 
clean_import_gene_List <- import.gene.list[ ,-c(1:5, 7:11, 13:26)]

merged_current_gene_list_result <- merge(clean_import_gene_List,filtered_current_ptsd_results, by = "gene_id")

merged_never_gene_list_result <- merge(clean_import_gene_List,filtered_never_ptsd_results, by = "gene_id")

merged_past_gene_list_result <- merge(clean_import_gene_List,filtered_past_ptsd_results, by = "gene_id")

ggplot(merged_current_gene_list_result, aes(x = NES, y = reorder(gene_name, NES), size = -log10(p_value_bonf), color = NES)) +
  geom_point(alpha = .7) + 
  labs (x = "Normalized Enrichment Score", y = "Genes", title = "Genes Associated with Current PTSD") + 
  scale_color_gradient(low = "tomato1", high = "cornflowerblue") +
  theme(axis.text.y = element_text(size = 6)) + 
  theme_minimal() + 
  theme(plot.title = element_text(face = "bold")) + 
  guides(
    color = guide_colorbar(order = 1),
    size = guide_legend(order = 2)), 
  theme(plot.title = element_text(hjust = 0.5))

ggsave('current_genes_bubble.png', 
       height = 9, width = 7)

ggplot(merged_never_gene_list_result, aes(x = NES, y = reorder(gene_name, NES), size = -log10(p_value_bonf), color = NES)) +
  geom_point(alpha = .7) + 
  labs (x = "Normalized Enrichment Score", y = "Genes", title = "Genes Associated with Never PTSD") + 
  scale_color_gradient(low = "tomato1", high = "cornflowerblue") +
  theme(axis.text.y = element_text(size = 6))  +
  theme_minimal() + 
  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave('never_genes_bubble.png', 
       height = 9, width = 7)

ggplot(merged_past_gene_list_result, aes(x = NES, y = reorder(gene_name, NES), size = -log10(p_value_bonf), color = NES)) +
  geom_point(alpha = .7) + 
  labs (x = "Normalized Enrichment Score", y = "Genes", title = "Genes Associated with Past PTSD") + 
  scale_color_gradient(low = "tomato1", high = "cornflowerblue") +
  theme(axis.text.y = element_text(size = 6)) + 
  theme_minimal() + 
  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5))
#theme(text = element_text(family = "Arial"))

ggsave('past_genes_bubble.png', 
       height = 6, width = 7)


#------------------------------------------------------------------------------------------
#----------------------------------GSEA----------------------------------------------------
#------------------------------------------------------------------------------------------

library(msigdbr)
library(msigdbdf)
library(clusterProfiler)
library(org.Hs.eg.db)

#Making a decreasing sorted vector based on p value and L2FC 

colnames(import.ptsd.results)[2] <- "gene_id"

gsea_starting_df <- merge(filtered.gene.list, import.ptsd.results, by = "gene_id")

gsea.never.sorted <- gsea_starting_df %>% 
  filter(boolean_attribute == "Never") 

gsea.past.sorted <- gsea_starting_df %>% 
  filter(boolean_attribute == "Past") 

gsea.current.sorted <- gsea_starting_df %>% 
  filter(boolean_attribute == "Current") 

rankeddf_never <- tibble(gene=gsea.never.sorted$gene_name,
                         rnk = -log(gsea.never.sorted$p_value_bonf)*sign(gsea.never.sorted$NES)) %>% 
  arrange(desc(rnk)) %>% 
  drop_na() 


rankeddf_past <- tibble(gene=gsea.past.sorted$gene_name,
                        rnk = -log(gsea.past.sorted$p_value_bonf)*sign(gsea.past.sorted$NES)) %>% 
  arrange(desc(rnk)) %>% 
  drop_na()

rankeddf_current <- tibble(gene=gsea.current.sorted$gene_name,
                           rnk = -log(gsea.current.sorted$p_value_bonf)*sign(gsea.current.sorted$NES)) %>% 
  arrange(desc(rnk)) %>% 
  drop_na()

# turns out there's LOTS of duplicates (prolly isoforms, so let's take the highest one and call it a day)
rankeddf_never_unique <- rankeddf_never %>%
  group_by(gene) %>%
  slice_max(order_by = abs(rnk), n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  arrange(desc(rnk))

write.table(rankeddf_never_unique, '25_11_02_rankeddf_never_unique.txt' )

rankeddf_past_unique <- rankeddf_past %>%
  group_by(gene) %>%
  slice_max(order_by = abs(rnk), n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  arrange(desc(rnk))

write_csv(rankeddf_past_unique, '25_11_02_rankeddf_past_unique.csv' )

rankeddf_current_unique <- rankeddf_current %>%
  group_by(gene) %>%
  slice_max(order_by = abs(rnk), n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  arrange(desc(rnk))

write_csv(rankeddf_current_unique, '25_11_02_rankeddf_current_unique.csv' )

# turning it into a vector for GSEA
never_genes_vector <-rankeddf_never_unique$rnk
names(never_genes_vector) <-rankeddf_never_unique$gene

past_genes_vector <-rankeddf_past_unique$rnk
names(past_genes_vector) <-rankeddf_past_unique$gene

current_genes_vector <-rankeddf_current_unique$rnk
names(current_genes_vector) <-rankeddf_current_unique$gene


# Load in gene sets 

hallmark_geneset <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

c2_geneset <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, gene_symbol)

cp_geneset <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP") %>% 
  dplyr::select(gs_name, gene_symbol)

c5_geneset <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, gene_symbol)

# Run GSEA 

never_geneset <- GSEA(never_genes_vector, 
                      TERM2GENE = c2_geneset, 
                      pAdjustMethod = 'BH', 
                      minGSSize = 10, 
                      maxGSSize = 500, 
                      eps = 1e-20, 
                      pvalueCutoff = 0.05,
                      scoreType = "pos")

current_geneset <- GSEA(current_genes_vector, 
                        TERM2GENE = hallmark_geneset, 
                        pAdjustMethod = 'BH', 
                        minGSSize = 10, 
                        maxGSSize = 500, 
                        eps = 1e-20, 
                        pvalueCutoff = 0.05,
                        scoreType = "pos")

current_geneset_c5 <- GSEA(current_genes_vector, 
                           TERM2GENE = c5_geneset, 
                           pAdjustMethod = 'BH', 
                           minGSSize = 10, 
                           maxGSSize = 500, 
                           eps = 1e-20, 
                           pvalueCutoff = 0.05,
                           scoreType = "pos")

past_geneset <- GSEA(past_genes_vector, 
                     TERM2GENE = hallmark_geneset, 
                     pAdjustMethod = 'BH', 
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     eps = 1e-20, 
                     pvalueCutoff = 0.05,
                     scoreType = "pos")

past_geneset_c5 <- GSEA(past_genes_vector, 
                        TERM2GENE = c5_geneset, 
                        pAdjustMethod = 'BH', 
                        minGSSize = 10, 
                        maxGSSize = 500, 
                        eps = 1e-20, 
                        pvalueCutoff = 0.05,
                        scoreType = "pos")

# turn results into a tibble so that we can graph it 

never_geneset_tibble <- as_tibble(never_geneset)

current_geneset_tibble <- as_tibble(current_geneset)

current_geneset_c5_tibble <- as_tibble(current_geneset_c5)

past_geneset_tibble <- as_tibble(past_geneset)

past_geneset_c5_tibble <- as_tibble(past_geneset_c5)

# filter tibbles so that its legible on the graphs (do this step after ggplotting to see which need it)
current_geneset_c5_sorted <- current_geneset_c5_tibble %>% 
  arrange(abs(NES)) %>% 
  slice_head(n = 30)

past_geneset_c5_sorted <- past_geneset_c5_tibble %>% 
  arrange(abs(NES)) %>% 
  slice_head(n = 30)

# plotting 
ggplot(never_geneset_tibble, aes(x = NES, y = reorder(ID, NES), fill = pvalue)) + 
  geom_col() + 
  scale_fill_gradient(low="tomato1", high="cornflowerblue") + 
  labs(title = "Never PTSD GSEA, C2 Geneset")

ggplot(current_geneset_tibble, aes(x = NES, y = reorder(ID, NES), fill = pvalue )) + 
  geom_col() + 
  scale_fill_gradient(low="tomato1", high="cornflowerblue") + 
  labs(title = "Current PTSD GSEA, Hallmark Geneset")

ggplot(current_geneset_c5_sorted, aes(x = NES, y = reorder(ID, NES), fill = pvalue )) + 
  geom_col() + 
  scale_fill_gradient(low="tomato1", high="cornflowerblue") + 
  labs(title = "Current PTSD GSEA, C5 Geneset")

ggplot(past_geneset_tibble, aes(x = NES, y = reorder(ID, NES), fill = pvalue)) + 
  geom_col() + 
  scale_fill_gradient(low="tomato1", high="cornflowerblue") + 
  labs(title = "Past PTSD GSEA, Hallmark")

ggplot(past_geneset_c5_sorted, aes(x = NES, y = reorder(ID, NES), fill = pvalue)) + 
  geom_col() + 
  scale_fill_gradient(low="tomato1", high="cornflowerblue") + 
  labs(title = "Past PTSD GSEA, C5 Geneset")

# Heatmap

heatmap_df <- merge (import.ptsd.results, filtered.gene.list, by = "gene_id" )

library(pheatmap)

heatmap_df <- heatmap_df %>%  
  dplyr::select(boolean_attribute, gene_name, NES, everything()) %>% 
  dplyr::mutate(rnk = -log(p_value_bonf) * sign(NES)) %>%
  group_by(gene_name) %>%
  slice_max(order_by = abs(rnk), n = 1, with_ties = FALSE) %>%
  ungroup()

heatmap_matrix <- heatmap_df %>%
  dplyr::select(boolean_attribute, gene_name, NES) %>%
  tidyr::pivot_wider(names_from = boolean_attribute, values_from = NES) 

rownames_saved <- heatmap_matrix$gene_name

heatmap_matrix <- as.matrix(heatmap_matrix[,-1])
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)

rownames(heatmap_matrix) <- rownames_saved

heatmap_matrix[is.na(heatmap_matrix)] <- 0

heatmap_matrix[is.na(heatmap_matrix)] <- 0

annotation_row <- heatmap_df %>%
  group_by(gene_name) %>%
  summarize(p_value_bonf = if_else(any(p_value_bonf < 0.05), "<.05", "No")) %>%
  column_to_rownames("gene_name")

top_genes <- apply(heatmap_matrix, 1, function(x) max(abs(x), na.rm = TRUE))

top_genes <- names(sort(top_genes, decreasing = TRUE))[1:50]

heatmap_matrix_top <- heatmap_matrix[top_genes, ]


pheatmap(
  heatmap_matrix_top,
  scale = "row",
  annotation_row = annotation_row,
  color = colorRampPalette(c("cornflowerblue", "white", "tomato1"))(100),
  # clustering_distance_rows = "euclidean",
  # clustering_method = "complete",
  main = "Top 50 Genes by NES across PTSD States"
)

ggsave('heatmap_all_states.png', 
       height = 9, width = 7)

