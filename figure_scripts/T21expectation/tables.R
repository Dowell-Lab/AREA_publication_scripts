library(tidyverse)
library(ggplot2)

# Input directories
indir <- "/Shares/down/public/INLCUDE_2024/kallisto_20241030/selfannoated/"
geneindir <- "/Shares/down/public/INLCUDE_2024/kallisto_20241030/metadata/"
outdir="/Users/allenma/hypergeomad/"

# Files
count_file <- paste0(indir, "kallisto_200401lines_participants_normcounts.csv")
genefile <- paste0(geneindir, "gencode.v27.chr_patch_hapl_scaff.annotation.gtf.df")
disorders_file <- paste0(indir, "full_HP_binary_attribute.csv")
whichT21file <- paste0(indir, "whichT21_participant.csv")
# disorders <- paste0(indir, "full_MONDO_binary_attribute.csv")
disorders <- paste0(indir, "full_HP_binary_attribute.csv")




#read in files
genes_df <- read.csv(genefile)
disorders_df <- read.csv(disorders_file, row.names=1)
whichT21_df <- read.csv(whichT21file)
counts_df <- read.csv(count_file, row.names=1)

#transform dataframes
row.names(counts_df) <- counts_df$Participant
counts_df <- counts_df %>% select(-Participant)
counts_dfT <- t(counts_df)
counts_dfT <- as.data.frame(counts_dfT)

disorders_df_withRNA <- disorders_df %>% filter(Participant %in% colnames(counts_dfT))

#subselect participant types
D21_Participants <- whichT21_df %>% filter(Disomic=="True") %>% filter(Participant %in% colnames(counts_dfT))
print("D21")
print(dim(D21_Participants))
T21_Participants <- whichT21_df %>% filter(MONDO_complete_trisomy_21=="True") %>% filter(Participant %in% colnames(counts_dfT))
print("T21")
print(dim(T21_Participants))
mosaic_T21_Participants <- whichT21_df %>% filter(MONDO_mosaic_trisomy_21=="True") %>% filter(Participant %in% colnames(counts_dfT))
print("mosaic_T21_Participants")
print(dim(mosaic_T21_Participants))
translocation_T21_Participants <- whichT21_df %>% filter(MONDO_translocation_Down_syndrome=="True") %>% filter(Participant %in% colnames(counts_dfT))
print("translocation_T21_Participants")
print(dim(translocation_T21_Participants))
#mosaic_translocation_T21_Participants <- whichT21_df %>% filter(MONDO_mosaic_translocation_Down_syndrome=="True") %>% filter(Participant %in% colnames(counts_dfT)) 
#print("mosaic_translocation_T21_Participants")
#print(dim(mosaic_translocation_T21_Participants))
DS_T21_Participants <- whichT21_df %>% filter(MONDO_Down_syndrome=="True") %>% filter(Participant %in% colnames(counts_dfT)) 
print("DS_T21_Participants")
print(dim(DS_T21_Participants))
# Add a column to each participant dataframe indicating genotype
D21_Participants$genotype <- "D21"
T21_Participants$genotype <- "T21"
mosaic_T21_Participants$genotype <- "mosaic_T21"
translocation_T21_Participants$genotype <- "translocation_T21"
#mosaic_translocation_T21_Participants$genotype <- "mosaic_translocation_T21"
#DS_T21_Participants$genotype <- "DS_T21"

lookup <- bind_rows(
  D21_Participants,
  T21_Participants,
  mosaic_T21_Participants,
  translocation_T21_Participants,
  #  mosaic_translocation_T21_Participants,
  DS_T21_Participants
)

#find the chromsome 21 genes

chr21genesdf <- genes_df %>% filter(chrom=="chr21")
chr21genes <- unique(chr21genesdf$gene_id)

#count dataframes for people types
D21_counts_df <-  counts_dfT %>% select(D21_Participants$Participant)
T21_counts_df <- counts_dfT %>% select(T21_Participants$Participant)
T21sim_counts_df <- D21_counts_df
T21sim_counts_df[rownames(T21sim_counts_df) %in% chr21genes, ] <- 
  T21sim_counts_df[rownames(T21sim_counts_df) %in% chr21genes, ] * 1.5
dim(T21sim_counts_df)
colnames(T21sim_counts_df) <- paste0("sim-", colnames(T21sim_counts_df))


#create and save D21 summary

D21summary <- data.frame(
  gene_id   = rownames(D21_counts_df),
  mean   = rowMeans(D21_counts_df, na.rm = TRUE),
  median = apply(D21_counts_df, 1, median, na.rm = TRUE),
  sd     = apply(D21_counts_df, 1, sd, na.rm = TRUE),
  mad    = apply(D21_counts_df, 1, mad, na.rm = TRUE),
  n_zeros = rowSums(D21_counts_df == 0, na.rm = TRUE)
)

D21summary$percent_people_0 <- D21summary$n_zeros/ncol(D21_counts_df)

filename = paste0(outdir, "Eupliod_summary_per_gene.csv")
write.csv(D21summary, filename)

#add Zscores to the count dataframes
add_zscore <- function(Zscoredf, D21summary) {
  Zscoredf <- Zscoredf %>%
    left_join(
      D21summary %>% select(gene_id, mean, sd),
      by = c("gene_id" = "gene_id")
    ) %>%
    mutate(zscore = (normcount - mean) / sd) 
  Zscoredf
}

add_modzscore <- function(Zscoredf, D21summary) {
  Zscoredf <- Zscoredf %>%
    left_join(
      D21summary %>% select(gene_id, median, mad),
      by = c("gene_id" = "gene_id")
    ) %>%
    mutate(modzscore = (0.6745*(normcount - median) )/ mad) 
  Zscoredf
}

Zscoredf <- counts_dfT %>%
  rownames_to_column(var = "gene_id") %>%   # move rownames into a column
  pivot_longer(
    cols = -gene_id,                        # all except Gene column
    names_to = "Participant", 
    values_to = "normcount"
  )

Zscoredfsim <- T21sim_counts_df %>%
  rownames_to_column(var = "gene_id") %>%   # move rownames into a column
  pivot_longer(
    cols = -gene_id,                        # all except Gene column
    names_to = "Participant", 
    values_to = "normcount"
  )

Zscoredf <- as.data.frame(rbind(Zscoredf,Zscoredfsim))
Zscoredf <- add_zscore(Zscoredf, D21summary)
Zscoredf <- add_modzscore(Zscoredf, D21summary)
Zscoredf$onchr21 <- Zscoredf$gene_id %in% chr21genes
Zscoredf <- left_join(Zscoredf, lookup %>% select(Participant, genotype), by = c("Participant" = "Participant"))
Zscoredf$genotype[is.na(Zscoredf$genotype)] <- "Unknown"

filename = paste0(outdir, "Zscoredf_per_gene_per_person.csv")
write.csv(Zscoredf, filename)

# count for hypergeometric 
count_bool_zscore <- function(wide_df, long_df, bool_cols, id_col = "Participant", zscore_col = "Zscore", zscore_threshold = 2) {
  # Precompute which IDs are above zscore threshold
  gene_names <- unique(as.character(long_df$gene_id))
  result <- data.frame(
    gene_id = character(0),
    column = character(0),
    with_condition = integer(0),
    count_FALSE = integer(0),
    count_FALSE_above_z = integer(0),
    count_TRUE_above_z = integer(0),
    count_FALSE_below_z = integer(0),
    count_TRUE_below_z = integer(0),
    stringsAsFactors = FALSE
  )
  for (gn in gene_names){
    long_df_onegene <- subset(long_df, gene_id == gn)
    id_above_z <- unique(long_df_onegene[[id_col]][long_df_onegene[[zscore_col]] > zscore_threshold])
    id_below_z <- unique(long_df_onegene[[id_col]][long_df_onegene[[zscore_col]] <= zscore_threshold])
    for (col in bool_cols) {
      bool_vals <- wide_df[[col]]
      ids_FALSE <- wide_df[[id_col]][bool_vals == "False"]
      ids_TRUE <- wide_df[[id_col]][bool_vals == "True"]
      count_FALSE <- length(ids_FALSE)
      count_TRUE <- length(ids_TRUE)
      count_FALSE_above_z <- sum(ids_FALSE %in% id_above_z)
      count_TRUE_above_z <- sum(ids_TRUE %in% id_above_z)
      count_FALSE_below_z <- sum(ids_FALSE %in% id_below_z)
      count_TRUE_below_z <- sum(ids_TRUE %in% id_below_z)
      result <- rbind(result, data.frame(gene_id = gn,
                                         column = col,
                                         count_TRUE = count_TRUE,
                                         count_FALSE = count_FALSE,
                                         count_FALSE_above_z = count_FALSE_above_z,
                                         count_TRUE_above_z = count_TRUE_above_z,
                                         count_FALSE_below_z = count_FALSE_below_z,
                                         count_TRUE_below_z = count_TRUE_below_z,
                                         stringsAsFactors = FALSE
      ))
    }}
  return(result)
}


T21_disorders_df <- disorders_df %>% filter(Participant %in% T21_Participants$Participant)
T21_disorders_df_just_disorders <- T21_disorders_df %>% select(-Participant)
# Convert character columns to logical (TRUE/FALSE)
T21_disorders_df_just_disorders[] <- lapply(T21_disorders_df_just_disorders, function(x) x == "True")

# Count TRUE values per column and only do hypergeometric on disorders with more than a min number of people
n_people_with = colSums(T21_disorders_df_just_disorders)
min_people = 10
n_people_with_filtered <- n_people_with[n_people_with > min_people]
bool_list <- names(n_people_with_filtered)


#smallset <- c(unique(T21_Zscoredf$gene_id)[1:4], c("ENSG00000159200.17"))
#T21_Zscoredf <- Zscoredf %>% filter(Participant %in% T21_Participants$Participant) %>% filter(onchr21==TRUE) %>% filter(gene_id %in% smallset)
T21_Zscoredf <- Zscoredf %>% filter(Participant %in% T21_Participants$Participant) %>% filter(onchr21==TRUE) 


zscore_threshold = 2
zscore_col = "zscore"
for_hypergeometirc_above <- count_bool_zscore(T21_disorders_df, T21_Zscoredf, bool_list, id_col = "Participant", zscore_col = "zscore", zscore_threshold = zscore_threshold)
for_hypergeometirc_above$n_world <- for_hypergeometirc_above$count_TRUE+for_hypergeometirc_above$count_FALSE
for_hypergeometirc_above$above_z <- for_hypergeometirc_above$count_TRUE_above_z+for_hypergeometirc_above$count_FALSE_above_z
for_hypergeometirc_above$hypergpval = phyper(for_hypergeometirc_above$count_TRUE_above_z-1, for_hypergeometirc_above$count_TRUE, for_hypergeometirc_above$n_world-for_hypergeometirc_above$count_TRUE, for_hypergeometirc_above$above_z , lower.tail = FALSE, log.p = FALSE)
for_hypergeometirc_above$expextation = for_hypergeometirc_above$n_world*(for_hypergeometirc_above$count_TRUE/for_hypergeometirc_above$n_world)*(for_hypergeometirc_above$above_z/for_hypergeometirc_above$n_world)
for_hypergeometirc_above<- for_hypergeometirc_above %>% arrange(hypergpval)

filename = paste0(outdir, "T21_4hypergeometric_zscoretype_",zscore_col,"thres",zscore_threshold,".csv")
write.csv(for_hypergeometirc_above, filename)

zscore_threshold = 3
zscore_col = "zscore"
for_hypergeometirc_above <- count_bool_zscore(T21_disorders_df, T21_Zscoredf, bool_list, id_col = "Participant", zscore_col = "zscore", zscore_threshold = zscore_threshold)
for_hypergeometirc_above$n_world <- for_hypergeometirc_above$count_TRUE+for_hypergeometirc_above$count_FALSE
for_hypergeometirc_above$above_z <- for_hypergeometirc_above$count_TRUE_above_z+for_hypergeometirc_above$count_FALSE_above_z
for_hypergeometirc_above$hypergpval = phyper(for_hypergeometirc_above$count_TRUE_above_z-1, for_hypergeometirc_above$count_TRUE, for_hypergeometirc_above$n_world-for_hypergeometirc_above$count_TRUE, for_hypergeometirc_above$above_z , lower.tail = FALSE, log.p = FALSE)
for_hypergeometirc_above$expextation = for_hypergeometirc_above$n_world*(for_hypergeometirc_above$count_TRUE/for_hypergeometirc_above$n_world)*(for_hypergeometirc_above$above_z/for_hypergeometirc_above$n_world)
for_hypergeometirc_above<- for_hypergeometirc_above %>% arrange(hypergpval)

filename = paste0(outdir, "T21_4hypergeometric_zscoretype_",zscore_col,"thres",zscore_threshold,".csv")
write.csv(for_hypergeometirc_above, filename)

zscore_threshold = 2
zscore_col = "modzscore"
for_hypergeometirc_above <- count_bool_zscore(T21_disorders_df, T21_Zscoredf, bool_list, id_col = "Participant", zscore_col = "zscore", zscore_threshold = zscore_threshold)
for_hypergeometirc_above$n_world <- for_hypergeometirc_above$count_TRUE+for_hypergeometirc_above$count_FALSE
for_hypergeometirc_above$above_z <- for_hypergeometirc_above$count_TRUE_above_z+for_hypergeometirc_above$count_FALSE_above_z
for_hypergeometirc_above$hypergpval = phyper(for_hypergeometirc_above$count_TRUE_above_z-1, for_hypergeometirc_above$count_TRUE, for_hypergeometirc_above$n_world-for_hypergeometirc_above$count_TRUE, for_hypergeometirc_above$above_z , lower.tail = FALSE, log.p = FALSE)
for_hypergeometirc_above$expextation = for_hypergeometirc_above$n_world*(for_hypergeometirc_above$count_TRUE/for_hypergeometirc_above$n_world)*(for_hypergeometirc_above$above_z/for_hypergeometirc_above$n_world)
for_hypergeometirc_above<- for_hypergeometirc_above %>% arrange(hypergpval)

filename = paste0(outdir, "T21_4hypergeometric_zscoretype_",zscore_col,"thres",zscore_threshold,".csv")
write.csv(for_hypergeometirc_above, filename)

zscore_threshold = 3
zscore_col = "modzscore"
for_hypergeometirc_above <- count_bool_zscore(T21_disorders_df, T21_Zscoredf, bool_list, id_col = "Participant", zscore_col = "zscore", zscore_threshold = zscore_threshold)
for_hypergeometirc_above$n_world <- for_hypergeometirc_above$count_TRUE+for_hypergeometirc_above$count_FALSE
for_hypergeometirc_above$above_z <- for_hypergeometirc_above$count_TRUE_above_z+for_hypergeometirc_above$count_FALSE_above_z
for_hypergeometirc_above$hypergpval = phyper(for_hypergeometirc_above$count_TRUE_above_z-1, for_hypergeometirc_above$count_TRUE, for_hypergeometirc_above$n_world-for_hypergeometirc_above$count_TRUE, for_hypergeometirc_above$above_z , lower.tail = FALSE, log.p = FALSE)
for_hypergeometirc_above$expextation = for_hypergeometirc_above$n_world*(for_hypergeometirc_above$count_TRUE/for_hypergeometirc_above$n_world)*(for_hypergeometirc_above$above_z/for_hypergeometirc_above$n_world)
for_hypergeometirc_above<- for_hypergeometirc_above %>% arrange(hypergpval)

filename = paste0(outdir, "T21_4hypergeometric_zscoretype_",zscore_col,"thres",zscore_threshold,".csv")
write.csv(for_hypergeometirc_above, filename)
