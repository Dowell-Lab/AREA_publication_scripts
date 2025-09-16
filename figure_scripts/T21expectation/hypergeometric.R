library(tidyverse)
library(ggplot2)

indir="/Users/allenma/hypergeomad/"
countsdir <- "/Shares/down/public/INLCUDE_2024/kallisto_20241030/selfannoated/"
geneindir <- "/Shares/down/public/INLCUDE_2024/kallisto_20241030/metadata/"

# Files
count_file <- paste0(countsdir, "kallisto_200401lines_participants_normcounts.csv")
genefile <- paste0(geneindir, "gencode.v27.chr_patch_hapl_scaff.annotation.gtf.df")
disorders_file <- paste0(countsdir, "full_HP_binary_attribute.csv")
whichT21file <- paste0(countsdir, "whichT21_participant.csv")
# disorders <- paste0(countsdir, "full_MONDO_binary_attribute.csv")
disorders <- paste0(countsdir, "full_HP_binary_attribute.csv")

outdir="/Users/allenma/hypergeomad/"

#read in files
genes_df <- read.csv(genefile)
disorders_df <- read.csv(disorders_file, row.names=1)
whichT21_df <- read.csv(whichT21file)
counts_df <- read.csv(count_file, row.names=1)

#subselect participant types
D21_Participants <- whichT21_df %>% filter(Disomic=="True") %>% filter(Participant %in% counts_df$Participant)
print("D21")
print(dim(D21_Participants))
T21_Participants <- whichT21_df %>% filter(MONDO_complete_trisomy_21=="True") %>% filter(Participant %in% counts_df$Participant)
print("T21")
print(dim(T21_Participants))
mosaic_T21_Participants <- whichT21_df %>% filter(MONDO_mosaic_trisomy_21=="True") %>% filter(Participant %in% counts_df$Participant)
print("mosaic_T21_Participants")
print(dim(mosaic_T21_Participants))
translocation_T21_Participants <- whichT21_df %>% filter(MONDO_translocation_Down_syndrome=="True") %>% filter(Participant %in% counts_df$Participant)
print("translocation_T21_Participants")
print(dim(translocation_T21_Participants))
#mosaic_translocation_T21_Participants <- whichT21_df %>% filter(MONDO_mosaic_translocation_Down_syndrome=="True") %>% filter(Participant %in% counts_df$Participant) 
#print("mosaic_translocation_T21_Participants")
#print(dim(mosaic_translocation_T21_Participants))
DS_T21_Participants <- whichT21_df %>% filter(MONDO_Down_syndrome=="True") %>% filter(Participant %in% counts_df$Participant) 
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


#Eupliod_summary_per_gene.csv  Zscoredf_per_gene_per_person.csv
D21summary = read.csv(paste0(indir, "Eupliod_summary_per_gene.csv"))
Zscoredf = read.csv(paste0(indir, "Zscoredf_per_gene_per_person.csv"))
