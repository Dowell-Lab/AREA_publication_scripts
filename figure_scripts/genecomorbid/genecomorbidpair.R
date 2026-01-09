library(tidyverse)
library(ggplot2)
library(cowplot)

# Input directories
indir <- "/Shares/down/public/INLCUDE_2024/kallisto_20241030/selfannoated/"
geneindir <- "/Shares/down/public/INLCUDE_2024/kallisto_20241030/metadata/"
outdir="/Users/allenma/graphsfortalk/"

# Files
count_file <- paste0(indir, "kallisto_200401lines_participants_normcounts.csv")
genefile <- paste0(geneindir, "gencode.v27.chr_patch_hapl_scaff.annotation.gtf.df")
disorders_file <- paste0(indir, "full_HP_binary_attribute.csv")
whichT21file <- paste0(indir, "whichT21_participant.csv")
# disorders <- paste0(indir, "full_MONDO_binary_attribute.csv")
disorders <- paste0(indir, "full_HP_binary_attribute.csv")


# Input directories
indir <- "/Shares/down/public/INLCUDE_2024/kallisto_20241030/selfannoated/"
geneindir <- "/Shares/down/public/INLCUDE_2024/kallisto_20241030/metadata/"
outdir="/Users/allenma/hypergeomad/"



#read in files
genes_df <- read.csv(genefile)
disorders_df <- read.csv(disorders_file, row.names=1)
whichT21_df <- read.csv(whichT21file)
counts_df_ori <- read.csv(count_file, row.names=1)

#transform dataframes
row.names(counts_df_ori) <- counts_df_ori$Participant
counts_df <- counts_df_ori %>% select(-Participant)
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
DS_T21_Participants <- whichT21_df %>% filter(MONDO_Down_syndrome=="True")

print("DS_T21_Participants")
print(dim(DS_T21_Participants))
# Add a column to each participant dataframe indicating genotype
D21_Participants$genotype <- "D21"
T21_Participants$genotype <- "T21"
mosaic_T21_Participants$genotype <- "mosaic_T21"
translocation_T21_Participants$genotype <- "translocation_T21"
#mosaic_translocation_T21_Participants$genotype <- "mosaic_translocation_T21"
DS_T21_Participants$genotype <- "DS_T21"

lookup <- bind_rows(
  D21_Participants,
  T21_Participants,
  mosaic_T21_Participants,
  translocation_T21_Participants,
  #  mosaic_translocation_T21_Participants,
  DS_T21_Participants
)

simple<- lookup %>% select(Participant, genotype)
write.csv(simple, paste0(indir,"simple_which_T21.csv"))
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


lookup_limin <- lookup %>% filter(genotype %in% c("D21", "T21")) %>% select(Participant, genotype)
lookup_liminT21 <- lookup %>% filter(genotype %in% c("T21")) %>% select(Participant, genotype)


oneset <- function(ensmbl, boolattribute){
  #counts_df, lookup_limin, disorders
  if (length(ensmbl)>1){ensmbl=ensmbl[1]}
  onegenecounts <- counts_df_ori %>% select(all_of(c("Participant", ensmbl)))
  colnames(onegenecounts) <- c("Participant", "normcounts") 
  one_disorder_df <- disorders_df %>% select(Participant, all_of(boolattribute))
  colnames(one_disorder_df) <- c("Participant", "disorder") 
  one_disorder_gene_counts <- merge(onegenecounts, one_disorder_df, by="Participant")
  one_disorder_gene_counts <- merge(one_disorder_gene_counts, lookup_limin, on="Participant")
  return(one_disorder_gene_counts)
}


#onegene = "RCAN1"
#boolattribute = "HP_Obstructive_sleep_apnea"

#onegene = "C2"
#boolattribute = "HP_Celiac_disease"

onegene = "CDKN2D"
boolattribute = "HP_Recurrent_otitis_media"


#onegene = "PRMT2"
#T21 alone
#boolattribute = "HP_Patent_foramen_ovale"
#boolattribute = "HP_Abnormality_of_the_eye"

#all people
#boolattribute = "HP_Abnormal_heart_morphology"
#boolattribute = "HP_Hypothyroidism"
#boolattribute = "HP_Patent_foramen_ovale"
#boolattribute = "HP_Hypotonia"
#boolattribute = "HP_Obstructive_sleep_apnea"
#boolattribute = "HP_Global_developmental_delay"
#boolattribute = "HP_Sleep_apnea"

#negative controls
#boolattribute = "HP_Abnormality_of_the_skeletal_system"
#boolattribute = "HP_Autoimmunity"
#boolattribute = "HP_Conductive_hearing_impairment"
#boolattribute = "HP_Abnormality_of_the_respiratory_system"

#onegene = "ATG9B"
#boolattribute = "HP_Hypothyroidism"


#onegene = "BRWD1"
#boolattribute = "HP_Hypotonia"


#onegene = "P2RX7"
#boolattribute = "HP_Gastroesophageal_reflux"

#onegene = "PRAG1"
#boolattribute = "HP_Ventricular_septal_defect"

#onegene = "PCSK7"
#boolattribute = "HP_Recurrent_otitis_media"

#onegene = "LINC01547"
#boolattribute = "HP_Hypothyroidism"


#onegene = "USP25"
#boolattribute = "HP_Hypothyroidism"


#onegene = "LSS"
#boolattribute = "HP_Abnormality_of_the_Eustachian_tube"

#onegene = "CU639417.2"
#boolattribute = "HP_Depression"

#onegene = "AL161785.1"
#onegene = "BRWD1"
#onegene = "PAXBP1"
#onegene = "SON"
#boolattribute = "HP_Abnormality_of_the_skeletal_system"


#onegene = "SAMSN1"
#boolattribute = "HP_Obesity"


onegenedf <-   genes_df %>% filter(gene_name==onegene) %>% select(chrom, gene_id, gene_name) %>% unique()
ensmbl <- onegenedf$gene_id

gene_disease_pair <- oneset(ensmbl, boolattribute)

gene_disease_pair_limit <- gene_disease_pair %>% filter(Participant %in% lookup_limin$Participant)
gene_disease_pair_limitT21 <- gene_disease_pair %>% filter(Participant %in% lookup_liminT21$Participant)



str(gene_disease_pair_limit)


plot2 <- ggplot(gene_disease_pair_limitT21, aes(y=normcounts, x=disorder, fill=disorder))+geom_violin()+geom_boxplot(width=0.1, color="black", alpha=0.2) +
  scale_fill_manual(values=c("False"="pink","True"="red"))+
  theme_bw()+
  labs(y = " ")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab(boolattribute)+ theme(legend.position="none")+geom_jitter(alpha=0.1)+
  theme(text = element_text(size = 20))

#filename <- paste0(outdir,"graphs/", onegene, "_and_", boolattribute, ".pdf")
#print(filename)
#ggsave(filename, device = "pdf")

plot1 <- ggplot(gene_disease_pair_limit, aes(y=normcounts, x=genotype, fill=genotype))+geom_violin()+geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_manual(values=c("D21"="black","T21"="darkred"))+
  labs(y = paste(onegene, " normcounts"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("genotype")+ theme(legend.position="none")+geom_jitter(alpha=0.05)+
  theme(text = element_text(size = 30))


#filename <- paste0(outdir,"graphs/", onegene, "_D21_T21",".pdf")
#ggsave(filename, device = "pdf")

# optional: add a spacer above plot2 if you still want it slightly higher

left_side <- plot_grid(
  NULL,   # invisible spacer
  plot1, 
  ncol = 1, 
  rel_heights = c(0.3, 1) # 0.0 means almost no spacer, or increase to push up
)

right_side <- plot_grid(
  NULL,   # invisible spacer
  plot2, 
  ncol = 1, 
  rel_heights = c(0.0, 1)  # 0.0 means almost no spacer, or increase to push up
)

# combine left (smaller plot1) and right (normal plot2)
combined <- plot_grid(
  left_side,
  right_side,
  nrow = 1, 
  rel_widths=c(0.6,1)
)


combined

# save
filename <- paste0(outdir,"graphs/","combined_", onegene, "_", boolattribute, ".pdf")
ggsave(filename, combined, width = 10, height = 5)


