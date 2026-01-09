library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(matrixStats)
library(ggrepel)  # for better label placement
library(readr)
library(stringr)

#step0 set the filters I'm going to use
min_disomic_mean_expression=500 #what the the minimun mean expression of disomic? normalized counts I'm willing to keep in any gene in this assay
min_T21_people_with_comorbid=50.8 #20% of the people with T21 #what is the minimum number of people with the comorbidity I'm willing to keep
howmanystd <- 2 #how many std is to high? Mean+x*std == to high
min_n_gene_to_high=50.8 # 20% # we will only calculate hypergenomtirc if there are more than this minimum number of people that are too high
numberofbootstraps = 100


#step 1 read in the variables
cfg <- read.csv("/Users/allenma/AREA_publication_scripts/figure_scripts/file_locations.csv",
                header = FALSE,
                stringsAsFactors = FALSE)
colnames(cfg) <- c("var", "path")  # first = variable name, second = file path


for (i in seq_len(nrow(cfg))) {
  nm  <- cfg$var[i]
  pth <- cfg$path[i]
  
  # skip rows with no name
  if (is.na(nm) || nm == "") next
  
  # allow missing path entries (e.g., Kalisto_tpm_df) to become NA
  assign(nm, if (is.na(pth) || pth == "") NA_character_ else pth, envir = .GlobalEnv)
}


#step 2 read in the comorbidity file, simplewhich, whichT21 and genes file

outdir="/Shares/down/public/INLCUDE_2024/AREA_paper_2025/output_data/hypergeo/"
simplewhichdf <- read.csv(simplewhich, row.names = 1)
comorbiddf = read.csv(HTP_Hakonarson_Disease_MONDO_df, row.names = 1) #n=996
whichT21df <- read.csv(whichT21_participant) #n=996
gtf <- read_tsv(gtf_first_400,
  comment = "#",
  col_names = c("seqname","source","feature","start","end",
                "score","strand","frame","attributes"),
  col_types = "ccccicccc"
)

gtf <- gtf %>%
  mutate(
    gene_id   = str_match(attributes, "gene_id \"([^\"]+)\"")[, 2],
    gene_name = str_match(attributes, "gene_name \"([^\"]+)\"")[, 2],
    gene_type = str_match(attributes, "gene_type \"([^\"]+)\"")[, 2]
  ) %>% filter(gene_type=="protein_coding")


switch_names <- gtf %>% select(gene_id, gene_name) %>% unique()

gtf_chr21 <- gtf %>% filter(seqname=="chr21")

#step 3 read in the normalized counts

normcountsdf <- read.csv(
  HTP_normcounts_patientid,
  row.names = 1,        # use first column as rownames
) #rows are people columns are genenames, need to flip
rownames(normcountsdf) <- normcountsdf$Participant
normcountsdf <- normcountsdf %>% select(-Participant)


normcountsdf <- normcountsdf[ , colMeans(normcountsdf) > min_mean_expression ]


m  <- as.matrix(normcountsdf)
ncdf <- t(m)

# after this:
rownames(ncdf) == colnames(normcountsdf)
colnames(ncdf) == rownames(normcountsdf)

ncdf <- as.data.frame(ncdf)

#step 4 filter so dataframes have only people with RNA
rownames(comorbiddf)<- comorbiddf$Participant
comorbiddf <- comorbiddf %>% select(-Participant) %>% filter(rownames(comorbiddf) %in% colnames(ncdf))
comorbiddf[] <- lapply(comorbiddf, as.logical)
comorbiddf[] <- as.data.frame(lapply(comorbiddf, as.integer))
DisomicParticipantsdf <- simplewhichdf %>% filter(genotype=="D21")
DisomicParticipantsdf_withRNA <- DisomicParticipantsdf %>% filter(Participant %in% colnames(ncdf))
CompleteT21Participantsdf <- simplewhichdf %>% filter(genotype=="T21")
CompleteT21Participantsdf_withRNA <- CompleteT21Participantsdf %>% filter(Participant %in% colnames(ncdf))
comorbiddf_T21 <- comorbiddf %>% filter(rownames(comorbiddf) %in% CompleteT21Participantsdf$Participant)
dim(comorbiddf_T21)
comorbiddf_T21 <- comorbiddf_T21[ , colSums(comorbiddf_T21) >= min_T21_people_with_comorbid]
dim(comorbiddf_T21)
ncdf_D21 <- ncdf %>% select(DisomicParticipantsdf_withRNA$Participant)
ncdf_T21 <- ncdf %>% select(CompleteT21Participantsdf_withRNA$Participant)


ncdf_T21sim_not21 <- ncdf_D21 %>%
  filter(!rownames(ncdf_D21) %in% gtf_chr21$gene_id)
ncdf_T21sim_21 <- ncdf_D21 %>%
  filter(rownames(ncdf_D21) %in% gtf_chr21$gene_id)
ncdf_T21sim_21 <- ncdf_T21sim_21 * 1.5

ncdf_T21sim <- rbind(ncdf_T21sim_not21, ncdf_T21sim_21)
colnames(ncdf_T21sim)<-paste0("sim",colnames(ncdf_T21sim))

#step 5 find the mean and std for each gene in disomic individals

D21mean <- rowMeans(ncdf_D21) #apply(ncdf_D21, 1, mean, trim=0.1)
D21st <- rowSds(as.matrix(ncdf_D21))
D21expectation <- cbind(rownames(ncdf_D21), D21mean, D21st)
colnames(D21expectation)<-c("ENSEMBL", "D21mean", "D21st")
D21expectation <- as.data.frame(D21expectation)
D21expectation$D21mean <- as.numeric(as.character(D21expectation$D21mean))
D21expectation$D21st <- as.numeric(as.character(D21expectation$D21st))
D21expectation$tohigh <- D21expectation$D21mean+(howmanystd*D21expectation$D21st)
D21expectation$tolow <- D21expectation$D21mean-(howmanystd*D21expectation$D21st)
D21expectation$CV <- D21expectation$D21st/D21expectation$D21mean

#I'm here
#step 6 use those to make a Zscore expression dataframe where each gene is Z scored by the disomic mean and std

ncdf_all <- ncdf_D21  |> rownames_to_column("id") |>
  full_join(ncdf_T21    |> rownames_to_column("id"), by = "id") |>
  full_join(ncdf_T21sim |> rownames_to_column("id"), by = "id") |>
  column_to_rownames("id")

zscoreallfromdisomic <- ncdf_all
zscoreallfromdisomic$ENSEMBL <- rownames(zscoreallfromdisomic)
zscoreallfromdisomic <- merge(zscoreallfromdisomic, D21expectation, on="ENSEMBL")
zscoreallfromdisomic <- zscoreallfromdisomic %>%
  mutate_at(vars(-D21mean,-D21st, -tohigh, -tolow, -CV, -ENSEMBL), funs((.-D21mean)/D21st ))
rownames(zscoreallfromdisomic)<-zscoreallfromdisomic$ENSEMBL
zscoreallfromdisomic <- zscoreallfromdisomic %>% dplyr::select(-D21mean,-D21st, -tohigh, -tolow, -CV, -ENSEMBL)
zscoreallfromdisomic <- as.data.frame(t(zscoreallfromdisomic))
zscoreallfromdisomic$Participant <- rownames(zscoreallfromdisomic)
#rows are Participants, columns are genes excpet the Participant column



# step 7 How many people are to high for each group?


peoplegenestohighdf <- as.data.frame((zscoreallfromdisomic > howmanystd) * 1L)

gene_ids_chr21 <- unique(gtf_chr21$gene_id)
peoplegenestohighdf_chr21 <- peoplegenestohighdf %>% select(any_of(gene_ids_chr21))

peoplegenestohighdf_D21 <- peoplegenestohighdf_chr21 %>% filter(rownames(peoplegenestohighdf) %in% colnames(ncdf_D21))
peoplegenestohighdf_T21  <- peoplegenestohighdf_chr21 %>% filter(rownames(peoplegenestohighdf) %in% colnames(ncdf_T21))
peoplegenestohighdf_simT21  <- peoplegenestohighdf_chr21 %>% filter(rownames(peoplegenestohighdf) %in% colnames(ncdf_T21sim))

peopletohighcounts <- as.data.frame(cbind(colSums(peoplegenestohighdf_D21)/nrow(peoplegenestohighdf_D21), colSums(peoplegenestohighdf_T21)/nrow(peoplegenestohighdf_T21), colSums(peoplegenestohighdf_simT21)/nrow(peoplegenestohighdf_simT21)))
colnames(peopletohighcounts) <- c("D21", "T21", "simT21")
rownames(peopletohighcounts) <-colnames(peoplegenestohighdf_chr21)

#the real T21 have more people than the real D21 (which is also the same as the number of sim T21 since that's how I made them)
#So I can calculate a range for the percentage of people we think will be to high


# 10 replicates of: sample 96 rows, take colSums
sample_sums_mat <- replicate(
  numberofbootstraps,
  {
    sampled <- peoplegenestohighdf_T21[sample(nrow(peoplegenestohighdf_T21), nrow(peoplegenestohighdf_D21)), , drop = FALSE]
    colSums(sampled)/nrow(sampled)
  }
)

# turn into data.frame and give informative column names
sample_sums_df <- as.data.frame(sample_sums_mat)
colnames(sample_sums_df) <- paste0("T21_96personsample", seq_len(ncol(sample_sums_df)))

# bind to existing counts (rows must correspond to same participants/genes as colSums)
peopletohighcounts <- cbind(peopletohighcounts, sample_sums_df)

peopletohighcounts$gene_id <- rownames(peopletohighcounts)
peopletohighcounts <- merge(peopletohighcounts, switch_names, on="gene_id")

peopletohighcounts <- peopletohighcounts %>%
  rowwise() %>%
  mutate(
    T21_sampling_mean = mean(c_across(starts_with("T21_96personsample"))),
    T21_sampling_sd   = sd(c_across(starts_with("T21_96personsample"))),
    T21_sampling_min   = min(c_across(starts_with("T21_96personsample"))),
    T21_sampling_max   = max(c_across(starts_with("T21_96personsample")))
  ) %>%
  ungroup()

D21expectation$gene_id <- rownames(D21expectation)

peopletohighcounts <- merge(peopletohighcounts, D21expectation, on="gene_id")


x_min <- peopletohighcounts$T21_sampling_min
x_max <- peopletohighcounts$T21_sampling_max

x_min

ggplot(peopletohighcounts) +
  # vertical band (still based on T21 sampling stats)
  annotate(
    "rect",
    xmin = x_min,
    xmax = x_max,
    ymin = peopletohighcounts$simT21 - 0.01,
    ymax = peopletohighcounts$simT21 + 0.01,
    fill = "lightblue", alpha = 0.2
  ) +
  # T21 points (red)
  geom_point(aes(x = T21, y = simT21), color = "red") +
  # D21 points (blue)
  geom_point(aes(x = D21, y = simT21), color = "blue" ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_bw() +
  labs(
    x = "percentage of individuals with T21 or D21 too high for this gene",
    y = "percentage of simulated T21 individuals too high for this gene"
  )


peopletohighcounts$log_D21mean <- log(peopletohighcounts$D21mean)

p<-ggplot(peopletohighcounts) +
  # vertical band (still based on T21 sampling stats)
  annotate(
    "rect",
    xmin = x_min,
    xmax = x_max,
    ymin = peopletohighcounts$simT21 - 0.01,
    ymax = peopletohighcounts$simT21 + 0.01,
    fill = "lightblue", alpha = 0.2
  ) +
  # T21 points (red)
  geom_point(aes(x = T21, y = simT21, alpha = log_D21mean), color = "red") +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", size = 1) + # Adds the 1:1 line
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_bw() +
  labs(
    x = "percentage of individuals with T21 too high for this gene",
    y = "percentage of simulated T21 individuals too high for this gene"
  )

# save as SVG
ggplot2::ggsave(
  filename = paste0(outdir,"percentofpeople_",howmanystd,"tohigh_T21_simT21.svg"),
  plot     = p,
  width    = 6,  # in inches; adjust as needed
  height   = 6
)

p


# step 8 count stuff for hypergeometric
comorbiddf_T21[1:4,1:6]
peoplegenestohighdf_T21[1:4,1:6]

dim(comorbiddf_T21)
dim(peoplegenestohighdf_T21)

# sums of each comorbidity column
comorb_sums   <- as.data.frame(colSums(comorbiddf_T21))
colnames(comorb_sums)<-c("people_with_comorbid")
comorb_sums$comorbidity <- rownames(comorb_sums)

# sums of each gene column
peopletohigh_sums   <- as.data.frame(colSums(peoplegenestohighdf_T21))
colnames(peopletohigh_sums)<-c("people_to_high_for_gene")
peopletohigh_sums$gene <- rownames(peopletohigh_sums)

# ensure same row order in both data frames
stopifnot(identical(rownames(comorbiddf_T21),
                    rownames(peoplegenestohighdf_T21)))

# matrix of counts: rows = comorbidities, cols = genes
overlap_counts <- t(as.matrix(comorbiddf_T21)) %*% as.matrix(peoplegenestohighdf_T21)

overlap_long <- as.data.frame.table(overlap_counts,
                                    responseName = "overlap")

colnames(overlap_long) <- c("comorbidity", "gene", "overlap")

overlap_long <- merge(overlap_long,peopletohigh_sums, on="gene")
overlap_long <- merge(overlap_long,comorb_sums, on="comorbidity")


# step 9 hypergeometic

N_total = nrow(peoplegenestohighdf_T21)
  
overlap_long$p_hyper <- with(
  overlap_long,
  phyper(
    q  = overlap - 1,
    m  = people_to_high_for_gene,
    n  = N_total - people_to_high_for_gene,
    k  = people_with_comorbid,
    lower.tail = FALSE
  )
)



overlap_long$padj_BH <- p.adjust(overlap_long$p_hyper, method = "BH")




table(overlap_long$padj_BH)


Rcan1="ENSG00000159200.17"

overlap_long %>% filter(gene==Rcan1) %>% filter(padj_BH!=1.0)

dim(overlap_long)
length(unique(overlap_long$gene))
length(unique(overlap_long$comorbidity))

overlap_long %>% filter(padj_BH<0.2)
