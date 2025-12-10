#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
renv::init()
#renv::install("bioc::recount")
#renv::install("bioc::recount3")
library(recount3)
library(tidyverse)

#this is getting all the human projects
#human_projects <- available_projects()

#this fileters by a minimum number of samples
#human_projects_largen <- human_projects %>% filter(n_samples>300)
#I ended up using the online version to download metadata for all samples
#so many of the projects with more than 300 samples were single cell samples and I Couldn't
#find a way to get rid of them

#found a project with lots of samples and marked for ptsd
#proj_info <- subset(
#  human_projects,
#  project == "SRP102999" & project_type == "data_sources"
#)


#this pulls down all the data
#rse <- create_rse(proj_info) #63856 genes
rse <- create_rse(proj_info, type="gene", annotation="refseq") #54042 genes

#saveRDS(rse, file = "my_rse_object.rds")
rse <- readRDS("my_rse_object.rds")

#the only kinds of counts are the raw counts per SRR
assayNames(rse)

#Lots of metatdata, but I can't find the ptsd column
recount3_cols <- colnames(colData(rse))

recount3_cols_groups <- table(gsub("\\..*", "", recount3_cols))
recount3_cols[grepl("ptsd", recount3_cols, ignore.case = TRUE)]

#I just pulled the metadata off of 
#https://www.ncbi.nlm.nih.gov/Traces/study/
metadatafile = "/Users/maryallen/Downloads/SRP102999_metadata.csv"

metadf <- read.csv(metadatafile)
metadf$aprox_seq_total_reads <- metadf$Bases/metadf$AvgSpotLen
metadf_mini <- metadf %>% select(Run, ptsd, aprox_seq_total_reads)

raw_counts <-assay(rse, "raw_counts")

lib_sizes <- colSums(raw_counts) / 1e6


gene_namedf <- rowData(rse)

gene_length <- rowData(rse)$bp_length
gene_length_kb <- gene_length / 1000

#reads in genes per sample, and samples marked with ptsd to make sure lib sizes are about the same
reads_in_genes <- as.data.frame(colSums(raw_counts))
colnames(reads_in_genes) <- c("total_reads_in_gene_counts")
reads_in_genes$Run <- rownames(reads_in_genes)
reads_in_genes <- merge(metadf_mini,reads_in_genes, by="Run")

ggplot(reads_in_genes, aes(y=total_reads_in_gene_counts, x=ptsd))+geom_violin()+geom_jitter()
ggplot(reads_in_genes, aes(y=aprox_seq_total_reads, x=ptsd))+geom_violin()+geom_jitter()


#how many samples of each type are there
table(reads_in_genes$ptsd)

#make Rpkm and tpm

#sweep is dividing all rows by the gene length in kb
rate <- sweep(raw_counts, 1, gene_length_kb, "/")

# Compute total reads per sample in matrix of counts that have been devided by gene length
scaling_factors <- colSums(rate)

# Normalize to TPM
#sweep is dividing all columns by the total counts in that column
tpm <- sweep(rate, 2, scaling_factors, "/") * 1e6

# Step 6: Store TPM in the RSE object
assays(rse)$TPM <- tpm

## Compute rpkm
rpkm <- t(t(raw_counts) / lib_sizes)  # divide each column by its lib size
rpkm <- rpkm / gene_length        # divide each row by gene length
rpkm <- rpkm * 1e9                # scale to RPKM

assays(rse)$RPKM <- rpkm

reads_in_genes$Past <- as.integer(reads_in_genes$ptsd == "Past")
reads_in_genes$Never <- as.integer(reads_in_genes$ptsd == "Never")
reads_in_genes$Current <- as.integer(reads_in_genes$ptsd == "Current")

boolean_file <- reads_in_genes %>% select(Run, Past, Never, Current)

outdir="/Users/maryallen/forarea/ptsd/"
write.csv(boolean_file, file=paste(outdir, "ptsd.boolean.csv", sep=""))

rank_file = t(tpm)
rank_file <- as.data.frame(rank_file)
rank_file$Run <- rownames(rank_file)


write.csv(rank_file, file=paste(outdir, "geneexptpm.rank.csv", sep=""))


write.csv(gene_namedf , file=paste(outdir, "ptsd.genedf.csv", sep=""))


#How many genes have 0 reads in anything
total_gene_tpm = as.data.frame(colSums(rank_file))
colnames(total_gene_tpm) <- c("sum_tpm")
table(total_gene_tpm$sum_tpm)
#not that many ~8000 our of 60000, AREA will filter that



min(metadf$Bases)/100
median(metadf$Bases)/100
max(metadf$Bases)/100
#> min(metadf$Bases)/100
#[1] 4595253
#> median(metadf$Bases)/100
#[1] 15737647
#> max(metadf$Bases)/100
#[1] 52447557

#> min(reads_in_genes$total_reads)
#[1] 183549808
#> median(reads_in_genes$total_reads)
#[1] 702280466
#> max(reads_in_genes$total_reads)
#[1] 2224928083

tpm_Medians_per_gene <- as.data.frame(rowMedians(tpm))
colnames(tpm_Medians_per_gene) <- c("tpm_per_gene")


ggplot(tpm_Medians_per_gene, aes(x=tpm_per_gene))+geom_histogram(bins=100)+ scale_x_log10(
  breaks = scales::trans_breaks("log10", function(x) 10^x, n = 10),
  labels = scales::trans_format("log10", scales::math_format(10^.x))
) +
  theme_minimal()

tpm_means_per_gene <- as.data.frame(rowMeans(tpm))
colnames(tpm_means_per_gene) <- c("tpm_per_gene")



ggplot(tpm_means_per_gene, aes(x=tpm_per_gene))+geom_histogram(bins=100)+ scale_x_log10()

