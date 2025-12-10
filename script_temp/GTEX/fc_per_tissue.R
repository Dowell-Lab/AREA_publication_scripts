library(dplyr)
library(tidyverse)
library("DESeq2")



print("Running collect a tissue script")

args <- commandArgs(trailingOnly = TRUE)

# Ensure there are exactly four arguments
if(length(args) != 4) {
  print(args)
  stop("Please provide exactly four arguments.")
}

print(args)

# Assign the arguments to variables
tissueprefix <- args[1]
tpmorreads <- args[2]
indir=args[3]
outdir=args[4]

#tissueprefix<- "Whole_Blood"
#tpmorreads <- "tpm"
#indir="/Shares/down/public/omics_databases/GTEX_hg38_RNA/RNA_per_tissue/"
#outdir="/Shares/down/public/omics_databases/GTEX_hg38_RNA/RNAfc_per_tissue/"
print("R thinks the varriables are:")
print(tpmorreads)
print(tissueprefix)
print(indir)
print(outdir)


infile <- paste(indir,tissueprefix , "_", tpmorreads, "_countstable.csv", sep="")

normcounts <- read.csv(infile,row.names = 1)
#I am here


row_stats <- t(apply(normcounts, 1, function(x) {
  # Calculate the min, max, 10th and 90th percentiles, 25th, 50th, and 75th percentiles,
  # count of 0s, and the standard deviation
  c(
    min = min(x),
    max = max(x),
    mean = mean(x),
    median = median(x),
    `10th_percentile` = quantile(x, probs = 0.1),
    `90th_percentile` = quantile(x, probs = 0.9),
    `25th_percentile` = quantile(x, probs = 0.25),
    `50th_percentile` = quantile(x, probs = 0.50),
    `75th_percentile` = quantile(x, probs = 0.75),
    count_zeros = sum(x == 0),
    std_dev = sd(x)  # Add standard deviation
  )
}))

row_stats <- as.data.frame(row_stats)
row_stats$percent_sample_0 <- row_stats$count_zeros/ncol(normcounts)
row_stats$CV <- row_stats$std_dev/row_stats$mean

outfilename = paste(outdir,tissueprefix , "_", tpmorreads, "_gene_population_quantiles.csv", sep="")

write.csv(row_stats, file=outfilename)
