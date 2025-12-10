library(dplyr)
library(tidyverse)

indir="/Shares/down/public/omics_databases/GTEX_hg38_RNA/RNA_original_files/"
outdir="/Shares/down/public/omics_databases/GTEX_hg38_RNA/metadata/"



#metadata
#SampleAttributes_metadatafilename = "GTEx_Analysis_v10_Annotations_SampleAttributesDD.xlsx" #explains col names
SampleAttributes_datafilename = paste(indir,"GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt", sep="")
#SubjectPhenotypes_metadatafilename = "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDD.xlsx"#explains col names
SubjectPhenotypes_datafilename = paste(indir,"GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt", sep="")
metatdataPhenotypes <- read.csv(SubjectPhenotypes_datafilename, sep="\t") # more in the protected data
metatdataAttributes <- read.csv(SampleAttributes_datafilename, sep="\t")

# create a column in metatdataAttributes that is the name of the data for that sample in in the GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz or GTEx_Analysis_v10_RNASeQCv2.4.2_gene_reads.gct.gz
metatdataAttributes$tpmcolnames <- gsub("-", ".", metatdataAttributes$SAMPID)

#pull out the donor idea so we can link the RNA-expression back to the doner
metatdataAttributes <- metatdataAttributes %>%
  mutate(
    consortium_donor = sub("^(.*?)-(.*?)-.*", "\\1-\\2", SAMPID)
  )

#pull out the aliquot_id
metatdataAttributes <- metatdataAttributes %>%
  mutate(
    aliquot_id = sub("^.*?-SM-(.*)$", "\\1", SAMPID)  # Capture everything after 'SM-'
  )

#pull out the tissue_replicate (unfortunately not every sample id has this)
metatdataAttributes <- metatdataAttributes %>%
  mutate(
    tissue_replicate = ifelse(grepl("-SM-", SAMPID),
                              sub("^.*?-(.*?)-SM-.*$", "\\1", SAMPID),
                              NA)
  )

#gather the samples tissue by removing sampid then the donor name from the consortium_donor id
metatdataAttributes <- metatdataAttributes %>%
  mutate(consortium_donor_tissue = str_remove(SAMPID, "-SM-.*"))
metatdataAttributes <- metatdataAttributes %>%
  mutate(tissue = str_remove(consortium_donor_tissue, paste(consortium_donor,"*")))
metatdataAttributes$tissue <- gsub("-", "", metatdataAttributes$tissue)


#for now we only want the GTEX samples to there own tissue files
dim(metatdataAttributes)
filtered_metatdataAttributes <- metatdataAttributes %>%
  filter(str_detect(SAMPID, "GTEX"))
dim(filtered_metatdataAttributes)
#now we need to filter to the RNA samples
colnames_donor <- filtered_metatdataAttributes %>% select(tpmcolnames, consortium_donor)
unique(filtered_metatdataAttributes$ANALYTE_TYPE)
filtered_metatdataAttributes <- filtered_metatdataAttributes %>% filter(ANALYTE_TYPE=="RNA:Total RNA")
dim(filtered_metatdataAttributes)
unique(filtered_metatdataAttributes$SMAFRZE)
filtered_metatdataAttributes <- filtered_metatdataAttributes %>% filter(SMAFRZE=="RNASEQ")
dim(filtered_metatdataAttributes)

tissues <- unique(filtered_metatdataAttributes$SMTSD)

n_donor_tissue <- numeric(length(tissues))

for (i in seq_along(tissues)){
  atissue <- tissues[i]
  atissue_filtered_metatdataAttributes <- filtered_metatdataAttributes %>% filter(SMTSD==atissue)
  donor_n <- length(unique(atissue_filtered_metatdataAttributes$consortium_donor))
  n_donor_tissue[i]<-donor_n
}

n_donor_tissue_df <- data.frame(tissues, n_donor_tissue)

filename_without_specialcharacters <- function(atissue) {
  atissuefilename <- atissue
  if (grepl(" ", atissuefilename)) {
    atissuefilename <- gsub(" ", "_", atissuefilename)
  }
  if (grepl("-", atissuefilename)) {
    atissuefilename <- gsub("-", "_", atissuefilename)
  }
  if (grepl("\\(", atissuefilename)) {
    atissuefilename <- gsub("\\(", "_", atissuefilename)
  }
  if (grepl("\\)", atissuefilename)) {
    atissuefilename <- gsub("\\)", "_", atissuefilename)
  }
  return(atissuefilename)
}

write.csv(n_donor_tissue_df, file=paste(outdir, "tissuesdf.csv", sep=""))


n_donor_tissue_df <- n_donor_tissue_df %>%
  mutate(prefixfilename = sapply(tissues, filename_without_specialcharacters))

n_donor_tissue_df$n_donor_tissue <- as.numeric(n_donor_tissue_df$n_donor_tissue )
n_donor_tissue_df_min_people <- n_donor_tissue_df %>% filter(n_donor_tissue>9) %>% arrange(tissues)
write.csv(n_donor_tissue_df_min_people,file=paste(outdir, "tissuesdf_mingreaterthan9samples.csv", sep=""))


