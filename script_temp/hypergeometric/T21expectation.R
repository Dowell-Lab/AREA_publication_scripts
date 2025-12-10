library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(matrixStats)
library(ggrepel)  # for better label placement




indir = "/Shares/down/public/HTP//RNAseq/"
#whichoutput="T21vsD21_DNAdosagecorrection_match"
#whichoutput="T21vsD21_DNAdosagecorrection"
whichoutput="T21vsD21_noDNAdosagecorrection"
comorbidfileroot = "Patient_mondo_logical.csv"
#comorbidfileroot = "Patient_pheno_logical.csv"
comorbidfile=paste0("/Shares/down/public/HTP/RNAseq/outputdata/", comorbidfileroot, sep="")
allgenesfile = paste(indir, "outputdata/", whichoutput, "/allgeneswithgenenames_res_not_unique.csv", sep="")
normcountsfile =paste(indir, "outputdata/",whichoutput,"/normcounts.csv", sep="")
outdir=paste(indir, "outputdata/",whichoutput,"/furtherprocessing/", sep="")
dir.create(outdir, showWarnings = FALSE)
min_n_gene_to_high=10
min_count_sum=100
min_people_with_comorbid=10

#read in who has which comorbidities
comorbidf <-  read.csv(comorbidfile, row.names = 1)
comorbidf
colnames(comorbidf) <- paste('comorbid', colnames(comorbidf), sep = '_')
comorbidf[] <- lapply(comorbidf, as.logical)
comorbidf$Patient <- gsub("-",".",rownames(comorbidf))

complete_trisomy_21 <- comorbidf %>% dplyr::filter(comorbid_complete_trisomy_21==1)
dim(complete_trisomy_21)
mosaic_trisomy_21 <-  comorbidf %>% dplyr::filter(comorbid_mosaic_trisomy_21==1)
dim(mosaic_trisomy_21)
translocation_Down_syndrome <-  comorbidf %>% dplyr::filter(comorbid_translocation_Down_syndrome==1)
dim(translocation_Down_syndrome)
mosaic_translocation_Down_syndrome <-  comorbidf %>% dplyr::filter(comorbid_mosaic_translocation_Down_syndrome==1)


#read in genes and their chromosomes
geneinfo = "/Shares/down/public/HTP/RNAseq/selfannotated/genes.csv"
genedf = read.csv(geneinfo)
genedf$ENSEMBL = genedf$gene_id 
genedf$SYMBOL = genedf$gene_name 
genedf$CHR = genedf$seqnames 

genesymbolchr <- genedf %>% dplyr::select(ENSEMBL, CHR, SYMBOL) %>% unique()
genesymbolchr_just21 <- genesymbolchr %>% filter(CHR==21)
genesymbolchr_just22 <- genesymbolchr %>% filter(CHR==22)
#ensb gene to symbol not 1 to 1
dim(genesymbolchr)
length(unique(genesymbolchr$ENSEMBL))
dim(genesymbolchr_just21)
length(unique(genesymbolchr_just21$ENSEMBL))


#read in which patients are at all T21 
mani <- read.csv(paste(indir, "inputdata/manifest_20230829_152940.csv", sep=""))
manimini <- mani %>% dplyr::select(Patient, sample_type)

manimini <- manimini %>%
  dplyr::mutate(Patient = Patient %>% 
                  str_replace_all("-", "."))

dim(manimini)
head(manimini)

maniminiD21 <- manimini %>% filter(sample_type=="D21")
maniminiT21 <- manimini %>% filter(sample_type=="T21")


complete_trisomy_21_with_RNA <- complete_trisomy_21$Patient[complete_trisomy_21$Patient %in% maniminiT21$Patient]
mosaic_trisomy_21_with_RNA <- mosaic_trisomy_21$Patient[mosaic_trisomy_21$Patient %in% maniminiT21$Patient]
translocation_Down_syndrome_with_RNA <- translocation_Down_syndrome$Patient[translocation_Down_syndrome$Patient %in% maniminiT21$Patient]
mosaic_translocation_Down_syndrome_with_RNA <- mosaic_translocation_Down_syndrome$Patient[mosaic_translocation_Down_syndrome$Patient %in% maniminiT21$Patient]


#read in conts per gene and make a D21/T21 count data set
ncdf <- read.csv(normcountsfile)
head(ncdf)
rownames(ncdf)<-ncdf$X
ncdf <- ncdf %>% dplyr::select(-X)
ncdfchr21<- ncdf %>% dplyr::filter(rownames(ncdf) %in% genesymbolchr_just21$ENSEMBL)
ncdfchr21_D21 <- ncdf %>% dplyr::filter(rownames(ncdf) %in% genesymbolchr_just21$ENSEMBL) %>% dplyr::select(maniminiD21$Patient)
ncdfnotchr21_D21 <- ncdf %>% dplyr::filter(!rownames(ncdf) %in% genesymbolchr_just21$ENSEMBL) %>% dplyr::select(maniminiD21$Patient)
ncdf_T21<-ncdf %>% dplyr::select(any_of(complete_trisomy_21_with_RNA))
ncdf_D21<-ncdf %>% dplyr::select(any_of(maniminiD21$Patient))
ncdfchr21_T21sim <- ncdfchr21_D21*1.5
ncdf_T21sim <- rbind(ncdfchr21_T21sim, ncdfnotchr21_D21)




D21mean <- rowMeans(ncdf_D21) #apply(ncdf_D21, 1, mean, trim=0.1)
D21st <- rowSds(as.matrix(ncdf_D21))
D21expectation <- cbind(rownames(ncdf_D21), D21mean, D21st)
colnames(D21expectation)<-c("ENSEMBL", "D21mean", "D21st")
D21expectation <- as.data.frame(D21expectation)
D21expectation$D21mean <- as.numeric(as.character(D21expectation$D21mean))
D21expectation$D21st <- as.numeric(as.character(D21expectation$D21st))
#%>% mutate_at(c('D21mean', 'D21st'), as.numeric)
howmanystd <- 2
D21expectation$tohigh <- D21expectation$D21mean+(howmanystd*D21expectation$D21st)
D21expectation$tolow <- D21expectation$D21mean-(howmanystd*D21expectation$D21st)
D21expectation$CV <- D21expectation$D21st/D21expectation$D21mean

zscoreallfromdisomic <- ncdf
zscoreallfromdisomic$ENSEMBL <- rownames(zscoreallfromdisomic)
zscoreallfromdisomic <- merge(zscoreallfromdisomic, D21expectation, on="ENSEMBL")
zscoreallfromdisomic <- zscoreallfromdisomic %>%
  mutate_at(vars(-D21mean,-D21st, -tohigh, -tolow, -CV, -ENSEMBL), funs((.-D21mean)/D21st ))
rownames(zscoreallfromdisomic)<-zscoreallfromdisomic$ENSEMBL
zscoreallfromdisomic <- zscoreallfromdisomic %>% dplyr::select(-D21mean,-D21st, -tohigh, -tolow, -CV, -ENSEMBL)
zscoreallfromdisomic <- as.data.frame(t(zscoreallfromdisomic))
zscoreallfromdisomic$Patient <- rownames(zscoreallfromdisomic)

gene_exp_long <- ncdf 
gene_exp_long$ENSEMBL <- rownames(gene_exp_long)
gene_exp_long <- gene_exp_long %>% pivot_longer(!ENSEMBL, names_to = "Patient", values_to = "normcounts")
gene_exp_long = gene_exp_long %>% filter(ENSEMBL %in% genesymbolchr_just21$ENSEMBL)
gene_mean <- gene_exp_long %>% group_by(ENSEMBL) %>% dplyr::summarize(Mean = mean(normcounts, na.rm=TRUE), Sum=sum(normcounts, na.rm=TRUE))
dim(gene_mean)
gene_mean <- gene_mean %>% filter(gene_mean$Sum>min_count_sum)
dim(gene_mean)



gene_exp_long_exp <- merge(gene_exp_long, D21expectation, on="ENSEMBL")
gene_exp_long_exp$genetohigh <- gene_exp_long_exp$normcounts>gene_exp_long_exp$tohigh
gene_exp_long_exp$genetolow <- gene_exp_long_exp$normcounts<gene_exp_long_exp$tolow

gene_exp_long_exp <- gene_exp_long_exp %>% mutate(sample_type =   case_when(
  Patient %in%  complete_trisomy_21$Patient ~ "complete_T21",
  Patient %in%  mosaic_trisomy_21$Patient ~ "mosaic_trisomy_21",
  Patient %in%  translocation_Down_syndrome$Patient ~ "translocation_Down_syndrome",
  Patient %in%  maniminiT21$Patient ~ "other_T21",
  TRUE ~ "D21"
))

gene_exp_long_exp_complete_T21= gene_exp_long_exp %>% filter(sample_type=="complete_T21") 
gene_exp_long_exp_complete_D21= gene_exp_long_exp %>% filter(sample_type=="D21") 

comorbidf_T21 <- comorbidf %>% filter(Patient %in% colnames(ncdf)) %>% filter(Patient %in% complete_trisomy_21$Patient)

comorbidf_noPatient <- comorbidf %>% filter(Patient %in% colnames(ncdf)) %>% filter(Patient %in% complete_trisomy_21$Patient) %>% select(-Patient)

comorbidcounts = as.data.frame(colSums(comorbidf_noPatient))
colnames(comorbidcounts)<- c("people_with_comorbid")

comorbidcounts <- comorbidcounts %>% filter(people_with_comorbid>=min_people_with_comorbid)
dim(comorbidcounts)

comorbidcounts$comorbid_comorbid <- rownames(comorbidcounts)
comorbidcounts_fixed <-comorbidcounts %>% separate_wider_delim(comorbid_comorbid, delim = "bid_", names = c("junk", "comorbid")) %>% select(people_with_comorbid, comorbid)

n_T21_people_to_high_pergene <- gene_exp_long_exp_complete_T21 %>% 
  group_by(
    ENSEMBL, genetohigh
  ) %>% 
  summarize(
    count = n()
  ) %>% 
  spread(
    key = genetohigh, 
    value = count, 
    fill = 0
  )

n_T21_people_to_high_pergene$people_count <- n_T21_people_to_high_pergene$`TRUE`+n_T21_people_to_high_pergene$`FALSE`

n_T21_people_to_high_pergene$T21_pecent_people_to_high <- n_T21_people_to_high_pergene$`TRUE`/n_T21_people_to_high_pergene$people_count



n_D21_people_to_high_pergene <- gene_exp_long_exp_complete_D21 %>% 
  group_by(
    ENSEMBL, genetohigh
  ) %>% 
  summarize(
    count = n()
  ) %>% 
  spread(
    key = genetohigh, 
    value = count, 
    fill = 0
  )

n_D21_people_to_high_pergene$people_count <- n_D21_people_to_high_pergene$`TRUE`+n_D21_people_to_high_pergene$`FALSE`

n_D21_people_to_high_pergene$D21_pecent_people_to_high <- n_D21_people_to_high_pergene$`TRUE`/n_D21_people_to_high_pergene$people_count


n_D21_people_to_high_pergene <- n_D21_people_to_high_pergene %>% select(ENSEMBL, D21_pecent_people_to_high)
n_T21_people_to_high_pergene <- n_T21_people_to_high_pergene %>% select(ENSEMBL, T21_pecent_people_to_high)

people_to_high_pergene <- merge(n_D21_people_to_high_pergene, n_T21_people_to_high_pergene, on="ENSEMBL")


people_to_high_pergene <- merge(people_to_high_pergene, D21expectation, on="ENSEMBL")

people_to_high_pergene$logmean <- log(people_to_high_pergene$D21mean)

ggplot(people_to_high_pergene, aes(x=D21_pecent_people_to_high, y=T21_pecent_people_to_high, size=logmean)) + geom_point()

ggplot(people_to_high_pergene, aes(x=exp_level, y=T21_pecent_people_to_high)) + geom_boxplot()

ggplot(people_to_high_pergene, aes(x=exp_level, y=T21_pecent_people_to_high, size=CV)) + geom_point()
#Turn the T21_pecent_people_to_high to a percent rather than a decimal
#theme_minimal
#x axis labels need to rotate 90 degrees
#x axis label should be "Mean expression level in disomic indivuals (n=96)
# title should be "Percentage of T21 indivuals that are higher than disomic range (mean+3std)
# size label should be "Coeffecent of variation in Disomics"
#Y axis label should be "Percentage of indivuals with T21 that are outside disomic range"


library(ggplot2)
library(scales)  # for percent formatting on y-axis

# Assuming people_to_high_pergene is your data frame

library(ggplot2)
library(scales)  # for percent formatting on y-axis

# Assuming people_to_high_pergene is your data frame

addnamesdf <-genedf %>% select(gene_name, ENSEMBL)

people_to_high_pergene <- merge(people_to_high_pergene,  addnamesdf, on=ENSEMBL)

p <- ggplot(people_to_high_pergene, aes(x = exp_level, y = T21_pecent_people_to_high, color = CV)) +
  geom_point() +
  scale_color_gradient(trans = "log",  # Apply log transformation to color scale
                       low = "blue", high = "red") +
  scale_y_continuous(labels = scales::percent) +  # Convert decimal to percent on y-axis
  theme_minimal() +  # Use minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels 90 degrees
  labs(
    x = "Gene expression D21 indivuals (n=96)",
    title = "Outside D21 range (mean+3std)",
    color = "Coeffecent of variation in D21",
    y = "Percentage T21 individuals "
  )

print(p)

ggsave(filename = paste(outdir,"Percentageoutsiderange.png", sep=""), plot = p, width = 8, height = 6, dpi = 300)


#label_data <- people_to_high_pergene %>%
#  filter(startsWith(gene_name, "IFN"))
#Final graph for grant

people_to_high_pergene <- people_to_high_pergene %>%
  mutate(CVround = round(CV, 2))

label_data <- people_to_high_pergene %>%
  filter(grepl("^(IFN|RUNX1|DYRK1A)", gene_name))

base_text_size <- 14  # Adjust this number as needed for overall size


text_theme <- theme_minimal(base_size = base_text_size) + 
  theme(
    axis.title = element_text(size = base_text_size + 2),    # axis titles larger
    axis.text = element_text(size = base_text_size),         # axis numbers/text bigger
    plot.title = element_text(size = base_text_size + 4, face = "bold"),  # bigger, bold title
    legend.text = element_text(size = base_text_size),       # if you have legends
    legend.title = element_text(size = base_text_size + 2)
  )

p <- ggplot(people_to_high_pergene, aes(x = exp_level, y = T21_pecent_people_to_high, color = CVround)) +
  geom_point() +
  geom_text_repel(data = label_data, aes(label = gene_name), size = 5, color = "black", nudge_x=1.3) +  # Add labels for IFN genes
  scale_color_gradient(
    trans = "log",
    low = "blue",
    high = "red",
    labels = scales::number_format(accuracy = 0.1)
  ) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(
    x = "Gene expression D21 indivuals (n=96)",
    color = "CV",
    y = "% T21 individuals outside D21 range"
  )+text_theme+theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme(legend.position = "inside",legend.position.inside = c(0.2, 0.7)) 


print(p)

ggsave(filename = paste(outdir,"Percentageoutsiderange2_withIFN_bigger.png", sep=""), plot = p, width = 8, height = 6, dpi = 300)


###finshpergrant
people_to_high_pergene <- people_to_high_pergene %>% 
  mutate(exp_level = cut(D21mean, 
                   c(-Inf, quantile(D21mean, c(0.001,.25, 0.50, .75, 0.90)), Inf), 
                   labels = c("Very Low", "Low", "Medium Low", "Medium High", "High", "Very High")))
table(people_to_high_pergene$exp_level)

ggplot(people_to_high_pergene,aes(x=exp_level, y=T21_pecent_people_to_high))+geom_violin()


ggplot(people_to_high_pergene,aes(x=T21_pecent_people_to_high))+geom_histogram()

#ggplot(people_to_high_pergene,aes(x=cat_CV, y=T21_pecent_people_to_high))+geom_violin()

people_to_high_pergene[people_to_high_pergene$ENSEMBL=="ENSG00000159200"]

mylist <- list() #create an empty list
i=1
for (thisENSEMBL in unique(gene_mean$ENSEMBL)){
  for (thiscomorbid in unique(rownames(comorbidcounts))){
  this_gene_exp_long_exp_complete_T21 <- gene_exp_long_exp_complete_T21 %>% filter(ENSEMBL==thisENSEMBL)
  peopleT21__gene_to_high <- this_gene_exp_long_exp_complete_T21 %>% filter(genetohigh==TRUE)
  peopleT21__gene_not_to_high <- this_gene_exp_long_exp_complete_T21 %>% filter(genetohigh==FALSE)
  peopleT21__thiscomorbid <- comorbidf_T21[comorbidf_T21[,thiscomorbid] == TRUE, ]
  peopleT21__notcomorbid <- comorbidf_T21[comorbidf_T21[,thiscomorbid] == FALSE, ]
  peopleT21__thiscomorbid_gene_to_high <- peopleT21__gene_to_high %>% filter(Patient %in% peopleT21__thiscomorbid$Patient)
  peopleT21__thiscomorbid_genenormal <- peopleT21__gene_not_to_high %>% filter(Patient %in% peopleT21__thiscomorbid$Patient)
  peopleT21__notcomorbid_gene_to_high <- peopleT21__gene_to_high %>% filter(Patient %in% peopleT21__notcomorbid$Patient)
  peopleT21__notcomorbid_genenormal <- peopleT21__gene_not_to_high %>% filter(Patient %in% peopleT21__notcomorbid$Patient)
  n_gene_to_high = nrow(peopleT21__gene_to_high)
  n_gene_normal = nrow(peopleT21__gene_not_to_high)
  n_comorbid = nrow(peopleT21__thiscomorbid)
  n_not_comorbid = nrow(peopleT21__notcomorbid)
  n_gene_to_high_and_comorbid <- nrow(peopleT21__thiscomorbid_gene_to_high)
  n_gene_to_high_and_notcomorbid <- nrow(peopleT21__notcomorbid_gene_to_high)
  n_genenormal_and_comorbid <- nrow(peopleT21__thiscomorbid_genenormal)
  n_ggenenormal_and__notcomorbid <- nrow(peopleT21__notcomorbid_genenormal)
  n_world = n_gene_to_high+n_gene_normal
  if (n_gene_to_high>min_n_gene_to_high){
    hypergpval = phyper(n_gene_to_high_and_comorbid-1, n_comorbid, n_world-n_comorbid, n_gene_to_high, lower.tail = FALSE, log.p = FALSE)
    exp = n_world*(n_comorbid/n_world)*(n_gene_to_high/n_world)
    mylist[[i]] <- c(thisENSEMBL, thiscomorbid, n_world, n_gene_to_high, n_gene_normal, n_comorbid, n_not_comorbid, n_gene_to_high_and_comorbid, exp, hypergpval, n_gene_to_high_and_notcomorbid,n_genenormal_and_comorbid,n_ggenenormal_and__notcomorbid)
    i = i+1
    }
}}

pairs <- do.call("rbind",mylist)

pairsdf = as.data.frame(pairs)
colnames(pairsdf) <- c("ENSEMBL", "comorbid", "n_world", "n_gene_to_high", "n_gene_normal", "n_comorbid", "n_not_comorbid", "n_gene_to_high_and_comorbid", "exp_gene_to_high_and_comorbid", "hypergpval_lowertail","n_gene_to_high_and_notcomorbid","n_genenormal_and_comorbid","n_ggenenormal_and__notcomorbid")

pairsdf$hypergpval_lowertail <- as.numeric(as.character(pairsdf$hypergpval_lowertail))
pairsdf$n_gene_to_high_and_comorbid <- as.numeric(as.character(pairsdf$n_gene_to_high_and_comorbid))
pairsdf$n_ggenenormal_and__notcomorbid <- as.numeric(as.character(pairsdf$n_ggenenormal_and__notcomorbid))
pairsdf$n_genenormal_and_comorbid <- as.numeric(as.character(pairsdf$n_genenormal_and_comorbid))
pairsdf$n_gene_to_high_and_notcomorbid <- as.numeric(as.character(pairsdf$n_gene_to_high_and_notcomorbid))
pairsdf$n_world <- as.numeric(as.character(pairsdf$n_world))
pairsdf$n_gene_to_high <- as.numeric(as.character(pairsdf$n_gene_to_high))
pairsdf$n_gene_normal <- as.numeric(as.character(pairsdf$n_gene_normal))
pairsdf$n_comorbid <- as.numeric(as.character(pairsdf$n_comorbid))
pairsdf$n_not_comorbid <- as.numeric(as.character(pairsdf$n_not_comorbid))

pairsdf$odds_ratio <- (pairsdf$n_gene_to_high_and_comorbid*pairsdf$n_ggenenormal_and__notcomorbid)/(pairsdf$n_genenormal_and_comorbid*n_gene_to_high_and_notcomorbid)
pairsdf$per_to_high <- pairsdf$n_gene_to_high/pairsdf$n_world
pairsdf$per_comorbid <- pairsdf$n_comorbid/pairsdf$n_world
pairsdf$per_to_high_ofcormobid <- pairsdf$n_gene_to_high_and_comorbid/pairsdf$n_comorbid
pairsdf$per_to_high_ofnotcormobid <- pairsdf$n_gene_to_high_and_notcomorbid/pairsdf$n_not_comorbid
pairsdf$per_comorbid_ingenetohigh <-pairsdf$n_gene_to_high_and_comorbid/pairsdf$n_gene_to_high
pairsdf$per_comorbid_ingenenormal <- pairsdf$n_genenormal_and_comorbid/pairsdf$n_gene_normal



ggplot(pairsdf, aes(y=hypergpval_lowertail)) + 
  geom_histogram()



genedf_mininfo <- genedf %>% select(ENSEMBL,SYMBOL, SYMBOL)
pairsdf <- pairsdf %>% arrange(hypergpval_lowertail)
pairsdf2 <- merge(pairsdf, genedf_mininfo, on="ENSEMBL")
pairsdf2 <- pairsdf2 %>% arrange(hypergpval_lowertail)
pairsdf2$adj_pval = nrow(pairsdf2)*pairsdf2$hypergpval_lowertail


pairsdf2 %>% select(SYMBOL, comorbid, odds_ratio,per_comorbid, hypergpval_lowertail, per_comorbid_ingenetohigh, per_comorbid_ingenenormal,adj_pval) %>% filter(comorbid=="comorbid_astigmatism") %>% filter(SYMBOL=="LINC00114")
pairsdf2 %>%filter(comorbid=="comorbid_astigmatism") %>% filter(SYMBOL=="LINC00114")


outfilename = paste(outdir, "all_pairs_gene_to_high_comorbid_","_std", howmanystd , "min_n_gene_to_high_", min_n_gene_to_high, "min_people_with_comorbid",  min_people_with_comorbid,"min_count_sum", min_count_sum, comorbidfileroot, sep="")
write.csv(pairsdf2, outfilename)


gene_biotype_df <- genedf %>% select(SYMBOL, ENSEMBL, gene_biotype)

pseaindir="/Shares/down/public/HTP/RNAseq/scripts/pangea_data/"
pvaldf = read.csv(paste0(pseaindir,"filter_minpeopleconditions10_normcounts_0.8_RNA_chr21_completeT21_pangea_original_patient_mondo.csv.csv",sep="") )
#pvaldf = read.csv(paste0(pseaindir,"filter_minpeopleconditions10_normcounts_0.8_RNA_chr21_completeT21_pangea_original_patient.csv.csv",sep="") )
pvaldf$ENSEMBL <- pvaldf$gene
D21expectation$log_D21mean <- log(D21expectation$D21mean)
D21expectation_minexp <- D21expectation %>% filter(log_D21mean>0)

pvaldf_exp <- pvaldf %>% filter(ENSEMBL %in% D21expectation_minexp$ENSEMBL) %>% drop_na()
pvaldf_exp <- merge(pvaldf_exp, comorbidcounts_fixed,on="comorbid") 
pvaldf_exp <- pvaldf_exp %>% filter(people_with_comorbid>17) %>% drop_na()

pvaldf_exp$adj_pval = nrow(pvaldf_exp)*pvaldf_exp$pval

pvaldf_exp <- merge(gene_biotype_df, pvaldf_exp, on=ENSEMBL)

pvaldf_exp$sig <- pvaldf_exp$adj_pval<1

pvaldfsig <- pvaldf_exp %>% filter(sig==TRUE) 

pvaldfsig <- pvaldfsig %>% arrange(pval)


pvaldfsig_genes <- pvaldfsig %>% select(SYMBOL, gene_biotype) %>% unique()

as.data.frame(table(pvaldfsig_genes$gene_biotype))


mondo_lncRNAs_staples <- c("LINC00114", "LINC00310", "LINC01684")
mondo_staple_disorders <- c("astigmatism", "hypothyroidism", "obstructive_sleep_apnea_syndrome", "pulmonary_hypertension", "patent_foramen_ovale", "gastroesophageal_reflux_disease", "ventricular_septal_defect", "patent_ductus_arteriosus")

agene="ENSG00000159200"
acomorbid="obstructive_sleep_apnea_syndrome"

TFs=c("AIRE", "BACH1", "ERG", "ETS2", "GABPA", "OLIG1", "OLIG2", "PKNOX1", "PRDM15", "RUNX1", "SIM2", "SON", "ZBTB21")
lincRNAs_comorbid <- pvaldfsig %>% filter(gene_biotype=="lincRNA") %>% select(ENSEMBL, comorbid)
TFs_comorbid <- pvaldfsig %>% filter(SYMBOL %in% TFs)%>% select(ENSEMBL, comorbid)

plotapair <- function(thisENSEMBL, thiscomorbid){
  genename = paste0(unlist(genesymbolchr[genesymbolchr$ENSEMBL==thisENSEMBL,]["SYMBOL"]))
  peopleT21__thiscomorbid <- comorbidf_T21[comorbidf_T21[,paste0("comorbid_", acomorbid)] == TRUE, ]
  this_gene_exp_long <- gene_exp_long_exp %>% filter(ENSEMBL==thisENSEMBL)
  this_gene_exp_long_D21_T21 <- this_gene_exp_long %>% filter(sample_type %in% c("complete_T21", "D21"))
  this_gene_exp_long_D21_T21$sample_type <- factor(this_gene_exp_long_D21_T21$sample_type, levels=c("D21", "complete_T21"))
  this_gene_exp_long_T21_comorbid <- this_gene_exp_long_D21_T21 %>% filter(sample_type=="complete_T21")
  this_gene_exp_long_T21_comorbid$comorbid <- this_gene_exp_long_T21_comorbid$Patient %in% peopleT21__thiscomorbid$Patient
  this_gene_NES_pval <- pvaldf_exp %>% filter(ENSEMBL==agene) %>% arrange(desc(adj_pval)) %>% filter(comorbid %in% mondo_staple_disorders)
  plot2 <- ggplot(this_gene_exp_long_T21_comorbid, aes(y=normcounts, x=comorbid, fill=comorbid))+geom_violin()+geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    ggtitle(genename)+ scale_fill_manual(values=c("FALSE"="pink","TRUE"="red"))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab(thiscomorbid)
  filename <- paste0(outdir,"graphs/", genename, "and", thiscomorbid, ".pdf")
  print(filename)
  ggsave(filename, device = "pdf")
  plot1 <- ggplot(this_gene_exp_long_D21_T21, aes(y=normcounts, x=sample_type, fill=sample_type))+geom_violin()+geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    ggtitle(genename)+ scale_fill_manual(values=c("D21"="black","complete_T21"="darkred"))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Chr 21")
  filename <- paste0(outdir,"graphs/", genename, "_D21_T21",".pdf")
  ggsave(filename, device = "pdf")
  plot3 <- ggplot(this_gene_NES_pval, aes(x=norm_NES, y=comorbid, fill=sig))+
    geom_bar(stat="identity")+theme_minimal()+ scale_fill_manual(values=c("FALSE"="black","TRUE"="red"))
  filename <- paste0(outdir,"graphs/", genename, "norm_NES",".pdf")
  ggsave(filename, device = "pdf")
  howfar="done"
}


plot1 <-  plotapair(agene, acomorbid)

for(i in 1:nrow(TFs_comorbid)) {
  row <- TFs_comorbid[i,]
  agene=unlist(row[,"ENSEMBL"])
  acomorbid=unlist(row[,"comorbid"])
  print(agene)
  print(acomorbid)
  plotapair(agene, acomorbid)
}

for(i in 1:nrow(lincRNAs_comorbid)) {
  row <- lincRNAs_comorbid[i,]
  agene=unlist(row[,"ENSEMBL"])
  acomorbid=unlist(row[,"comorbid"])
  print(agene)
  print(acomorbid)
  plotapair(agene, acomorbid)
}

pairsdf2 %>%filter(comorbid=="comorbid_patent_foramen_ovale") %>% filter(SYMBOL=="FP236383.1") #OR 3.7
pvaldfsig %>%filter(comorbid=="patent_foramen_ovale") %>% filter(SYMBOL=="FP236383.1") # p-value 2.623501e-06, adj_pval = 0.04143032


pairsdf2 %>%filter(comorbid=="comorbid_pulmonary_hypertension") %>% filter(SYMBOL=="GABPA") #OR 22.55
pvaldfsig %>%filter(comorbid=="pulmonary_hypertension") %>% filter(SYMBOL=="GABPA") # p-value 1.696323e-06, adj_pval = 0.02678834


pairsdf2 %>%filter(comorbid=="comorbid_obstructive_sleep_apnea_syndrome") %>% filter(SYMBOL=="RCAN1") #OR 29.02041
pvaldfsig %>%filter(comorbid=="obstructive_sleep_apnea_syndrome") %>% filter(SYMBOL=="RCAN1") # p-value 6.280812e-05, adj_pval = 0.99


basicspairs <- pairsdf2 %>% select(ENSEMBL, comorbid, odds_ratio)


basicspairs$pairname <-paste0(basicspairs$ENSEMBL,basicspairs$comorbid)
basicspairs <- basicspairs %>% select(pairname, odds_ratio)

pvaldfsig$pairname <- paste0(pvaldfsig$ENSEMBL,"comorbid_",pvaldfsig$comorbid)

pvaldfsig3 <- merge(basicspairs,pvaldfsig, on="pairname")

pvaldfsig3 <- pvaldfsig3 %>% filter(norm_NES<0) 
