#Script for analyzing and visualizing gene correlation networks
library(DESeq2)
library(tidyverse)
library(WGCNA);options(stringsAsFactors = F);
#library(CorLevelPlot)
library(gridExtra)
#library(Rmisc)

options(bitmapType = "cairo")

setwd("/gpfs1/cl/pbio3990/Transcriptomics/")

#Step 1. Import counts data
countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header =T, row.names = 1)

tail(countsTable)


countsTableRound <- round(countsTable) # because DESeq doesn't like decimals

tail(countsTableRound)

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = T, stringsAsFactors = T, row.names = 1)

traitData = read.table("/gpfs1/cl/pbio3990/Trait_Data.txt", header = T, 
                       row.names = 1)
#filter the matrix to just base data because those are the data for which we have traits measured
filtered_count_matrix_BASEonly <- countsTable[,conds$FinalTemp == "BASE"]
filtered_sample_metadata_BASEonly <- conds[conds$FinalTemp == "BASE",]
rounded_filtered_count_maxtrix <- round(filtered_count_matrix_BASEonly)

#Step 2: Detecting outlier
#detecting outlier genes
gsg <- goodSamplesGenes(t(rounded_filtered_count_maxtrix))
summary(gsg)

table(gsg$goodGenes)
# FALSE  TRUE
# 37235 82203

table(gsg$goodSamples) #All TRUE
#filter out bad gene

data_WGCNA <- rounded_filtered_count_maxtrix[gsg$goodGenes == TRUE, ]
dim(data_WGCNA)

#use clustering with a tree demdrogram to identify outlier samples
htree <- hclust(dist(t(data_WGCNA)), method = "average")
plot(htree)

#PCA - outlier detection method
pca <- prcomp(t(data_WGCNA))
pca_data <- pca$x
#make a df
pca_data <- as.data.frame(pca_data)

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits =2)

ggplot(pca_data, aes(PC1,PC2))+
  geom_point()+
  geom_text(label = rownames(pca_data))+
  labs(x = paste0("PC1: ", pca.var.percent[1],"%"),
       y = paste0("PC2: ", pca.var.percent[2], "%"))

#Step 3: Noramlization
colData <- row.names(filtered_sample_metadata_BASEonly)

dds_WGCNA <- DESeqDataSetFromMatrix(countData = data_WGCNA,
                                    colData = filtered_sample_metadata_BASEonly,
                                    design = ~1) #there are no specified groups

dds_WGCNA_75 <- dds_WGCNA[rowSums(counts(dds_WGCNA)>= 15) >=6,]
nrow(dds_WGCNA_75)
dds_norm <- vst(dds_WGCNA_75) # variance stabilization

#save normalized counts to use below
norm.counts <- assay(dds_norm) %>% 
  t()

#Step 4 :Network construction

#Choose a set of soft-threshold powers
power <- c(c(1:10), seq(fro = 12, to = 50, by = 2))

#Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)
sft.data <- sft$fitIndices

#Plot to pick power
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power))+
  geom_point()+
  geom_text(nudge_y = 0.1)+
  geom_hline(yintercept = 0.8, color = "red")+
  labs(x = "Power", y = "Scale free topology model fit, signed R^2")+
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power))+
  geom_point()+
  geom_text(nudge_y = 0.1)+
  geom_hline(yintercept = 0.8, color = "red")+
  labs(x = "Power", y = "Mean Connectivity")+
  theme_classic()

grid.arrange(a1,a2, nrow =2)
