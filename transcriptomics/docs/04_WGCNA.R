#Script for analyzing and visualizing gene correlation networks
library(DESeq2)
library(tidyverse)
library(WGCNA);options(stringsAsFactors = F);
library(CorLevelPlot)
library(gridExtra)
library(Rmisc)

options(bitmapType = "cairo")

setwd("/gpfs1/cl/pbio3990/Transcriptomics/")

#Step 1. Import counts data--------------------------
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

#Step 2: Detecting outlier---------------------
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

#Step 3: Noramlization-------------------------------
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

#Step 4 :Network construction----------------------

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

soft_power <- 26
temp_cor <- cor
cor <- WGCNA::cor #set the temp_cor func


#to use WGCNA's correlation function

norm.counts[] <- sapply(norm.counts, as.numeric)

#blockwiseModules() creates the network and identifies modules based on parameters that we choose
bwnet26 <- blockwiseModules(norm.counts,
                            maxBlockSize = 30000,
                            TOMType = "signed",
                            power = soft_power,
                            mergeCutHeight = 0.25,
                            randomSeed = 1234,
                            verbose = 3) #TOMType controls  interpreting - and + correlation
#signed is + and unsigned is - and +

cor <- temp_cor #resets cor() to baseR's cor func instead of using WGCNA's cor func

saveRDS(bwnet26, file = "~/Projects/eco_genomics/transcriptomics/outputs/bwnet26.rds")

#load the bwnet file in later use:
#bwnet26 <- readRDS("outputs/bwnet26.rds")

#Step 5: Explore module Eigengenes-------------------------------

module_eigengenes <- bwnet26$MEs
head(module_eigengenes)
dim(module_eigengenes)

#get the number of genes for wach module
table(bwnet26$colors)

#Plot dendrogram and module colors(merge and unmerged)[based on similarity cut off]
plotDendroAndColors(bwnet26$dendrograms[[1]], 
                    cbind(bwnet26$unmergedColors, bwnet26$colors),
                    c("unmerged", "merged"),
                    dendroLabels = F,
                    addGuide = T,
                    hang = 0.03,
                    guideHang = 0.05)

#Step 6:Correlation of modules with traits-----------------------

#Define the numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

#Test for a correlation between module eigengenes and trait data
modules.trait.corr <- cor(module_eigengenes, traitData, use = "p")

#Calculate p-vals for eeach correlation
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

#Visualize module-trait association as a heatmap
heatmap.data <- merge(module_eigengenes, traitData, by = "row.names")
head(heatmap.data)

#Address error of row.names not being numeric
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = "Row.names")

names(heatmap.data)

#Make pretty heat map of correlations
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[42:44], #the values may chnage based on number of eigengenes
             y = names(heatmap.data)[1:41],
             col = c("blue2", "skyblue", "white", "pink", "red"))

