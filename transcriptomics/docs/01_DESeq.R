### Code for analyzing RNAseq data using DESeq2

library(DESeq2)
library(tidyverse)

options(bitmapType = "cairo")

setwd("/gpfs1/cl/pbio3990/Transcriptomics/")

# Import counts matrix
countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header =T, row.names = 1)

tail(countsTable)


countsTableRound <- round(countsTable) # because DESeq doesn't like decimals

tail(countsTableRound)

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = T, stringsAsFactors = T, row.names = 1)

############################################################################################################
#Explore counts matrix

#See how may read we have from each sample
colSums(countsTableRound)
mean(colSums(countsTableRound))

barplot(colSums(countsTableRound), names.arg = colnames(countsTableRound), 
        cex.names = 0.5, las = 2, ylim = c(0,30000000))

abline(h=mean(colSums(countsTableRound)), col = "blue4", lwd = 2)

#Average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) #3244.739
median(rowSums(countsTableRound)) #64

apply(countsTableRound, 2, mean) #gives a sense of variation in sequencing

############################################################################################################
#Start analysis in DESeq2

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds,
                              design = ~ DevTemp + FinalTemp)

dim(dds)

dds <- dds[rowSums(counts(dds) >= 10) >= 15,]
nrow(dds)#35527 = number of transcripts with more than 10 read and more than or equal to 15 reads

# Run the DESeq model to test for global differential gene expression
dds <-  DESeq(dds)

#List the results you've generated 
resultsNames(dds) #"Intercept"    "DevTemp_D22_vs_D18"    "FinalTemp_A33_vs_A28"  "FinalTemp_BASE_vs_A28"

#Visualize our global gene expression patterns using PCA
#First, transform the data for plotting using variance stabilization

vsd <- vst(dds, blind = F)

pcaData <-  plotPCA(vsd, intgroup = c("DevTemp", "FinalTemp"), returnData = T)

percentVar <- round(100*attr(pcaData,"percentVar"))

final_temp_colors <- c("BASE" = "grey",
                       "A28" = "hotpink",
                       "A33" = "red")

shapes_choose <- c("D18" = 16,
                   "D22" = 18)

p <- ggplot(pcaData, aes(PC1, PC2, color = FinalTemp, shape = DevTemp)) +
  geom_point(size = 5) +
  scale_shape_manual(values = shapes_choose) +
  scale_color_manual(values = final_temp_colors)+
  labs(x = paste0("PC1: ", percentVar[1], "%"), y = paste0("PC2: ", percentVar[2], "%")) +
  theme_bw(base_size = 16)

