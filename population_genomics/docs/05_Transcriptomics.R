### Code for analyzing RNAseq data using DESeq2

library(DESeq2)
library(tidyverse)

setwd("/gpfs1/cl/pbio3990/Transcriptomics/")

# Import counts matrix
countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header =T, row.names = 1)
countsTableRound <- round(countsTable) # because DESeq doesn't like decimals

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = T, stringsAsFactors = T, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds,
                              design = ~ DevTemp + FinalTemp)
