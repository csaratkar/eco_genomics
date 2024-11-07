# Are the genes that are differentially expressed in response to 28 °C and 33 °C relative to
# baseline conditions at 18 °C (and 22 °C) changing in expression in the same direction
# relative to baseline? 
# I.e., Do the shared differentially expressed genes have a positive or negative relationship? 
# To address this question, you would need to perform two contrasts in DESeq2 for each of the 
# developmental temperatures: 18_BASE versus 18_A28 and 18_BASE versus 18_A33. You could make 
# a Euler plot of these genes. 
# You can combine the results matrices and color code the genes significant in both contrasts,
# then make a scatterplot of their log2FoldChange expression responses  as we did in class.
# What is your answer? This can be repeated for contrasts for the 22 °C developmental 
# temperature

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(eulerr)


options(bitmapType = "cairo")

setwd("/gpfs1/cl/pbio3990/Transcriptomics/")

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                                                                      header =T, row.names = 1)
tail(countsTable)

countsTableRound <- round(countsTable) # because DESeq doesn't like decimals

tail(countsTableRound)

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = T, stringsAsFactors = T, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds,
                              design = ~ DevTemp + FinalTemp)

rownames(dds)

dds <- dds[rowSums(counts(dds) >= 10) >= 15,]
nrow(dds)

dds <-  DESeq(dds)

resultsNames(dds)


dds$group <-  factor(paste0(dds$DevTemp, dds$FinalTemp))
design(dds) <- ~ group
dds <- DESeq(dds)

dim(dds)#35527    21

resultsNames(dds)

#1. BASE and A28 to D18---------------------------

res_D18_BASE_D18_A28 <- results(dds, contrast = c("group", "D18BASE", "D18A28"), alpha = 0.05)
res_D18_BASE_D18_A28 <- res_D18_BASE_D18_A28[!is.na(res_D18_BASE_D18_A28$padj),]
res_D18_BASE_D18_A28 <- res_D18_BASE_D18_A28[order(res_D18_BASE_D18_A28$padj),]
head(res_D18_BASE_D18_A28)
summary(res_D18_BASE_D18_A28)

#make a list of which genes in our comparisons of interest are differentially expressed (list of DEGs)
degs_D18_BASE_D18_A28 <-  row.names(res_D18_BASE_D18_A28[res_D18_BASE_D18_A28$padj < 0.05,])

plotMA(res_D18_BASE_D18_A28, ylim = c(-4,4))

#2. BASE and A33 to D18---------------------------

res_D18_BASE_D18_A33 <- results(dds, contrast = c("group", "D18BASE", "D18A33"), alpha = 0.05)
res_D18_BASE_D18_A33 <- res_D18_BASE_D18_A33[!is.na(res_D18_BASE_D18_A33$padj),]
res_D18_BASE_D18_A33 <- res_D18_BASE_D18_A33[order(res_D18_BASE_D18_A33$padj),]
head(res_D18_BASE_D18_A33)
summary(res_D18_BASE_D18_A33)

#make a list of which genes in our comparisons of interest are differentially expressed (list of DEGs)
degs_D18_BASE_D18_A33 <-  row.names(res_D18_BASE_D18_A33[res_D18_BASE_D18_A33$padj < 0.05,])

plotMA(res_D18_BASE_D18_A33, ylim = c(-4,4))

#Euler plot-----------------------------------

#Numbers of DEGs per Final Temp
length(degs_D18_BASE_D18_A28)#41
length(degs_D18_BASE_D18_A33) #332

#Look at overlap in which genes are differentially expressed in multiple contrasts
length(intersect(degs_D18_BASE_D18_A28,degs_D18_BASE_D18_A33)) #34

41-34 #7 diff for A28
332-34 #298 diff for A33


myEuler18 <- euler(c("A28" = 7, 
                   "A33" = 298, 
                   "A28&A33" = 34))

plot(myEuler18, lty=1:2, quantities = T, fill = c("lightskyblue", "indianred1", "orchid2"))

title("Number of differential genes between FinalTemps at DevTemp 18 C")


#Scatterplot--------------------

#contrast D18_A28vsBASE
res_D18_BASEvsA28 <- as.data.frame(results(dds, contrast = c("group", "D18BASE", "D18A28"), alpha = 0.05))


#contrast D18_A33vsBASE
res_D18_BASEvsA33 <- as.data.frame(results(dds, contrast = c("group", "D18BASE", "D18A33"), alpha = 0.05))

#merge dataframes
res_df18 <- merge(res_D18_BASEvsA28, res_D18_BASEvsA33, by = "row.names", suffixes = c(".28", ".33"))
rownames(res_df18) <- res_df18$Row.names
res_df18 <- res_df18[,-1]

#define color mapping logic with mutate
res_df18 <- res_df18 %>% 
  mutate(fill=case_when(
    padj.28<0.05 & stat.28<0 ~ "skyblue",
    padj.28<0.05 & stat.28>0 ~ "dodgerblue",
    padj.33<0.05 & stat.33<0 ~ "salmon",
    padj.33<0.05 & stat.33>0 ~ "red"
  ))

#Count the number of points per fill color
color_counts18 <- res_df18 %>% 
  group_by(fill) %>% 
  summarise(count = n())

label_positions18 <- data.frame(fill = c("salmon", "dodgerblue", "red", "skyblue"),
                              x_pos=c(1,5,0,-7.5),
                              y_pos = c(-5,0,9,3))

label_data18 <- merge(color_counts18, label_positions18, by = "fill")
#Plot

scatter18 <- ggplot(res_df18, aes(x = log2FoldChange.28, y = log2FoldChange.33, color = fill)) +
  geom_point(alpha = 0.8) +
  scale_color_identity()+
  geom_text(data = label_data18, aes(x = x_pos, y = y_pos, label = count, color = fill), size = 5)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")+
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey")+
  xlim(-10,10) + ylim(-10,10)+
  labs(x = "Log2FoldChange 28 vs. BASE at 18",
       y = "Log2FoldChange 33 vs. BASE at 18",
       title = "How does response to DevTemp 18 vary by FinalTemp?")+
  theme_minimal()

ggsave("~/Projects/eco_genomics/transcriptomics/figures/scatter_dev18.png", scatter18, width = 6, height = 6)
#4. BASE and A28 to D22---------------------------

res_D22_BASE_D22_A28 <- results(dds, contrast = c("group", "D22BASE", "D22A28"), alpha = 0.05)
res_D22_BASE_D22_A28 <- res_D22_BASE_D22_A28[!is.na(res_D22_BASE_D22_A28$padj),]
res_D22_BASE_D22_A28 <- res_D22_BASE_D22_A28[order(res_D22_BASE_D22_A28$padj),]
head(res_D22_BASE_D22_A28)
summary(res_D22_BASE_D22_A28)

#make a list of which genes in our comparisons of interest are differentially expressed (list of DEGs)
degs_D22_BASE_D22_A28 <-  row.names(res_D22_BASE_D22_A28[res_D22_BASE_D22_A28$padj < 0.05,])

plotMA(res_D22_BASE_D22_A28, ylim = c(-4,4))

#5. BASE and A33 to D22---------------------------

res_D22_BASE_D22_A33 <- results(dds, contrast = c("group", "D22BASE", "D22A33"), alpha = 0.05)
res_D22_BASE_D22_A33 <- res_D22_BASE_D22_A33[!is.na(res_D22_BASE_D22_A33$padj),]
res_D22_BASE_D22_A33 <- res_D22_BASE_D22_A33[order(res_D22_BASE_D22_A33$padj),]
head(res_D22_BASE_D22_A33)
summary(res_D22_BASE_D22_A33)

#make a list of which genes in our comparisons of interest are differentially expressed (list of DEGs)
degs_D22_BASE_D22_A33 <-  row.names(res_D22_BASE_D22_A33[res_D22_BASE_D22_A33$padj < 0.05,])

plotMA(res_D22_BASE_D22_A33, ylim = c(-4,4))

#Euler plot------------------------------

#Numbers of DEGs per Final Temp
length(degs_D22_BASE_D22_A28)#289
length(degs_D22_BASE_D22_A33) #1564

#Look at overlap in which genes are differentially expressed in multiple contrasts
length(intersect(degs_D22_BASE_D22_A28,degs_D22_BASE_D22_A33)) #144

289 - 144 #145 for A28
1564-144 #1420 for A33

myEuler33 <- euler(c("A28" = 145, 
                   "A33" = 1420, 
                   "A28&A33" = 144))

plot(myEuler33, lty=1:2, quantities = T, fill = c("blue", "red", "purple"))
plot(myEuler33, lty=1:2, quantities = T, fill = c("dodgerblue3", "firebrick3", "darkorchid"))

#Scatterplot--------------------

#contrast D22_A28vsBASE
res_D22_BASEvsA28 <- as.data.frame(results(dds, contrast = c("group", "D22BASE", "D22A28"), alpha = 0.05))


#contrast D22_A33vsBASE
res_D22_BASEvsA33 <- as.data.frame(results(dds, contrast = c("group", "D22BASE", "D22A33"), alpha = 0.05))

#merge dataframes
res_df22 <- merge(res_D22_BASEvsA28, res_D22_BASEvsA33, by = "row.names", suffixes = c(".28", ".33"))
rownames(res_df22) <- res_df22$Row.names
res_df22 <- res_df22[,-1]

#define color mapping logic with mutate
res_df22 <- res_df22 %>% 
  mutate(fill=case_when(
    padj.28<0.05 & stat.28<0 ~ "slateblue1",
    padj.28<0.05 & stat.28>0 ~ "slateblue4",
    padj.33<0.05 & stat.33<0 ~ "tomato2",
    padj.33<0.05 & stat.33>0 ~ "red3"
  ))

#Count the number of points per fill color
color_counts22 <- res_df22 %>% 
  group_by(fill) %>% 
  summarise(count = n())

label_positions22 <- data.frame(fill = c("tomato2", "slateblue4", "red3", "slateblue1"),
                                x_pos=c(1,5,0,-7.5),
                                y_pos = c(-5,0,9,3))

label_data22 <- merge(color_counts22, label_positions22, by = "fill")
#Plot

scatter22 <- ggplot(res_df22, aes(x = log2FoldChange.28, y = log2FoldChange.33, color = fill)) +
  geom_point(alpha = 0.8) +
  scale_color_identity()+
  geom_text(data = label_data22, aes(x = x_pos, y = y_pos, label = count, color = fill), size = 5)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")+
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey")+
  xlim(-10,10) + ylim(-10,10)+
  labs(x = "Log2FoldChange 28 vs. BASE at 22",
       y = "Log2FoldChange 33 vs. BASE at 22",
       title = "How does response to DevTemp 22 vary by FinalTemp?")+
  theme_minimal()

ggsave("~/Projects/eco_genomics/transcriptomics/figures/scatter_dev22.png", scatter22, width = 6, height = 6)

