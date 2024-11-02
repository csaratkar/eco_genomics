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

# res_D22vsD18 <- results(dds, name = "DevTemp_D22_vs_D18", alpha = 0.5)
# 
# Order by significance
# res_D22vsD18 <- res_D22vsD18[order(res_D22vsD18$padj),]
# head(res_D22vsD18)

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

#3. A28 and A33 to D18---------------------------

res_D18_A28_D18_A33 <- results(dds, contrast = c("group", "D18A28", "D18A33"), alpha = 0.05)
res_D18_A28_D18_A33 <- res_D18_A28_D18_A33[!is.na(res_D18_A28_D18_A33$padj),]
res_D18_A28_D18_A33 <- res_D18_A28_D18_A33[order(res_D18_A28_D18_A33$padj),]
head(res_D18_A28_D18_A33)
summary(res_D18_A28_D18_A33)

#make a list of which genes in our comparisons of interest are differentially expressed (list of DEGs)
degs_D18_A28_D18_A33 <-  row.names(res_D18_A28_D18_A33[res_D18_A28_D18_A33$padj < 0.05,])

plotMA(res_D18_A28_D18_A33, ylim = c(-4,4))

#Numbers of DEGs per Final Temp
length(degs_D18_BASE_D18_A28)#41
length(degs_D18_BASE_D18_A33) #332
#length(degs_D18_A28_D18_A33) #234

#Look at overlap in which genes are differentially expressed in multiple contrasts
length(intersect(degs_D18_BASE_D18_A28,degs_D18_BASE_D18_A33)) #34
#length(intersect(degs_D18_BASE_D18_A28, degs_D18_A28_D18_A33)) #18
#length(intersect(degs_D18_BASE_D18_A33, degs_D18_A28_D18_A33)) #163
#length(intersect(degs_D18_BASE_D18_A28, intersect(degs_D18_BASE_D18_A33, degs_D18_A28_D18_A33))) #17

41-34-18+17 #6 diff expressed genes for A28 at 18
332-34-163+17 #152 diff expressed genes between BASE and A33 at 18
234-18-163+17 #70 diff expressed genes between A33 and A28 at 18..
34-17 #17 diff expressed genes between BASE and A28 at 18

41-34 #7 diff for A28
332-34 #298 diff for A33


myEuler18 <- euler(c("A28" = 7, 
                   "A33" = 298, 
                   "A28&A33" = 34))

plot(myEuler18, lty=1:2, quantities = T, fill = c("skyblue", "tomato", "lavender"))

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
    padj.28<0.05 & stat.28<0 ~ "turquoise2",
    padj.28<0.05 & stat.28>0 ~ "magenta1",
    padj.33<0.05 & stat.33<0 ~ "blue2",
    padj.33<0.05 & stat.33>0 ~ "red"
  ))

#Count the number of points per fill color
color_counts18 <- res_df18 %>% 
  group_by(fill) %>% 
  summarise(count = n())

label_positions18 <- data.frame(fill = c("blue2", "magenta1", "red", "turquoise2"),
                              x_pos=c(1,5,0,-7.5),
                              y_pos = c(-5,0,9,3))

label_data <- merge(color_counts18, label_positions18, by = "fill")
#Plot

plot18 <- ggplot(res_df18, aes(x = log2FoldChange.28, y = log2FoldChange.33, color = fill)) +
  geom_point(alpha = 0.8) +
  scale_color_identity()+
  geom_text(data = label_data, aes(x = x_pos, y = y_pos, label = count, color = fill), size = 5)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")+
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey")+
  xlim(-10,10) + ylim(-10,10)+
  labs(x = "Log2FoldChange 28 vs. BASE at 18",
       y = "Log2FoldChange 33 vs. BASE at 18",
       title = "How does response to DevTemp 18 vary by FinalTemp?")+
  theme_minimal()


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

#6. A28 and A33 to D22---------------------------

res_D22_A28_D22_A33 <- results(dds, contrast = c("group", "D22A28", "D22A33"), alpha = 0.05)
res_D22_A28_D22_A33 <- res_D22_A28_D22_A33[!is.na(res_D22_A28_D22_A33$padj),]
res_D22_A28_D22_A33 <- res_D22_A28_D22_A33[order(res_D22_A28_D22_A33$padj),]
head(res_D22_A28_D22_A33)
summary(res_D22_A28_D22_A33)

#make a list of which genes in our comparisons of interest are differentially expressed (list of DEGs)
degs_D22_A28_D22_A33 <-  row.names(res_D22_A28_D22_A33[res_D22_A28_D22_A33$padj < 0.05,])

plotMA(res_D22_A28_D22_A33, ylim = c(-4,4))

#Numbers of DEGs per Final Temp
length(degs_D22_BASE_D22_A28)#289
length(degs_D22_BASE_D22_A33) #1564
#length(degs_D22_A28_D22_A33) #200

#Look at overlap in which genes are differentially expressed in multiple contrasts
length(intersect(degs_D22_BASE_D22_A28,degs_D22_BASE_D22_A33)) #144
#length(intersect(degs_D22_BASE_D22_A28, degs_D22_A28_D22_A33)) #8
#length(intersect(degs_D22_BASE_D22_A33, degs_D22_A28_D22_A33)) #149
#length(intersect(degs_D22_BASE_D22_A28, intersect(degs_D22_BASE_D22_A33, degs_D22_A28_D22_A33))) #3

289 - 144 #145 for A28
1564-144 #1420 for A33

myEuler33 <- euler(c("A28" = 145, 
                   "A33" = 1420, 
                   "A28&A33" = 144))

plot(myEuler33, lty=1:2, quantities = T, fill = c("blue", "red", "purple"))

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
    padj.28<0.05 & stat.28<0 ~ "turquoise2",
    padj.28<0.05 & stat.28>0 ~ "magenta1",
    padj.33<0.05 & stat.33<0 ~ "blue2",
    padj.33<0.05 & stat.33>0 ~ "red"
  ))

#Count the number of points per fill color
color_counts22 <- res_df22 %>% 
  group_by(fill) %>% 
  summarise(count = n())

label_positions22 <- data.frame(fill = c("blue2", "magenta1", "red", "turquoise2"),
                                x_pos=c(1,5,0,-7.5),
                                y_pos = c(-5,0,9,3))

label_data <- merge(color_counts22, label_positions22, by = "fill")
#Plot

plot22 <- ggplot(res_df22, aes(x = log2FoldChange.28, y = log2FoldChange.33, color = fill)) +
  geom_point(alpha = 0.8) +
  scale_color_identity()+
  geom_text(data = label_data, aes(x = x_pos, y = y_pos, label = count, color = fill), size = 5)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")+
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey")+
  xlim(-10,10) + ylim(-10,10)+
  labs(x = "Log2FoldChange 28 vs. BASE at 22",
       y = "Log2FoldChange 33 vs. BASE at 22",
       title = "How does response to DevTemp 22 vary by FinalTemp?")+
  theme_minimal()

