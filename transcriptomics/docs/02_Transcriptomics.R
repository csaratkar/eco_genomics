#Load in 01_DESeq data and run DESeq object
library(pheatmap)

options(bitmapType = "cairo")

resultsNames(dds) #"Intercept"   "DevTemp_D22_vs_D18"    "FinalTemp_A33_vs_A28"  "FinalTemp_BASE_vs_A28"

#Creates results files for Developmental temperature 22 vs 18
res_D22vsD18 <- results(dds, name = "DevTemp_D22_vs_D18", alpha = 0.5)

#Order by significance
res_D22vsD18 <- res_D22vsD18[order(res_D22vsD18$padj),]
head(res_D22vsD18)
# log2FoldChange is the quantity of gene expression of D22 relative to D18

#look at counts of a specific top gene that we're interested in to validate that the model is working
d <- plotCounts(dds, gene = "TRINITY_DN140854_c0_g5_i2", int = (c("DevTemp", "FinalTemp")), returnData = T)
d

p <- ggplot(d, aes(x = DevTemp, y = count, color = DevTemp, shape = FinalTemp)) +
  theme_minimal() +
  theme(text = element_text(size=20), panel.grid.major=element_line(color = "grey")) +
  geom_point(position=position_jitter(w=0.2, h=0), size = 3)
p

#Log fold change average (MA) plot
plotMA(res_D22vsD18, ylim = c(-4, 4))

#Volcano plot
#convert our DESeq results object into na data frame to plot
res_df <- as.data.frame(res_D22vsD18)

#add a column to dataframe to denote if a gene is significantly differentially expressed
res_df$Significant <- ifelse(res_df$padj <0.05 & abs(res_df$log2FoldChange)>1, "Significant", "Not Significant")

#Plot
ggplot(res_df, aes(x=log2FoldChange, y = -log10(padj), color = Significant)) +
  geom_point(alpha=0.8)+
  scale_color_manual(values = c("slateblue", "tomato")) +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", title = "Volcano plot") +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_hline(yintercept=-log10(0.05), linetype = "dashed", color = "orange")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "orange")


#Heat map
vsd <- vst(dds, blind = F)
topgenes <- head(rownames(res_D22vsD18), 20)
mat <- assay(vsd)[topgenes,]
df <- as.data.frame(colData(dds)[,c("DevTemp", "FinalTemp")])
pheatmap(mat, annotation_col=df, show_rownames=F, cluster_cols=T, cluster_rows=T)
# columns are samples, rows are genes