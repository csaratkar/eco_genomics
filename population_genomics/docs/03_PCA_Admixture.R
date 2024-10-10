library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)


options(bitmapType = "cairo")

setwd("~/Projects/eco_genomics/population_genomics/")

vcf <- read.vcfR("outputs/vcf_final.filtered.vcf.gz")

#Need to thin in the SNPs for LD (Linkage Disequilibrium) before running PCA and Admixture analyses
# to account for the assumption of indepence for loci

vcf.thin <- distance_thin(vcf, min.distance = 500)
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
dim(meta)

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]) , ]
dim(meta2)

write.vcf(vcf.thin, "outputs/vcf_final.filtered.thinned.vcf.gz")
 # hide the uncompressed vcf because too big for github

system("gunzip -c ~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz > ~/vcf_final.filtered.thinned.vcf")

geno <- vcf2geno(input.file = "/gpfs1/home/c/s/csaratka/vcf_final.filtered.thinned.vcf",
                 output.file = "outputs/vcf_final.filtered.thinned.geno")

CentPCA <- LEA ::pca("~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.geno", scale = TRUE)

CentPCA <- load.pcaProject("vcf_final.filtered.thinned.pcaProject")

show(CentPCA)

plot(CentPCA) # shows eigenvalues(x), first dot is how useful PC1 is, second dot is how useful PC2 is etc.

# to find the percentages of how much a PC axis explains the variance it is eigenvalue of the axis over the number of eigenvalues
CentPCA$eigenvalues[1]/sum(CentPCA$eigenvalues) # the % variance explained by PC1

plot(CentPCA$projections,
     col=as.factor(meta2$region))

ggplot(as.data.frame(CentPCA$projections),
       aes(x = V1, y = V2, color = meta2$region, shape=meta2$continent)) +
       geom_point(alpha = .5) +
  labs(title = " Centaurea genetic PCA", x = "PC1(2.2%)", y = "PC2(1.1%)", color = "Region", shape = "Continents") 
# + lims(x = c(-10,10), y = c(-10,10))

ggsave("figures/CentPCA_PC1vPC2.pdf", width = 6, height = 6, units = "in")

ggplot(as.data.frame(CentPCA$projections),
       aes(x = V2, y = V3, color = meta2$region, shape=meta2$continent)) +
  geom_point(alpha = .5) +
  labs(title = " Centaurea genetic PCA", x = "PC2", y = "PC3", color = "Region", shape = "Continents") 
# + lims(x = c(-10,10), y = c(-10,10))

# Run Admixture analyses and create plots using the func snmf() in the LEA R package

CentAdmix <-snmf("outputs/vcf_final.filtered.thinned.geno",
                 K = 1:10,
                 entropy = TRUE,
                 repetitions = 3,
                 project = "new") #if adding to analysis, you can choose project = "continue"

par(mfrow = c(2,1))
plot(CentAdmix, main = "SNMF", col = "blue4")
#Pay attention to the "elbow" of the graph, the values you want to focus on 

plot(CentPCA$eigenvalues[1:10], ylab = "Eigen values", xlab = "Number of PCs", col = "blue4", main="PCA")
dev.off() #resets plotting window


myK = 3

CE = cross.entropy(CentAdmix, K = myK)
best = which.min(CE) #lowest croos entropy value

myKQ = Q(CentAdmix, K = myK, run = best) 
#dataset of each individuals' percentage ancestry from each K group

myKQmeta = cbind(myKQ, meta2)

my.colors = c("blue4", "gold", "tomato", "lightblue", "olivedrab")

myKQmeta = as_tibble(myKQmeta) %>% 
  group_by(continent) %>% 
  arrange(region, pop, .by_group = TRUE)


pdf("figures/Admixture_K4.pdf", width = 10, height = 5)
barplot(as.matrix(t(myKQmeta[ ,1:myK])),
        border = NA,
        space = 0,
        col = my.colors[1:myK],
        xlab = "Geographic regions",
        ylab = "Ancestry proportions",
        main = paste0("Ancestry matrix K =", myK))
axis(1,
     at = 1:length(myKQmeta$region),
     labels = myKQmeta$region,
     tick = F,
     cex.axis = 0.5,
     las = 3)
dev.off()
