library(vcfR)
library(tidyverse)
library(qqman)
library(SNPfiltR)
library(gt)
library(LEA)
library(pcadapt)

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")


vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")

dna <-  ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format = "fasta")

gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep = "\t", quote = "")

chr1 <- create.chromR(name = "Chromosome 1", vcf = vcf, seq = dna, ann = gff)


vcf.filt <- hard_filter(vcf, depth = 3) #Could do higher depth values for more accuracy
vcf.filt <- max_depth(vcf.filt, maxdepth = 60) # filters out genotypes with >60 reads/SNP

meta <- read.csv("metadata/meta4vcf.csv", header = TRUE)

meta2 <- meta[,c(1,4)]

names(meta2) <- c("id", "pop")

meta2$id = as.factor(meta2$id)
meta2$pop =  as.factor(meta2$pop)

vcf.filt.indMiss.75 <- missing_by_sample(vcf.filt, 
                                      popmap = meta2,
                                      cutoff = 0.75) 

vcf.filt.indMiss.25 <- missing_by_sample(vcf.filt, 
                                      popmap = meta2,
                                      cutoff = 0.25)

write.vcf(vcf.filt.indMiss.25,
          "~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.25.vcf.gz")


vcf.filt.indMiss.5 <- missing_by_sample(vcf.filt, 
                                      popmap = meta2,
                                      cutoff = 0.5)
write.vcf(vcf.filt.indMiss.5,
          "~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.5.vcf.gz")

X11.options(type="cairo")
options(bitmapType = "cairo")

#############################################################################################################################################
#2. a-d for vcf.filt.indMiss.25

vcf25 <- read.vcfR("~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.25.vcf.gz")

meta25<- meta[meta$id %in% colnames(vcf25@gt[,-1]),] 
dim(meta25)

vcf.div25 <- genetic_diff(vcf25,
                        pops = as.factor(meta25$region),
                        method = "nei")

str(vcf.div25)


chr.main25 <- unique(vcf.div25$CHROM)[1:8]

chrnum25 <- as.data.frame(cbind(chr.main25,seq(1,8,1)))

vcf.div.MHplot25 <- left_join(chrnum25, vcf.div25, join_by(chr.main25==CHROM))
dim(vcf.div.MHplot25)

vcf.div.MHplot25 <- vcf.div.MHplot25 %>% 
  filter(Gst>0) %>% 
  mutate(SNP = paste0(chr.main25, "_", POS))
str(vcf.div.MHplot25)

vcf.div.MHplot25$V2= as.numeric(vcf.div.MHplot25$V2)
vcf.div.MHplot25$POS= as.numeric(vcf.div.MHplot25$POS)

#manhattan(vcf.div.MHplot25,
          # chr = "V2",
          # bp = "POS",
          # p = "Gst",
          # col = c("blue4", "orange3"),
          # logp = FALSE,
          # ylab = "Fst among regions",
          # suggestiveline = quantile(vcf.div.MHplot25$Gst, 0.5))

#suggestline - suggests that the region above the line is experiencing a higher level of differentiation
#write.csv(vcf.div.MHplot, "~/Projects/eco_genomics/population_genomics/outputs/Genetic_Diff_byRegion",
          #quote = F,
          #row.names = F)

#a.
num_indSNP25 <- tibble(Individuals = nrow(vcf.div.MHplot25), Loci = n_distinct(vcf.div.MHplot25$SNP))

num_indSNP25|>
  gt() |>
  tab_header(title = "# of inds and loci (25% filter)")

names(vcf.div.MHplot25) 
vcf.div.MHplot25 %>% 
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  ggplot(aes(x = value, fill = name)) +
  geom_histogram(position = "identity",alpha = 0.5, bins = 50)+
  labs(x = "Gene diversity (Hs) withing Regions (25% filter)", y = "Counts of SNPs", title = "Genome-wide expected heterozygosity (Hs)", fill = "Regions")

#ggsave("Histogram_GenomDiversity_byRegion.25.pdf", 
#       path = "~/Projects/eco_genomics/population_genomics/figures")

#b.

filter_summary_25 <- vcf.div.MHplot25 %>% 
  as_tibble() %>%   
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>% 
  filter(value!=0 & value<0.5) %>% 
  summarise(Avg_Hs=mean(value), SD_Hs = sd(value), N_Hs=n())


Hs_SD25table <- filter_summary_25 |>
  gt() |>
  tab_header(
    title = "Hs and SD for each region (25% filter)"
  ) |>
  cols_label(name = "Region",
             Avg_Hs = "Average (Hs)",
             SD_Hs = "SD (HS)", 
             N_Hs = "N (Hs)")

#c.


Hs25  <- vcf.div.MHplot25[,4:9]

CEU25_0 <- sum(Hs25$Hs_CEU==0)
NE25_0 <- sum(Hs25$Hs_NE ==0)
NEU25_0 <- sum(Hs25$Hs_NEU ==0)
PNW25_0 <- sum(Hs25$Hs_PNW ==0)
SEU25_0 <- sum(Hs25$Hs_SEU ==0)
WEU25_0 <- sum(Hs25$Hs_WEU ==0)

CEU25 <- sum(Hs25$Hs_CEU>0)
NE25 <- sum(Hs25$Hs_NE >0)
NEU25 <- sum(Hs25$Hs_NEU >0)
PNW25 <- sum(Hs25$Hs_PNW >0)
SEU25 <- sum(Hs25$Hs_SEU >0)
WEU25 <- sum(Hs25$Hs_WEU >0)

Nloci_25_0 <-  tibble(Region = c("CEU", "NE", "NEU", "PNW", "SEU", "WEU"),
       Equal_zero = c(CEU25_0, NE25_0, NEU25_0, PNW25_0, SEU25_0, WEU25_0),
       Above_zero = c(CEU25, NE25, NEU25, PNW25, SEU25, WEU25))

Nloci_25table <- Nloci_25_0|>
  gt() |>
  tab_header(
    title = "Number of loci = 0 or > 0 (25% filter)"
  ) |>
  cols_label(Region = "Region",
             Equal_zero = "Equal 0",
             Above_zero = "Above 0")

#d.
#vcf <- read.vcfR("outputs/vcf_final.filtered.vcf.gz")

#Need to thin in the SNPs for LD (Linkage Disequilibrium) before running PCA and Admixture analyses
# to account for the assumption of indepence for loci

vcf.thin25 <- distance_thin(vcf25, min.distance = 500)
# meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
# dim(meta)
# 
# meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]) , ]
# dim(meta2)

write.vcf(vcf.thin25, "~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.25.vcf.gz")
# hide the uncompressed vcf because too big for github

setwd("~/Projects/eco_genomics/population_genomics/")

system("gunzip -c ~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.25.vcf.gz > ~/vcf_final.filtered.thinned.25.vcf")

geno25 <- vcf2geno(input.file = "/gpfs1/home/c/s/csaratka/vcf_final.filtered.thinned.25.vcf",
                 output.file = "/gpfs1/home/c/s/csaratka/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.25.geno")

CentPCA25 <- LEA::pca("~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.25.geno", scale = TRUE)

CentPCA25 <- load.pcaProject("vcf_final.filtered.thinned.25.pcaProject")

show(CentPCA)

plot(CentPCA)

vcf.pcadapt25 <- read.pcadapt("/gpfs1/home/c/s/csaratka/vcf_final.filtered.thinned.25.vcf", type = "vcf")

#vcfR <- read.vcfR("~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

vcfR25 <- read.vcfR("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.vcf")

# meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
# meta2 <- meta[meta$id %in% colnames(vcfR@gt[,-1]),]



pcadapt.pca25 <- pcadapt(vcf.pcadapt25,
                       K=2,
                       method = "componentwise",
                       min.maf = 0.01,
                       LD.clumping = list(size=500, thr=0.2)) #min.maf = minimum minor frequency allele

summary(pcadapt.pca25)
plot(pcadapt.pca25, options = "scores",
     pop = meta25$region,
     i = 1, j = 2, K=2)

view(head(vcfR25@fix))
vcfR.fix25 <- as.data.frame(vcfR25@fix[,1:2])


chr.main25.admix <-  unique(vcfR.fix25$CHROM)[1:8]
chrnum25.admix <- as.data.frame(cbind(chr.main25.admix, seq(1,8,1)))

Pval25 <- pcadapt.pca25$pvalues

pcadapt.MHplot25 <- cbind(vcfR.fix25,Pval25)
pcadapt.MHplot25 <- left_join(chrnum25.admix, pcadapt.MHplot25, join_by(chr.main25.admix==CHROM))

pcadapt.MHplot <- pcadapt.MHplot %>% 
  mutate(SNP=paste0(chr.main,"_",POS))

str(pcadapt.MHplot)

pcadapt.MHplot$V2 = as.numeric(pcadapt.MHplot$V2)
pcadapt.MHplot$POS = as.numeric(pcadapt.MHplot$POS)
pcadapt.MHplot$pPC1 = as.numeric(pcadapt.MHplot[,4])
pcadapt.MHplot$pPC2 = as.numeric(pcadapt.MHplot[,5])

pcadapt.MHplot <- pcadapt.MHplot %>% drop_na(pPC1)

manhattan(pcadapt.MHplot,
          chr = "V2",
          bp = "POS",
          p = "pPC1",
          col = c("blue4", "orange3"),
          logp = T,
          ylab = "-log10 p-value",
          genomewideline = F,
          main = "PCAdapt genome scan for selection (PC1)")

manhattan(pcadapt.MHplot,
          chr = "V2",
          bp = "POS",
          p = "pPC2",
          col = c("blue4", "orange3"),
          logp = T,
          ylab = "-log10 p-value",
          genomewideline = F,
          main = "PCAdapt genome scan for selection (PC2)")

View(pcadapt.MHplot %>% 
       filter(pPC1<quantile(pcadapt.MHplot$pPC1, 0.001))) %>% 
  select(chr.main, POS, pPC1)




CentAdmix25 <-snmf("outputs/vcf_final.filtered.thinned.25.geno",
                 K = 1:10,
                 entropy = TRUE,
                 repetitions = 3,
                 project = "new") #if adding to analysis, you can choose project = "continue"

par(mfrow = c(2,1))

plot(CentAdmix25, main = "SNMF", col = "blue4")
#Pay attention to the "elbow" of the graph, the values you want to focus on 
## Elbow is 2-4


#plot(CentPCA$eigenvalues[1:10], ylab = "Eigen values", xlab = "Number of PCs", col = "blue4", main="PCA")
#dev.off() #resets plotting window


myK25 = 3

CE25 = cross.entropy(CentAdmix25, K = myK25)
best25 = which.min(CE25) #lowest cross entropy value

myKQ25 = Q(CentAdmix25, K = myK25, run = best25) 
#dataset of each individuals' percentage ancestry from each K group

myKQmeta25 = cbind(myKQ25, meta25)

my.colors25 = c("blue4", "gold", "tomato", "lightblue", "olivedrab")

myKQmeta25 = as_tibble(myKQmeta25) %>% 
  group_by(continent) %>% 
  arrange(region, pop, .by_group = TRUE)


#pdf("figures/Admixture_K4.pdf", width = 10, height = 5)
barplot(as.matrix(t(myKQmeta25[ ,1:myK25])),
        border = NA,
        space = 0,
        col = my.colors25[1:myK25],
        xlab = "Geographic regions",
        ylab = "Ancestry proportions",
        main = paste0("Ancestry matrix K =", myK25))
axis(1,
     at = 1:length(myKQmeta25$region),
     labels = myKQmeta25$region,
     tick = F,
     cex.axis = 0.5,
     las = 3)
dev.off()


#############################################################################################################################################
#2. a-d for vcf.filt.indMiss.5

vcf50 <- read.vcfR("~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.5.vcf.gz")

meta50<- meta[meta$id %in% colnames(vcf50@gt[,-1]),] # @ is like $ by for this file class
dim(meta50)

vcf.div50 <- genetic_diff(vcf50,
                        pops = as.factor(meta50$region),
                        method = "nei")

str(vcf.div50)

chr.main50 <- unique(vcf.div50$CHROM)[1:8]

chrnum50 <- as.data.frame(cbind(chr.main50,seq(1,8,1)))

vcf.div.MHplot50 <- left_join(chrnum50, vcf.div50, join_by(chr.main50==CHROM))
dim(vcf.div.MHplot50)

vcf.div.MHplot50 <- vcf.div.MHplot50 %>% 
  filter(Gst>0) %>% 
  mutate(SNP = paste0(chr.main50, "_", POS))
str(vcf.div.MHplot50)

vcf.div.MHplot50$V2= as.numeric(vcf.div.MHplot50$V2)
vcf.div.MHplot50$POS= as.numeric(vcf.div.MHplot50$POS)

#manhattan(vcf.div.MHplot50,
#          chr = "V2",
          # bp = "POS",
          # p = "Gst",
          # col = c("blue4", "orange3"),
          # logp = FALSE,
          # ylab = "Fst among regions",
          # suggestiveline = quantile(vcf.div.MHplot50$Gst, 0.5))

#suggestline - suggests that the region above the line is experiencing a higher level of differentiation
#write.csv(vcf.div.MHplot, "~/Projects/eco_genomics/population_genomics/outputs/Genetic_Diff_byRegion",
#quote = F,
#row.names = F)

num_indSNP50 <- tibble(Individuals = nrow(vcf.div.MHplot50), Loci = n_distinct(vcf.div.MHplot50$SNP))

num_indSNP50|>
  gt() |>
  tab_header(title = "# of inds and loci (50% filter)")

names(vcf.div.MHplot50) 
vcf.div.MHplot50 %>% 
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  ggplot(aes(x = value, fill = name)) +
  geom_histogram(position = "identity",alpha = 0.5, bins = 50)+
  labs(x = "Gene diversity (Hs) withing Regions (50% filter)", y = "Counts of SNPs", title = "Genome-wide expected heterozygosity (Hs)", fill = "Regions")

#ggsave("Histogram_GenomDiversity_byRegion.5.pdf", 
 #      path = "~/Projects/eco_genomics/population_genomics/figures")

filter_summary_50 <- vcf.div.MHplot50 %>% 
  as_tibble() %>%   
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>% 
  filter(value!=0 & value<0.5) %>% 
  summarise(avg_Hs=mean(value), sd_Hs = sd(value), N_Hs=n())

Hs50  <- vcf.div.MHplot50[,4:9]

Hs_SD50table <- filter_summary_50 |>
  gt() |>
  tab_header(
    title = "Hs and SD for each region (50% filter)"
  ) |>
  cols_label(name = "Region",
             avg_Hs = "Average (Hs)",
             sd_Hs = "SD (HS)", 
             N_Hs = "N (Hs)")

#c.

CEU50_0 <- sum(Hs50$Hs_CEU==0)
NE50_0 <- sum(Hs50$Hs_NE ==0)
NEU50_0 <- sum(Hs50$Hs_NEU ==0)
PNW50_0 <- sum(Hs50$Hs_PNW ==0)
SEU50_0 <- sum(Hs50$Hs_SEU ==0)
WEU50_0 <- sum(Hs50$Hs_WEU ==0)

CEU50 <- sum(Hs50$Hs_CEU>0)
NE50 <- sum(Hs50$Hs_NE >0)
NEU50 <- sum(Hs50$Hs_NEU >0)
PNW50 <- sum(Hs50$Hs_PNW >0)
SEU50 <- sum(Hs50$Hs_SEU >0)
WEU50 <- sum(Hs50$Hs_WEU >0)

Nloci_50_0 <-  tibble(Region = c("CEU", "NE", "NEU", "PNW", "SEU", "WEU"),
                      Equal_zero = c(CEU50_0, NE50_0, NEU50_0, PNW50_0, SEU50_0, WEU50_0),
                      Above_zero = c(CEU50, NE50, NEU50, PNW50, SEU50, WEU50))

Nloci_50table <- Nloci_50_0|>
  gt() |>
  tab_header(
    title = "Number of loci = 0 or > 0 (50% filter)"
  ) |>
  cols_label(Region = "Region",
             Equal_zero = "Equal 0",
             Above_zero = "Above 0")

vcf.thin50 <- distance_thin(vcf50, min.distance = 500)
# meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
# dim(meta)
# 
# meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]) , ]
# dim(meta2)

write.vcf(vcf.thin50, "~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.5.vcf.gz")
# hide the uncompressed vcf because too big for github

system("gunzip -c ~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.5.vcf.gz > ~/vcf_final.filtered.thinned.5.vcf")

geno25 <- vcf2geno(input.file = "/gpfs1/home/c/s/csaratka/vcf_final.filtered.thinned.5.vcf",
                   output.file = "/gpfs1/home/c/s/csaratka/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.5.geno")


CentAdmix50 <-snmf("outputs/vcf_final.filtered.thinned.5.geno",
                   K = 1:10,
                   entropy = TRUE,
                   repetitions = 3,
                   project = "new") #if adding to analysis, you can choose project = "continue"

par(mfrow = c(2,1))
plot(CentAdmix50, main = "SNMF", col = "blue4")
#Pay attention to the "elbow" of the graph, the values you want to focus on 
## Elbow is 2-4


#plot(CentPCA$eigenvalues[1:10], ylab = "Eigen values", xlab = "Number of PCs", col = "blue4", main="PCA")
#dev.off() #resets plotting window


myK50 = 3

CE50 = cross.entropy(CentAdmix50, K = myK50)
best50 = which.min(CE50) #lowest cross entropy value

myKQ50 = Q(CentAdmix50, K = myK50, run = best50) 
#dataset of each individuals' percentage ancestry from each K group

myKQmeta50 = cbind(myKQ50, meta50)

my.colors50 = c("blue3", "goldenrod", "red", "cyan", "darkgreen")

myKQmeta50 = as_tibble(myKQmeta50) %>% 
  group_by(continent) %>% 
  arrange(region, pop, .by_group = TRUE)


#pdf("figures/Admixture_K4.pdf", width = 10, height = 5)
barplot(as.matrix(t(myKQmeta50[ ,1:myK50])),
        border = NA,
        space = 0,
        col = my.colors50[1:myK50],
        xlab = "Geographic regions",
        ylab = "Ancestry proportions",
        main = paste0("Ancestry matrix K =", myK50))
axis(1,
     at = 1:length(myKQmeta50$region),
     labels = myKQmeta50$region,
     tick = F,
     cex.axis = 0.5,
     las = 3)
dev.off()