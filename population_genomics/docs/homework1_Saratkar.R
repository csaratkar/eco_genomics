library(vcfR)
library(tidyverse)
library(qqman)
library(SNPfiltR)

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

list.files()

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
                                      cutoff = 0.9)

write.vcf(vcf.filt.indMiss.25,
          "~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.25.vcf.gz")


vcf.filt.indMiss.5 <- missing_by_sample(vcf.filt, 
                                      popmap = meta2,
                                      cutoff = 0.5)
write.vcf(vcf.filt.indMiss.5,
          "~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.5.vcf.gz")

#############################################################################################################################################
#2. a-d for vcf.filt.indMiss.25

vcf25 <- read.vcfR("~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.25.vcf.gz")

meta25<- meta[meta$id %in% colnames(vcf25@gt[,-1]),] # @ is like $ by for this file class
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

manhattan(vcf.div.MHplot25,
          chr = "V2",
          bp = "POS",
          p = "Gst",
          col = c("blue4", "orange3"),
          logp = FALSE,
          ylab = "Fst among regions",
          suggestiveline = quantile(vcf.div.MHplot25$Gst, 0.5))

#suggestline - suggests that the region above the line is experiencing a higher level of differentiation
#write.csv(vcf.div.MHplot, "~/Projects/eco_genomics/population_genomics/outputs/Genetic_Diff_byRegion",
          #quote = F,
          #row.names = F)

#1
num_indSNP25 <- tibble(Individuals = nrow(vcf.div.MHplot25), Loci = n_distinct(vcf.div.MHplot25$SNP))

names(vcf.div.MHplot25) 
vcf.div.MHplot25 %>% 
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  ggplot(aes(x = value, fill = name)) +
  geom_histogram(position = "identity",alpha = 0.5, bins = 50)+
  labs(x = "Gene diversity (Hs) withing Regions (25% filter)", y = "Counts of SNPs", title = "Genome-wide expected heterozygosity (Hs)", fill = "Regions")

ggsave("Histogram_GenomDiversity_byRegion.25.pdf", 
       path = "~/Projects/eco_genomics/population_genomics/figures")

filter_summary_25 <- vcf.div.MHplot25 %>% 
  as_tibble() %>%   
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>% 
  filter(value!=0 & value<0.5) %>% 
  summarise(avg_Hs=mean(value), sd_Hs = sd(value), N_Hs=n())

Hs_sum25 <- Hs25 %>% 
  as_tibble() %>%   
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>% 
  sum(Hs25) 

Hs25 <- vcf.div.MHplot25 %>% 
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  group_by(name)

Hs25 <- vcf.div.MHplot25[,4:9]

Hs_sum25 <- sum(Hs25 == 0)

CEU25 <- sum(Hs25$Hs_CEU==0)
NE25 <- sum(Hs25$Hs_NE ==0)
NEU25 <- sum(Hs25$Hs_NEU ==0)
PNW25 <- sum(Hs25$Hs_PNW ==0)

Hs_sum25 <- data_frame(matrix(nrow = 6, ncol = 3))
colnames(Hs_sum25) <- c("names","equal_zero","above_zero")

Hs_sum25[1,2] <- 3

for(i in 1:ncol(Hs_sum25)){
  Hs_sum25["names[,2]"] <- sum(Hs25[,i]==0)
}

#############################################################################################################################################
#2. a-d for vcf.filt.indMiss.5

vcf50 <- read.vcfR("~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.5.vcf.gz")

meta50<- meta[meta$id %in% colnames(vcf50@gt[,-1]),] # @ is like $ by for this file class
dim(meta50)

vcf.div50 <- genetic_diff(vcf50,
                        pops = as.factor(meta4$region),
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

manhattan(vcf.div.MHplot50,
          chr = "V2",
          bp = "POS",
          p = "Gst",
          col = c("blue4", "orange3"),
          logp = FALSE,
          ylab = "Fst among regions",
          suggestiveline = quantile(vcf.div.MHplot50$Gst, 0.5))

#suggestline - suggests that the region above the line is experiencing a higher level of differentiation
#write.csv(vcf.div.MHplot, "~/Projects/eco_genomics/population_genomics/outputs/Genetic_Diff_byRegion",
#quote = F,
#row.names = F)

num_indSNP50 <- tibble(Individuals = nrow(vcf.div.MHplot50), Loci = n_distinct(vcf.div.MHplot50$SNP))

names(vcf.div.MHplot50) 
vcf.div.MHplot50 %>% 
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  ggplot(aes(x = value, fill = name)) +
  geom_histogram(position = "identity",alpha = 0.5, bins = 50)+
  labs(x = "Gene diversity (Hs) withing Regions (50% filter)", y = "Counts of SNPs", title = "Genome-wide expected heterozygosity (Hs)", fill = "Regions")

ggsave("Histogram_GenomDiversity_byRegion.5.pdf", 
       path = "~/Projects/eco_genomics/population_genomics/figures")

filter_summary_50 <- vcf.div.MHplot50 %>% 
  as_tibble() %>%   
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>% 
  filter(value!=0 & value<0.5) %>% 
  summarise(avg_Hs=mean(value), sd_Hs = sd(value), N_Hs=n())
