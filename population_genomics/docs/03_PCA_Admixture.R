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

meta2 <- meta[meta$id %in% colnames(vcf.thin@gt[,-1]) , ]
dim(meta2)

write.vcf(vcf.thin, "outputs/vcf_final.filtered.thinned.vcf.gz")
 # hide the uncompressed vcf because too big for github

system("gunzip -c ~/Projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz > ~/vcf_final.filtered.thinned.vcf")

geno <- vcf2geno(input.file = "/gpfs1/home/c/s/csaratka/vcf_final.filtered.thinned.vcf",
                 output.file = "outputs/vcf_final.filtered.thinned.geno")

CentPCA <- LEA ::pca("outputs/vcf_final.filtered.thinned.geno", scale = TRUE)

plot(CentPCA$projections,
     col=as.factor(meta2$region))
