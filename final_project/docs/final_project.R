library(DESeq2)
library(tidyverse)
library(ape)
library(vcfR)

dna <-  ape::read.dna("/users/a/r/armccrac/data/AK_pycno_metagen/rawdata", format = "fasta")
