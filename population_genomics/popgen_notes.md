# Coding and data notes for the population genetics module

## Author: Chanchal Saratkar

### 09-10-2024 - Intro to Centaurea GBS data and working with VCF files

We'll be analyzing the GBS data from three regions (Eu, Northeast, Pacific Northwest) starting today with variant call format files (VCF's)

- .gz is g zip file
- ll is the shell command for long list
- zcat is a shell command to view a zipped file
- zcat file.gz | head is the shell command to read an unfiltered gz file
- fastQ files include a sequence line and a Q-score line
- The Q-scores in a fastQ file (.fq) are written as a varitety of symbols with letters being at the higher end

### 09-12-2024 Viewing VCF files and talking about 

- N50 - # of base pairs in 50% of the reads

Shell commanads

- spack load {package} - load a package
- samtools - indexes reads
- bcftools - genotypes reads

### 09-17-2024 Filtering VCF

- quantile() - makes quartiles
- heatmap.bp() - heat map
- SNPflitR - package that filters SNPS
- hard_filter() - filters based on depth
- max_depth() - establishes max depth
- missing_by_sample() - filters individuals with missing data
- missing_by_SNP() - filters loci
- write.vcf() - makes VCF 