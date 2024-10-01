# Coding and data notes for the population genetics module

## Author: Chanchal Saratkar

### 09-10-2024 - Intro to Centaurea GBS data and working with VCF files in 01 VCF filtering.r

We'll be analyzing the GBS data from three regions (Eu, Northeast, Pacific Northwest) starting today with variant call format files (VCF's)

- .gz is g zip file
- ll is the shell command for long list
- zcat is a shell command to view a zipped file
- zcat file.gz | head is the shell command to read an unfiltered gz file
- fastQ files include a sequence line and a Q-score line
- The Q-scores in a fastQ file (.fq) are written as a varitety of symbols with letters being at the higher end

### 09-12-2024 Viewing VCF files and talking about in 01 VCF filtering.r

- N50 - # of base pairs in 50% of the reads

Shell commands

- spack load {package} - load a package
- samtools - indexes reads
- bcftools - genotypes reads

### 09-17-2024 Filtering VCF in 01 VCF filtering.r

- quantile() - makes quartiles
- heatmap.bp() - heat map
- SNPflitR - package that filters SNPS
- hard_filter() - filters based on depth
- max_depth() - establishes max depth
- missing_by_sample() - filters individuals with missing data
- missing_by_SNP() - filters loci
- write.vcf() - makes VCF 

### 09-19-2024 Created Manhattan plot in 02_Diversity_Differentiation.R

- X11.options(type="cairo") - helps fix plotting issues
- @ is used for vcfs instead of $
- genetic_diff() - helps genetically differentiate
- manhattan() - creates Manhattan plots

### 09-24-2024 Created PCA plot

- LEA - package for creating PCAs
- PCAs are just based on math

### 09-26-2024 - Admixture Analysis

- Admixture Analysis (Structure) - has genetic model (H-W equilibrium)
- Helps find population structure
- Groups the population into clusters which is decided by K

1. K has to be chosen
2. Assign individuals to one of those K groups, can have prior knowledge or not
3. Calc allele frequencies
4. Calc 2pq -> observed frequency of heterozygosity
5. Go back to #2

- Q = fractional ancestry
- Creates matrix that gives percentages about what group an individual belongs to
- Cross-validation - train a model on 80% of the data and then give the 20% to analysis then move on to the next 20%
- The K value that has the lowest cross-validation value
- No linkage disequilibrium

### 10-01-2024 - Admixture Analysis cont.

- Pay attention to the values at the "elbow" of the cross-entropy score plot
- Goal is to reduce complexity while increasing 