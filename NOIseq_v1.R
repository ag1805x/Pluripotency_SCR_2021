#######################################################
#
# Script name: NOISeq_v1.R
# Language: R
# Program: NOISeq
# Purpose: Explorative analysis of count data
# Author: Arindam Ghosh
# Date: 4 September 2019
#
#######################################################

library(NOISeq)
counts <- as.data.frame(read.table("PRJEB4879_counts_clean_MultiOverlap.tsv"))
dim(counts)

factors <- as.data.frame(read.table("factors.txt", sep="\t", header = TRUE, row.names=c(1)))
geneLength <- as.data.frame(read.csv("/home/dell/Documents/Arindam/Work/ReferenceGenome/Human_84/NoiSeq/GRCh38_GeneLenth.csv", header=T, row.names=c(1)))
geneGC <- as.data.frame(read.table("/home/dell/Documents/Arindam/Work/ReferenceGenome/Human_84/NoiSeq/GRCh38_GeneGC.txt", sep="\t", header=T, row.names=c(1)))
geneGC <- geneGC[,-1, FALSE]
geneBioTypes <- as.data.frame(read.csv("/home/dell/Documents/Arindam/Work/ReferenceGenome/Human_84/NoiSeq/GRCh38_GeneBioType.csv", header=T, row.names=c(1)))
geneChrom <- as.data.frame(read.csv("/home/dell/Documents/Arindam/Work/ReferenceGenome/Human_84/NoiSeq/GRCh38_Chrom.csv", header=T, row.names=c(1)))


data <- readData(data = counts, length = geneLength, gc = geneGC, biotype= geneBioTypes, chromosome = geneChrom, factors = factors)
data


### Biotype detection plot (per sample)
biodetection <- dat(data, k = 0, type = "biodetection", factor = NULL)
par(mfrow = c(1, 2))
explo.plot(biodetection, samples = c(1, 2), plottype = "persample")


### Biotype detection plot (compare biotype detection in two samples)
par(mfrow = c(1, 2))
explo.plot(biodetection, samples = c(1, 2), toplot = "protein_coding", plottype = "comparison")

### Count distribution per sample
countsbio = dat(data, factor = NULL, type = "countsbio")
explo.plot(countsbio, toplot = "protein_coding", samples = NULL, plottype = "boxplot")

### PCA
pca <- dat(data, type="PCA")

explo.plot(pca, factor = "SampleType")

### QC reports 
QCreport(data, samples = NULL, factor = "SampleType", norm = FALSE)


