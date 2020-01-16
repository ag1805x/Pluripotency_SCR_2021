#######################################################
#
# Script name: FeatureCounts_v1.R
# Language: R
# Program: FeatureCounts
# Purpose: Estimation of gene counts from bam file
# Author: Arindam Ghosh
# Date: 4 September 2019
#
#######################################################

library(Rsubread)

bam.list <- dir(path="../Align/", pattern=".bam$", recursive=TRUE, full.names=TRUE)


fc <- featureCounts(bam.list, annot.ext = "/home/dell/Documents/Arindam/Work/ReferenceGenome/Human_84/Annotation/gtf/Homo_sapiens.GRCh38.84.gtf", isGTFAnnotationFile = T, GTF.featureType = 'gene', useMetaFeatures=T, allowMultiOverlap = T, countMultiMappingReads = T, nthreads = 10, isPairedEnd=F, strandSpecific=0)

write.table(fc$counts, "PRJNA319819_counts_Gene_MultimapOverlap.tsv", sep="\t", col.names=NA, quote=F)
write.table(fc$annotation, "PRJNA319819_annotation_Gene_MultimapOverlap.tsv", sep="\t", row.names=F, quote=F)
write.table(fc$stat, "PRJNA319819_stats_Gene_MultimapOverlap.tsv", sep="\t", col.names=NA, quote=F)

temp <- fc$counts
colnames(temp)
colnames(temp) <- substr(colnames(temp), 11, 20)
colnames(temp)
write.table(temp, "PRJNA319819_counts_clean_Gene_MultimapOverlap.tsv", sep="\t", col.names=NA, quote=F)












