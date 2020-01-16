#######################################################
#
# Script name: DESeq2_v1
# Language: R
# Program: DESeq2
# Purpose: Perform differential gene expression analysis with DESeq2 and generate associated graphs
# Author: Arindam Ghosh
# Date: 11 September 2019
#
#######################################################


library(DESeq2)
library(edgeR)
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(geneplotter)



# Project ID for which analysis is being performed
ProjectID <- c("PRJNA286204")

# Read all gene counts
counts <- read.table("/home/dell/Documents/Arindam/WorkB/PRJNA286204/FeatureCounts_Gene_MultimapOverlap/PRJNA286204_counts_clean_Gene_MultimapOverlap.tsv", row.names=c(1), header=T)
sprintf ("Total Genes")
dim(counts)

# Read list of protein coding genes Ensembl IDs
protein_coding <- as.matrix(read.table("/home/dell/Documents/Arindam/ReferenceGenome/Human_Ensembl84/Annotation/gff/protein_coding"))
sprintf ("Total Protein Coding Genes")
dim(protein_coding)

# Extract counts for protein coding genes
counts_pc <- counts[protein_coding,]
dim(counts_pc)

# Filter genes having low counts
keep <- rowSums(cpm(counts_pc)>1) >1
counts_pc_filt <- counts_pc[keep,]
sprintf ("Number of genes after filtering low count genes")
dim(counts_pc_filt)

# Read metadata about samples
factors <- read.table("/home/dell/Documents/Arindam/WorkB/CombinedStudy/Data_Gene_MultimapOverlap/factors.txt", sep="\t", header=T, row.names=c(1))
factors <- factors[factors$StudyAccession == ProjectID,]

# DESeq2 differential expression analysis
dds <- DESeqDataSetFromMatrix(countData=counts_pc_filt, colData=factors, design=~SampleType)
dds <- DESeq(dds)
res <- results(dds)

sprintf ("DE test info")
res

write.table(res, file=paste(ProjectID, "_DESeq2_result.tsv", sep=""), sep="\t", col.names=NA) 

res_clean = na.exclude(as.data.frame(res))
sprintf ("Number of genes after NA exclude (DESeq2)")
dim(res_clean)

upreg_deseq <- res_clean[(res_clean$log2FoldChange>1 & res_clean$padj<0.05),]
sprintf ("Number of upregulated genes (DESeq2)")
dim(upreg_deseq)
write.table(upreg_deseq, file=paste(ProjectID, "_DESeq2_upreg.tsv", sep=""), sep="\t", col.names=NA)

downreg_deseq <- res_clean[(res_clean$log2FoldChange<(-1) & res_clean$padj<0.05),]
sprintf ("Number of downregulated genes (DESeq2)")
dim(downreg_deseq)
write.table(downreg_deseq, file=paste(ProjectID, "_DESeq2_downreg.tsv", sep=""), sep="\t", col.names=NA)



# PLOT SAMPLE TO SAMPLE DISTANCE
rld <- rlogTransformation(dds, blind=T)
tiff(filename=paste(ProjectID, "_Sample_heatmap1.tiff", sep=""), height=10, width=10, units='in', res=300)
distRL = dist(t(assay(rld)))
mat=as.matrix(distRL)
hmcol = colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13,13))
dev.off()
tiff(filename=paste(ProjectID, "_Sample_heatmap2.tiff", sep=""), height=10, width=10, units='in', res=300)
pheatmap(mat, clustering_distance_rows=distRL, clustering_distance_cols=distRL, col=rev(hmcol))
dev.off()



#PCA
tiff(filename=paste(ProjectID, "_Sample_PCA.tiff", sep=""), height=10, width=10, units='in', res=300)
pca <- DESeq2::plotPCA(rld, intgroup=c("SampleType"), returnData=TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
ggplot(pca,aes(x=PC1,y=PC2,color=SampleType, label=row.names(pca) )) +
geom_point(size = 3) +
xlab(paste0("PC1: ", percentVar[1], "% variance")) +
ylab(paste0("PC2: ", percentVar[2], "% variance"))
dev.off()

# Plot Dispersion Estimates
tiff(filename=paste(ProjectID, "_Dispersion.tiff", sep=""), height=10, width=10, units='in', res=300)
plotDispEsts(dds)
dev.off()

# MA Plot
tiff(filename=paste(ProjectID, "_MA.tiff", sep=""), height=10, width=10, units='in', res=300)
plotMA(res)
dev.off()

# Volcano Plot
tiff(filename=paste(ProjectID, "_Volcano.tiff", sep=""), height=10, width=10, units='in', res=300)
with(res_clean, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano Plot", col="grey", xlim=c(-10,10)))
with(subset(res_clean,padj<0.05 & log2FoldChange>1),points(log2FoldChange, -log10(padj), pch=20, col="green"))
with(subset(res_clean,padj<0.05 & log2FoldChange<(-1)),points(log2FoldChange, -log10(padj), pch=20, col="red"))
abline(h=2, lty=2)
abline(v=-1, lty=2)
abline(v=1, lty=2)
dev.off()

# ECDF Plot
tiff(filename=paste(ProjectID, "_ECDF.tiff", sep=""), height=10, width=10, units='in', res=300)
multiecdf(counts(dds, normalized=TRUE)[,], xlab="MeanCounts", xlim=c(0,1000))
dev.off()

# Plot Density
tiff(filename=paste(ProjectID, "_Density.tiff", sep=""), height=10, width=10, units='in', res=300)
multidensity(counts(dds, normalized=TRUE)[,], xlab="MeanCounts", xlim=c(0,1000))
dev.off()

















































