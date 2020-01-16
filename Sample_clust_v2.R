#######################################################
#
# Script name: Sample_clust_v2
# Language: R
# Program: hclust + PCA
# Purpose: Hierarchial clustering and PCA of samples based on expression
# Author: Arindam Ghosh
# Date: 24 September 2019
#
#######################################################

library(dendextend)
library(edgeR)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(pheatmap)
# ****************************************************************************
countsPC <- read.table("CountsPC_Gene_MultimapOverlap.tsv", header=T, row.names=c(1))
factors <- as.data.frame(read.table("factors.txt", sep="\t", header = TRUE, row.names=c(1)))

SamplesUsed <- "NP123456789"
OutputNumber <- 11

keep <- rowSums(cpm(countsPC)>1) >1
countsPC_filt <- countsPC[keep,]


# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

geneLength <- read.csv("/home/bioinfo/Documents/Data/ReferenceGenome/Human_Ensembl84/NoiSeq/GRCh38_GeneLenth.csv", row.names=c(1))
geneLength_filt <- geneLength[row.names(countsPC_filt),, drop=FALSE]
DGElist <- DGEList(counts = countsPC_filt, group = factors$SampleType, genes = data.frame(Length = geneLength_filt$Length))
DGElist <- calcNormFactors(DGElist, method="TMM")
countsPC_norm <- rpkm(DGElist, log=TRUE)
#countsPC_norm <- cpm(countsPC_filt, log=TRUE)


# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


dist <- dist(t(countsPC_norm), method="euclidean")
cluster <- hclust(dist)

dend <- as.dendrogram(cluster)

#Add colors to label
label_colors_to_use <- as.numeric(factors$SampleType) #Input groups of sample type for label color
label_colors_to_use <- label_colors_to_use[order.dendrogram(dend)] #Order group annotation as per in dendrogram

#Use to add custom label color
for(i in 1:length(label_colors_to_use)){ 
 label_colors_to_use[i] <- switch(as.numeric(label_colors_to_use[i]), "#e6194B", "#000000") 
}

#Add colored bars showing study accession
StudyAccession_colors_to_use <- as.numeric(factors$StudyAccession) #Input groups of sample type for label color
StudyAccession_colors_to_use <- StudyAccession_colors_to_use[order.dendrogram(dend)] #Order group annotation as per in dendrogram

#Use to add custom color to bar
for(i in 1:length(StudyAccession_colors_to_use)){ 
 StudyAccession_colors_to_use[i] <- switch(as.numeric(StudyAccession_colors_to_use[i]), "#e6194B", "#098f19", "#42d4f4", "#000075", "#f58231", "#911eb4", "#ffff00", "#f032e6", "#aaffc3")
}


dend <- set(dend, "labels_colors", label_colors_to_use) #Apply colors to labels
dend <- hang.dendrogram(dend, hang=-1) #Leaves at same height
dend <- set(dend, "branches_lwd", 2) #Thick branches
dend <- set(dend, "labels_cex", 1) #label font size

#Plot the dendrogram
tiff(filename=paste0("Sample_hclust", OutputNumber, ".tiff"), height=10, width=20, units='in', res=300)
par(mar = c(10,5,3,1))
plot(dend, main = paste0("Sample Clustering\nDist.method=Eucleadian Samples=", SamplesUsed, " Normalization=RPKM/DGEList", ylab = "Height")) 
colored_bars(colors = StudyAccession_colors_to_use, dend = dend, rowLabels = "Project", sort_by_labels_order = FALSE)
dev.off()


#PCA analysis
pca <- prcomp(t(countsPC_norm))
percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
pca_data <- data.frame(pca$x, SampleType=factors$SampleType, StudyAccession=factors$StudyAccession)
tiff(filename=paste0("Sample_PCA", OutputNumber, ".tiff"), height=10, width=10, units='in', res=300)
ggplot(pca_data,aes(x=PC1,y=PC2, shape=SampleType, col=StudyAccession )) +
 geom_point(size = 4) +
 labs(title="Sample PCA", subtitle=paste0("Samples = ", SamplesUsed, " Normalization=RPKM/DGEList"))+
 xlab(paste0("PC1: ", percentage[1], "% variance")) +
 ylab(paste0("PC2: ", percentage[2], "% variance")) +
 theme(
       axis.text = element_text(size=15), 
       axis.title = element_text(size=16), 
       plot.title = element_text(size=18, hjust=0.5), 
       plot.subtitle = element_text(size=14, hjust=0.5),
       legend.text = element_text(size=14, hjust=0.5),
       legend.title = element_text(size=16, hjust=0.5)     
      )
dev.off()

#Sample correlation
tiff(filename="Sample_heatmap1.tiff", height=10, width=10, units='in', res=300)
distRL = dist(t(assay(vsd)))
mat=as.matrix(distRL)
hmcol = colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13,13))
dev.off()
tiff(filename="Sample_heatmap2.tiff", height=10, width=10, units='in', res=300)
pheatmap(mat, clustering_distance_rows=distRL, clustering_distance_cols=distRL, col=rev(hmcol))
dev.off()
