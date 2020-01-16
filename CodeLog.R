#### Code Log

#### 4Sep2019 # DESeq2 ###########################################################################################################################################################################


library(edgeR)
library(DESeq2)
library(geneplotter)
library(gplots)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

counts <- read.table("AllCounts_MultiOverlap.tsv", row.names=c(1), header=T)
dim(counts)
protein_coding = as.matrix(read.table("/home/dell/Documents/Arindam/Work/ReferenceGenome/Human_84/Annotation/gff/protein_coding"))
counts_pc <- counts[protein_coding,]
keep <- rowSums(cpm(counts_pc)>1) >1
counts_pc_filt <- counts_pc[keep, ] #16311


factors <- as.data.frame(read.table("factors.txt", sep="\t", header = TRUE, row.names=c(1)))
dds = DESeqDataSetFromMatrix(countData=counts_pc_filt, colData=factors, design=~SampleType)
dds = DESeq(dds)
res = results(dds)
res_clean = na.exclude(as.data.frame(res))


rld = rlogTransformation(dds, blind=T)
pca = DESeq2::plotPCA(rld, intgroup=c("SampleType"), returnData=TRUE)


#####

> protein_coding = as.matrix(read.table("/home/dell/Documents/Arindam/Work/ReferenceGenome/Human_84/Annotation/gff/protein_coding"))
> counts_pc <- counts[protein_coding,]
> keep <- rowSums(cpm(counts_pc)>1) >1


##### 6Sep2019 Merge counts ###########################################################################################################################################################################
NP1 <- read.table("/home/dell/Documents/Arindam/WorkB/PRJNA286204/FeatureCounts_Gene_MultimapOverlap/PRJNA286204_counts_clean_Gene_MultimapOverlap.tsv")

NP2 <- read.table("/home/dell/Documents/Arindam/WorkB/PRJNA360413/FeatureCounts_Gene_MultimapOverlap/PRJNA360413_counts_clean_Gene_MultimapOverlap.tsv")

NP3 <- read.table("/home/dell/Documents/Arindam/WorkB/PRJNA268222/FeatureCounts_Gene_MultimapOverlap/PRJNA268222_counts_clean_Gene_MultimapOverlap.tsv")

NP4 <- read.table("/home/dell/Documents/Arindam/WorkB/PRJNA339108/FeatureCounts_Gene_MultimapOverlap/PRJNA339108_counts_clean_Gene_MultimapOverlap.tsv")

NP5 <- read.table("/home/dell/Documents/Arindam/WorkB/PRJEB7132/FeatureCounts_Gene_MultimapOverlap/PRJEB7132_counts_clean_Gene_MultimapOverlap.tsv")

NP6 <- read.table("/home/dell/Documents/Arindam/WorkB/PRJEB20388/FeatureCounts_Gene_MultimapOverlap/PRJEB20388_counts_clean_Gene_MultimapOverlap.tsv")

NP7 <- read.table("/home/dell/Documents/Arindam/WorkB/PRJEB4879/FeatureCounts_Gene_MultimapOverlap/PRJEB4879_counts_clean_Gene_MultimapOverlap.tsv")

NP8 <- read.table("/home/dell/Documents/Arindam/WorkB/PRJNA326944/FeatureCounts_Gene_MultimapOverlap/PRJNA326944_counts_clean_Gene_MultimapOverlap.tsv")

NP9 <- read.table("/home/dell/Documents/Arindam/WorkB/PRJNA319819/FeatureCounts_Gene_MultimapOverlap/PRJNA319819_counts_clean_Gene_MultimapOverlap.tsv")


##### 9 September 2019 ###########################################################################################################################################################################


countsPC <- read.table("CountsPC_Gene_MultimapOverlap.tsv", header=T, row.names=c(1))
factors <- as.data.frame(read.table("factors.txt", sep="\t", header = TRUE, row.names=c(1)))
> ddsMat <- DESeqDataSetFromMatrix(countData=countsPC, colData=factors, design=~SampleType)
> nrow(ddsMat) #19826
> ddsMat <- ddsMat[ rowSums(counts(ddsMat))>1,]
> nrow(ddsMat) #19491
> vsd <- vst(ddsMat, blind = FALSE)
> rld <- rlog(ddsMat, blind = FALSE)
> sampleDists <- dist(t(assay(vsd)))

> library("pheatmap")
> library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )

rownames(sampleDistMatrix) <- paste( vsd$SampleType, rownames(sampleDistMatrix), sep = " - " )

colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col = colors)
plotPCA(vsd, intgroup = c("SampleType", "StudyAccession"))
pcaData <- plotPCA(vsd, intgroup = c( "SampleType", "StudyAccession"), returnData = TRUE)

> library(ggplot2)

ggplot(pcaData, aes(x = PC1, y = PC2, color = SampleType, shape = SampleAccession)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() #error

ddsMat <- DESeq(ddsMat)

res <- results(ddsMat)


#### 13Sep2019 Clustering ###########################################################################################################################################################################

#library(ggdendro)
#library(ggplot2)
library(dendextend)

countsPC <- read.table("CountsPC_Gene_MultimapOverlap.tsv", header=T, row.names=c(1))
factors <- as.data.frame(read.table("factors.txt", sep="\t", header = TRUE, row.names=c(1)))

#colnames(countsPC) <- paste(factors$SampleType, colnames(countsPC), sep = "-") 
dist <- dist(t(countsPC), method="euclidean")
cluster <- hclust(dist)


# main, sub, xlab, ylab


dend <- as.dendrogram(cluster)


#Add colors to label
label_colors_to_use <- as.numeric(factors$SampleType) #Input groups of sample type for label color
label_colors_to_use <- label_colors_to_use[order.dendrogram(dend)] #Order group annotation as per in dendrogram

#Use to add custom label color
for(i in 1:length(label_colors_to_use)){ + 
  if(label_colors_to_use[i] == 1){ +
       label_colors_to_use[i] = "#00022e" +
  }else if(label_colors_to_use[i] == 2){ +
       label_colors_to_use[i] = "#004400" +
  } +
 }

#Add colored bars showing study accession
StudyAccession_colors_to_use <- as.numeric(factors$StudyAccession) #Input groups of sample type for label color
StudyAccession_colors_to_use <- StudyAccession_colors_to_use[order.dendrogram(dend)] #Order group annotation as per in dendrogram

dend <- set(dend, "labels_colors", label_colors_to_use) #Apply colors to labels
dend <- hang.dendrogram(dend, hang=-1) #Leaves at same height
dend <- set(dend, "branches_lwd", 2) #Thick branches
dend <- set(dend, "labels_cex", 1) #label font size


par(mar = c(10,2,1,1))
colored_bars(colors = StudyAccession_colors_to_use, dend = dend, rowLabels = "Project", sort_by_labels_order = FALSE))


plot(dend) #Plot the dendrogram


##### 19 September 2019 ###################################################################################################################################################################################

# Remove certain projects from sample matrix

list <- c("PRJNA319819", "PRJNA326944", "PRJNA339108")
remove <- vector()

for(i in 1:length(countsPC)){
+ if(as.character(factors$StudyAccession[i]) %in% list){
+ remove <- c(remove, i)
+ }
+ }

countsPC_new <- countsPC[, - remove]
factors_new <- factors[- remove,]   

# PCA plot



> pca <- prcomp(t(countsPC))
> plot(pca$x)

> percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)

> pca_data <- as.data.frame(pca$x)
> ggplot(pca_data,aes(x=PC1,y=PC2 )) +
+ geom_point(size = 3) +
+ xlab(paste0("PC1: ", percentage[1], "% variance")) +
+ ylab(paste0("PC2: ", percentage[2], "% variance"))




##### 23 October 2019 ###########################################################################################################################################################################
# Silhouettes with hclust
library(cluster)
hsil <- silhouette(cutree(cluster, k=2), dist)
summary(hsil)$avg.width


##### 19 November 2019 ###########################################################################################################################################################################
# Heatmap of gene expression



tiff("Heatmap_DEG.tif", height=10, width=20, units='in', res=300)
heatmap.2(countsPC_batch, trace="none", col=redgreen(75), scale="row")






# Box plot of normalized and unnormalized counts
tiff("Boxplot.tif", height=10, width=20, units='in', res=300)
par(mfrow=c(2,1))
boxplot(log2(countsPC_filt+1), outline=FALSE, col=box_color, main="Unnormalized counts", xlab="Samples", ylab="counts")
boxplot(countsPC_batch, outline=FALSE, col=box_color, main="Normalized counts", xlab="Samples", ylab="counts")
dev.off()



##### 22 November 2019 ###########################################################################################################################################################################
# WGCNA scale free plot


k <- softConnectivity(datExpr, power=9, blockSize=20000)
scaleFreePlot(k)


##### 9 December 2019 ###########################################################################################################################################################################
# Sample to Sample heatmap

tiff(filename="Sample_heatmap1.tiff", height=10, width=10, units='in', res=300)
distRL = dist(t(assay(vsd)))
mat=as.matrix(distRL)
hmcol = colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13,13))
dev.off()
tiff(filename="Sample_heatmap2.tiff", height=10, width=10, units='in', res=300)
pheatmap(mat, clustering_distance_rows=distRL, clustering_distance_cols=distRL, col=rev(hmcol))
dev.off()


##### 10 December 2019
# TOM plot color reverse

library(gplots)
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
tiff("Network heatmap plot2.tif", width=14, height=10, units="in", res=100)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes", col=myheatcol)
dev.off()

# Module Sample type relation table

Data_ModuleTraitAssociation_SampleType <- data.frame(corr=as.data.frame(moduleTraitCor)$SampleType, p.val=as.data.frame(moduleTraitPvalue)$SampleType, row.names=row.names(as.data.frame(moduleTraitCor)))
write.table(Data_ModuleTraitAssociation_SampleType, "Data_ModuleTraitAssociation_SampleType.tsv", sep="\t", col.names=NA)



####################


intModules <- c("blue", "turquoise")
for (module in intModules)
{
modGenes <- (moduleColors==module)
modGeneID <- connectivity[modGenes, ]
fileName <- paste("Connectivity-", module, ".txt", sep="");
write.table(as.data.frame(modGeneID), file = fileName, sep="\t", row.names = TRUE, col.names = NA)
}




##### 19 December 2019 ###########################################################################################################################################################################
# WGCNA test

library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
library(DDPNA)
WNGCAadjust <- wgcnatest(datExpr, maxBlockSize=20000, TOMType="unsigned")



# Clustering of raw counts
distRL_raw = dist(t(log2(countsPC_filt+1)))
mat_raw =as.matrix(distRL_raw)
pheatmap(mat_raw, clustering_distance_rows=distRL_raw, clustering_distance_cols=distRL_raw, col=rev(hmcol))




##### 23 December 2019 ###########################################################################################################################################################################
# Names of genes
gene.annotation0 <- read.csv("/home/dell/Documents/Arindam/ReferenceGenome/Human_Ensembl84/Annotation/gtf/GeneAnnotation.csv", row.names=c(1))
upreg.annot <- cbind(upreg, gene.annotation0[row.names(upreg),]$AssociatedGeneName)
colnames(upreg.annot)[7] <- "AssociatedGeneName"




# HEatmap of top 50 preferentially expressed genes
row.names(counts.int) <- paste(gene.annotation0[row.names(counts.int),]$AssociatedGeneName, row.names(counts.int), sep =" | ")
ColumnColors <- ifelse(factors_new$SampleType=="Naive", "green", "yellow")
tiff("Heatmap_DEG2.tif", height=12, width=8, units='in', res=300)
heatmap.2(counts.int, trace="none", col=redgreen(75), scale="row", dendrogram="column", keysize=0.4, density.info="none", margins=c(7,11), key.title=NA, key.xlab=NA, key.ylab=NA, ColSideColors=ColumnColors, scale="row")
dev.off()




##### 26 December 2019 ###########################################################################################################################################################################
# Bar graph of gene counts in each module 
library(ggplot2)
tiff("BarPlot.tif", height=10, width=14, units='in', res=300)
axis.title <- element_text(face="bold", size=16)
axis.text.x <- element_text(face="bold", size=13, angle=45,  vjust =1, hjust=1)
axis.text.y <- element_text(face="bold", size=13)
ggplot(box.data, aes(x=row.names(box.data), y=box.data$Freq, fill=row.names(box.data))) + geom_bar(stat="identity") + geom_text(aes(label=box.data$Freq), size=5, vjust=-0.1) + scale_fill_manual(values=row.names(box.data)) + labs(x="Modules", y="Number of genes") + theme(legend.position="none", axis.title=axis.title, axis.text.x=axis.text.x, axis.text.y=axis.text.y, panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()








##### 27 December 2019 ###########################################################################################################################################################################
# Heapmap of module gene epression

geneInfo <- read.csv("/home/dell/Documents/Arindam/WorkB/CombinedStudy/WGCNA/19_SignedHybrid/geneInfo.csv")
genes.int <- geneInfo[geneInfo$moduleColor %in% c("blue", "turquoise"), c(2,8)]
countsPC.int <- countsPC_batch[as.matrix(genes.int$EnsemblGeneID),]
annotation_col <- data.frame(row.names=row.names(factors_new), factors_new$SampleType)
colnames(annotation_col) <- "SampleType"
annotation_row <- data.frame(genes.int, row.names=c(1))
colnames(annotation_row) <- "Module"
annotation_colors <- list(Module = c(turquoise="turquoise", blue="blue"))
tiff("Heatmap_ModuleGenes.tif", height=12, width=8, units='in', res=300)
pheatmap(counts.int.ordered, color=rev(redgreen(75)), scale="row", show_rownames=F, cluster_rows=F, legend=F, annotation_col=annotation_col, annotation_row=annotation_row, annotation_legend=T, annotation_colors=annotation_colors, fontsize=12)
dev.off()





TOM <- TOMsimilarityFromExpr(datExpr, power=11, networkType="signed hybrid")


modules <- c("turquoise")
inModule <- is.finite(match(moduleColors, modules))
modGeneID <- annot.keep[inModule]
modGeneNames <- annot$AssociatedGeneName[match(modGeneID, annot$EnsemblGeneID)]
mod.adj <- adj[inModule, inModule]
dimnames(mod.adj) <- list(modGeneID, modGeneID)

cyt = exportNetworkToCytoscape(mod.adj, 
  edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""), 
  nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""), 
  weighted = TRUE, 
  threshold = 0.313, 
  nodeNames = modGeneID, 
  altNodeNames = modGeneNames, 
  nodeAttr = moduleColors[inModule])





/home/dell/Documents/Arindam/ReferenceGenome/Human_Ensembl84/Annotation/gtf








