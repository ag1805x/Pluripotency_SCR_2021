#######################################################
#
# Script name: WGCNA_v1
# Language: R
# Program: WGCNA
# Purpose: Weighted gene co-expression network analysis
# Author: Arindam Ghosh
# Date: 19 November 2019
#
#######################################################


library(edgeR)
library(DESeq2)
library(WGCNA)


# Load expression data and data annotations
countsPC <- read.table("/home/dell/Documents/Arindam/WorkB/CombinedStudy/Data_Gene_MultimapOverlap/CountsPC_Gene_MultimapOverlap.tsv", header=T, row.names=c(1))
factors <- as.data.frame(read.table("/home/dell/Documents/Arindam/WorkB/CombinedStudy/Data_Gene_MultimapOverlap/factors.txt", sep="\t", header = TRUE, row.names=c(1)))

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Use this to select subset of data on which to perform analysis
#
#
## Remove data from selected project IDS
   list <- c("PRJNA268222", "PRJNA339108", "PRJEB4879", "PRJNA319819", "PRJNA326944")
   ProjectRemove <- vector()
   for(i in 1:length(countsPC)){
     if(as.character(factors$StudyAccession[i]) %in% list){
     ProjectRemove <- c(ProjectRemove, i)}
    }

## Extract counts from projectIDs for which to perform analysis
   countsPC_new <- countsPC[, - ProjectRemove]
   dim(countsPC_new)
   factors_new <- factors[- ProjectRemove,] 
   factors_new <- droplevels(factors_new)
   dim(factors_new)

## Low count genes filter
   keep <- rowSums(cpm(countsPC_new)>1)>1
   countsPC_filt <- countsPC_new[keep,]
   dim(countsPC_filt)
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Low count genes filter
#keep <- rowSums(cpm(countsPC)>1) >1
#countsPC_filt <- countsPC[keep,]
#dim(countsPC_filt) 

# Normalization
dds <- DESeqDataSetFromMatrix(countData=countsPC_filt, colData=factors_new, design = ~ StudyAccession + SampleType)
dds <- DESeq(dds)

vsd <- vst(dds, blind=FALSE)
mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat, vsd$StudyAccession)
assay(vsd) <- mat
countsPC_batch <- assay(vsd)

####################################


# Data input and cleaning
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
datExpr <- as.data.frame(t(countsPC_batch))
as.numeric.factor <- function(x){as.numeric(x)-1} #Converting to numeric starting 0
trait <- data.frame("StudyAccession" = as.numeric.factor(factors_new$StudyAccession), "SampleType" = as.numeric.factor(factors_new$SampleType), "InstrumentModel" = as.numeric.factor(factors_new$InstrumentModel), "LibraryLayout" = as.numeric.factor(factors_new$LibraryLayout)) #Data annotation required as numeric
rownames(trait) <- row.names(factors_new)

# Dendrogram of sample and trait relation
sampleTree <- hclust(dist(datExpr, method="euclidean"))
traitColors <- numbers2colors(trait, signed=FALSE, colors=colors()) #Use colors for different color instead of variant colors
tiff("SampleTraitHeatmap.tif", width=14, height=10, units="in", res=300)
plotDendroAndColors(sampleTree, traitColors, groupLabels=names(trait), main="Sample dendrogram and trait heatmap")
dev.off()

# Automatic construction of the gene network and identification of modules
# Choosing the soft-thresholding power: analysis of network topology
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, blockSize=20000, networkType="signed hybrid")
tiff("SoftThreshold.tif", width=14, height=10, units="in", res=100, bg='white')
#sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
abline(h=0.80,col="blue")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# One-step network construction and module detection
net <- blockwiseModules(datExpr,   maxBlockSize = 20000, power = 14, networkType="signed hybrid", TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "NP_TOM", verbose = 3)
table(net$colors)
#sizeGrWindow(12, 9)
mergedColors <- labels2colors(net$colors)
tiff("ClusterDendrogram.tif", width=14, height=10, units="in", res=300)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

# Relating modules to external information and identifying important genes
# Quantifying module–trait associations
nSamples <- nrow(datExpr)
nGenes <- ncol(datExpr)
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, trait, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
#sizeGrWindow(10,6)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor)

tiff("ModuleTraitAssociation.tif", width=6, height=10, units="in", res=300)
par(mar = c(5, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(trait), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = greenWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.6, zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()

# Export Modute Trait (SampleType) corrlation and p-value to a
Data_ModuleTraitAssociation_SampleType <- data.frame(corr=as.data.frame(moduleTraitCor)$SampleType, p.val=as.data.frame(moduleTraitPvalue)$SampleType, row.names=row.names(as.data.frame(moduleTraitCor)))
write.table(Data_ModuleTraitAssociation_SampleType, "Data_ModuleTraitAssociation_SampleType.tsv", sep="\t", col.names=NA)

#  Gene relationship to trait and important modules: Gene Significance and Module Membership
SampleType <- as.data.frame(trait$SampleType)
names(SampleType) <- "SampleType"
modNames <- substring(names(MEs),3)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) <- paste("MM", modNames, sep="");
names(MMPvalue) <- paste("p.MM", modNames, sep="");
geneTraitSignificance <- as.data.frame(cor(datExpr, SampleType, use = "p"));
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) <- paste("GS.", names(SampleType), sep="");
names(GSPvalue) <- paste("p.GS.", names(SampleType), sep="");

#  Intramodular analysis: identifying genes with high GS and MM
module <- "blue"
column <- match(module, modNames)
moduleGenes <- moduleColors == module
tiff(filename=paste0("ModuleTraitAssociation-", module, ".tif"), width=8, height=8, units="in", res=300)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]), xlab = paste("Module Membership in", module, "module"), ylab = "Gene significance for SampleType", main = paste("Module membership vs. gene significance\n"), cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

# Summary output of network analysis results
#names(datExpr) [moduleColors=="red"] #will return probe IDs belonging to the red module
annot <- read.csv("/home/dell/Documents/Arindam/ReferenceGenome/Human_Ensembl84/Annotation/gtf/GeneAnnotation.csv", header=T)
dim(annot)
annot.keep <- names(datExpr)
genes2annot <- match(annot.keep, annot$EnsemblGeneID)
sum(is.na(genes2annot))
geneInfo0 <- data.frame(EnsemblGeneID = annot.keep, Description=annot$Description[genes2annot], ChromosomeName=annot$ChromosomeName[genes2annot], AssociatedGeneName=annot$AssociatedGeneName[genes2annot], GCcontent=annot$GCcontent[genes2annot], GeneType=annot$GeneType[genes2annot], moduleColor=moduleColors, geneTraitSignificance, GSPvalue)
modOrder <- order(-abs(cor(MEs, SampleType, use="p")))

for (mod in 1:ncol(geneModuleMembership))
{
oldNames <- names(geneInfo0)
geneInfo0 <- data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], MMPvalue[, modOrder[mod]]);
names(geneInfo0) <- c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""), paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.SampleType));
geneInfo <- geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo.csv")


# Interfacing network analysis with other data such as functional annotation and gene ontology

allGeneID <- annot$EnsemblGeneID[genes2annot]
intModules <- c("blue", "brown", "lightgreen", "pink", "yellow", "greenyellow", "black", "turquoise")
for (module in intModules)
{
modGenes <- (moduleColors==module)
modGeneID <- allGeneID[modGenes]
fileName <- paste("EnsemblGeneIDs-", module, ".txt", sep="");
write.table(as.data.frame(modGeneID), file = fileName, row.names = FALSE, col.names = FALSE)
}

# Scale-free plot
tiff("scaleFreePlot.tif", res=100)
k <- softConnectivity(datExpr, power=11, blockSize=20000)

sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k, col="darkgrey")
scaleFreePlot(k)
dev.off()

# IntraModular connectivity
adj <- adjacency(datExpr, type="signed hybrid", power=11)
connectivity <- intramodularConnectivity(adj, moduleColors)
connectivity.annot <- merge(connectivity, annot[, c("EnsemblGeneID", "AssociatedGeneName")], by.x="row.names", by.y="EnsemblGeneID")
colnames(connectivity.annot)[1] <- "EnsemblGeneID"

intModules <- c("blue", "turquoise")
for (module in intModules)
{
modGenes <- (moduleColors==module)
connectivity.modGenes <- connectivity.annot[modGenes, ]
fileName <- paste("Connectivity-", module, ".txt", sep="");
write.table(as.data.frame(connectivity.modGenes), file = fileName, sep="\t", row.names = FALSE)

# Network visualization using WGCNA functions
# Visualizing the gene network
dissTOM <- 1 - TOMsimilarityFromExpr(datExpr, power=14, networkType="signed hybrid")
plotTOM <- dissTOM^7 #How to decide
diag(plotTOM) <- NA
#sizeGrWindow(9,9)
tiff("Network heatmap plot.tif", width=14, height=10, units="in", res=100)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()

# Visualizing the network of eigengenes
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MET <- orderMEs(cbind(MEs, SampleType))
#sizeGrWindow(5,7.5);
tiff("EigengeneNetworks.tif", width=14, height=10, units="in", res=100)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)
dev.off()


# Exporting a gene network to external visualization software
# Exporting to Cytoscape
# TOM <- TOMsimilarityFromExpr(datExpr, power=11, networkType="signed hybrid")
# modules <- c("blue")

#inModule <- is.finite(match(moduleColors, modules))
#modGeneID <- annot.keep[inModule]
#modGeneNames <- annot$AssociatedGeneName[match(modGeneID, annot$EnsemblGeneID)]
#modTOM <- TOM[inModule, inModule]
#dimnames(modTOM) <- list(modGeneID, modGeneID)

modules <- c("turquoise")
inModule <- is.finite(match(moduleColors, modules))
modGeneID <- annot.keep[inModule]
modGeneNames <- annot$AssociatedGeneName[match(modGeneID, annot$EnsemblGeneID)]
mod.adj <- adj[inModule, inModule]
dimnames(mod.adj) <- list(modGeneID, modGeneID)

cyt = exportNetworkToCytoscape(modTOM, 
  edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""), 
  nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""), 
  weighted = TRUE, 
  threshold = 0.313, 
  nodeNames = modGeneID, 
  altNodeNames = modGeneNames, 
  nodeAttr = moduleColors[inModule])
}
