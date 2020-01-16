#######################################################
#
# Script name: UpSetR_Plot_v1
# Language: R
# Program: UpSetR
# Purpose: Check intersection between groups of data
# Author: Arindam Ghosh
# Date: 18 September 2019
#
#######################################################


library(UpSetR)
NP1_upreg <- read.table("/home/dell/Documents/Arindam/WorkB/PRJNA286204/DESeq2_Gene_MultimapOverlap/PRJNA286204_DESeq2_upreg.tsv", header=T, row.names=c(1))
NP2_upreg <- read.table("/home/dell/Documents/Arindam/WorkB/PRJNA360413/DESeq2_Gene_MultimapOverlap/PRJNA360413_DESeq2_upreg.tsv", header=T, row.names=c(1))
NP3_upreg <- read.table("/home/dell/Documents/Arindam/WorkB/PRJNA268222/DESeq2_Gene_MultimapOverlap/PRJNA268222_DESeq2_upreg.tsv", header=T, row.names=c(1))
NP4_upreg <- read.table("/home/dell/Documents/Arindam/WorkB/PRJNA339108/DESeq2_Gene_MultimapOverlap/PRJNA339108_DESeq2_upreg.tsv", header=T, row.names=c(1))
NP5_upreg <- read.table("/home/dell/Documents/Arindam/WorkB/PRJEB7132/DESeq2_Gene_MultimapOverlap/PRJEB7132_DESeq2_upreg.tsv", header=T, row.names=c(1))
NP6_upreg <- read.table("/home/dell/Documents/Arindam/WorkB/PRJEB20388/DESeq2_Gene_MultimapOverlap/PRJEB20388_DESeq2_upreg.tsv", header=T, row.names=c(1))
NP7_upreg <- read.table("/home/dell/Documents/Arindam/WorkB/PRJEB4879/DESeq2_Gene_MultimapOverlap/PRJEB4879_DESeq2_upreg.tsv", header=T, row.names=c(1))

upreg.list <- list(NP1_upreg = rownames(NP1_upreg), NP2_upreg = rownames(NP2_upreg), NP3_upreg = rownames(NP3_upreg), NP4_upreg = rownames(NP4_upreg), NP5_upreg = rownames(NP5_upreg), NP6_upreg = rownames(NP6_upreg), NP7_upreg = rownames(NP7_upreg))
names(upreg.list)
 
length(upreg.list$NP1_upreg)
length(upreg.list$NP2_upreg)
length(upreg.list$NP3_upreg)
length(upreg.list$NP4_upreg)
length(upreg.list$NP5_upreg)
length(upreg.list$NP6_upreg)
length(upreg.list$NP7_upreg)

upreg.genes <- fromList(upreg.list)

tiff(filename="UpSetPlot_Deg_Upreg.tiff", height=10, width=50, units='in', res=300)
upset(upreg.genes, order.by ="degree", nset=7, nintersects=NA, text.scale=c(1.75, 1.5, 1.5, 1.5, 2, 1.75), empty.intersections = "off", mainbar.y.label = "Number of genes", sets.x.label = "Number of genes per set", shade.alpha=1, queries = list(list(query = intersects, params=list("NP1_upreg", "NP2_upreg", "NP3_upreg", "NP4_upreg", "NP5_upreg", "NP6_upreg", "NP7_upreg"), color="red", active=T)))
dev.off()

##########################################################

NP1_downreg <- read.table("/home/dell/Documents/Arindam/WorkB/PRJNA286204/DESeq2_Gene_MultimapOverlap/PRJNA286204_DESeq2_downreg.tsv", header=T, row.names=c(1))
NP2_downreg <- read.table("/home/dell/Documents/Arindam/WorkB/PRJNA360413/DESeq2_Gene_MultimapOverlap/PRJNA360413_DESeq2_downreg.tsv", header=T, row.names=c(1))
NP3_downreg <- read.table("/home/dell/Documents/Arindam/WorkB/PRJNA268222/DESeq2_Gene_MultimapOverlap/PRJNA268222_DESeq2_downreg.tsv", header=T, row.names=c(1))
NP4_downreg <- read.table("/home/dell/Documents/Arindam/WorkB/PRJNA339108/DESeq2_Gene_MultimapOverlap/PRJNA339108_DESeq2_downreg.tsv", header=T, row.names=c(1))
NP5_downreg <- read.table("/home/dell/Documents/Arindam/WorkB/PRJEB7132/DESeq2_Gene_MultimapOverlap/PRJEB7132_DESeq2_downreg.tsv", header=T, row.names=c(1))
NP6_downreg <- read.table("/home/dell/Documents/Arindam/WorkB/PRJEB20388/DESeq2_Gene_MultimapOverlap/PRJEB20388_DESeq2_downreg.tsv", header=T, row.names=c(1))
NP7_downreg <- read.table("/home/dell/Documents/Arindam/WorkB/PRJEB4879/DESeq2_Gene_MultimapOverlap/PRJEB4879_DESeq2_downreg.tsv", header=T, row.names=c(1))

downreg.list <- list(NP1_downreg = rownames(NP1_downreg), NP2_downreg = rownames(NP2_downreg), NP3_downreg = rownames(NP3_downreg), NP4_downreg = rownames(NP4_downreg), NP5_downreg = rownames(NP5_downreg), NP6_downreg = rownames(NP6_downreg), NP7_downreg = rownames(NP7_downreg))
names(downreg.list)
 
length(downreg.list$NP1_downreg)
length(downreg.list$NP2_downreg)
length(downreg.list$NP3_downreg)
length(downreg.list$NP4_downreg)
length(downreg.list$NP5_downreg)
length(downreg.list$NP6_downreg)
length(downreg.list$NP7_downreg)

downreg.genes <- fromList(downreg.list)

tiff(filename="UpSetPlot_Deg_Downreg.tiff", height=10, width=50, units='in', res=300)
upset(downreg.genes, order.by ="degree", nset=7, nintersects=NA, text.scale=c(1.75, 1.5, 1.5, 1.5, 2, 1.75), empty.intersections = "off", mainbar.y.label = "Number of genes", sets.x.label = "Number of genes per set", shade.alpha=1, queries = list(list(query = intersects, params=list("NP1_downreg", "NP2_downreg", "NP3_downreg", "NP4_downreg", "NP5_downreg", "NP6_downreg", "NP7_downreg"), color="red", active=T)))
dev.off()

##############################################################

list <- list(NP1_upreg = rownames(NP1_upreg), NP2_upreg = rownames(NP2_upreg), NP3_upreg = rownames(NP3_upreg), NP4_upreg = rownames(NP4_upreg), NP5_upreg = rownames(NP5_upreg), NP6_upreg = rownames(NP6_upreg), NP7_upreg = rownames(NP7_upreg), NP1_downreg = rownames(NP1_downreg), NP2_downreg = rownames(NP2_downreg), NP3_downreg = rownames(NP3_downreg), NP4_downreg = rownames(NP4_downreg), NP5_downreg = rownames(NP5_downreg), NP6_downreg = rownames(NP6_downreg), NP7_downreg = rownames(NP7_downreg))

all.genes <- fromList(list)
metadata <- read.table("/home/dell/Desktop/metedata.txt", header=T)


tiff(filename="UpSetPlot_Deg_all.tiff", height=10, width=250, units='in', res=72)

upset(all.genes, order.by ="degree", nset=14, nintersects=NA, sets = c("NP1_upreg", "NP2_upreg", "NP3_upreg", "NP4_upreg", "NP5_upreg", "NP6_upreg", "NP7_upreg", "NP1_downreg", "NP2_downreg", "NP3_downreg", "NP4_downreg", "NP5_downreg", "NP6_downreg", "NP7_downreg"), text.scale=c(1.75, 1.5, 1.5, 1.5, 2, 1.75), mainbar.y.label = "Number of genes", sets.x.label = "Number of genes per set", keep.order=T, set.metadata=list(data=metadata, plots=list(list(type='matrix_rows', column='Dataset', assign=2, alpha=1, colors=c(NP1="#ff0101", NP2="#ffd8b1", NP3="#ffe119", NP4="#aaffc3", NP5="#42d4f4", NP6="#e6beff", NP7="#a9a9a9")))))

dev.off()
































