#######################################################
#
# Script name: Sample_hclust_v1
# Language: R
# Program: hclust
# Purpose: Create hierarchial clustering of samples based on expression and plot dendrogram
# Author: Arindam Ghosh
# Date: 19 September 2019
#
#######################################################

library(dendextend)

countsPC <- read.table("CountsPC_Gene_MultimapOverlap.tsv", header=T, row.names=c(1))
factors <- as.data.frame(read.table("factors.txt", sep="\t", header = TRUE, row.names=c(1)))

dist <- dist(t(countsPC), method="euclidean")
cluster <- hclust(dist)


# main, sub, xlab, ylab


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
 StudyAccession_colors_to_use[i] <- switch(as.numeric(StudyAccession_colors_to_use[i]), "#e6194B", "#098f19", "#42d4f4", "#000075", "#f58231", "#911eb4", "#42d4f4", "#f032e6", "#aaffc3")
}


dend <- set(dend, "labels_colors", label_colors_to_use) #Apply colors to labels
dend <- hang.dendrogram(dend, hang=-1) #Leaves at same height
dend <- set(dend, "branches_lwd", 2) #Thick branches
dend <- set(dend, "labels_cex", 1) #label font size


tiff(filename="Sample_hclust.tiff", height=10, width=20, units='in', res=300)

par(mar = c(10,5,3,1))
plot(dend, main = "Sample Clustering\nDist.method=Eucleadian Samples=NP123456789 Normalization=None", ylab = "Height") #Plot the dendrogram
colored_bars(colors = StudyAccession_colors_to_use, dend = dend, rowLabels = "Project", sort_by_labels_order = FALSE)
dev.off()



