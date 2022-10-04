library(ggh4x)
library(reshape2)
library(ggplot2)
library(data.table)
library(tidyr)
library(rlist)
library(stringr)
library(readxl)
library(forcats)
library(tidyverse)
library(WGCNA)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dirname(rstudioapi::getActiveDocumentContext()$path)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
lnames = load(file = "mrna-01-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
allowWGCNAThreads()
enableWGCNAThreads()
#Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 10, to=30, by=1))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers,networkType = "signed", 
                        corFnc="bicor",
                        corOptions=list(maxPOutliers=0.1),
                        verbose = 5)
power<-  sft$fitIndices %>% filter(SFT.R.sq >= 0.9) %>%  
        filter(row_number() == 1)
power <- power$Power
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

collectGarbage()

#?blockwiseModules
net <- blockwiseModules(datExpr, power = power, corType = "bicor",
                       TOMType = "signed", networkType = "signed",
                       minModuleSize = 15,
                       reassignThreshold = 0, mergeCutHeight = 0.20,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM_signed_2022",
                       verbose = 3,maxBlockSize = 20000)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;

geneTree = net$dendrograms[[1]]

save(MEs, moduleLabels, moduleColors, geneTree,
     file = "mrna-02-networkConstruction-auto_signed_20220927.RData")
