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
options(scipen = 999)
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
options(stringsAsFactors = FALSE)

traitData = read.csv("meta_combined_allTimes.csv",row.names = 1)
ISI <- "ISI_400"
traitData <- subset(traitData,traitData$Time =="4h")#only 4h data
traitData <- subset(traitData,!(traitData$ID %in% ISI))#remove faulty probe (1783B	2Gy	4h, resequenced probes on all doses too different, not used)
rm(ISI)
mrna <- read.csv("MRNA_log_counts_DESEq_normalized_normalized_ENSEMBL_proteincoding_only.csv",row.names = 1)
rownames(mrna) <- mrna$ENSEMBL
mrna <- mrna[colnames(mrna) %in% traitData$ID]
traitData <- traitData[traitData$ID %in% colnames(mrna),]
#write.csv(traitData,"Meta_50T.csv")
lncrna <- read.csv("ISI_ID_lncRNA_4h_log_deseq2.csv",row.names = 1)
lncrna <- lncrna[colnames(lncrna) %in% traitData$ID]
sum(is.na(lncrna))
lncrna <- lncrna[,colnames(mrna)]
mrna <- rbind(mrna,lncrna)
test <- mrna
test$gene <- row.names(test)
test <- melt(test, id.vars= c("gene"))

mrna <- data.frame(t(mrna))
mrna_gsg = goodSamplesGenes(mrna, verbose = 3)#check "verbose"
mrna_gsg$allOK
rm(mrna_gsg)
#Next sample clustering to check for obvious outliers
sampleTree = hclust(dist(mrna), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 50, height = 20);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

dev.off()
# Determine cluster under the line
?hclust
clust = cutreeStatic(sampleTree, cutHeight = 150, minSize = 10)

#sampleTree$labels[sampleTree$height >= 110]
#outlier <- c("ISI_522","ISI_521","ISI_523","ISI_428","ISI_429","ISI_430","ISI_464","ISI_465","ISI_466")
#outlier_meta <- subset(traitData, traitData$ID %in% outlier)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = mrna[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
traitData$LDIR <- ifelse(traitData$Dose == "0.05 Gy",1,0)
traitData$HDIR <- ifelse(traitData$Dose == "2 Gy",1,0)
traitData$SI <- ifelse(traitData$Dose == "0 Gy",1,0)
traitData$CO <- ifelse(traitData$Group == "CO",1,0)
traitData$FPN <- ifelse(traitData$Group == "FPN",1,0)
traitData$SPN <- ifelse(traitData$Group == "SPN",1,0)

allTraits = traitData[, c(4,8,16,18:23)]
allTraits$Age_Sampling <-as.numeric(allTraits$Age_Sampling)
#allTraits$RekID <- as.numeric(as.factor(allTraits$RekID))
allTraits$Sex <- as.numeric(as.factor(allTraits$Sex))
colnames(allTraits)[1] <- "ID"
colnames(allTraits)[3] <- "Age"
allTraits <- allTraits[,c(1,2,3,6,4,5,7,8,9)]
dim(allTraits)
names(allTraits)
# Form a data frame analogous to expression data that will hold the clinical traits.
Samples = rownames(datExpr)
traitRows = match(Samples, allTraits$ID)
datTraits = allTraits[traitRows, -1]
rownames(datTraits) = allTraits[traitRows, 1]
collectGarbage()

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)

# Plot the sample dendrogram and the colors underneath.
#png("Dendrogram_Heatmap.png",height=12000,width=30000,res=600)  
pdf(file = "sampleClustering_traits.pdf", width = 100, height = 50);
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()
save(datExpr, datTraits, file = "mrna-01-dataInput.RData")

write.csv(rownames(df),"background_lncrna_mrna.csv",row.names = F)
