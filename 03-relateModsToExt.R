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

allowWGCNAThreads()
enableWGCNAThreads()
lnames = load(file = "mrna-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
  lnames = load(file = "mrna-02-networkConstruction-auto_signed_20220927.RData");

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
#moduleTraitCor = cor(MEs, datTraits[,c(1,2,18:35)], use = "p")
moduleTraitCor = cor(MEs, datTraits, use = "p")#student asymptotic p-value for given correlations.
#?cor
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#??corPvalueStudent
colnames <- colnames(moduleTraitPvalue)
rownames <- rownames(moduleTraitPvalue)
moduleTraitPvalue <- matrix(p.adjust(as.vector(as.matrix(moduleTraitPvalue)), method='fdr'),ncol=8)
colnames(moduleTraitPvalue) <- colnames
rownames(moduleTraitPvalue) <- rownames
moduleTraitPvalue <- -log10(moduleTraitPvalue)


####Submit to plotting####
