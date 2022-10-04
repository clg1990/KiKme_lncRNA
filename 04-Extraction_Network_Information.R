setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dirname(rstudioapi::getActiveDocumentContext()$path)
library(WGCNA)
library(ggplot2)
library(org.Hs.eg.db)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "mrna-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
  lnames = load(file = "mrna-02-networkConstruction-auto_signed_20220927.RData");
lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
moduleTraitPvalue <- matrix(p.adjust(as.vector(as.matrix(moduleTraitPvalue)), method='fdr'),ncol=8)
# names (colors) of the modules
HDIR = as.data.frame(datTraits$HDIR);
names(HDIR)<- "HDIR"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificanceHDIR = as.data.frame(cor(datExpr, HDIR, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificanceHDIR), nSamples));
#names(geneTraitSignificance) = paste("GS.", names(datExpr), sep="");
#names(GSPvalue) = paste("p.GS.", names(HDIR), sep="");
adjacency = adjacency(datExpr, power = 24,
                      type = "signed")
mart <-read.csv("mart_data.csv")
####correct some faulty/missing notation####
mart$hgnc_symbol <- ifelse(mart$ensembl_gene_id == "ENSG00000260912", "AL158206.1",
                                        ifelse(mart$ensembl_gene_id =="ENSG00000277287","AL109976.1",
                                               ifelse(mart$ensembl_gene_id == "ENSG00000273033","LINC02035",
                                                      ifelse(mart$ensembl_gene_id == "ENSG00000257167","TMPO-AS1",
                                                             ifelse(mart$ensembl_gene_id == "ENSG00000259768"," AC004943.2",
                                                                    ifelse(mart$ensembl_gene_id =="ENSG00000261061","AC092718.4",
                                                                           ifelse(mart$ensembl_gene_id =="ENSG00000271936","AC012073.1",
                                                                                  ifelse(mart$ensembl_gene_id =="ENSG00000273084","AC092171.5",
                                                                                         ifelse(mart$ensembl_gene_id =="ENSG00000280206","AC026401.3",
                                                                                                mart$hgnc_symbol)))))))))
####lightgreen####
#Intramodular analysis: identifying genes with high GS and MM
module = "lightgreen"
column = match(module, modNames);
moduleGenes = moduleColors==module;
#identifying genes with high GS and MM
intra_modular_analysis=data.frame(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificanceHDIR[moduleGenes, 1]))
rownames(intra_modular_analysis) = colnames(datExpr)[moduleColors==module] #only the pink module
intra_modular_analysis$bt <- mart$gene_biotype[match(row.names(intra_modular_analysis), mart$ensembl_gene_id)]
intra_modular_analysis$symbol <- mart$hgnc_symbol[match(row.names(intra_modular_analysis), mart$ensembl_gene_id)]

#high intramodular connectivity ~ high kwithin => hub genes (kwithin: connectivity of each driver gene in the module to all other genes in the module)
connectivity=intramodularConnectivity(adjacency, moduleColors)
connectivity = connectivity[colnames(datExpr)[moduleColors==module],] 
order.kWithin = order(connectivity$kWithin, decreasing = TRUE)
connectivity = connectivity[order.kWithin,] #order rows following kWithin
connectivity$bt <- mart$gene_biotype[match(row.names(connectivity), mart$ensembl_gene_id)]
connectivity$symbol <- mart$hgnc_symbol[match(row.names(connectivity), mart$ensembl_gene_id)]
#save connectivity information
write.csv(connectivity,"connectivity_lightgreen.csv")

#merge gene significance and connectivity
intra_modular_analysis$kwithin <- connectivity$kWithin[match(rownames(intra_modular_analysis),rownames(connectivity))]
intra_modular_analysis$kout <- connectivity$kOut[match(rownames(intra_modular_analysis),rownames(connectivity))]
intra_modular_analysis$ktotal <- connectivity$kTotal[match(rownames(intra_modular_analysis),rownames(connectivity))]
intra_modular_analysis$Module <- "Lightgreen Module"
ima_lightgreen <- intra_modular_analysis
#save data for plotting#
save(ima_lightgreen,file = "ima_lightgreen.rds")

####salmon####
#Intramodular analysis: identifying genes with high GS and MM
module = "salmon"
column = match(module, modNames);
moduleGenes = moduleColors==module;
#identifying genes with high GS and MM
intra_modular_analysis=data.frame(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificanceHDIR[moduleGenes, 1]))
rownames(intra_modular_analysis) = colnames(datExpr)[moduleColors==module] #only the pink module
intra_modular_analysis$bt <- mart$gene_biotype[match(row.names(intra_modular_analysis), mart$ensembl_gene_id)]
intra_modular_analysis$symbol <- mart$hgnc_symbol[match(row.names(intra_modular_analysis), mart$ensembl_gene_id)]

#high intramodular connectivity ~ high kwithin => hub genes (kwithin: connectivity of the each driver gene in the module to all other genes in the module)
connectivity=intramodularConnectivity(adjacency, moduleColors)
connectivity = connectivity[colnames(datExpr)[moduleColors==module],] 
order.kWithin = order(connectivity$kWithin, decreasing = TRUE)
connectivity = connectivity[order.kWithin,] #order rows following kWithin
connectivity$bt <- mart$gene_biotype[match(row.names(connectivity), mart$ensembl_gene_id)]
connectivity$symbol <- mart$hgnc_symbol[match(row.names(connectivity), mart$ensembl_gene_id)]
#save connectivity information
write.csv(connectivity,"connectivity_salmon.csv")

#merge gene significance and connectivity
intra_modular_analysis$kwithin <- connectivity$kWithin[match(rownames(intra_modular_analysis),rownames(connectivity))]
intra_modular_analysis$kout <- connectivity$kOut[match(rownames(intra_modular_analysis),rownames(connectivity))]
intra_modular_analysis$ktotal <- connectivity$kTotal[match(rownames(intra_modular_analysis),rownames(connectivity))]
intra_modular_analysis$Module <- "Salmon Module"
ima_salmon <- intra_modular_analysis
#save for plotting#
save(ima_salmon,file = "ima_salmon.rds")

####red####
#Intramodular analysis: identifying genes with high GS and MM
module = "red"
column = match(module, modNames);
moduleGenes = moduleColors==module
#identifying genes with high GS and MM
intra_modular_analysis=data.frame(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificanceHDIR[moduleGenes, 1]))
rownames(intra_modular_analysis) = colnames(datExpr)[moduleColors==module] #only the red module
intra_modular_analysis$bt <- mart$gene_biotype[match(row.names(intra_modular_analysis), mart$ensembl_gene_id)]
intra_modular_analysis$symbol <- mart$hgnc_symbol[match(row.names(intra_modular_analysis), mart$ensembl_gene_id)]

#high intramodular connectivity ~ high kwithin => hub genes (kwithin: connectivity of the each driver gene in the module to all other genes in the module)
connectivity=intramodularConnectivity(adjacency, moduleColors)
connectivity = connectivity[colnames(datExpr)[moduleColors==module],] 
order.kWithin = order(connectivity$kWithin, decreasing = TRUE)
connectivity = connectivity[order.kWithin,] #order rows following kWithin
connectivity$bt <- mart$gene_biotype[match(row.names(connectivity), mart$ensembl_gene_id)]
connectivity$symbol <- mart$hgnc_symbol[match(row.names(connectivity), mart$ensembl_gene_id)]
#save connectivity information
write.csv(connectivity,"connectivity_red.csv")

#merge gene significance and connectivity
intra_modular_analysis$kwithin <- connectivity$kWithin[match(rownames(intra_modular_analysis),rownames(connectivity))]
intra_modular_analysis$kout <- connectivity$kOut[match(rownames(intra_modular_analysis),rownames(connectivity))]
intra_modular_analysis$ktotal <- connectivity$kTotal[match(rownames(intra_modular_analysis),rownames(connectivity))]

intra_modular_analysis$Module <- "Red Module"
ima_red <- intra_modular_analysis
save(ima_red,file = "ima_red.rds")

####Loop through data to get adjacency data for all lncRNAs in the radiation-response modules####
list <- rownames(intra_modular_analysis)[intra_modular_analysis$bt == "lncRNA"]
list
#names(list) <- mapIds(org.Hs.eg.db, list, 'SYMBOL', 'ENSEMBL')
#names(list)
for (i in 1:length(list)){
  filter <- list[i]
  adjacency = data.frame(adjacency(datExpr, selectCols =filter,
                                                     power = 24,
                                                     type = "signed"))  
  adjacency$bt <- mart$gene_biotype[match(row.names(adjacency), mart$ensembl_gene_id)]
  adjacency$symbol <- mart$hgnc_symbol[match(row.names(adjacency), mart$ensembl_gene_id)]
  colnames(adjacency)[1] <- "Adjacency"
  adjacency <- adjacency[order(-adjacency$Adjacency),]
  adjacency <- adjacency[,c(3,2,1)]
  write.csv(adjacency,file = paste("adjacency_",adjacency$symbol[adjacency$Adjacency ==1],".csv",sep=""))
  }