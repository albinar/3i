# Written by Albina Rahim
# Date: March 21, 2016

remove(list=ls())

setwd("/home/rstudio/code/Projects/3i/Panel_M-cell/")


library("png")

######################################################################

outliersJM <- function(x, SDs = 3){
  mean <- mean(x)
  sd <- sd(x)
  indices <- union(which(x >= (mean + SDs*sd)),which(x <= (mean - SDs*sd)))
  outliers <- x[indices]
  return(list(out=outliers, ind=indices))
}

######################################################################

#results.dir <-  "/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/SPLEEN/Results"
#results.dir <- "/home/rstudio/results/IMPC/3i/Panel_T-cell/SPLEEN/Results-Original"
# load(paste0(results.dir,"/store.allFCS.Rdata"))
# load(paste0(results.dir,"/ind.marg.neg.clean.all.Rdata"))
# load(paste0(results.dir, "/res.clean.Rdata"))
# load(paste0(results.dir, "/all.gthres.Rdata"))
# load(paste0(results.dir, "/failedGating.files.Rdata"))
# allProportions_Table <- read.csv(paste0(results.dir, "/allProportions_Table.csv"))
# flaggedFiles <- read.csv(paste0(results.dir,"/flagged.FCS.TcellSPLEEN.Feedback.csv"))
# ## Reading the TvsF flagged csv file which was sent to Adam for feedback
# flagged.FCS <- as.matrix(read.csv(paste0(results.dir, "/flagged.FCS.TcellSPLEEN.Feedback.csv"), sep = ","))
# ## List the indices of the flagged FCS files based on TvsF algorithm
# flagged.FCS.index <- which(res.clean[c('Has the file passed'),] == "F")
# 
# rownames(flagged.FCS) <- flagged.FCS.index
# 
# ## Manually including some flagged files for further analysis
# manualExclude.flagged.FCS <- flagged.FCS[-which(flagged.FCS[,c('Final.Action')] == "Keep"),]
# manualExclude.flagged.FCS.index <- as.integer(rownames(manualExclude.flagged.FCS))
# 
# ## Removing flagged files based on Adam's feedback from further analysis after recieving confirmation from Adam (Jan 04, 2017)
# if ( length(manualExclude.flagged.FCS.index) != 0){
#   store.allFCS <- store.allFCS[-manualExclude.flagged.FCS.index,]
#   res.clean <- res.clean[,-manualExclude.flagged.FCS.index]
#   ind.marg.neg.clean.all[manualExclude.flagged.FCS.index] <- NULL
# }
# ## There are files which failed the Gating in our previous step
# ## So we exclude the thresholds of these files
# if (nrow(failedGating.files) != 0){
#   failedGating.files.index <- as.integer(rownames(failedGating.files))
#   store.allFCS <- store.allFCS[-failedGating.files.index,]
#   res.clean <- res.clean[,-failedGating.files.index]
#   ind.marg.neg.clean.all[failedGating.files.index] <- NULL
#   all.gthres <- all.gthres[-failedGating.files.index,]
#   allProportions_Table_Initial <- allProportions_Table_Initial[-failedGating.files.index,]
# }
# 
# 
# ## I am manually removing "SPLN_L000061297_A01_001.labelled.fcs" file based on my observation during the first phase of the run
# ## and also based on Ania's feedback (January, 2017), confirming this WT file to be an outlier.
# manual.remove.Ind <- which(store.allFCS[,c('Barcodes')] == "L000061297")
# store.allFCS <- store.allFCS[-manual.remove.Ind,]
# res.clean <- res.clean[,-manual.remove.Ind]
# ind.marg.neg.clean.all[manual.remove.Ind] <- NULL
# all.gthres <- all.gthres[-manual.remove.Ind,]
# allProportions_Table_Initial <-  allProportions_Table_Initial[-manual.remove.Ind,]
# 
# rownames(allProportions_Table_Initial) <- 1:nrow(allProportions_Table_Initial)
# rownames(store.allFCS) <- 1:nrow(store.allFCS)


# boxplot(all.gthres, pch=19, cex=0.2)
results.dir <- "/home/rstudio/results/IMPC/3i/Panel_M-cell/SPLEEN/Outliers"
data.dir <- "/home/rstudio/data/IMPC/3i/Panel_M-cell/SPLEEN/Results"
allProportions_Table.temp <- read.csv(paste0(data.dir,"/Events_Proportions_Table_20170802_1113.csv"))
allProportions_Table <- allProportions_Table.temp[,c(1:29)]
populationsNames <- colnames(allProportions_Table)


#############################################################################################################

## Finding outliers for 'Lymphocytes'
boxplot(allProportions_Table[,13], pch=19, cex=0.2, main = "Lymphocytes", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,13], SDs = 4.5)$ind, outliersJM(allProportions_Table[,13], SDs = 4.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/home/rstudio/data/IMPC/3i/Panel_M-cell/Panel_M-cell_scatterplots/SPLEEN/Results/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/Total_", allProportions_Table[ind.with.store[k],"FCS.files"], ".png")
  file.copy(name, paste0(results.dir,"/Lymphocytes/"))
}

outliers.Lymphocytes <- as.matrix(allProportions_Table[ind.with.store,c(1:8,13)])
save(outliers.Lymphocytes, file = paste0(results.dir,"/Lymphocytes/outliers.Lymphocytes.Rdata"))


###########################################################################################################

## Finding outliers for 'CD45+'
boxplot(allProportions_Table[,14], pch=19, cex=0.2, main = "CD45+", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,14], SDs = 5)$ind, outliersJM(allProportions_Table[,14], SDs = 5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/home/rstudio/data/IMPC/3i/Panel_M-cell/Panel_M-cell_scatterplots/SPLEEN/Results/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/Total_", allProportions_Table[ind.with.store[k],"FCS.files"], ".png")
  file.copy(name, paste0(results.dir,"/CD45/"))
  
  # x11(width = 16, height = 9)
  # img <- readPNG(name)
  # grid::grid.raster(img)
}

outliers.CD45 <- as.matrix(allProportions_Table[ind.with.store,c(1:8,14)])
save(outliers.CD45, file = paste0(results.dir,"/CD45/outliers.CD45.Rdata"))

###########################################################################################################

## Finding outliers for Lineage neg
boxplot(allProportions_Table[,15], pch=19, cex=0.2, main = "Lineage neg", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,15], SDs = 5.5)$ind, outliersJM(allProportions_Table[,15], SDs = 5.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/home/rstudio/data/IMPC/3i/Panel_M-cell/Panel_M-cell_scatterplots/SPLEEN/Results/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/Total_", allProportions_Table[ind.with.store[k],"FCS.files"], ".png")
  file.copy(name, paste0(results.dir,"/LineageNeg/"))
}

outliers.LineageNeg <- as.matrix(allProportions_Table[ind.with.store,c(1:8,15)])
save(outliers.LineageNeg, file = paste0(results.dir,"/LineageNeg/outliers.LineageNeg.Rdata"))

###########################################################################################################


## Finding outliers for LinMac
boxplot(allProportions_Table[,16], pch=19, cex=0.2, main = "Lin- Mac-", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,16], SDs = 3.75)$ind, outliersJM(allProportions_Table[,16], SDs = 3.75)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/home/rstudio/data/IMPC/3i/Panel_M-cell/Panel_M-cell_scatterplots/SPLEEN/Results/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/Total_", allProportions_Table[ind.with.store[k],"FCS.files"], ".png")
  file.copy(name, paste0(results.dir,"/LinMac/"))
}

outliers.LinMac <- as.matrix(allProportions_Table[ind.with.store,c(1:8,16)])
save(outliers.LinMac, file = paste0(results.dir,"/LinMac/outliers.LinMac.Rdata"))


#############################################################################################################

## Finding outliers for RP Machrophages
boxplot(allProportions_Table[,17], pch=19, cex=0.2, main = "RP Machrophages", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,17], SDs = 3.75)$ind, outliersJM(allProportions_Table[,17], SDs = 3.75)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/home/rstudio/data/IMPC/3i/Panel_M-cell/Panel_M-cell_scatterplots/SPLEEN/Results/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/Total_", allProportions_Table[ind.with.store[k],"FCS.files"], ".png")
  file.copy(name, paste0(results.dir,"/RPmacrophages/"))
}

outliers.RPmac <- as.matrix(allProportions_Table[ind.with.store,c(1:8,17)])
save(outliers.RPmac, file = paste0(results.dir,"/RPmacrophages/outliers.RPmac.Rdata"))

#############################################################################################################

## Finding outliers for Ly6G-
boxplot(allProportions_Table[,18], pch=19, cex=0.2, main = "Ly6G-", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,18], SDs = 3.75)$ind, outliersJM(allProportions_Table[,18], SDs = 3.75)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/home/rstudio/data/IMPC/3i/Panel_M-cell/Panel_M-cell_scatterplots/SPLEEN/Results/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/Total_", allProportions_Table[ind.with.store[k],"FCS.files"], ".png")
  file.copy(name, paste0(results.dir,"/Ly6G/"))
}

outliers.Ly6G <- as.matrix(allProportions_Table[ind.with.store,c(1:8,18)])
save(outliers.Ly6G, file = paste0(results.dir,"/Ly6G/outliers.Ly6G.Rdata"))

#############################################################################################################

## Finding outliers for Neutrophils
boxplot(allProportions_Table[,19], pch=19, cex=0.2, main = "Neutrophils", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,19], SDs = 4.5)$ind, outliersJM(allProportions_Table[,19], SDs = 4.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/home/rstudio/data/IMPC/3i/Panel_M-cell/Panel_M-cell_scatterplots/SPLEEN/Results/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/Total_", allProportions_Table[ind.with.store[k],"FCS.files"], ".png")
  file.copy(name, paste0(results.dir,"/Neutrophils/"))
}

outliers.Neutrophils <- as.matrix(allProportions_Table[ind.with.store,c(1:8,19)])
save(outliers.Neutrophils, file = paste0(results.dir,"/Neutrophils/outliers.Neutrophils.Rdata"))

#############################################################################################################

## Finding outliers for Eosinophils
boxplot(allProportions_Table[,20], pch=19, cex=0.2, main = "Eosinophils", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,20], SDs = 3.75)$ind, outliersJM(allProportions_Table[,20], SDs = 3.75)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/home/rstudio/data/IMPC/3i/Panel_M-cell/Panel_M-cell_scatterplots/SPLEEN/Results/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/Total_", allProportions_Table[ind.with.store[k],"FCS.files"], ".png")
  file.copy(name, paste0(results.dir,"/Eosinophils/"))
}

outliers.Eosinophils <- as.matrix(allProportions_Table[ind.with.store,c(1:8,20)])
save(outliers.Eosinophils, file = paste0(results.dir,"/Eosinophils/outliers.Eosinophils.Rdata"))

#############################################################################################################

## Finding outliers for Monocytes Ly6c hi and NOT(Monocytes Ly6c hi)
boxplot(allProportions_Table[,21], pch=19, cex=0.2, main = " Monocytes Ly6c hi", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,21], SDs = 4)$ind, outliersJM(allProportions_Table[,21], SDs = 4)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/home/rstudio/data/IMPC/3i/Panel_M-cell/Panel_M-cell_scatterplots/SPLEEN/Results/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/Total_", allProportions_Table[ind.with.store[k],"FCS.files"], ".png")
  file.copy(name, paste0(results.dir,"/MonocytesLy6c/"))
}

outliers.MonocytesLy6c <- as.matrix(allProportions_Table[ind.with.store,c(1:8,21,22)])
save(outliers.MonocytesLy6c, file = paste0(results.dir,"/MonocytesLy6c/outliers.MonocytesLy6c.Rdata"))

#############################################################################################################

## Finding outliers for pDC and NOT(pDC)
boxplot(allProportions_Table[,23], pch=19, cex=0.2, main = "pDC", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,23], SDs = 5)$ind, outliersJM(allProportions_Table[,23], SDs = 5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/home/rstudio/data/IMPC/3i/Panel_M-cell/Panel_M-cell_scatterplots/SPLEEN/Results/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/Total_", allProportions_Table[ind.with.store[k],"FCS.files"], ".png")
  file.copy(name, paste0(results.dir,"/pDC/"))
}

outliers.pDC <- as.matrix(allProportions_Table[ind.with.store,c(1:8,23,24)])
save(outliers.pDC, file = paste0(results.dir,"/pDC/outliers.pDC.Rdata"))

#############################################################################################################

## Finding outliers for cDC and Misc
boxplot(allProportions_Table[,25], pch=19, cex=0.2, main = "cDC", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,25], SDs = 4)$ind, outliersJM(allProportions_Table[,25], SDs = 4)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/home/rstudio/data/IMPC/3i/Panel_M-cell/Panel_M-cell_scatterplots/SPLEEN/Results/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/Total_", allProportions_Table[ind.with.store[k],"FCS.files"], ".png")
  file.copy(name, paste0(results.dir,"/cDC/"))
}

outliers.cDC <- as.matrix(allProportions_Table[ind.with.store,c(1:8,25,26)])
save(outliers.cDC, file = paste0(results.dir,"/cDC/outliers.cDC.Rdata"))

#############################################################################################################

## Finding outliers for CD8ATypeDC and CD11B+CD86Lo
boxplot(allProportions_Table[,27], pch=19, cex=0.2, main = "CD8A Type DC", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,27], SDs = 4)$ind, outliersJM(allProportions_Table[,27], SDs = 4)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/home/rstudio/data/IMPC/3i/Panel_M-cell/Panel_M-cell_scatterplots/SPLEEN/Results/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/Total_", allProportions_Table[ind.with.store[k],"FCS.files"], ".png")
  file.copy(name, paste0(results.dir,"/CD8ATypeDC/"))
}

outliers.CD8ATypeDC <- as.matrix(allProportions_Table[ind.with.store,c(1:8,27,28)])
save(outliers.CD8ATypeDC, file = paste0(results.dir,"/CD8ATypeDC/outliers.CD8ATypeDC.Rdata"))

#############################################################################################################

## Finding outliers for CD103+ DC
boxplot(allProportions_Table[,29], pch=19, cex=0.2, main = "CD103+ DC", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,29], SDs = 4)$ind, outliersJM(allProportions_Table[,29], SDs = 4)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/home/rstudio/data/IMPC/3i/Panel_M-cell/Panel_M-cell_scatterplots/SPLEEN/Results/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/Total_", allProportions_Table[ind.with.store[k],"FCS.files"], ".png")
  file.copy(name, paste0(results.dir,"/CD103posDC/"))
}

outliers.CD103posDC <- as.matrix(allProportions_Table[ind.with.store,c(1:8,29)])
save(outliers.CD103posDC, file = paste0(results.dir,"/CD103posDC/outliers.CD103posDC.Rdata"))
