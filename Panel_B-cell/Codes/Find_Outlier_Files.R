# Written by Albina Rahim
# Date: March 21, 2016

remove(list=ls())

setwd("/home/rstudio/code/Projects/3i/Panel_B-cell/")


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
results.dir <- "/home/rstudio/results/IMPC/3i/Panel_B-cell/SPLEEN/Outliers"
data.dir <- "/home/rstudio/data/IMPC/3i/Panel_B-cell/SPLEEN/Results"
allProportions_Table <- read.csv(paste0(data.dir,"/allProportions_Table.csv"))

populationsNames <- colnames(allProportions_Table)


#############################################################################################################

## Finding outliers for 'Lymphocytes'
boxplot(allProportions_Table[,11], pch=19, cex=0.2, main = "Lymphocytes", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,11], SDs = 4.5)$ind, outliersJM(allProportions_Table[,11], SDs = 4.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0(data.dir,"/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/ScatterPlot_", allProportions_Table[ind.with.store[k],"Barcodes"], ".png")
  file.copy(name, paste0(results.dir,"/Lymphocytes/"))
}

outliers.Lymphocytes <- as.matrix(allProportions_Table[ind.with.store,c(1:8,11)])
save(outliers.Lymphocytes, file = paste0(results.dir,"/Lymphocytes/outliers.Lymphocytes.Rdata"))


###########################################################################################################

## Finding outliers for 'CD45+'
boxplot(allProportions_Table[,14], pch=19, cex=0.2, main = "CD45+", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,14], SDs = 5)$ind, outliersJM(allProportions_Table[,14], SDs = 5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0(data.dir,"/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/ScatterPlot_", allProportions_Table[ind.with.store[k],"Barcodes"], ".png")
  file.copy(name, paste0(results.dir,"/CD45+/"))
  
  # x11(width = 16, height = 9)
  # img <- readPNG(name)
  # grid::grid.raster(img)
}

outliers.CD45 <- as.matrix(allProportions_Table[ind.with.store,c(1:8,14)])
save(outliers.CD45, file = paste0(results.dir,"/CD45+/outliers.CD45.Rdata"))

###########################################################################################################

## Finding outliers for B-cells
boxplot(allProportions_Table[,15], pch=19, cex=0.2, main = "B-cells", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,15], SDs = 4.5)$ind, outliersJM(allProportions_Table[,15], SDs = 4.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0(data.dir,"/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/ScatterPlot_", allProportions_Table[ind.with.store[k],"Barcodes"], ".png")
  file.copy(name, paste0(results.dir,"/Bcells/"))
}

outliers.Bcells <- as.matrix(allProportions_Table[ind.with.store,c(1:8,15)])
save(outliers.Bcells, file = paste0(results.dir,"/Bcells/outliers.Bcells.Rdata"))

###########################################################################################################


## Finding outliers for Plasma
boxplot(allProportions_Table[,16], pch=19, cex=0.2, main = "Plasma", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,16], SDs = 4)$ind, outliersJM(allProportions_Table[,16], SDs = 4)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0(data.dir,"/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/ScatterPlot_", allProportions_Table[ind.with.store[k],"Barcodes"], ".png")
  file.copy(name, paste0(results.dir,"/Plasma/"))
}

outliers.Plasma <- as.matrix(allProportions_Table[ind.with.store,c(1:8,16)])
save(outliers.Plasma, file = paste0(results.dir,"/Plasma/outliers.Plasma.Rdata"))


#############################################################################################################

## Finding outliers for B1 and B2 B-cells
boxplot(allProportions_Table[,17], pch=19, cex=0.2, main = "B1 B-cells", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,17], SDs = 4.5)$ind, outliersJM(allProportions_Table[,17], SDs = 4.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0(data.dir,"/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/ScatterPlot_", allProportions_Table[ind.with.store[k],"Barcodes"], ".png")
  file.copy(name, paste0(results.dir,"/B1-Bcells/"))
}

outliers.B1Bcells <- as.matrix(allProportions_Table[ind.with.store,c(1:8,17,18)])
save(outliers.B1Bcells, file = paste0(results.dir,"/B1-Bcells/outliers.B1Bcells.Rdata"))

#############################################################################################################

## Finding outliers for GC and NOT(GC)
boxplot(allProportions_Table[,19], pch=19, cex=0.2, main = "GC", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,19], SDs = 4.5)$ind, outliersJM(allProportions_Table[,19], SDs = 4.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0(data.dir,"/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/ScatterPlot_", allProportions_Table[ind.with.store[k],"Barcodes"], ".png")
  file.copy(name, paste0(results.dir,"/GC/"))
}

outliers.GC <- as.matrix(allProportions_Table[ind.with.store,c(1:8,19,20)])
save(outliers.GC, file = paste0(results.dir,"/GC/outliers.GC.Rdata"))

#############################################################################################################


## Finding outliers for Late GC and Early GC
boxplot(allProportions_Table[,21], pch=19, cex=0.2, main = "Late GC", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,21], SDs = 4.5)$ind, outliersJM(allProportions_Table[,21], SDs = 4.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0(data.dir,"/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/ScatterPlot_", allProportions_Table[ind.with.store[k],"Barcodes"], ".png")
  file.copy(name, paste0(results.dir,"/lateGC/"))
}

outliers.lateGC <- as.matrix(allProportions_Table[ind.with.store,c(1:8,21,22)])
save(outliers.lateGC, file = paste0(results.dir,"/lateGC/outliers.lateGC.Rdata"))

#############################################################################################################

## Finding outliers for MZ+MZP B-cells
boxplot(allProportions_Table[,23], pch=19, cex=0.2, main = "MZ+MZP B-cells", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,23], SDs = 3)$ind, outliersJM(allProportions_Table[,23], SDs = 3)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0(data.dir,"/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/ScatterPlot_", allProportions_Table[ind.with.store[k],"Barcodes"], ".png")
  file.copy(name, paste0(results.dir,"/MZ+MZP-Bcells/"))
}

outliers.MZposMZP.Bcells <- as.matrix(allProportions_Table[ind.with.store,c(1:8,23)])
save(outliers.MZposMZP.Bcells, file = paste0(results.dir,"/MZ+MZP-Bcells/outliers.MZposMZP.Bcells.Rdata"))

#############################################################################################################

## Finding outliers for Transitional B-cells
boxplot(allProportions_Table[,24], pch=19, cex=0.2, main = "Transitional B-cells", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,24], SDs = 3.5)$ind, outliersJM(allProportions_Table[,24], SDs = 3.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0(data.dir,"/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/ScatterPlot_", allProportions_Table[ind.with.store[k],"Barcodes"], ".png")
  file.copy(name, paste0(results.dir,"/Transitional-Bcells/"))
}

outliers.Transitional.Bcells <- as.matrix(allProportions_Table[ind.with.store,c(1:8,24)])
save(outliers.Transitional.Bcells, file = paste0(results.dir,"/Transitional-Bcells/outliers.Transitional.Bcells.Rdata"))

#############################################################################################################


## Finding outliers for Fo B-cells
boxplot(allProportions_Table[,25], pch=19, cex=0.2, main = "Fo B-cells", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,25], SDs = 3.5)$ind, outliersJM(allProportions_Table[,25], SDs = 3.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0(data.dir,"/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/ScatterPlot_", allProportions_Table[ind.with.store[k],"Barcodes"], ".png")
  file.copy(name, paste0(results.dir,"/Fo-Bcells/"))
}

outliers.Fo.Bcells <- as.matrix(allProportions_Table[ind.with.store,c(1:8,25)])
save(outliers.Fo.Bcells, file = paste0(results.dir,"/Fo-Bcells/outliers.Fo.Bcells.Rdata"))

#############################################################################################################

## Finding outliers for IgG+ Memory B-cells
boxplot(allProportions_Table[,26], pch=19, cex=0.2, main = "IgG+ Memory B-cells", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,26], SDs = 3.5)$ind, outliersJM(allProportions_Table[,26], SDs = 3.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0(data.dir,"/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/ScatterPlot_", allProportions_Table[ind.with.store[k],"Barcodes"], ".png")
  file.copy(name, paste0(results.dir,"/IgG+Memory-Bcells/"))
}

outliers.IgGposMemory.Bcells <- as.matrix(allProportions_Table[ind.with.store,c(1:8,26)])
save(outliers.IgGposMemory.Bcells, file = paste0(results.dir,"/IgG+Memory-Bcells/outliers.IgGposMemory.Bcells.Rdata"))

#############################################################################################################

## Finding outliers for MZP B cell and MZ B cells
boxplot(allProportions_Table[,27], pch=19, cex=0.2, main = "MZP B-cells", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,27], SDs = 3.5)$ind, outliersJM(allProportions_Table[,27], SDs = 3.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0(data.dir,"/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/ScatterPlot_", allProportions_Table[ind.with.store[k],"Barcodes"], ".png")
  file.copy(name, paste0(results.dir,"/MZP-Bcells/"))
}

outliers.MZP.Bcells <- as.matrix(allProportions_Table[ind.with.store,c(1:8,27,28)])
save(outliers.MZP.Bcells, file = paste0(results.dir,"/MZP-Bcells/outliers.MZP.Bcells.Rdata"))

#############################################################################################################


## Finding outliers for  Transitional-2 and Transitional-1
boxplot(allProportions_Table[,29], pch=19, cex=0.2, main = "Transitional-2", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,29], SDs = 4)$ind, outliersJM(allProportions_Table[,29], SDs = 4)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0(data.dir,"/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/ScatterPlot_", allProportions_Table[ind.with.store[k],"Barcodes"], ".png")
  file.copy(name, paste0(results.dir,"/Transitional-2/"))
}

outliers.Transitional2 <- as.matrix(allProportions_Table[ind.with.store,c(1:8,29,30)])
save(outliers.Transitional2, file = paste0(results.dir,"/Transitional-2/outliers.Transitional2.Rdata"))

#############################################################################################################

## Finding outliers for  FoI and FoII
boxplot(allProportions_Table[,31], pch=19, cex=0.2, main = "FoI", ylab="Proportion")
indices.check <- cbind(outliersJM(allProportions_Table[,31], SDs = 4)$ind, outliersJM(allProportions_Table[,31], SDs = 4)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0(data.dir,"/Figures/ScatterPlots/",
                 allProportions_Table[ind.with.store[k],"Genotype"], "/ScatterPlot_", allProportions_Table[ind.with.store[k],"Barcodes"], ".png")
  file.copy(name, paste0(results.dir,"/FoI/"))
}

outliers.FoI <- as.matrix(allProportions_Table[ind.with.store,c(1:8,31,32)])
save(outliers.FoI, file = paste0(results.dir,"/FoI/outliers.FoI.Rdata"))

#############################################################################################################







