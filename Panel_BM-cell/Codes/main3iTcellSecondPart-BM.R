# Developed by Albina Rahim
# Date: August 02, 2016

remove(list=ls())

# This code creates Matrix3iTcell.csv, which is a matrix of all the proportions. Each file vertically, each phenotype horizontally.

setwd("/code/Projects/3i/Panel_BM-cell/")

library("flowCore")
library("flowBin")
library("flowDensity")
library("flowType")
library("stringr")

source("Codes/3iTcellfunctions-BM.R")
source("Codes/rotate.data-BM.R")

load("Results/uniqueGT.Rdata")
load("Results/channels.ind.NoLiveNoCD45.Rdata")
## Loading matrix which contains information for unique FCS files (no duplicates based on barcodes)
load("Results/store.allFCS.unique.Rdata")
load (file =  paste("Results/After_Clean/", store.allFCS.unique[1,1],"/AftClean_", store.allFCS.unique[1,2],".Rdata",sep="") )

start <- Sys.time()

# path to the flowType files
pathFT <- paste("/code/Projects/3i/Panel_BM-cell/Results/FlowType")
# Reads all folders and files in current path folder and makes a list of all of their paths
allFT  <- dir(pathFT,  full.names=T, recursive = F) 

wildTypes <- c(allFT[which(regexec("[+]_[+]", allFT) !=-1)], allFT[which(regexec("[+]_Y"  , allFT) !=-1)])

allFT <- allFT[-which(regexec("[+]_[+]", allFT) !=-1)]
allFT <- allFT[-which(regexec("[+]_Y"  , allFT) !=-1)]

## Contains information for duplicated FCS WT files (duplicate Barcodes)
allWildTypes <- dir(wildTypes,  full.names=T, recursive = T) 
## Removing WTs files with duplicated barcodes
allWildTypes <- cbind(allWildTypes, sapply(1:length(allWildTypes), function(x){str_extract(strsplit(allWildTypes[x],split="/")[[1]][9],"L[0-9]+")})) 
allWildTypes <- allWildTypes[!duplicated(allWildTypes[,2]),]

## Contains information for duplicated FCS KO files (duplicate Barcodes)
allKnockOuts  <- dir(allFT,  full.names=T, recursive = T)
## Removing KOs with duplicated barcodes
allKnockOuts <- cbind(allKnockOuts, sapply(1:length(allKnockOuts), function(x){str_extract(strsplit(allKnockOuts[x],split="/")[[1]][9],"L[0-9]+")})) 
allKnockOuts <- allKnockOuts[!duplicated(allKnockOuts[,2]),]


WildTypeTypes <- unlist(lapply ( 1:nrow(allWildTypes), function(x){
  strsplit(allWildTypes[x,1],split="/")[[1]][8]
})) # Extracting the Wild type names
KnockOutTypes <- unlist(lapply ( 1:nrow(allKnockOuts), function(x){
  strsplit(allKnockOuts[x,1],split="/")[[1]][8]
})) # Extracting the Knock-out names

Genotypes <- unlist(c(WildTypeTypes, KnockOutTypes)) # Contains the names of Wild Types and Knock-outs
#GenotypesLong <- c(allWildTyes, allKnockOuts) # Contains the whole path of all the FT files

# save( Genotypes, file =  paste("../3iTcell_Results/Genotypes.Rdata",sep=""))
# save( GenotypesLong, file =  paste("../3iTcell_Results/GenotypesLong.Rdata",sep=""))

allFT <- c(allWildTypes[,1], allKnockOuts[,1])  # Contains the whole path of all the FT files

save( allFT, file =  paste("Results/allFT.Rdata",sep=""))
load (file =  allFT[1] )
save( flowType.res, file =  paste("Results/flowType.res_Sample.Rdata",sep=""))

## Matrix where each column is an individual file, each row phenotype. Contains the proportions of each phenotype for each file
Matrix3iTcell <- matrix(0, calcNumPops(c(rep(2,10)), 8), length(allFT) )

## This is the simplest way of Normalization. But we may need to change this part based on Finak's paper.
for ( q1 in 1:length(allFT) ) {
    print(paste(q1, ": Starting ", allFT[q1], sep="" ))
    load (file =  allFT[q1] )
    PropMarkers <- flowType.res@PropMarkers
    Proportions <- flowType.res@CellFreqs; ## Containing the cell frequencies measured for each phenotype
    Proportions <- Proportions / max(Proportions) # Scaling between the values of 0 to 1. This is a Normalization (See the RchyOptimyx Vigenette)
    Matrix3iTcell[,q1] <- Proportions
}

print("Completed inserting Proportions into Matrix3iTcell matrix")

Phenotypes <- unlist(lapply(flowType.res@PhenoCodes, function(x){return(decodePhenotype(
                          x,as.vector(f@parameters@data$desc[PropMarkers]), flowType.res@PartitionsPerMarker))}))

rownames(Matrix3iTcell) <- Phenotypes


save( Phenotypes, file =  paste("Results/Phenotypes.Rdata",sep=""))
save( Genotypes, file =  paste("Results/Genotypes.Rdata",sep=""))
save( allWildTypes, file =  paste("Results/allWildTypes.Rdata",sep=""))
save( allKnockOuts, file =  paste("Results/allKnockOuts.Rdata",sep=""))
write.table(Matrix3iTcell, "Results/Matrix3iTcell.csv", sep=",", row.names=F)

cat("Time is: ",TimeOutput(start),sep="")





