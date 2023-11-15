# Developed by Albina Rahim
# Date: April 19, 2016

remove(list=ls())

# This code creates Matrix3iTcell.csv, which is a matrix of all the proportions. Each file vertically, each phenotype horizontally.

#setwd("/code/Projects/KCL_WTSI_2015-08-27/")
setwd("/data/Panel_M-cell")


library("flowCore")
library("flowBin")
library("flowDensity")
library("flowType")

source("Codes/3iTcellfunctions-M.R")
source("Codes/rotate.data-M.R")

load("Results/uniqueGT.Rdata")
load("Results/channels.ind.Rdata")
load("Results/store.allFCS.Rdata")
load (file =  paste("Results/After_Clean/", store.allFCS[1,1],"/AftClean_", store.allFCS[1,2],".Rdata",sep="") )

start <- Sys.time()

# path to the flowType files
pathFT <- "Results/FlowType/"
# Reads all folders and files in current path folder and makes a list of all of their paths
allFT  <- dir(pathFT,  full.names=T, recursive = F) 

wildTypes <- c(allFT[which(regexec("[+]_[+]", allFT) !=-1)], allFT[which(regexec("[+]_Y"  , allFT) !=-1)])

allFT <- allFT[-which(regexec("[+]_[+]", allFT) !=-1)]
allFT <- allFT[-which(regexec("[+]_Y"  , allFT) !=-1)]

allWildTypes <- dir(wildTypes,  full.names=T, recursive = T)
allKnockOuts  <- dir(allFT,  full.names=T, recursive = T)

# # Remove files with too few cells and which have an abnormal distn
# files.to.remove.barcodes <- c(43170:43182, 43215:43330, 44310:44323, 88659, 88655, 90915, 95648, 72301, 72306)
# FT.remove.idx <- lapply(files.to.remove.barcodes, function(x){ grep(x, allWildTypes)})
# FT.remove.idx <- unlist(FT.remove.idx[lapply(FT.remove.idx, length)>0])
# allWildTypes <- allWildTypes[-FT.remove.idx]
# FT.remove.idx <- lapply(files.to.remove.barcodes, function(x){ grep(x, allKnockOuts)})
# FT.remove.idx <- unlist(FT.remove.idx[lapply(FT.remove.idx, length)>0])
# allKnockOuts <- allKnockOuts[-FT.remove.idx]

KnockOutTypes <- lapply ( 1:length(allKnockOuts), function(x){
                strsplit(allKnockOuts[x],split="/")[[1]][4]
}) # Extracting the Knock-out names
WildTypeTypes <- lapply ( 1:length(allWildTypes), function(x){
                strsplit(allWildTypes[x],split="/")[[1]][4]
}) # Extracting the Wild type names

Genotypes <- unlist(c(WildTypeTypes, KnockOutTypes)) # Contains the names of Wild Types and Knock-outs
GenotypesLong <- c(allWildTypes, allKnockOuts) # Contains the whole path of all the FT files

allFT <- c(allWildTypes, allKnockOuts)
save( allFT, file =  paste("Results/allFT.Rdata",sep=""))

load (file =  allFT[1] )
save( flowType.res, file =  paste("Results/flowType.res_Sample.Rdata",sep=""))

PropMarkers <- flowType.res@PropMarkers
Phenotypes <- unlist(lapply(flowType.res@PhenoCodes, function(x){return(decodePhenotype(
  x,as.vector(f@parameters@data$desc[PropMarkers]), flowType.res@PartitionsPerMarker))}))

## Matrix where each column is an individual file, each row phenotype. Contains the proportions of each phenotype for each file
Matrix3iTcell <- matrix(0, calcNumPops(c(2,2,2,2,2,2,2,2,2), 9), length(allFT) )

## This is the simplest way of Normalization. But we may need to change this part based on Finak's paper.
for ( q1 in 1:length(allFT) ) {
    print(paste(q1, ": Starting ", allFT[q1], sep="" ))
    load (file =  allFT[q1] )
    Proportions <- flowType.res@CellFreqs; ## Containing the cell frequencies measured for each phenotype
    Proportions <- Proportions / max(Proportions) # Scaling between the values of 0 to 1. This is a Normalization (See the RchyOptimyx Vigenette)
    Matrix3iTcell[,q1] <- Proportions
}
rownames(Matrix3iTcell) <- Phenotypes

save( Phenotypes, file =  paste("Results/Phenotypes.Rdata",sep=""))
write.table(Matrix3iTcell, "Results/Matrix3iTcell.csv", sep=",", row.names=F)
save( Genotypes, file =  paste("Results/Genotypes.Rdata",sep=""))
save( GenotypesLong, file =  paste("Results/GenotypesLong.Rdata",sep=""))
save(Matrix3iTcell, file =  paste("Results/Matrix3iTcell.Rdata",sep=""))

cat("Time is: ",TimeOutput(start),sep="")



