# Developed by Albina Rahim
# Date: July 11, 2016

remove(list=ls())

# This code creates RchyOptimyx plots from Matrix3iTcell.csv.

setwd("/code/Projects/3i/Panel_T-cell_MLN/")

library("flowCore")
library("flowBin")
library("flowDensity")
library("flowType")
library("flowClean")
library("RchyOptimyx")
source("Codes/3iTcellfunctions.R")

start <- Sys.time()

load( file =  paste("Results/Phenotypes.Rdata",sep=""))
load( file =  paste("Results/Genotypes.Rdata",sep=""))
#load( file =  paste("Results/GenotypesLong.Rdata",sep=""))
load( file =  paste("Results/allFT.Rdata",sep=""))
load( file =  paste("Results/flowType.res_Sample.Rdata",sep=""))
load("Results/GlobalFrame657.Rdata")

dir.create ( "Results/WilcoxTest")
#dir.create ( "Results/T_Test")
dir.create ( "Results/RchyOptimyx")
dir.create ( "Results/RchyOptimyx_Non_Adj")


Matrix3iTcell <- read.csv("Results/Matrix3iTcell.csv")

TotalAdjPvals <- NULL

uniqueWT <- unique(Genotypes)

# nameGenotypes <- (Genotypes[duplicated(Genotypes)])
# nameGenotypes <- nameGenotypes[which(nameGenotypes == uniqueWT[3])[1]:length(nameGenotypes)]
# nameGenotypes <- unique(nameGenotypes)

WTindices <- c(which(Genotypes == uniqueWT[1]), which(Genotypes == uniqueWT[2]) )

flowType.resWT <- list()
# flowType.resWT_Avg <- matrix(0, calcNumPops(c(3,2,2,2,2,2,2,2,3,2), 10), 1 )
flowType.resWT_Avg <- matrix(0, calcNumPops(c(rep(2,10)), 8), 1 )

# Average Cell Frequency of the WTs
for ( r1 in 1:length(WTindices)){
    load (file =  allFT[WTindices[r1]] )
    flowType.resWT[[r1]] <- flowType.res@CellFreqs  
    flowType.resWT_Avg <- flowType.resWT_Avg + flowType.res@CellFreqs
}
flowType.resWT_Avg <- flowType.resWT_Avg / length(WTindices)


ForLoopIndices <- 3:length(uniqueWT)

for ( k1 in ForLoopIndices){
    print(paste(k1,": Starting", uniqueWT[k1]))
    KOindices <- c(which(Genotypes == uniqueWT[k1]))

    if(length(KOindices) <= 2){
        print("Skipped, 2 or less files")
        next
    }
    
    Labels <- c(WTindices, KOindices) # Contains the indices of both the WildTypes and KnockOutTypes
    
    # Pvals       <- vector()
    WCpValues   <- vector()
    print("Starting Wilcoxon test")
    for (r1 in 1:dim(Matrix3iTcell)[1]){
        # for the row with all proportions of 1 
        if (r1 == 1){
          print("Row with value 1")
            # Pvals[r1] <- 1
            WCpValues[r1] <- 1
            next;
        }
        # other proportions
        # print("Starting t.test")
        # temp <- t.test(Matrix3iTcell[r1, WTindices], Matrix3iTcell[r1, KOindices])
        # Pvals[r1] <- temp$p.value
       
        temp <- wilcox.test( as.numeric(Matrix3iTcell[r1, WTindices]), as.numeric(Matrix3iTcell[r1, KOindices]) )
        WCpValues[r1] <- temp$p.value
    }
    
    print("End of Wilcoxon test")
    # Pvals    [is.nan(Pvals)]     <- 1
    WCpValues[is.nan(WCpValues)] <- 1

    # names(Pvals) <- Phenotypes
    names(WCpValues) <- Phenotypes
    
    print("Saving Wilcoxon values")
    # save(Pvals, file=paste0("Results/T_Test/", uniqueWT[k1], ".Rdata", sep=""))
    save(WCpValues, file=paste0("Results/WilcoxTest/", uniqueWT[k1], ".Rdata", sep=""))
}

write.table(flowType.resWT_Avg, "Results/flowType.resWT_Avg.csv", sep=",", row.names=F)

