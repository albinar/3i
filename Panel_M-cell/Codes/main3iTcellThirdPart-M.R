# Developed by Albina Rahim

# This code creates RchyOptimyx plots from Matrix3iTcell.csv.

#setwd("/home/arahim/Desktop/3iTcell_AFC_Oct_2014/")
setwd("/data/Panel_M-cell")

library("flowCore")
library("flowBin")
library("flowDensity")
library("flowType")
library("flowClean")
library("RchyOptimyx")
source("Codes/3iTcellfunctions-M.R")

start <- Sys.time()

load( file =  paste("Results/Phenotypes.Rdata",sep=""))
load( file =  paste("Results/Genotypes.Rdata",sep=""))
load( file =  paste("Results/GenotypesLong.Rdata",sep=""))
load( file =  paste("Results/allFT.Rdata",sep=""))
load( file =  paste("Results/flowType.res_Sample.Rdata",sep=""))
load("Results/GlobalFrame657.Rdata")

dir.create ( "Results/WilcoxTest")
dir.create ( "Results/T_Test")

#Matrix3iTcell <- read.csv("Results/Matrix3iTcell.csv")
load(file =  paste("Results/Matrix3iTcell.Rdata",sep=""))

TotalAdjPvals <- NULL

uniqueWT <- unique(Genotypes)

nameGenotypes <- (Genotypes[duplicated(Genotypes)])
#nameGenotypes <- nameGenotypes[which(nameGenotypes == uniqueWT[3])[1]:length(nameGenotypes)]    # uniqueWT[3] only contains 1 file for me
nameGenotypes <- nameGenotypes[which(nameGenotypes == uniqueWT[4])[1]:length(nameGenotypes)]
nameGenotypes <- unique(nameGenotypes)

WTindices <- c(which(Genotypes == uniqueWT[1]), which(Genotypes == uniqueWT[2]) )

# ########################################################################################################################
# # Shapiro-Wilk test for normality of WT data
# library(dplyr)
# library(broom)
# c.idx <- intersect(which(Matrix3iTcell[,1] == Matrix3iTcell[,2]), which(Matrix3iTcell[,2] == Matrix3iTcell[,3])) # rows of constants
# shap.wilk <- apply(Matrix3iTcell[-c.idx, WTindices], 1, function(x) tidy(shapiro.test(x)))
# shap.wilk <- do.call('rbind',shap.wilk)
# x11()
# plot(rowSums(Matrix3iTcell[-c.idx, WTindices]), t(shap.wilk['p.value']))
# x11()
# # plotting histograns for some of the rows (phenotypes) with the highest p values
# par(mfrow=c(2, 2)) 
# temp <- Matrix3iTcell[-c.idx, WTindices]
# hist(temp[440,], breaks = 100)
# hist(temp[1313,], breaks = 100)
# hist(temp[3191,], breaks = 100)
# hist(temp[8450,], breaks = 100)
# ########################################################################################################################

ForLoopIndices <- 3:length(uniqueWT)
l <- nrow(Matrix3iTcell)

for ( k1 in ForLoopIndices){

    print(paste(k1,": Starting", uniqueWT[k1]))
    KOindices <- c(which(Genotypes == uniqueWT[k1]))
    startP <- Sys.time()

    if(length(KOindices) <= 2){
        print("Skipped, 2 or less files")
        next
    }
    
    Labels <- c(WTindices, KOindices)

    WCpValues[2:l] <- apply(Matrix3iTcell[2:l,], 1, function(x) { (wilcox.test(x[WTindices], x[KOindices]))$p.value })
    Pvals[2:l] <- apply(Matrix3iTcell[2:l,], 1, function(x) { (t.test(x[WTindices], x[KOindices]))$p.value })
    WCpValues[1] <- 1
    Pvals[1] <- 1
    
    Pvals    [is.nan(Pvals)]     <- 1
    WCpValues[is.nan(WCpValues)] <- 1

    names(Pvals) <- Phenotypes
    names(WCpValues) <- Phenotypes
    
    cat("Time is: ",TimeOutput(startP),sep="")

    save(Pvals, file=paste0("Results/T_Test/", uniqueWT[k1], ".Rdata", sep=""))
    save(WCpValues, file=paste0("Results/WilcoxTest/", uniqueWT[k1], ".Rdata", sep=""))
}



#     Pvals       <- vector()
#     WCpValues   <- vector()
#     
#     for (r1 in 1:dim(Matrix3iTcell)[1]){
# 
#       
#         # for the column with all proportions of 1 
#         if (r1 == 1){
#             Pvals[r1] <- 1
#             WCpValues[r1] <- 1
#             next;
#         }
#         # other proportions
#         temp <- t.test(Matrix3iTcell[r1, WTindices], Matrix3iTcell[r1, KOindices])
#         Pvals[r1] <- temp$p.value
#         temp <- wilcox.test(Matrix3iTcell[r1, WTindices], Matrix3iTcell[r1, KOindices] )
#         WCpValues[r1] <- temp$p.value
# 
#         
#     }
