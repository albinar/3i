# Developed by Albina Rahim
# Date: April 13, 2016

# This code finds if any gates are outliers and updated them to an average of all gates.
remove(list=ls())

setwd("/code/Projects/3i/Panel_T-cell_MLN/")

library("flowCore")
library("flowBin")
library("flowDensity")
library("flowType")
library("flowClean")
library("arrayMvout")

source("Codes/3iTcellfunctions.R")
source("Codes/rotate.data.R")

# Creating directories
suppressWarnings ( dir.create ( "Results/Gating_Thresholds_Updated") )
suppressWarnings ( dir.create ( "Results/Gating_Thresholds_CD45_Updated") )

load("Results/uniqueGT.Rdata")
#load("Results/channels.ind.Rdata")
load("Results/store.allFCS.Rdata")

start <- Sys.time()

# path to the FCS_Groups files
pathFCS <- paste("/code/Projects/3i/Panel_T-cell_MLN/FCS_Groups")

# Reads all folders and files in current path folder and makes a list of all of their paths
allFCS <- dir(pathFCS, full.names=T, recursive=T) 


# Create sub-directories
invisible(sapply(1:length(uniqueGT), function(x){
  suppressWarnings(dir.create ( paste("Results/Gating_Thresholds_Updated/", uniqueGT[x], sep="") ))
  suppressWarnings(dir.create ( paste("Results/Gating_Thresholds_CD45_Updated/", uniqueGT[x], sep="") ))
}))

# load all gthres into ThreshTable
ThreshTable <- NULL

for(i in 1:nrow(store.allFCS)){
  load(file =  paste("Results/Gating_Thresholds/", store.allFCS[i,1],"/Gthres_", store.allFCS[i,2],".Rdata",sep=""))
  ThreshTable <- rbind(ThreshTable, gthres)
}

rownames(ThreshTable) <- 1:length(ThreshTable[,1])
colnames(ThreshTable) <- c("singlets.gate.h", "singlets.gate.l", "live.gate", "fsc.a.gate.high", "fsc.a.gate.low", "ssc.a.gate", "cd45.gate.low", "cd45.gate.high", "cd4.gate.slant", "tcrd.gate.slanta", "tcrd.gate.slantb", "tcrd.gate.low", "tcrd.gate.high", "cd4.gate",
                           "cd62.gate", "cd44.gate", "klrg1.gate", "gitr.gate", "cd5.gate", "cd161.gate.low", "cd161.gate.high", "cd8.gate", "cd44.gate.high", "cd25.gate.low", "cd25.gate.high")

save (ThreshTable, file =  paste("Results/ThreshTable.Rdata",sep="") )  

ThreshTableOld <- ThreshTable

if ( !file.exists(file =  paste("Results/GoodFakeData.Rdata",sep=""))){
    # Create fake data that has a a small standard deviation and no outliers to use to trick ArrayOutliers to find outliers from a vector.
    GoodFakeData <- NULL
    for ( x in 1:4 ) {
        GoodFakeData <- cbind(GoodFakeData, runif(nrow(ThreshTable),5,5.5))
    }
    save(GoodFakeData, file =  paste("Results/GoodFakeData.Rdata",sep=""))
} else {
    load (file =  paste("Results/GoodFakeData.Rdata",sep=""))
}


Average <- NULL
sink(file="Results/Threshold_Change.txt", append = T)
for ( x in 1:ncol(ThreshTable)) { 
    Average[x] <- mean(ThreshTable[ which(!is.na(ThreshTable[,x])) ,x])
    ThreshTable[ which(is.na(ThreshTable[,x])) ,x] <- Average[x]
    
    Res <- ArrayOutliers(data.frame(cbind(ThreshTable[,c(x)],GoodFakeData)), alpha=0.5)$outl # find all outliers for each marker
    temp <- as.numeric(rownames(Res))

    if ( length( temp )  != 0 ) {
        print(paste("Changing ", store.allFCS[temp,1], " / ", store.allFCS[temp,2], " from ", ThreshTable[temp, x],
                    " to ", Average[x], "; Marker: ", x, "; Marker Name: ", colnames(ThreshTable)[x], sep=""))
        ThreshTable[temp, x] <-  Average[x]
    }
}
sink()

save (ThreshTable, file =  paste("Results/ThreshTableUpdated.Rdata",sep="") )  

# save all updated gthres & gates.cd45
for(i in 1:nrow(store.allFCS)){
    gthres <- ThreshTable[i,]
    gates.cd45 <- ThreshTable[i,c("cd44.gate", "cd62.gate", "cd25.gate.high", "cd8.gate", "live.gate", "klrg1.gate", "cd5.gate", "cd45.gate.low", "cd161.gate.low", "cd4.gate", "gitr.gate", "tcrd.gate.low")]
    save ( gthres, file =  paste("Results/Gating_Thresholds_Updated/", store.allFCS[i,1],"/Gthres_", store.allFCS[i,2],".Rdata",sep="") )
    save ( gates.cd45, file =  paste("Results/Gating_Thresholds_CD45_Updated/", store.allFCS[i,1],"/Gthres_", store.allFCS[i,2],".Rdata",sep="") )
}

## Optional Part: Plotting threshold value plots for the markers for all the FCS files (WTs+ KOs) (This was done for Adam)
suppressWarnings ( dir.create ( "Results/Figures/Gating_Threshold_Plots") )
for(i in 1:ncol(ThreshTable)){
  png(file = paste("Results/Figures/Gating_Threshold_Plots/", "Plot_", colnames(ThreshTable)[i], ".png", sep = ""))
  #max.x <- max(ThreshTable[,i]) + 0.5
  #min.x <- min(ThreshTable[,i]) -0.5
    plot(ThreshTable[,i], main = colnames(ThreshTable)[i], devn = T, ylab= "Threshold Values")
    dev.off()
}

cat("Total time is: ",TimeOutput(start),"\n",sep="")

