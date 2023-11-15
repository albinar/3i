# Developed by Albina Rahim
# Date: June 09, 2016

# This code finds if any gates are outliers and updated them to an average of all gates for the Bone Marrow Panel.
remove(list=ls())

setwd("/code/Projects/3i/Panel_BM-cell/")

library("flowCore")
library("arrayMvout")

source("Codes/3iTcellfunctions-BM.R")
source("Codes/rotate.data-BM.R")

# Creating directories
suppressWarnings ( dir.create ( "Results/Gating_Thresholds_Updated") )
suppressWarnings ( dir.create ( "Results/Gating_Thresholds_CD45_Updated") )

load("Results/uniqueGT.Rdata")
#load("Results/channels.ind.Rdata")
load("Results/store.allFCS.Rdata")

start <- Sys.time()

# path to the FCS_Groups files for the Bone Marrow Panel
pathFCS <- paste("/code/Projects/3i/Panel_BM-cell/FCS_Groups")

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

colnames(ThreshTable) <- c("singlets.gate.h", "singlets.gate.l", "live.gate", "fsc.a.gate.high", "fsc.a.gate.low", "ssc.a.gate", "cd45.gate.low", "cd45.gate.high", "cd43.gate.high",
                           "gr1.gate", "b220.gate", "cd3.gate", "cd138.gate", "cd11b.gate", "cd43.gate", "cd24.gate", "igm.gate", "igd.gate", "bp1.gate")

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


# ## Optional: The following was done because we observed that for the channel cd24.gate, the outliers for some reason were not fixed.
# ## We have to find a good solution for this. But for now, anything that is less than 1, we replace it with the average
# cd24.mean <- mean(ThreshTable[,c("cd24.gate")])
# cd24.indexToChange <- which(ThreshTable[,c("cd24.gate")] < 1.0)
# for(i in 1:length(cd24.indexToChange)){
#   ThreshTable[cd24.indexToChange[i], c("cd24.gate")] <- cd24.mean
# }
save (ThreshTable, file =  paste("Results/ThreshTableUpdated.Rdata",sep="") )  

# save all updated gthres
for(i in 1:nrow(store.allFCS)){
    gthres <- ThreshTable[i,]
    gates.cd45 <- ThreshTable[i,c("igd.gate", "cd43.gate", "cd24.gate", "gr1.gate", "live.gate", "igm.gate", "cd11b.gate", "cd45.gate.low", "cd138.gate", "cd3.gate", "bp1.gate", "b220.gate")]
    #colnames(gates.cd45) <- c("IgD", "CD43", "CD24", "GR1", "Live", "IgM", "CD11b", "CD45", "CD138", "CD3", "BP1", "B220")
    
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

