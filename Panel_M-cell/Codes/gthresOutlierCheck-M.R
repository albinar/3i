# Developed by Albina Rahim
# Date: June 09, 2016

# This code finds if any gates are outliers and updated them to an average of all gates.
remove(list=ls())

#setwd("/data/Panel_M-cell")
setwd("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_M-cell/MLN/")

library("arrayMvout")
library("flowCore")
# library("foreach")

source("/data/Panel_M-cell/Codes/3iTcellfunctions-M.R")
source("/data/Panel_M-cell/Codes/rotate.data-M.R")

results.dir <- "Results"


# Creating directories
suppressWarnings(dir.create(paste0(results.dir, "/Gating_Thresholds_Updated")))

#load(paste0(results.dir, "/store.allFCS.Rdata"))

start <- Sys.time()

# load all gthres into ThreshTable
load(file = paste0(results.dir, "/ThreshTable.Rdata"))  
ThreshTable <- gthres
ThreshTableOld <- ThreshTable

save (gthres, file = paste0(results.dir, "/ThreshTable.Rdata")) 


if(!file.exists(file = paste(results.dir, "/GoodFakeData.Rdata"))){
    # Create fake data that has a a small standard deviation and no outliers to use to trick ArrayOutliers to find outliers from a vector.
    GoodFakeData <- NULL
    for ( x in 1:4 ) {
        GoodFakeData <- cbind(GoodFakeData, runif(nrow(ThreshTable),5,5.5))
    }
    save(GoodFakeData, file = paste0(results.dir, "/GoodFakeData.Rdata"))
} else {
    load(file = paste(results.dir, "/GoodFakeData.Rdata"))
}

Average <- NULL
sink(file = paste0(results.dir, "/Threshold_Change.txt"), append = T)

thres.to.correct <- 6:ncol(ThreshTable)
# I am not correcting the linneg.gate.slant threshold because all the corrections here were including cells that they should not have been including
# For the MLN data set cannot it doensn't make sense to find outliers for linneg.macneg.gate.slanta and linneg.macneg.gate.slantb because 
# there are two different ways that files can be gated
# See main3iTcell-M_MLN.R for more details
# indices.to.remove <- c(which(colnames(ThreshTable) == 'linneg.gate.slant'), which(colnames(ThreshTable) == 'lin.gate.slant'),
#                        which(colnames(ThreshTable) == 'fsc.a.gate.low'), which(colnames(ThreshTable) == 'lincd3cd19cd161.gate.high')) # SPLEEN
# indices.to.remove <- c(which(colnames(ThreshTable) == 'linneg.gate.slant'), which(colnames(ThreshTable) == 'lin.gate.slant'), # MLN
#                        which(colnames(ThreshTable) == 'linneg.macneg.gate.slanta'), which(colnames(ThreshTable) == 'linneg.macneg.gate.slantb'),
#                        which(colnames(ThreshTable) == 'fsc.a.gate.low'))
#thres.to.correct <- thres.to.correct[-which(thres.to.correct %in% indices.to.remove)]

for ( x in thres.to.correct) { 
    Average[x] <- mean(ThreshTable[ which(!is.na(ThreshTable[,x])) ,x])
    ThreshTable[ which(is.na(ThreshTable[,x])) ,x] <- Average[x]
  
    Res <- ArrayOutliers(data.frame(cbind(ThreshTable[,c(x)],GoodFakeData)), alpha=0.5)$outl # find all outliers for each marker
    temp <- as.numeric(rownames(Res))

    if ( length( temp )  != 0 ) {
        print(paste("Changing ", ThreshTable[temp, 4], " from ", ThreshTable[temp, x],
                    " to ", Average[x], "; Marker: ", x, "; Marker Name: ", colnames(ThreshTable)[x], sep=""))
        ThreshTable[temp, x] <-  Average[x]
    }
}
sink()
closeAllConnections() 

gthres <- ThreshTable
save(gthres, file =  paste0(results.dir, "/ThreshTableUpdated.Rdata"))  

## Optional Part: Plotting threshold value plots for the markers for all the FCS files (WTs+ KOs) (This was done for Adam)
suppressWarnings(dir.create(paste0(results.dir, "/Figures/Gating_Threshold_Plots")))
for(i in 6:ncol(ThreshTable)){
    png(file = paste0(results.dir, "/Figures/Gating_Threshold_Plots/", "Plot_", colnames(ThreshTable)[i], ".png"))
    plot(ThreshTable[,i], main = colnames(ThreshTable)[i], ylab= "Threshold Values") 
    dev.off()
}

cat("Total time is: ",TimeOutput(start),"\n",sep="")

