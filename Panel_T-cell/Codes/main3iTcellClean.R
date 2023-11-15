# Originally written by Justin Meskas 
# Modified & Updated by Albina Rahim
# Date: March 16, 2016
remove(list=ls())
setwd("/code/Projects/3i/Panel_T-cell_MLN/")

load("Results/lgl657.Rdata")
load("Results/Genotype.Rdata")
load("Results/uniqueGT.Rdata")
load("Results/channels.ind.Rdata")
load("Results/store.allFCS.Rdata")


library("flowCore")
library("flowBin")
library("flowDensity")
library("flowType")
library("flowClean")


source("Codes/3iTcellfunctions.R")


suppressWarnings ( dir.create ( "Results/After_Clean") )
#suppressWarnings ( dir.create ( "Results/Before_Transform") )
suppressWarnings ( dir.create ( "Results/FlowType") )
suppressWarnings ( dir.create ( "Results/Figures") )
suppressWarnings ( dir.create ( "Results/Figures/Clean") )
suppressWarnings ( dir.create ( "Results/Figures/ScatterPlots") )

verbose_debris <- T
verbose_margin <- T
verbose_flowClean <- T


start <- Sys.time()

# path to the FCS_Groups files
pathFCS <- paste("/code/Projects/3i/Panel_T-cell_MLN/FCS_Groups")

# Reads all folders and files in current path folder and makes a list of all of their paths
allFCS <- dir(pathFCS, full.names=T, recursive=T) 

invisible(sapply(1:length(uniqueGT), function(x){
  suppressWarnings(dir.create ( paste("Results/FlowType/", uniqueGT[x], sep="") ))
  suppressWarnings ( dir.create ( paste("Results/After_Clean/", uniqueGT[x], sep="") ) )
  #suppressWarnings ( dir.create ( paste("Results/Before_Transform/", uniqueGT[x], sep="") ) )
  suppressWarnings ( dir.create ( paste("Results/Figures/ScatterPlots/", uniqueGT[x], sep="") ) )}))

#run flowClean only for 5-6 markers, based on what Herve suggested
#that if there's a clog in one channel, it affects other channels as well
clean.chan <- as.vector(sample(channels.ind,6))
index.Clean.Error <- 0
error.Message <- NA

for (i in 1:length(allFCS)){
      
        start2 <- Sys.time()
        print(paste(i, ": Starting", store.allFCS[i,1],"/",store.allFCS[i,2], sep = ""))
        
        # read in one flowFrame
        f <- read.FCS(allFCS[i])
        
        ## Removing Margin events------------------------------------------------------------------------------------------
        scat.chans <- c(grep (colnames(f),pattern = "FSC*"),grep (colnames(f),pattern = "SSC*"))
        ## Removing margin events in Scatter channels
        f <- removeMargins(f, chans=scat.chans, verbose= verbose_margin)
        #Removing negative values in scatter channels
        f <- removeMargins(f, chans=scat.chans, debris=T, neg=T, verbose= verbose_debris) 
        
        
        ## Compensation
        print("Starting Compensation")
        if( det(f@description$SPILL)==1 ){
          print("Check the spillover matrix, it's probably an identity matrix!")
          return(NULL)
        } else{
          f <- compensate(f, f@description$SPILL)}
        
        # ## Saving the data before Transformation (will be needed for the MFI calculation)
        # f.beforeTransform <- f
        
        ## Transformation
        print("Starting Transformation")
        f <- transform(f, lgl657)

        # flowClean--------------------------------------------------------------------------------------------------------
        #set.seed(5678)
          print("Starting flowClean")
          removed.cell<- c()
          f.clean <- try(clean(f, vectMarkers=clean.chan, diagnostic = T,binSize = 0.005,nCellCutoff = 300, announce = verbose_flowClean,
                               ext="fcs",filePrefixWithDir=paste("Results/Figures/Clean/", store.allFCS[i,2], sep="")))
          to.be.removed <-c()
          if(class(f.clean)=="try-error")
          {
            print("flowClean with the current values doesn't work")
            save ( f , file =  paste("Results/After_Clean/", store.allFCS[i,1],"/AftClean_", store.allFCS[i,2],".Rdata",sep="") )
            index.Clean.Error <- c(index.Clean.Error,i)
            error.Message <-c(error.Message, f.clean[1])
            next
            #return(list(frame=f,inds= to.be.removed))
          }
          rg <- rectangleGate(filterId="gvb", list("GoodVsBad"=c(0, 9999)))
          idx <- filter(f.clean, rg)
          if (length(which(idx@subSet==F))>nrow(f)*.08){
            print( warning("More than 8% of cells would be removed: flowClean is ignored, check the data"))
            removed.cell <- c(removed.cell,identifier(f))
            idx@subSet=T
          }
          if(length(which(idx@subSet==F))!=0)
          {
            to.be.removed <- which(idx@subSet==F)
            exprs(f.clean) <- exprs(f.clean)[-which(idx@subSet==F),]
            f <- f.clean
            clean.ind <- which (f@parameters@data[,2]=="GoodVsBad")
            f<- f[,-clean.ind]
          }else{
            f <- f
          }

          print("flowClean Completed")

        save ( f , file =  paste("Results/After_Clean/", store.allFCS[i,1],"/AftClean_", store.allFCS[i,2],".Rdata",sep="") )
        cat("Time is: ",TimeOutput(start2),"\n",sep="")

        # save ( f.beforeTransform , file =  paste("Results/Before_Transform/", store.allFCS[i,1],"/BfrTrans_", store.allFCS[i,2],".Rdata",sep="") )
        # cat("Time is: ",TimeOutput(start2),"\n",sep="")
        
}
#remove(f.clean)
index.Clean.Error <-index.Clean.Error[index.Clean.Error !=0]
error.Message <- error.Message[!is.na(error.Message)]
if(length(index.Clean.Error) != 0 && length(error.Message)!=0){
  store.Clean.Error <- cbind(index.Clean.Error, store.allFCS[index.Clean.Error, 2], error.Message)
}

save ( store.Clean.Error , file =  paste("Results/store.Clean.Error.Rdata",sep="") )
cat("Total time is: ",TimeOutput(start),"\n",sep="")




# for(i in 1:length(uniqueGT)){
#   pathFrom <- paste0("/code/Projects/3i/Panel_T-cell/August_Dataset_Results/Results/After_Clean/", uniqueGT[1])
#   allFrom <- dir(pathFrom, full.names=T, recursive=T) 
#   pathTo <- paste0("/code/Projects/3i/Panel_T-cell/Results/After_Clean/", uniqueGT[1])
#   allTo <- dir(pathTo, full.names=T, recursive=T) 
#   if (!file.exists(paste(FCS.Groups, "/", Genotype2[temp], "/", store.allMice2[x,3], sep="")) ) {
# }