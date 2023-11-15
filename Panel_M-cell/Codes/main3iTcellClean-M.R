# Originally written by Justin Meskas 
# Modified & Updated by Albina Rahim
# Date: April 28, 2016

remove(list=ls())
#setwd("/code/Projects/3i/Panel_BM-cell/")
setwd("/data/Panel_M-cell")


load("Results/lgl657.Rdata")
load("Results/Genotype.Rdata")
load("Results/uniqueGT.Rdata")
load("Results/channels.ind.Rdata")
load("Results/store.allFCS.Rdata")

#lgl <- lgl657

library("flowCore")
library("flowBin")
# library("flowDensity")
# library("flowType")
library("flowClean")
# library("snowfall")
# library("foreach")

source("Codes/3iTcellfunctions-M.R")


suppressWarnings ( dir.create ( "Results/After_Clean") )
suppressWarnings ( dir.create ( "Results/Before_Transform") )
suppressWarnings ( dir.create ( "Results/FlowType") )
suppressWarnings ( dir.create ( "Results/Figures") )
suppressWarnings ( dir.create ( "Results/Figures/Clean") )
suppressWarnings ( dir.create ( "Results/Figures/ScatterPlots") )
suppressWarnings ( dir.create ( "Results/Figures/CleaningProblems") )

verbose_debris <- F
verbose_margin <- F
verbose_flowClean <- F

tissue <- "MLN" 

start <- Sys.time()

# path to the FCS_Groups files
pathFCS <- paste("/data/Panel_M-cell/FCS_Groups")

# Reads all folders and files in current path folder and makes a list of all of their paths
allFCS <- dir(pathFCS, full.names=T, recursive=T)

invisible(sapply(1:length(uniqueGT), function(x){
  suppressWarnings(dir.create ( paste("Results/FlowType/", uniqueGT[x], sep="") ))
  suppressWarnings ( dir.create ( paste("Results/After_Clean/", uniqueGT[x], sep="") ) )
  suppressWarnings ( dir.create ( paste("Results/Before_Transform/", uniqueGT[x], sep="") ) )
  suppressWarnings ( dir.create ( paste("Results/Figures/ScatterPlots/", uniqueGT[x], sep="") ) )}))

#run flowClean only for 5-6 markers, based on what Herve suggested
#that if there's a clog in one channel, it affects other channels as well
clean.chan <- as.vector(sample(channels.ind,6))
index.Clean.Error <- 0
error.Message <- NA
error.File.Index <- NULL
flowClean.ignored.idx <- NULL
low.cell.count.idx <- NULL

for (i in 1:length(allFCS)){
  
  possibleError <- tryCatch({
    
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
    
    # Check if have enough cells
    if(nrow(f) < 20000){
      low.cell.count.idx <- c(low.cell.count.idx,i)
      next
    }
    ## Compensation
    print("Starting Compensation")
    if( det(f@description$SPILL)==1 ){
      print("Check the spillover matrix, it's probably an identity matrix!")
      return(NULL)
    } else{
      f <- compensate(f, f@description$SPILL)}
    
    ## Saving the data before Transformation (will be needed for the MFI calculation)
    f.beforeTransform <- f
    
    # ## Transformation
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
      save ( f.beforeTransform , file =  paste("Results/Before_Transform/", store.allFCS[i,1],"/BfrTrans_", store.allFCS[i,2],".Rdata",sep="") )
      index.Clean.Error <- c(index.Clean.Error,i)
      error.Message <-c(error.Message, f.clean[1]) 
      png(file = paste0( "Results/Figures/CleaningProblems/", store.allFCS[i,2], '.png'))
      plotDens(f, channels = c(19,14), main = store.allFCS[i,2])
      dev.off()
      next
      #return(list(frame=f,inds= to.be.removed))
    }
    rg <- rectangleGate(filterId="gvb", list("GoodVsBad"=c(0, 9999)))
    idx <- filter(f.clean, rg)
    if (length(which(idx@subSet==F))>nrow(f)*.08){
      print( warning("More than 8% of cells would be removed: flowClean is ignored, check the data"))
      removed.cell <- c(removed.cell,identifier(f))
      flowClean.ignored.idx <- c(flowClean.ignored.idx, i)
      idx@subSet=T
      png(file = paste0( "Results/Figures/CleaningProblems/", store.allFCS[i,2], ".png"))
      plotDens(f, channels = c(19,14), main = store.allFCS[i,2])
      dev.off()
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
    save ( f.beforeTransform , file =  paste("Results/Before_Transform/", store.allFCS[i,1],"/BfrTrans_", store.allFCS[i,2],".Rdata",sep="") )
    
    cat("Time is: ",TimeOutput(start2),"\n",sep="")
    
  },error = function(err) {
    err$message <- paste0(err$message, " (in FCS file index ", i, ")")
    return(err)
  }) # end of tryCatch
  
  # if there is an error in Cleaning
  if(inherits(possibleError, "error")){
    error.File.Index <- c(error.File.Index, i)
    print(paste0("Error in cleaning: ", possibleError))
  }
        
}
# #remove(f.clean)
# index.Clean.Error <-index.Clean.Error[index.Clean.Error !=0]
# error.Message <- error.Message[!is.na(error.Message)]
# store.Clean.Error <- cbind(index.Clean.Error, error.Message)
# save( store.Clean.Error,    file =  paste("Results/store.Clean.Error.Rdata",sep=""))
# cat("Total time is: ",TimeOutput(start),"\n",sep="")
# 
# ##################################################################
# # Print out summary of files where cleaning failed or was not done
# ##################################################################
# 
# if(length(flowClean.ignored.idx)>0){
#   print("More than 8% of cells would be removed: flowClean is ignored for the following files")
#   print(flowClean.ignored.idx)
#   print(store.allFCS[flowClean.ignored.idx,2])
# }
# save(flowClean.ignored.idx, file = paste("Results/store.Ignore.Clean.Rdata",sep=""))
# if(length(error.File.Index)>0){
#   print('The following files failed cleaning')
#   print(error.File.Index)
#   print(store.allFCS[error.File.Index,2])
# }
# if(length(low.cell.count.idx)>0){
#   print("< 20,000 files after margin event removal. Cleaning not performed.")
#   print(low.cell.count.idx)
#   print(store.allFCS[low.cell.count.idx,2])
# }
# 
##############################################
# Files to fail - ie. remove from store.allFCS
##############################################
if(tissue == "SPLEEN"){

  load("Results/store.allFCS.Rdata")

  # Fail because < 20000 cells after margin event removal
  files.to.fail.barcodes <- c(88659, 88655)
  files.to.fail.idx <- unlist(lapply(files.to.fail.barcodes, function(x){ grep(x, store.allFCS[,2])}))
  if(length(files.to.fail.idx) > 0){
    store.allFCS <- store.allFCS[-files.to.fail.idx,]
  }
  # I confirmed with Adam Laing: These files should be disregarded/failed as one of the lineage antibodies wasn’t added
  #to the panel so this has a major impact on the downstream data.
  files.to.fail.barcodes <- c(43170:43182, 43215:43226, 43228:43231, 44310:44323)
  files.to.fail.idx <- unlist(lapply(files.to.fail.barcodes, function(x){ grep(x, store.allFCS[,2])}))
  if(length(files.to.fail.idx) > 0){
    store.allFCS <- store.allFCS[-files.to.fail.idx,]
  }

  # Save store.allFCS so that I don't have to remove again in main3iTcell-M.R and main3iTcellFT-M.R
  save(store.allFCS, file = "Results/store.allFCS.Rdata")
}else if(tissue == "MLN"){

  save(store.allFCS, file = "Results/store.allFCS.original.Rdata") # Back-up in case I need it

  load("Results/store.allFCS.Rdata")
  # Fail because < 20000 cells after margin event removal
  files.to.fail.barcodes <- c(43326, 43318, 43324, 108024, 43331, 43321, 54190, 43322)
  # These are the files from flowClean.ignored.idx and error.File.Index - I double checked the plots here and it looks like there is a fluidics issue for these files
  files.to.fail.barcodes <- c(files.to.fail.barcodes,108025, 108022, 108027, 110170, 108021, 88678, 73866, 88665, 88666, 71749,87201,71752,108026,106663,57391,71739,71740,71744,108029,93645,90935,93635,90462,52423,087203,103184,7398,43323)
  # Fail these files
  files.to.fail.idx <- unlist(lapply(files.to.fail.barcodes, function(x){ grep(x, store.allFCS[,2])}))
  if(length(files.to.fail.idx) > 0){
    store.allFCS <- store.allFCS[-files.to.fail.idx,]
  }
  # I confirmed with Adam Laing:  These files should be disregarded/failed as one of the lineage antibodies wasn’t added
  #to the panel so this has a major impact on the downstream data.
  files.to.fail.barcodes <- c(43315:43320,43325:43330,44324:44335)
  # I confirmed with Adam Laing:  These files should also be failed because there are not enough cells (e.g. live, Lin- ect). For most, the gating fails.
  files.to.fail.barcodes <- c(files.to.fail.barcodes,95655,102527,102536,83706,96545,108028,108031,93645,90935,93635,108024,87201,90932,101043,77890, 97258, 90926, 85042, 113732)
  # I confirmed with Adam Laing: I should fail these files because there is an obvious fluidics issue
  files.to.fail.barcodes <- c(files.to.fail.barcodes,71750,71751,108744,47020,47026)
  # Some additional similar files I found
  files.to.fail.barcodes <- c(files.to.fail.barcodes,108030, 96529, 71745, 96552)
  files.to.fail.idx <- unlist(lapply(files.to.fail.barcodes, function(x){ grep(x, store.allFCS[,2])}))
  if(length(files.to.fail.idx) > 0){
    store.allFCS <- store.allFCS[-files.to.fail.idx,]
  }

  # Save store.allFCS so that I don't have to remove again in main3iTcell-M.R and main3iTcellFT-M.R
  save(store.allFCS, file = "Results/store.allFCS.Rdata")
}
