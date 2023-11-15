# Originally written by Justin Meskas 
# Modified & Updated by Albina Rahim
# Date: April 28, 2016

remove(list=ls())
#setwd("/code/Projects/3i/Panel_BM-cell/")
setwd("/data/Panel_M-cell")


load("Results_SPLEEN/lgl657.Rdata")
load("Results_SPLEEN/Genotype.Rdata")
load("Results_SPLEEN/uniqueGT.Rdata")
load("Results_SPLEEN/channels.ind.Rdata")
load("Results_SPLEEN/store.allFCS.Rdata")

#lgl <- lgl657

library("flowCore")
library("flowBin")
# library("flowDensity")
# library("flowType")
#library("flowClean")
# library("snowfall")
library("foreach")

source("Codes/3iTcellfunctions-M.R")
source("Codes/timeVsFluorescence.R")

save.dir_TvsF <- "Results_SPLEEN/Figures/Clean/TvsF"

suppressWarnings ( dir.create ( "Results_SPLEEN/After_Clean") )
suppressWarnings ( dir.create ( "Results_SPLEEN/Before_Transform") )
suppressWarnings ( dir.create ( "Results_SPLEEN/FlowType") )
suppressWarnings ( dir.create ( "Results_SPLEEN/FlowType_LineagePositive") )
suppressWarnings ( dir.create ( "Results_SPLEEN/FlowType_LineageNegative") )
suppressWarnings ( dir.create ( "Results_SPLEEN/Figures") )
suppressWarnings ( dir.create ( "Results_SPLEEN/Figures/Clean") )
suppressWarnings ( dir.create ( "Results_SPLEEN/Figures/ScatterPlots") )
suppressWarnings ( dir.create ( save.dir_TvsF ) )

verbose_debris <- F
verbose_margin <- F
verbose_flowClean <- F


start <- Sys.time()

# path to the FCS_Groups files
pathFCS <- paste("/data/Panel_M-cell/FCS_Groups_SPLEEN")

# Reads all folders and files in current path folder and makes a list of all of their paths
allFCS <- dir(pathFCS, full.names=T, recursive=T)

invisible(sapply(1:length(uniqueGT), function(x){
  suppressWarnings(dir.create ( paste("Results_SPLEEN/FlowType/", uniqueGT[x], sep="") ))
  suppressWarnings(dir.create ( paste("Results_SPLEEN/FlowType_LineagePositive/", uniqueGT[x], sep="") ))
  suppressWarnings(dir.create ( paste("Results_SPLEEN/FlowType_LineageNegative/", uniqueGT[x], sep="") ))
  suppressWarnings ( dir.create ( paste("Results_SPLEEN/After_Clean/", uniqueGT[x], sep="") ) )
  suppressWarnings ( dir.create ( paste("Results_SPLEEN/Before_Transform/", uniqueGT[x], sep="") ) )
  suppressWarnings ( dir.create ( paste("Results_SPLEEN/Figures/ScatterPlots_LinearGates/", uniqueGT[x], sep="") ) )
  suppressWarnings ( dir.create ( paste("Results_SPLEEN/Figures/ScatterPlots/", uniqueGT[x], sep="") ) )}))

#run flowClean only for 5-6 markers, based on what Herve suggested
#that if there's a clog in one channel, it affects other channels as well
#clean.chan <- as.vector(sample(channels.ind,6))

clean.chan <- c(7:10, 12:14)

library(doMC)
no_cores <- detectCores() - 1
registerDoMC(no_cores)

start2 <- Sys.time()

results <- foreach (i = 1:length(allFCS)) %dopar% { 
  
  res <- tryCatch({
  
    print(paste(i, ": Starting", store.allFCS[i,1],"/",store.allFCS[i,2], sep = ""))
    
    # read in one flowFrame - !Maybe it'd make sense to double check that store.allFCS[i] really corresponds to allFCS[i]
    f <- read.FCS(allFCS[i])
    
    ## Removing Margin events------------------------------------------------------------------------------------------
    scat.chans <- c(grep (colnames(f),pattern = "FSC*"),grep (colnames(f),pattern = "SSC*"))
    ## Removing margin events in Scatter channels
    f <- removeMargins(f, chans=scat.chans, verbose= verbose_margin)
    #Removing negative values in scatter channels
    f <- removeMargins(f, chans=scat.chans, debris=T, neg=T, verbose= verbose_debris) 
    
    # Check if have enough cells - if not, skip cleaning
    if(nrow(f) > 20000){
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
      
      # Clean--------------------------------------------------------------------------------------------------------
      
      f.cleaned <- try(flowAI:::flow_auto_qc(f, remove_from="FR", second_fractionFR = 0.1))
      if(class(f.cleaned)=="try-error"){ # use TvsF
        cat(paste0("flowAI failed, using TvsF", "\n"))
        f.clean.rem <- timeVsFluorescence(f, 500, 1*10^(-4), save.dir_TvsF, name1="Name 1 for figures", name2="Name 2 for figures", Plt=F)
        to.be.removed <- f.clean.rem$inds
        f.clean.rem <- f.clean.rem$frame
        
        png(file = paste0("Results_SPLEEN/Figures/Clean/TvsF/", store.allFCS[i,1],"_",store.allFCS[i,2] , '.png'),
            width=1800, height=1800)
        par(mfrow=c(4,1),mar=(c(5, 5, 4, 2) + 0.1))
        plotDens(f, channels = c(ncol(f), 12))
        points(exprs(f)[to.be.removed, c(ncol(f), 12)],col=1,pch=".")
        plotDens(f, channels = c(ncol(f), 14))
        points(exprs(f)[to.be.removed, c(ncol(f), 14)],col=1,pch=".")
        plotDens(f.clean.rem, channels = c(ncol(f), 12))
        plotDens(f.clean.rem, channels = c(ncol(f), 14))
        dev.off()
        
      } else {
        #print("flowAI worked")
        f.clean.rem <- f.cleaned
        # rem.ind <- which(f.cleaned@parameters@data$name == "remove_from_all")
        rem.ind <- which(f.cleaned@parameters@data$name == "remove_from_FR")
        to.be.removed <- which(f.cleaned@exprs[,rem.ind]>=9999)
        if (length(to.be.removed) >= 0.25*nrow(f.cleaned) ){ # use TvsF
          cat(paste0("flowAI removed too much, using TvsF"))
          f.clean.rem <- timeVsFluorescence(f, 500, 1*10^(-4), save.dir_TvsF, name1="Name 1 for figures", name2="Name 2 for figures", Plt=F)
          to.be.removed <- f.clean.rem$inds
          f.clean.rem <- f.clean.rem$frame
          
          png(file = paste0("Results_MLN/Figures/Clean/TvsF/", store.allFCS[i,1],"_",store.allFCS[i,2] , '.png'),
              width=1800, height=1800)
          par(mfrow=c(4,1),mar=(c(5, 5, 4, 2) + 0.1))
          plotDens(f, channels = c(ncol(f), 12))
          points(exprs(f)[to.be.removed, c(ncol(f), 12)],col=1,pch=".")
          plotDens(f, channels = c(ncol(f), 14))
          points(exprs(f)[to.be.removed, c(ncol(f), 14)],col=1,pch=".")
          plotDens(f.clean.rem, channels = c(ncol(f), 12))
          plotDens(f.clean.rem, channels = c(ncol(f), 14))
          dev.off()
          
        } else {
          
          png(file = paste0("Results_SPLEEN/Figures/Clean/", store.allFCS[i,1],"_",store.allFCS[i,2] , '.png'),
              width=1800, height=1800)
          par(mfrow=c(4,1),mar=(c(5, 5, 4, 2) + 0.1))

          if( length(which(f.cleaned@exprs[,rem.ind]>=9999)) >= 1 ){
            f.clean.rem@exprs<-f.clean.rem@exprs[-to.be.removed,]
            plotDens(f, channels = c(ncol(f), 12))
            points(exprs(f)[to.be.removed, c(ncol(f), 12)],col=1,pch=".")
            plotDens(f, channels = c(ncol(f), 14))
            points(exprs(f)[to.be.removed, c(ncol(f), 14)],col=1,pch=".")
            plotDens(f.clean.rem, channels = c(ncol(f), 12))
            plotDens(f.clean.rem, channels = c(ncol(f), 14))
          }else{
            plotDens(f, channels = c(ncol(f), 12))
            plotDens(f, channels = c(ncol(f), 14))
          }
          
          dev.off()
        }
      }
      
      f <- f.clean.rem
      
      save ( f , file =  paste("Results_SPLEEN/After_Clean/", store.allFCS[i,1],"/AftClean_", store.allFCS[i,2],".Rdata",sep="") )
      save ( f.beforeTransform , file =  paste("Results_SPLEEN/Before_Transform/", store.allFCS[i,1],"/BfrTrans_", store.allFCS[i,2],".Rdata",sep="") )
      
    }
    
  },error = function(e) {
    return(e)
  }) # end of tryCatch
  
  if(exists("f")){
    if(nrow(f) < 20001){
      res <- c(res, paste0(store.allFCS[i,2], " has less than 20000 cells"))
    }
  }
        
}

cat("Time is: ",TimeOutput(start2),"\n",sep="")

results <-unlist(results)
#save(error.files, file = "Results_SPLEEN/notEnoughCells.Rdata")

fail.files <- F

if(fail.files == T){
  save(store.allFCS, file = "Results_SPLEEN/store.allFCS.original.Rdata") # Back-up in case I need it
  
  # I confirmed with Adam Laing: These files should be disregarded/failed as one of the lineage antibodies wasn’t added
  #to the panel so this has a major impact on the downstream data.
  files.to.fail.barcodes <- c(43170:43182, 43215:43226, 43228:43231, 44310:44323)
  # Add the files which had < 20,000 cells
  files.to.fail.barcodes <- c(files.to.fail.barcodes, 88659, 88655)
  files.to.fail.idx <- unlist(lapply(files.to.fail.barcodes, function(x){ grep(x, store.allFCS[,2])}))
  if(length(files.to.fail.idx) > 0){
    store.allFCS <- store.allFCS[-files.to.fail.idx,]
    #!Maybe it'd make sense to delecte allFCS[files.to.fail.idx]?
  }
  
  # Save store.allFCS so that I don't have to remove again in main3iTcell-M.R and main3iTcellFT-M.R
  save(store.allFCS, file = "Results_SPLEEN/store.allFCS.Rdata")
  
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
# ##############################################
# # Files to fail - ie. remove from store.allFCS
# ##############################################
# if(tissue == "SPLEEN"){
# 
#   load("Results/store.allFCS.Rdata")
# 
#   # Fail because < 20000 cells after margin event removal
#   files.to.fail.barcodes <- c(88659, 88655)
#   files.to.fail.idx <- unlist(lapply(files.to.fail.barcodes, function(x){ grep(x, store.allFCS[,2])}))
#   if(length(files.to.fail.idx) > 0){
#     store.allFCS <- store.allFCS[-files.to.fail.idx,]
#   }
#   # I confirmed with Adam Laing: These files should be disregarded/failed as one of the lineage antibodies wasn’t added
#   #to the panel so this has a major impact on the downstream data.
#   files.to.fail.barcodes <- c(43170:43182, 43215:43226, 43228:43231, 44310:44323)
#   files.to.fail.idx <- unlist(lapply(files.to.fail.barcodes, function(x){ grep(x, store.allFCS[,2])}))
#   if(length(files.to.fail.idx) > 0){
#     store.allFCS <- store.allFCS[-files.to.fail.idx,]
#   }
# 
#   # Save store.allFCS so that I don't have to remove again in main3iTcell-M.R and main3iTcellFT-M.R
#   save(store.allFCS, file = "Results/store.allFCS.Rdata")
# }else if(tissue == "MLN"){
# 
#   save(store.allFCS, file = "Results/store.allFCS.original.Rdata") # Back-up in case I need it
# 
#   load("Results/store.allFCS.Rdata")
#   # Fail because < 20000 cells after margin event removal
#   files.to.fail.barcodes <- c(43326, 43318, 43324, 108024, 43331, 43321, 54190, 43322)
#   # These are the files from flowClean.ignored.idx and error.File.Index - I double checked the plots here and it looks like there is a fluidics issue for these files
#   files.to.fail.barcodes <- c(files.to.fail.barcodes,108025, 108022, 108027, 110170, 108021, 88678, 73866, 88665, 88666, 71749,87201,71752,108026,106663,57391,71739,71740,71744,108029,93645,90935,93635,90462,52423,087203,103184,7398,43323)
#   # Fail these files
#   files.to.fail.idx <- unlist(lapply(files.to.fail.barcodes, function(x){ grep(x, store.allFCS[,2])}))
#   if(length(files.to.fail.idx) > 0){
#     store.allFCS <- store.allFCS[-files.to.fail.idx,]
#   }
#   # I confirmed with Adam Laing:  These files should be disregarded/failed as one of the lineage antibodies wasn’t added
#   #to the panel so this has a major impact on the downstream data.
#   files.to.fail.barcodes <- c(43315:43320,43325:43330,44324:44335)
#   # I confirmed with Adam Laing:  These files should also be failed because there are not enough cells (e.g. live, Lin- ect). For most, the gating fails.
#   files.to.fail.barcodes <- c(files.to.fail.barcodes,95655,102527,102536,83706,96545,108028,108031,93645,90935,93635,108024,87201,90932,101043,77890, 97258, 90926, 85042, 113732)
#   # I confirmed with Adam Laing: I should fail these files because there is an obvious fluidics issue
#   files.to.fail.barcodes <- c(files.to.fail.barcodes,71750,71751,108744,47020,47026)
#   # Some additional similar files I found
#   files.to.fail.barcodes <- c(files.to.fail.barcodes,108030, 96529, 71745, 96552)
#   files.to.fail.idx <- unlist(lapply(files.to.fail.barcodes, function(x){ grep(x, store.allFCS[,2])}))
#   if(length(files.to.fail.idx) > 0){
#     store.allFCS <- store.allFCS[-files.to.fail.idx,]
#   }
# 
#   # Save store.allFCS so that I don't have to remove again in main3iTcell-M.R and main3iTcellFT-M.R
#   save(store.allFCS, file = "Results/store.allFCS.Rdata")
# }
