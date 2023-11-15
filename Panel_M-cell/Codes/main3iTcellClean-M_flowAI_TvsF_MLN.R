# Originally written by Justin Meskas 
# Modified & Updated by Albina Rahim
# Date: April 28, 2016

remove(list=ls())
#setwd("/code/Projects/3i/Panel_BM-cell/")
setwd("/data/Panel_M-cell")


load("Results_MLN/lgl657.Rdata")
load("Results_MLN/Genotype.Rdata")
load("Results_MLN/uniqueGT.Rdata")
load("Results_MLN/channels.ind.Rdata")
load("Results_MLN/store.allFCS.Rdata")

#lgl <- lgl657

library("flowCore")
library("flowBin")
library("flowDensity")
# library("flowType")
#library("flowClean")
# library("snowfall")
library("foreach")
library('arrayMvout')

source("Codes/3iTcellfunctions-M.R")
source("Codes/timeVsFluorescence.R")

save.dir_TvsF <- "Results_MLN/Figures/Clean/TvsF"

suppressWarnings ( dir.create ( "Results_MLN/After_Clean") )
suppressWarnings ( dir.create ( "Results_MLN/Before_Transform") )
suppressWarnings ( dir.create ( "Results_MLN/FlowType") )
suppressWarnings ( dir.create ( "Results_MLN/Figures") )
suppressWarnings ( dir.create ( "Results_MLN/Figures/Clean") )
suppressWarnings ( dir.create ( "Results_MLN/Figures/ScatterPlots") )
suppressWarnings ( dir.create ( save.dir_TvsF ) )

verbose_debris <- F
verbose_margin <- F
verbose_flowClean <- F


start <- Sys.time()

# path to the FCS_Groups files
pathFCS <- paste("/data/Panel_M-cell/FCS_Groups_MLN")

# Reads all folders and files in current path folder and makes a list of all of their paths
allFCS <- dir(pathFCS, full.names=T, recursive=T)

invisible(sapply(1:length(uniqueGT), function(x){
  suppressWarnings(dir.create ( paste("Results_MLN/FlowType/", uniqueGT[x], sep="") ))
  suppressWarnings(dir.create ( paste("Results_MLN/FlowType_LineagePositive/", uniqueGT[x], sep="") ))
  suppressWarnings(dir.create ( paste("Results_MLN/FlowType_LineageNegative/", uniqueGT[x], sep="") ))
  suppressWarnings ( dir.create ( paste("Results_MLN/After_Clean/", uniqueGT[x], sep="") ) )
  suppressWarnings ( dir.create ( paste("Results_MLN/Before_Transform/", uniqueGT[x], sep="") ) )
  suppressWarnings ( dir.create ( paste("Results_MLN/Figures/ScatterPlots/", uniqueGT[x], sep="") ) )}))

#run flowClean only for 5-6 markers, based on what Herve suggested
#that if there's a clog in one channel, it affects other channels as well
#clean.chan <- as.vector(sample(channels.ind,6))

clean.chan <- c(7:10, 12:14)

library(doMC)
#no_cores <- detectCores() - 1
no_cores <- detectCores() - 2
registerDoMC(no_cores)

start2 <- Sys.time()

results <- foreach (i = 1:length(allFCS)) %dopar% { 
  
  res <- tryCatch({
  
    print(paste(i, ": Starting", store.allFCS[i,1],"/",store.allFCS[i,2], sep = ""))
    
    # read in one flowFrame
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
        print("flowAI failed, using TvsF")
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
        #print("flowAI worked")
        f.clean.rem <- f.cleaned
        # rem.ind <- which(f.cleaned@parameters@data$name == "remove_from_all")
        rem.ind <- which(f.cleaned@parameters@data$name == "remove_from_FR")
        to.be.removed <- which(f.cleaned@exprs[,rem.ind]>=9999)
        if (length(to.be.removed) >= 0.25*nrow(f.cleaned) ){ # use TvsF
          print("flowAI removed too much, using TvsF")
          
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
          
          png(file = paste0("Results_MLN/Figures/Clean/", store.allFCS[i,1],"_",store.allFCS[i,2] , '.png'),
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
      
      save ( f , file =  paste("Results_MLN/After_Clean/", store.allFCS[i,1],"/AftClean_", store.allFCS[i,2],".Rdata",sep="") )
      save ( f.beforeTransform , file =  paste("Results_MLN/Before_Transform/", store.allFCS[i,1],"/BfrTrans_", store.allFCS[i,2],".Rdata",sep="") )
      
    }
    
  },error = function(e) {
    return(e)
  }) # end of tryCatch
  
  if(nrow(f) < 20001){
    res <- c(res, "Less than 20000 cells")
  }
        
}

cat("Time is: ",TimeOutput(start2),"\n",sep="")

unlist(results)
error.files <- which(results == "Less than 20000 cells")
error.files <- store.allFCS[error.files,]
save(error.files, file = "Results_MLN/notEnoughCells.Rdata")

# For the MLN organ I need to fail some files that get through cleaning, but where the time domain data is obviuosly niO
# when I did flowClean, I ran some of these files by Adam Liang and he confirmed that I should fail them.

fail.files <- F

if(fail.files == T){
  
  save(store.allFCS, file = "Results_MLN/store.allFCS.original.Rdata") # Back-up in case I need it
  files.to.fail.barcodes <- c(43315:43322, 43324:43331, 44324:44335, 
                              108024, 54190,  # in addition to some of the above barcodes, these have < 20,000 cells
                              71737, 71739, 71740, 71744:71746, 71748:71752,
                              90935, 93635, 93645, 108024, #not enough live cells
                              95655, 102527, 102536, # the CD45/Lin(cd3,cd19,cd103) distribution has low values for both markers, and only a very small number of cells get recognized as cd45+.
                              88666, 180744, 
                              47020, 47026)
  files.to.fail.idx <- unlist(lapply(files.to.fail.barcodes, function(x){ grep(x, store.allFCS[,2])}))
  if(length(files.to.fail.idx) > 0){
    store.allFCS <- store.allFCS[-files.to.fail.idx,]
  }
  
  # Save store.allFCS so that I don't have to remove again in main3iTcell-M.R and main3iTcellFT-M.R
  save(store.allFCS, file = "Results_MLN/store.allFCS.Rdata")
  
}

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

