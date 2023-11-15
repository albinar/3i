# Developed by Albina Rahim
# Date: April 13, 2016

remove(list=ls())

Do_flowType <- F #For doing flowType at the end of the code

#setwd("/code/Projects/3i/Panel_T-cell/")
setwd("/data/Panel_M-cell")

library("flowCore")
library("flowBin")
library("flowDensity")
library("flowType")
library("flowClean")
library("stringr")
library('MASS')

source("Codes/3iTcellfunctions-M.R")
source("Codes/rotate.data-M.R")

suppressWarnings ( dir.create ( "Results/FlowType") )
suppressWarnings ( dir.create ( "Results/Figures/ScatterPlots_Updated") )

load("Results/uniqueGT.Rdata")
load("Results/channels.ind.Rdata")
load("Results/store.allFCS.Rdata")
#load("Results/channels.ind.Rdata")

start <- Sys.time()

gating.error.indices <- NULL

# Create directories
invisible(sapply(1:length(uniqueGT), function(x){
  suppressWarnings(dir.create ( paste("Results/FlowType/", uniqueGT[x], sep="") ))
  suppressWarnings(dir.create ( paste("Results/Figures/ScatterPlots_Updated/", uniqueGT[x], sep="") ))}))

all.props.updated <- matrix(nrow = nrow(store.allFCS), ncol = 20)
all.events.updated <- matrix(nrow = nrow(store.allFCS), ncol = 20)
MFI.mean.OriginalData <- matrix(nrow =nrow(store.allFCS), ncol = 342)
MFI.median.OriginalData <- matrix(nrow = nrow(store.allFCS), ncol = 342)
MFI.mean.TransformedData <- matrix(nrow = nrow(store.allFCS), ncol = 342)
MFI.median.TransformedData <- matrix(nrow = nrow(store.allFCS), ncol = 342)

load(file =  paste("Results/ThreshTableUpdated.Rdata",sep="") )  

write.data <- T

loop.idx <- 1:nrow(store.allFCS)
#loop.idx <- c(929)
#loop.idx <- c(1685)

for(i in loop.idx){
  start2 <- Sys.time()
  print(paste(i, ": Starting ", store.allFCS[i,1], " / ", store.allFCS[i,2]), sep="" )
  
  
  load (file =  paste("Results/After_Clean/", store.allFCS[i,1],"/AftClean_", store.allFCS[i,2],".Rdata",sep="") )
  # Reading FCS file before Transform. Will need this for MFI calculation
  load(file =  paste("Results/Before_Transform/", store.allFCS[i,1],"/BfrTrans_", store.allFCS[i,2],".Rdata",sep="") )
  
  possibleError <- tryCatch({
    
    all.events.updated[i,1] <- nrow(f)
    all.props.updated[i,1] <- 1
    
    ## 1st method of plotting FSC-A_SSC-W
    singlets.flowD.h <- flowDensity(f, channels = c("SSC-A", "SSC-W"), position = c(NA, F), gates = c(NA, gthres[i,1]))
    singlets.flowD.l <- flowDensity(singlets.flowD.h, channels = c("SSC-A", "SSC-W"), position = c(NA, T), gates = c(NA, gthres[i,2]))
    singlets.flowD.l@proportion <- (singlets.flowD.l@cell.count/nrow(f))*100
    singlets <- singlets.flowD.l@flow.frame
    
    all.props.updated[i,2] <- singlets.flowD.l@proportion
    all.events.updated[i,2] <- singlets.flowD.l@cell.count
    
    ###############################################################################################
    ## Gating Singlets. Plotting Live/Dead_SSC-A
    live.flowD <- flowDensity(singlets, channels = c(11, 4), position = c(F, NA), gates = c(gthres[i,3], NA))
    live <- live.flowD@flow.frame
    
    all.props.updated[i,3] <- live.flowD@proportion
    all.events.updated[i,3] <- live.flowD@cell.count
    
    ###############################################################################################
    ## Gating Live. Plotting FSC-A_SSC-A
    temp <- flowDensity(live, channels = c(1, 4), position = c(T, F), gates = c(gthres[i,5], gthres[i,6]))
    lymph.flowD <- flowDensity(temp, channels = c(1, 4), position = c(F, F), gates = c(gthres[i,4], gthres[i,6]))
    lymph.flowD@proportion <- (lymph.flowD@cell.count/live.flowD@cell.count)*100
    lymph <- lymph.flowD@flow.frame
    
    all.props.updated[i,4] <- lymph.flowD@proportion
    all.events.updated[i,4] <- lymph.flowD@cell.count
    
    ###############################################################################################
    ## Gating Lymphocytes. Plotting CD45_Lin(CD3, CD109, CD161)
    temp <- flowDensity(lymph, channels = c(14, 12), position = c(T, F), gates = c(gthres[i,7], gthres[i,9]))
    cd45.flowD <- flowDensity(temp, channels = c(14, 12), position = c(F, F), gates = c(gthres[i,8], gthres[i,9]))
    cd45 <- cd45.flowD@flow.frame
    cd45.flowD@proportion <- (cd45.flowD@cell.count/lymph.flowD@cell.count)*100
    
    all.props.updated[i,5] <- cd45.flowD@proportion
    all.events.updated[i,5] <- cd45.flowD@cell.count
    
    ###############################################################################################
    ## Gating CD45 population. Plotting CD45_Lin(CD3, CD109, CD161) to get Lineage negative cells
    
    theta0 = -pi*25/180
    cd45s <- cd45
    cd45s <- rotate.data(cd45, c(14, 12), theta = -theta0)$data
    R <- matrix(c(cos(theta0), sin(theta0), -sin(theta0), cos(theta0)), 2, 2)
    
    cd45sTempb <- flowDensity(cd45s, channels=c(14, 12),position=c(NA, F), gates = c(NA, gthres[i,10]))
    cd45sTempb.ind <- cd45sTempb@index
    lin.neg <- cd45sTempb@flow.frame
    lin.neg@exprs[, c(14, 12)] <- t(t(R) %*% t(lin.neg@exprs[, c(14, 12)]))
    
    cd45sTempb@filter <- rotate.data(cd45sTempb@filter, c(14, 12),theta = theta0)$data
    cd45sTempb@flow.frame <- rotate.data(cd45sTempb@flow.frame, c(14, 12), theta = theta0)$data
    cd45sTempb@proportion <- (cd45sTempb@cell.count/cd45.flowD@cell.count)*100
    lin.neg <- cd45sTempb@flow.frame

    all.props.updated[i,6] <- cd45sTempb@proportion
    all.events.updated[i,6] <- cd45sTempb@cell.count
    
    ###############################################################################################
    ## Gating Lin-. Plotting F4/80_CD11b
    
    theta0 = pi*25/180
    lin.neg.s <- lin.neg
    lin.neg.s <- rotate.data(lin.neg, c(8, 13), theta = -theta0)$data
    R <- matrix(c(cos(theta0), sin(theta0), -sin(theta0), cos(theta0)), 2 , 2)

    linneg.Tempb <- flowDensity(lin.neg.s, channels=c(8, 13),position=c(F, NA), gates = c(gthres[i,11], NA))
    linneg.Tempb.ind <- linneg.Tempb@index
    linneg.macneg <- linneg.Tempb@flow.frame
    linneg.macneg@exprs[, c(8, 13)] <- t(t(R) %*% t(linneg.macneg@exprs[, c(8, 13)]))
    
    linneg.Tempb@filter <- rotate.data(linneg.Tempb@filter, c(8, 13), theta = theta0)$data
    linneg.Tempb@flow.frame <- rotate.data(linneg.Tempb@flow.frame, c(8, 13),theta = theta0)$data
    linneg.Tempb@proportion <- (linneg.Tempb@cell.count/cd45sTempb@cell.count)*100
    linneg.macneg <- linneg.Tempb@flow.frame
    
    all.props.updated[i,7] <- linneg.Tempb@proportion
    all.events.updated[i,7] <- linneg.Tempb@cell.count
    
    # Take complement of the Lin-Mac- cells to obtain RP Macrophages
    # Cell count of RP Macrophages = Cell count of Lin- minus cell count of Lin-Mac- 
    RP.mac.cell.count <- cd45sTempb@cell.count -linneg.Tempb@cell.count
    all.props.updated[i,8] <- RP.mac.cell.count/cd45sTempb@cell.count*100
    all.events.updated[i,8] <- RP.mac.cell.count
    
    # Indices of the RP Mac Population. Required for the MFI calculation.
    RP.mac.ind <- setdiff(cd45sTempb@index,linneg.Tempb@index)
    
#     # Optional Plot: Done to observe if the correct cells are getting saved wrt f
#     plotDens(lin.neg, c(8,13), main='Lineage neg')
#     points(cd45sTempb@flow.frame@exprs[linneg.Tempb@index, c(8,13)], col=2, cex=0.2)
#     points(cd45sTempb@flow.frame@exprs[RP.mac.ind, c(8,13)], col=3, cex=0.2)
    
    ############################################################################################
    ## Gating Lin-Mac-. Plotting Ly6G_SSC-A to obtain Neutrophils2, Ly6G-, Eosinohils3  
    
    theta0 = atan(-150000)
    linneg.macneg.s <- linneg.macneg
    linneg.macneg.s <- rotate.data(linneg.macneg,c(9,4),theta = -theta0)$data
    R <- matrix(c(cos(theta0), sin(theta0), -sin(theta0), cos(theta0)), 2 , 2)

    eosinTempb0 <- flowDensity(linneg.macneg.s, channels =c(9, 4), position = c(NA, T), gates = c(NA, gthres[i,12]))
    eosinTempb <-flowDensity(eosinTempb0, channels =c(9, 4), position = c(T, NA), gates = c(gthres[i,13], NA)) 
    eosinTempb.ind <- eosinTempb@index
    eosin3 <- eosinTempb@flow.frame
    eosin3@exprs[, c(9, 4)] <- t(t(R) %*% t(eosin3@exprs[, c(9, 4)]))
    eosinTempb@filter <- rotate.data(eosinTempb@filter, c(9, 4),theta = theta0)$data
    eosinTempb@flow.frame <- rotate.data(eosinTempb@flow.frame, c(9, 4), theta = theta0)$data
    eosinTempb@proportion <- (eosinTempb@cell.count/linneg.Tempb@cell.count)*100
    eosin3 <- eosinTempb@flow.frame
    
    Ly6GnegTempb <-flowDensity(eosinTempb0, channels =c(9, 4), position = c(F, NA), gates = c(gthres[i,13], NA)) 
    Ly6GnegTempb.ind <- Ly6GnegTempb@index
    Ly6Gneg <- Ly6GnegTempb@flow.frame
    Ly6Gneg@exprs[, c(9, 4)] <- t(t(R) %*% t(Ly6Gneg@exprs[, c(9, 4)]))
    Ly6GnegTempb@filter <- rotate.data(Ly6GnegTempb@filter, c(9, 4),theta = theta0)$data
    Ly6GnegTempb@flow.frame <- rotate.data(Ly6GnegTempb@flow.frame, c(9, 4), theta = theta0)$data
    Ly6GnegTempb@proportion <- (Ly6GnegTempb@cell.count/linneg.Tempb@cell.count)*100
    Ly6Gneg <- Ly6GnegTempb@flow.frame
    
    ncells.neu2 <- linneg.Tempb@cell.count - (Ly6GnegTempb@cell.count + eosinTempb@cell.count)
    if(ncells.neu2 > 1){ 
      neu2Tempb <- flowDensity(linneg.macneg.s, channels =c(9, 4), position = c(NA, F), gates = c(NA, gthres[i,12]))
      neu2Tempb.ind <- neu2Tempb@index
      neu2 <- neu2Tempb@flow.frame
      neu2@exprs[, c(9, 4)] <- t(t(R) %*% t(neu2@exprs[, c(9, 4)]))
      neu2Tempb@filter <- rotate.data(neu2Tempb@filter, c(9, 4),theta = theta0)$data
      neu2Tempb@flow.frame <- rotate.data(neu2Tempb@flow.frame, c(9, 4), theta = theta0)$data
      neu2Tempb@proportion <- (neu2Tempb@cell.count/linneg.Tempb@cell.count)*100
      neu2 <- neu2Tempb@flow.frame
      
      all.events.updated[i,10] <- neu2Tempb@cell.count
      all.props.updated[i,10] <- neu2Tempb@cell.count/linneg.Tempb@cell.count*100
    }else{
      all.events.updated[i,10] <- ncells.neu2
      all.props.updated[i,10] <- ncells.neu2/linneg.Tempb@cell.count*100
    }

    all.props.updated[i,9] <- Ly6GnegTempb@cell.count/linneg.Tempb@cell.count*100
    all.props.updated[i,11] <- eosinTempb@cell.count/linneg.Tempb@cell.count*100
    all.events.updated[i,9] <- Ly6GnegTempb@cell.count 
    all.events.updated[i,11] <- eosinTempb@cell.count
    
#     # Optional Plot: Done to observe if the correct cells are geing saved wrt f
#     plotDens(linneg.macneg, c(9,4), main='lin-mac-')
#     points(cd45sTempb@flow.frame@exprs[neu2Tempb@index, c(9,4)], col=2, cex=0.2)
#     points(cd45sTempb@flow.frame@exprs[eosinTempb@index, c(9,4)], col=3, cex=0.2)
#     points(cd45sTempb@flow.frame@exprs[Ly6GnegTempb@index, c(9,4)], col=4, cex=0.2)
     
    ############################################################################################
    ## Gating Ly6G-. Plotting Ly6C_CD11b to obtain Monocytes Ly6c hi and NOT(Monocytes Ly6c hi)
    
    mono.Ly6c.hi.flowD <- flowDensity(Ly6Gneg, channels = c(10, 13), position = c(T, T), gates = c(gthres[i,14], gthres[i,15]))
    mono.Ly6c.hi <- mono.Ly6c.hi.flowD@flow.frame
    
    all.props.updated[i,12] <- mono.Ly6c.hi.flowD@proportion
    all.events.updated[i,12] <- mono.Ly6c.hi.flowD@cell.count
    
    # Grab the NOT(Monocytes Ly6c hi) using notSubFrame() function
    not.mono.Ly6c.hi.flowD  <- notSubFrame(Ly6Gneg, channels=c(10, 13), position= "logical", gates = "missing", mono.Ly6c.hi.flowD@filter)
    not.mono.Ly6c.hi <- not.mono.Ly6c.hi.flowD@flow.frame
    
    all.props.updated[i,13] <- not.mono.Ly6c.hi.flowD@cell.count/Ly6GnegTempb@cell.count*100
    all.events.updated[i,13] <- not.mono.Ly6c.hi.flowD@cell.count
    
    ############################################################################################
    ## Gating NOT(Monocytes Ly6c hi). Plotting Ly6C_CD3
    
    pDC.flowD <- flowDensity(not.mono.Ly6c.hi, channels = c(10, 15), position = c(T, T), gates = c(gthres[i,16], gthres[i,17]))
    pDC <- pDC.flowD@flow.frame
    
    all.props.updated[i,14] <- pDC.flowD@proportion
    all.events.updated[i,14] <- pDC.flowD@cell.count
    
    # Grab NOT(pDC) using notSubFrame() function
    if(pDC.flowD@cell.count >0){
      not.pDC.flowD  <- notSubFrame(not.mono.Ly6c.hi, channels=c(10, 15), position= "logical", gates = "missing", pDC.flowD@filter)
      not.pDC <- not.pDC.flowD@flow.frame
    }else{
      not.pDC.flowD <- not.mono.Ly6c.hi.flowD
      not.pDC <- not.mono.Ly6c.hi
    }
    
    all.props.updated[i,15] <-  not.pDC.flowD@cell.count/not.mono.Ly6c.hi.flowD@cell.count*100
    all.events.updated[i,15] <- not.pDC.flowD@cell.count
    
#     # Optional Plot: Done to observe if the correct cells are geing saved wrt f
#     plotDens(not.mono.Ly6c.hi, c(10,15), main='NOT(monocytes Ly6c hi')
#     points(not.mono.Ly6c.hi.flowD@flow.frame@exprs[not.pDC.flowD@index, c(10,15)], col=2, cex=0.2)
    
    ############################################################################################ 
    ## Gating NOT(pDC). Plotting Ly6C_CD317
    
    cDC.flowD <- flowDensity(not.pDC, channels = c(16, 7), position = c(T, T), gates = c(gthres[i,19],  gthres[i,18]))
    cDC <- cDC.flowD@flow.frame
    
    all.props.updated[i,16] <- cDC.flowD@proportion
    all.events.updated[i,16] <- cDC.flowD@cell.count
    
#     # Take complement of the cDC data to obtain Misc
#     misc <- not.pDC
#     misc@exprs <- misc@exprs[-cDC.flowD@index, ]
    
    # Grab misc using notSubFrame() function
    misc.flowD  <- notSubFrame(not.pDC, channels=c(16, 7), position= "logical", gates = "missing", cDC.flowD@filter)
    
    all.props.updated[i,17] <- misc.flowD@cell.count/not.pDC.flowD@cell.count*100
    all.events.updated[i,17] <- misc.flowD@cell.count
    
    ############################################################################################
    ## Gating cDC. Plotting CD11b_CD86 
    
    theta0 = -pi*45/180
    cDC.s <- cDC
    cDC.s <- rotate.data(cDC, c(13, 18), theta = -theta0)$data
    R <- matrix(c(cos(theta0), sin(theta0), -sin(theta0), cos(theta0)), 2 , 2)
    
    returnTry <- tryCatch({
      cd8a.typedc.Tempb <- flowDensity(cDC.s, channels=c(13, 18),position=c(NA, T), gates = c(NA, gthres[i,20]))
      cd8a.typedc.Tempb.ind <- cd8a.typedc.Tempb@index
      cd8a.typedc <- cd8a.typedc.Tempb@flow.frame
      cd8a.typedc@exprs[, c(13, 18)] <- t(t(R) %*% t(cd8a.typedc@exprs[, c(13, 18)]))
      
      cd8a.typedc.Tempb@filter <- rotate.data(cd8a.typedc.Tempb@filter, c(13, 18), theta = theta0)$data
      cd8a.typedc.Tempb@flow.frame <- rotate.data(cd8a.typedc.Tempb@flow.frame, c(13, 18), theta = theta0)$data
      cd8a.typedc.Tempb@proportion <- (cd8a.typedc.Tempb@cell.count/nrow(cDC))*100
      cd8a.typedc <-cd8a.typedc.Tempb@flow.frame
    },error = function(err) {
      print("Error with Normal"); return(0)
    })
    
    all.props.updated[i,18] <- cd8a.typedc.Tempb@proportion
    all.events.updated[i,18] <- cd8a.typedc.Tempb@cell.count
    
#     # Take complement of the  CD8A Type DC to obtain CD11B+ CD86Lo cells
#     cd11b.plus.cd86.lo <- cDC
#     cd11b.plus.cd86.lo@exprs <- cd11b.plus.cd86.lo@exprs[-cd8a.typedc.Tempb@index, ]
    
    # Grab CD11B+ CD86Lo cells using notSubFrame() function
    cd11b.plus.cd86.lo.flowD <- notSubFrame(cDC, channels=c(13, 18), position= "logical", gates = "missing", cd8a.typedc.Tempb@filter)
    
    all.props.updated[i,19] <- cd11b.plus.cd86.lo.flowD@cell.count/cDC.flowD@cell.count*100
    all.events.updated[i,19] <- cd11b.plus.cd86.lo.flowD@cell.count
    
    ############################################################################################
    ## Gating CD8A Type DC. Plotting CD103_CD11b to obtain CD103+ DC
    
    cd103plus.dc.flowD <- flowDensity(cd8a.typedc, channels = c(17, 13), position = c(T, NA), gates = c(gthres[i,21], NA))
    cd103plus.dc <- cd103plus.dc.flowD@flow.frame
    
    all.props.updated[i,20] <- cd103plus.dc.flowD@proportion
    all.events.updated[i,20] <- cd103plus.dc.flowD@cell.count
    
    ####################################################################################################

    # Saving the MFIs (mean) of Original Data (Data before Transformed)
    singlets.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[singlets.flowD.l@index, x], na.rm = T)})
    live.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[live.flowD@index, x], na.rm = T)})
    lymph.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[lymph.flowD@index, x], na.rm = T)})
    cd45.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd45.flowD@index, x], na.rm = T)})
    linneg.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd45sTempb@index, x], na.rm = T)})
    linneg.macneg.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[linneg.Tempb@index, x], na.rm = T)})
    RPmac.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[RP.mac.ind, x], na.rm = T)})
    Ly6Gneg.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[Ly6GnegTempb@index, x], na.rm = T)})
    neu2.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[neu2Tempb@index, x], na.rm = T)})
    eosin3.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[eosinTempb@index, x], na.rm = T)})
    mono.Ly6c.hi.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[mono.Ly6c.hi.flowD@index, x], na.rm = T)})
    not.mono.Ly6c.hi.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[not.mono.Ly6c.hi.flowD@index, x], na.rm = T)})
    pDC.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[pDC.flowD@index, x], na.rm = T)})
    not.pDC.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[not.pDC.flowD@index, x], na.rm = T)})
    cDC.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cDC.flowD@index, x], na.rm = T)})
    misc.MFI <-  sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[misc.flowD@index, x], na.rm = T)})
    CD8A.TypeDC.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd8a.typedc.Tempb@index, x], na.rm = T)})
    CD11bplus.CD86lo.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd11b.plus.cd86.lo.flowD@index, x], na.rm = T)}) 
    CD103plusDC.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd103plus.dc.flowD@index, x], na.rm = T)}) 
    
    MFI.mean.OriginalData[i, ] <- c(singlets.MFI, live.MFI, lymph.MFI, cd45.MFI, linneg.MFI, linneg.macneg.MFI, RPmac.MFI, Ly6Gneg.MFI, neu2.MFI,
                                    eosin3.MFI, mono.Ly6c.hi.MFI, not.mono.Ly6c.hi.MFI, pDC.MFI, not.pDC.MFI, cDC.MFI, misc.MFI, CD8A.TypeDC.MFI, 
                                    CD11bplus.CD86lo.MFI, CD103plusDC.MFI)

    # Saving the MFIs (median) of Original Data (Data before Transformed)
    singlets.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[singlets.flowD.l@index, x], na.rm = T)})
    live.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[live.flowD@index, x], na.rm = T)})
    lymph.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[lymph.flowD@index, x], na.rm = T)})
    cd45.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd45.flowD@index, x], na.rm = T)})
    linneg.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd45sTempb@index, x], na.rm = T)})
    linneg.macneg.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[linneg.Tempb@index, x], na.rm = T)})
    RPmac.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[RP.mac.ind, x], na.rm = T)})
    Ly6Gneg.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[Ly6GnegTempb@index, x], na.rm = T)})
    neu2.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[neu2Tempb@index, x], na.rm = T)})
    eosin3.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[eosinTempb@index, x], na.rm = T)})
    mono.Ly6c.hi.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[mono.Ly6c.hi.flowD@index, x], na.rm = T)})
    not.mono.Ly6c.hi.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[not.mono.Ly6c.hi.flowD@index, x], na.rm = T)})
    pDC.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[pDC.flowD@index, x], na.rm = T)})
    not.pDC.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[not.pDC.flowD@index, x], na.rm = T)})
    cDC.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cDC.flowD@index, x], na.rm = T)})
    misc.MFI <-  sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[misc.flowD@index, x], na.rm = T)})
    CD8A.TypeDC.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd8a.typedc.Tempb@index, x], na.rm = T)})
    CD11bplus.CD86lo.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd11b.plus.cd86.lo.flowD@index, x], na.rm = T)}) 
    CD103plusDC.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd103plus.dc.flowD@index, x], na.rm = T)}) 
    
    MFI.median.OriginalData[i, ] <- c(singlets.MFI, live.MFI, lymph.MFI, cd45.MFI, linneg.MFI, linneg.macneg.MFI, RPmac.MFI, Ly6Gneg.MFI, neu2.MFI,
                                    eosin3.MFI, mono.Ly6c.hi.MFI, not.mono.Ly6c.hi.MFI, pDC.MFI, not.pDC.MFI, cDC.MFI, misc.MFI, CD8A.TypeDC.MFI, 
                                    CD11bplus.CD86lo.MFI, CD103plusDC.MFI)   
    
    # Saving the MFIs (mean) of Transformed (Data after Transformation)
    singlets.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[singlets.flowD.l@index, x], na.rm = T)})
    live.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[live.flowD@index, x], na.rm = T)})
    lymph.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[lymph.flowD@index, x], na.rm = T)})
    cd45.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd45.flowD@index, x], na.rm = T)})
    linneg.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd45sTempb@index, x], na.rm = T)})
    linneg.macneg.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[linneg.Tempb@index, x], na.rm = T)})
    RPmac.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[RP.mac.ind, x], na.rm = T)})
    Ly6Gneg.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[Ly6GnegTempb@index, x], na.rm = T)})
    neu2.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[neu2Tempb@index, x], na.rm = T)})
    eosin3.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[eosinTempb@index, x], na.rm = T)})
    mono.Ly6c.hi.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[mono.Ly6c.hi.flowD@index, x], na.rm = T)})
    not.mono.Ly6c.hi.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[not.mono.Ly6c.hi.flowD@index, x], na.rm = T)})
    pDC.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[pDC.flowD@index, x], na.rm = T)})
    not.pDC.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[not.pDC.flowD@index, x], na.rm = T)})
    cDC.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cDC.flowD@index, x], na.rm = T)})
    misc.MFI <-  sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[misc.flowD@index, x], na.rm = T)})
    CD8A.TypeDC.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd8a.typedc.Tempb@index, x], na.rm = T)})
    CD11bplus.CD86lo.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd11b.plus.cd86.lo.flowD@index, x], na.rm = T)}) 
    CD103plusDC.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd103plus.dc.flowD@index, x], na.rm = T)}) 
    
    MFI.mean.TransformedData[i, ] <- c(singlets.MFI, live.MFI, lymph.MFI, cd45.MFI, linneg.MFI, linneg.macneg.MFI, RPmac.MFI, Ly6Gneg.MFI, neu2.MFI,
                                    eosin3.MFI, mono.Ly6c.hi.MFI, not.mono.Ly6c.hi.MFI, pDC.MFI, not.pDC.MFI, cDC.MFI, misc.MFI, CD8A.TypeDC.MFI, 
                                    CD11bplus.CD86lo.MFI, CD103plusDC.MFI)
    
    # Saving the MFIs (median) of Transformed (Data after Transformation)
    singlets.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[singlets.flowD.l@index, x], na.rm = T)})
    live.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[live.flowD@index, x], na.rm = T)})
    lymph.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[lymph.flowD@index, x], na.rm = T)})
    cd45.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd45.flowD@index, x], na.rm = T)})
    linneg.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd45sTempb@index, x], na.rm = T)})
    linneg.macneg.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[linneg.Tempb@index, x], na.rm = T)})
    RPmac.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[RP.mac.ind, x], na.rm = T)})
    Ly6Gneg.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[Ly6GnegTempb@index, x], na.rm = T)})
    neu2.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[neu2Tempb@index, x], na.rm = T)})
    eosin3.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[eosinTempb@index, x], na.rm = T)})
    mono.Ly6c.hi.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[mono.Ly6c.hi.flowD@index, x], na.rm = T)})
    not.mono.Ly6c.hi.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[not.mono.Ly6c.hi.flowD@index, x], na.rm = T)})
    pDC.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[pDC.flowD@index, x], na.rm = T)})
    not.pDC.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[not.pDC.flowD@index, x], na.rm = T)})
    cDC.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cDC.flowD@index, x], na.rm = T)})
    misc.MFI <-  sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[misc.flowD@index, x], na.rm = T)})
    CD8A.TypeDC.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd8a.typedc.Tempb@index, x], na.rm = T)})
    CD11bplus.CD86lo.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd11b.plus.cd86.lo.flowD@index, x], na.rm = T)}) 
    CD103plusDC.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd103plus.dc.flowD@index, x], na.rm = T)}) 
    
    MFI.median.TransformedData[i, ] <- c(singlets.MFI, live.MFI, lymph.MFI, cd45.MFI, linneg.MFI, linneg.macneg.MFI, RPmac.MFI, Ly6Gneg.MFI, neu2.MFI,
                                      eosin3.MFI, mono.Ly6c.hi.MFI, not.mono.Ly6c.hi.MFI, pDC.MFI, not.pDC.MFI, cDC.MFI, misc.MFI, CD8A.TypeDC.MFI, 
                                      CD11bplus.CD86lo.MFI, CD103plusDC.MFI)   
    
    ####################################################################################################
    
  },error = function(err) {
    err
    })
  
  if(inherits(possibleError, "error")){
    gating.error.indices <- c(gating.error.indices, i)
    print(paste0("Error in gating: ", possibleError))
    #next
  }
  
  
  tryCatch({
    #--------Start Big Png------------
    png ( file = paste("Results/Figures/ScatterPlots_Updated/", store.allFCS[i,1], "/", "Total_", store.allFCS[i,2], ".png", sep = "" ), width=2500, height=2500*3/5)
    par(mfrow=c(3,5),mar=(c(5, 5, 4, 2) + 0.1))
    
    plotDens(f,c("SSC-A","SSC-W"), main="Ungated", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(h = gthres[i,2], lwd=1)
    abline(h = gthres[i,1], lwd=1)
    
    # Plotting Live/Dead_SSC-A 
    plotDens(singlets,  c(11,4), main= "Singlets", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v=gthres[i,3], lwd=1); 
    
    # Plotting FSC-A_SSC-A
    plotDens(live, c(1,4), main = "Live", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v=gthres[i,5], lwd=1); abline(v=gthres[i,4], lwd=1)
    abline(h=gthres[i,6], lwd=1)
    
    # Plotting CD45_Lin(CD3,CD19,CD161)
    plotDens(lymph, channels = c(14,12), main='Lymphocytes', cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v=gthres[i,7], lwd=1); abline(v=gthres[i,8], lwd=1)
    abline(h=gthres[i,9], lwd=1)
    
    plotDens(cd45, channels = c(14,12), main="CD45", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(cd45sTempb@filter, lwd=1)
    
    # Plotting F4/80_CD11b
    plotDens(lin.neg, channels = c(8,13), main="Lineage neg", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(linneg.Tempb@filter,lwd=1)
    
    plotDens(linneg.macneg, channels = c(9,4), main="Lin-Mac-", cex.lab = 2, cex.axis = 2, cex.main=2)
    if(ncells.neu2 >1){
      lines(neu2Tempb@filter)
    }
    lines(Ly6GnegTempb@filter)
    lines(eosinTempb@filter)
    
    # Plotting Ly6C_CD11b
    plotDens(Ly6Gneg, channels = c(10,13), main="Ly6G-", cex.lab = 2, cex.axis = 2, cex.main=2)
    segments(gthres[i,14], gthres[i,15], gthres[i,14], par('usr')[4], lwd=1)
    segments(gthres[i,14], gthres[i,15], par('usr')[2], gthres[i,15], lwd=1)
    
    # Plotting Ly6C_CD317 
    plotDens(not.mono.Ly6c.hi, channels = c(10,15), main="NOT(Monocytes Ly6c hi)", cex.lab = 2, cex.axis = 2, cex.main=2)
    #lines(pDC.flowD@filter,lwd=1)
    segments(gthres[i,16], gthres[i,17], gthres[i,16], par('usr')[4], lwd=1)
    segments(gthres[i,16], gthres[i,17], par('usr')[2], gthres[i,17], lwd=1)
    
    # Plotting CD11b_MHCII
    plotDens(not.pDC, channels = c(16, 7), main="NOT pDC", cex.lab = 2, cex.axis = 2, cex.main=2)
    segments(gthres[i,19], gthres[i,18], gthres[i,19], par('usr')[4], lwd=1)
    segments(gthres[i,19], gthres[i,18], par('usr')[2], gthres[i,18], lwd=1)
    
#     # Plotting CD11b_CD86
#     plotDens(cDC, channels = c(13, 18), main="cDC", cex.lab = 2, cex.axis = 2, cex.main=2)
#     lines(cd8a.typedc.Tempb@filter, lwd=1)
#     
    # Plotting CD11b_CD86 with contour
    plotDens(cDC, channels = c(13, 18), main="cDC", cex.lab = 2, cex.axis = 2, cex.main=2)
    data.new <- exprs(cDC)[,c(13,18)]
    data.new <- data.new[complete.cases(data.new),]
    z <- kde2d(data.new[, 1], data.new[, 2])
    contour(z, drawlabels = FALSE, add = TRUE, lty = 2)
    lines(cd8a.typedc.Tempb@filter, lwd=1)
    
#     # Plotting CD103_CD11b
#     plotDens(cd8a.typedc, channels = c(17, 13), main="CD8A Type DC", cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(-1,4), ylim = c(-1,4))
#     abline(v=gthres[i,21], lwd=1)
#     
    # Plotting CD103_CD11b with contour
    plotDens(cd8a.typedc, channels = c(17, 13), main="CD8A Type DC", cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(-1,4), ylim = c(-1,4))
    data.new <- exprs(cd8a.typedc)[,c(17,13)]
    data.new <- data.new[complete.cases(data.new),]
    z <- kde2d(data.new[, 1], data.new[, 2])
    contour(z, drawlabels = FALSE, add = TRUE, lty = 2)
    abline(v=gthres[i,21], lwd=1)
    
    dev.off()
    #par(mfrow=c(1,1))
  },error = function(err) {
    print("Error in plotting")
    })
  
  
  if (Do_flowType == T) {
    print("Start FlowType")
    
    theta0 = -pi*45/180
    linneg.macneg.ss <- linneg.macneg.s
    linneg.macneg.ss  <- rotate.data(linneg.macneg.s, c(13, 18), theta = -theta0)$data
    
    listGthres <- list(gthres[i,'linneg.macneg.gate.slanta'], gthres[i,'linneg.macneg.gate.slantb'], cd11b.gate.s,  
                    min(gthres[i,'Ly6c.gate'],gthres[i,'Ly6c.gate2']), cd317.gate2, gthres[i,'mhcII.gate'],
                    gthres[i,'cd11c.gate'], gthres[i,'cDC.gate.slant'], gthres[i,'cd103.gate'])
    flowType.res <- flowType(Frame=linneg.macneg.ss, PropMarkers= c(4,9,13,10,15,7,16,18,17), MaxMarkersPerPop = NULL, PartitionsPerMarker=2,
                             Methods='Thresholds', Thresholds= listGthres, verbose=F, MemLimit=400)
    
    save ( flowType.res, file =  paste("Results/FlowType/", store.allFCS[i,1],"/FT_", store.allFCS[i,2],".Rdata",sep="") )
    
  }
  
  cat("Time is: ",TimeOutput(start2),"\n",sep="")
}

if(write.data){
  
  save(gating.error.indices, file = paste("Results/store.Gate.Error.Updated.Rdata",sep=""))
  colnames(all.props.updated) <- c("All Events-%Parent", "Singlets-%Parent", "Live-%Parent", "Lymphocytes-%Parent", "CD45-%Parent", "Lineage neg-%Parent", 
                                   "Lin-Mac--%Parent", "RP macrophages-%Parent", "Ly6G--%Parent", "Neutrophils 2-%Parent", "Eosinophils 3-%Parent", 
                                   "Monocytes Ly6c hi-%Parent","NOT(Monocytes Ly6c hi)-%Parent", "pDC-%Parent", "NOT pDC-%Parent", "cDC-%Parent", 
                                   "Misc-%Parent", "CD8A Type DC-%Parent","CD11B+ CD86Lo-%Parent", "CD103+ DC-%Parent")
  rownames(all.props.updated) <- 1:length(all.props.updated[,1])
  
  colnames(all.events.updated) <- c("All Events-#Events", "Singlets-#Events", "Live-#Events", "Lymphocytes-#Events", "CD45-#Events", "Lineage neg-#Events", 
                                    "Lin-Mac--#Events", "RP macrophages-#Events", "Ly6G--#Events", "Neutrophils 2-#Events", "Eosinophils 3-#Events", 
                                    "Monocytes Ly6c hi-#Events","NOT(Monocytes Ly6c hi)-#Events", "pDC-#Events", "NOT pDC-#Events", "cDC-%#Events", 
                                    "Misc-#Events", "CD8A Type DC-#Events","CD11B+ CD86Lo-#Events", "CD103+ DC-#Events")
  rownames(all.events.updated) <- 1:length(all.props.updated[,1])
  
  save(gating.error.indices, file = paste("Results/store.Gate.Error.Updated.Rdata",sep=""))
  
  Events_Proportions_Table <- NULL
  for(i in 1:ncol(all.props.updated)){
    Events_Proportions_Table <- cbind(Events_Proportions_Table,all.events.updated[,i],all.props.updated[,i])
  }
  colnames(Events_Proportions_Table) <- c(rbind(colnames(all.events.updated), colnames(all.props.updated)))
  rownames(Events_Proportions_Table) <- 1:length(Events_Proportions_Table[,1])
  
  # # Albina constructed CSVfile.Rdata by binding rows of the following .csv files in Separate_File_Folders-M.R
  # # /mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2015-08-27/attachments/3i_IMPC_Data_Genotypes.csv
  # # /mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2016-03-08/attachments/SPLN data RB.csv
  # # I diffed these files - most, but not all of the lines in te first file are duplicated in the second,
  # # Thus, CSVfile.Rdata has many duplicated rows. 
  # 
  load("Results/CSVfile.Rdata")
  Barcodes <- str_extract(store.allFCS[, 2],"L[0-9]+")
  Genotype <- str_replace(store.allFCS[, 1], '_', '/')
  Mouse_Label <-  CSVfile[, 'Label.Barcode']   
  # Thus when I get the indices for the barcodes I just arbitrarily pick the largest because as far as I can tell
  # the rows in CSVfile pertaining to one barcode are the same
  indices <- sapply(1:length(Barcodes), function(xdat){current.index <- max(grep(Barcodes[xdat], Mouse_Label))})
  AssayDate <- CSVfile[indices,'Assay.Date']
  Events_Proportions_Table <- cbind(store.allFCS[,1:2], Barcodes, AssayDate, store.allFCS[,3], Events_Proportions_Table)
  colnames(Events_Proportions_Table)[5] <- "Total Number of Cells"
  
  all.props.updated <- cbind(store.allFCS[,1:2], all.props.updated)
  all.events.updated <- cbind(store.allFCS[,1:3], all.events.updated)
  
  date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv") 
  write.csv(Events_Proportions_Table, file =  paste("Results/Events_Proportions_Table_Updated", date.time, sep=""))
  write.csv(all.props.updated, file =  paste("Results/allProportions_Table_Updated", date.time, sep=""))
  write.csv(all.events.updated, file =  paste("Results/allEvents_Table_Updated", date.time, sep=""))
  
  colnames.MFI <- c(sapply(1:(ncol(f@exprs)-1), function(x){paste("Singlets:", (f@parameters@data$desc[[x]]))}),
                    sapply(1:(ncol(f@exprs)-1), function(x){paste("Live:", (f@parameters@data$desc[[x]]))}),
                    sapply(1:(ncol(f@exprs)-1), function(x){paste("Lymphocytes:", (f@parameters@data$desc[[x]]))}),
                    sapply(1:(ncol(f@exprs)-1), function(x){paste("CD45+:", (f@parameters@data$desc[[x]]))}),
                    sapply(1:(ncol(f@exprs)-1), function(x){paste("Lineage neg:", (f@parameters@data$desc[[x]]))}),
                    sapply(1:(ncol(f@exprs)-1), function(x){paste("Lin-Mac-:", (f@parameters@data$desc[[x]]))}),
                    sapply(1:(ncol(f@exprs)-1), function(x){paste("RP Macrophages:", (f@parameters@data$desc[[x]]))}),
                    sapply(1:(ncol(f@exprs)-1), function(x){paste("Ly6G-:", (f@parameters@data$desc[[x]]))}),
                    sapply(1:(ncol(f@exprs)-1), function(x){paste("Neutrophils2:", (f@parameters@data$desc[[x]]))}),
                    sapply(1:(ncol(f@exprs)-1), function(x){paste("Eosinophils3:", (f@parameters@data$desc[[x]]))}),
                    sapply(1:(ncol(f@exprs)-1), function(x){paste("Monocytes Ly6c hi:", (f@parameters@data$desc[[x]]))}),
                    sapply(1:(ncol(f@exprs)-1), function(x){paste("NOT(Monocytes Ly6c hi):", (f@parameters@data$desc[[x]]))}),
                    sapply(1:(ncol(f@exprs)-1), function(x){paste("pDC:", (f@parameters@data$desc[[x]]))}),
                    sapply(1:(ncol(f@exprs)-1), function(x){paste("NOT(pDC):", (f@parameters@data$desc[[x]]))}),
                    sapply(1:(ncol(f@exprs)-1), function(x){paste("cDC:", (f@parameters@data$desc[[x]]))}),
                    sapply(1:(ncol(f@exprs)-1), function(x){paste("misc:", (f@parameters@data$desc[[x]]))}),
                    sapply(1:(ncol(f@exprs)-1), function(x){paste("CD8A Type DC:", (f@parameters@data$desc[[x]]))}),
                    sapply(1:(ncol(f@exprs)-1), function(x){paste("CD11B+ CD86lo:", (f@parameters@data$desc[[x]]))}),
                    sapply(1:(ncol(f@exprs)-1), function(x){paste("CD103+ DC:", (f@parameters@data$desc[[x]]))}))
  
  colnames(MFI.mean.OriginalData) <- colnames.MFI
  colnames(MFI.median.OriginalData) <- colnames.MFI
  colnames(MFI.mean.TransformedData) <- colnames.MFI
  colnames(MFI.median.TransformedData) <- colnames.MFI
  
  MFI.mean.OriginalData <- cbind(store.allFCS[,1:2], Barcodes, AssayDate, store.allFCS[,3], MFI.mean.OriginalData)
  MFI.median.OriginalData <- cbind(store.allFCS[,1:2], Barcodes, AssayDate, store.allFCS[,3], MFI.median.OriginalData)
  MFI.mean.TransformedData <- cbind(store.allFCS[,1:2], Barcodes, AssayDate, store.allFCS[,3], MFI.mean.TransformedData)
  MFI.median.TransformedData <- cbind(store.allFCS[,1:2], Barcodes, AssayDate, store.allFCS[,3], MFI.median.TransformedData)
  
  write.csv(MFI.mean.OriginalData, file =  paste("Results/MFI_mean_OriginalData", date.time, sep=""))
  write.csv(MFI.median.OriginalData, file =  paste("Results/MFI_median_OriginalData", date.time, sep=""))
  write.csv(MFI.mean.TransformedData, file =  paste("Results/MFI_mean_TransformedData", date.time, ep=""))
  write.csv(MFI.median.TransformedData, file =  paste("Results/MFI_median_TransformedData", date.time, sep=""))
}

cat("Total time is: ",TimeOutput(start),"\n",sep="")

