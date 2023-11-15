# Developed by Albina Rahim
# Date: April 13, 2016

remove(list=ls())

Do_flowType <- T #For doing flowType at the end of the code

setwd("/code/Projects/3i/Panel_T-cell_MLN/")

library("flowCore")
library("flowBin")
library("flowDensity")
library("flowType")
library("flowClean")


source("Codes/3iTcellfunctions.R")
source("Codes/rotate.data.R")

suppressWarnings ( dir.create ( "Results/FlowType") )
suppressWarnings ( dir.create ( "Results/Figures/ScatterPlots_Updated") )

load("Results/uniqueGT.Rdata")
load("Results/channels.ind.Rdata")
load("Results/store.allFCS.Rdata")
#load("Results/channels.ind.Rdata")

verbose_debris <- T
verbose_margin <- T

start <- Sys.time()

# path to the FCS_Groups files
pathFCS <- paste("/code/Projects/3i/Panel_T-cell_MLN/FCS_Groups")

# Reads all folders and files in current path folder and makes a list of all of their paths
allFCS <- dir(pathFCS, full.names=T, recursive=T) 

all.props <- matrix(nrow = nrow(store.allFCS), ncol = 45)
all.events <- matrix(nrow = nrow(store.allFCS), ncol = 45)
MFI.mean.OriginalData <- matrix(nrow =nrow(store.allFCS), ncol = 792)
MFI.median.OriginalData <- matrix(nrow = nrow(store.allFCS), ncol = 792)
MFI.mean.TransformedData <- matrix(nrow = nrow(store.allFCS), ncol = 792)
MFI.median.TransformedData <- matrix(nrow = nrow(store.allFCS), ncol = 792)

errorFileIndex <- NULL

# Create directories
invisible(sapply(1:length(uniqueGT), function(x){
  suppressWarnings(dir.create ( paste("Results/FlowType/", uniqueGT[x], sep="") ))
  suppressWarnings(dir.create ( paste("Results/Figures/ScatterPlots_Updated/", uniqueGT[x], sep="") ))
  }))

load(file =  paste("Results/After_Clean/", store.allFCS[1,1],"/AftClean_", store.allFCS[1,2],".Rdata",sep="") )

## We will not include the Live marker while doing the FlowType
markers <- c("CD5","CD161","CD62L","CD44","CD8","CD4","CD25","GITR","CD24","gdTCR|TCRd|Syto*","KLRG", "CD45")
channels.ind.NoLive <- Find.markers(frame=f,marker.list=markers)
names(channels.ind.NoLive)[grep(names(channels.ind.NoLive),pattern = "TCRd*")] <- "TCRd"
channels.ind.NoLive <- sort(channels.ind.NoLive)

## We will not include the Live marker and the CD45 marker while doing the FlowType
markers <- c("CD5","CD161","CD62L","CD44","CD8","CD4","CD25","GITR","CD24","gdTCR|TCRd|Syto*","KLRG")
channels.ind.NoLiveNoCD45 <- Find.markers(frame=f,marker.list=markers)
names(channels.ind.NoLiveNoCD45)[grep(names(channels.ind.NoLiveNoCD45),pattern = "TCRd*")] <- "TCRd"
channels.ind.NoLiveNoCD45 <- sort(channels.ind.NoLiveNoCD45)


for(q2 in 1:nrow(store.allFCS)){
        start2 <- Sys.time()
        
        print(paste(q2, ": Starting ", store.allFCS[q2,1], " / ", store.allFCS[q2,2]), sep="" )
        
        ## Reading FCS files before Transform. Will need this for MFI calculation
        f.beforeTransform <- read.FCS(filename = allFCS[q2])
        ## Removing Margin events------------------------------------------------------------------------------------------
        scat.chans <- c(grep (colnames(f.beforeTransform),pattern = "FSC*"),grep (colnames(f.beforeTransform),pattern = "SSC*"))
        ## Removing margin events in Scatter channels
        f.beforeTransform <- removeMargins(f.beforeTransform, chans=scat.chans, verbose= verbose_margin)
        #Removing negative values in scatter channels
        f.beforeTransform <- removeMargins(f.beforeTransform, chans=scat.chans, debris=T, neg=T, verbose= verbose_debris) 
        
        
        ## Compensation
        print("Starting Compensation")
        if( det(f.beforeTransform@description$SPILL)==1 ){
          print("Check the spillover matrix, it's probably an identity matrix!")
          return(NULL)
        } else{
          f.beforeTransform <- compensate(f.beforeTransform, f.beforeTransform@description$SPILL)}
        
        load(file =  paste("Results/After_Clean/", store.allFCS[q2,1],"/AftClean_", store.allFCS[q2,2],".Rdata",sep="") )
        load(file =  paste("Results/Gating_Thresholds_Updated/", store.allFCS[q2,1],"/Gthres_", store.allFCS[q2,2],".Rdata",sep="") )  
        load(file =  paste("Results/Gating_Thresholds_CD45_Updated/", store.allFCS[q2,1],"/Gthres_", store.allFCS[q2,2],".Rdata",sep="") )  
        
        # Indices in the original clean FCS file (After Clean)
        f.ind <- c(1:nrow(f@exprs))
        all.events[q2,1] <- nrow(f)
        all.props[q2,1] <- (nrow(f)/nrow(f))*100
        ## Gates used for gating the CD45+. These gates will be used for the flowType
        listGthres <- as.list(gates.cd45)
        listGthres.NoLive <- listGthres[which(listGthres != listGthres$live.gate)] ## Excluding the Live marker threholds
        listGthres.NoLiveNoCD45 <- listGthres[which(listGthres != listGthres$live.gate & listGthres != listGthres$cd45.gate.low)] ## Excluding the Live marker & CD45 marker thresholds
        
        tryCatch({
              ## Plotting FSC-A_SSC-W. Gating all the Events to get the Singlets
              singlets.flowD.h <- flowDensity(f, channels = c("FSC-A", "SSC-W"), position = c(NA, F), gates = c(NA, gthres[1]))
              singlets.flowD.h.ind <- singlets.flowD.h@index
              singlets.flowD.l <- flowDensity(singlets.flowD.h, channels = c("FSC-A", "SSC-W"), position = c(NA, T), gates = c(NA, gthres[2]))
              singlets.flowD.l@proportion <- (singlets.flowD.l@cell.count/nrow(f))*100
              singlets <- getflowFrame(singlets.flowD.l)
              
              # Indices in the Singlets Population with respect to f. Will be needed for the MFI calculation
              singlets.flowD.l.ind <- singlets.flowD.l@index
              # # Optional Plot: Done to observe if the correct cells of the singlets population are being saved with respect to f
              # plotDens(f,c("FSC-A","SSC-W"),main="Ungated", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(f@exprs[singlets.flowD.l.ind,c(1,6)], col=2, cex = 0.1)
              
              all.events[q2,2] <- singlets.flowD.l@cell.count
              all.props[q2,2] <- singlets.flowD.l@proportion 
              
              # # Indices in the Singlets Population. Will be needed for the MFI calculation
              # singlets.ind <- intersect(intersect(f.ind, singlets.flowD.h.ind), singlets.flowD.l.ind)
              
              ###############################################################################################
              ## Gating Singlets to get the Live population. Plotting Live/Dead_SSC-A
              #live.flowD <- flowDensity(singlets, channels = c(11,4), position = c(F,NA), gates = c(gthres[3],NA))
              live.flowD <- flowDensity(singlets.flowD.l, channels = c(11,4), position = c(F,NA), gates = c(gthres[3],NA))
              live.flowD@proportion <- (live.flowD@cell.count/singlets.flowD.l@cell.count)*100
              live <- getflowFrame(live.flowD)
              
              # Indices in the Live Population with respect to f. Will be needed for the MFI calculation
              live.flowD.ind <- live.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the live population are being saved with respect to f
              # plotDens(singlets,  c(11,4), main= "Singlets", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(singlets.flowD.l@flow.frame@exprs[live.flowD.ind,c(11,4)], col=2, cex = 0.1)
              
              all.events[q2,3] <- live.flowD@cell.count
              all.props[q2,3] <- live.flowD@proportion
             
              ###############################################################################################
              ## Gating Live to get the Lymphocytes population. Plotting FSC-A_SSC-A
              lymph.flowD.temp <- flowDensity(live.flowD, channels = c(1,4), position = c(T, F), gates = c(gthres[5], gthres[6]))
              lymph.flowD.temp.ind <- lymph.flowD.temp@index
              lymph.flowD <- flowDensity(lymph.flowD.temp, channels = c(1,4), position = c(F,F), gates = c(gthres[4], gthres[6]))
              lymph.flowD@proportion <- (lymph.flowD@cell.count/live.flowD@cell.count)*100
              lymph <- getflowFrame(lymph.flowD)
              
              # Indices in the Lymphocytes Population with respect to f. Will be needed for the MFI calculation
              lymph.flowD.ind <- lymph.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the lymph population are being saved with respect to f
              # plotDens(live, c(1,4), main = "Live", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(live.flowD@flow.frame@exprs[lymph.flowD.ind,c(1,4)], col=2, cex = 0.1)
            
              all.events[q2,4] <- lymph.flowD@cell.count
              all.props[q2,4] <- lymph.flowD@proportion
    
              ###############################################################################################
              ## Gating Lymphocytes to get the CD45 population. Plotting CD45_CD161
              cd45.flowD.temp <- flowDensity(lymph.flowD, channels = c(14,15), position = c(T,F), gates = c(gthres[7], gthres[21]))
              cd45.flowD.temp.ind <- cd45.flowD.temp@index
              cd45.flowD <- flowDensity(cd45.flowD.temp, channels = c(14,15), position = c(F,F), gates = c(gthres[8], gthres[21]))
              cd45.flowD@proportion <- (cd45.flowD@cell.count/lymph.flowD@cell.count)*100
              cd45 <- getflowFrame(cd45.flowD)
              
              # Indices in the CD45 Population with respect to f. Will be needed for the MFI calculation
              cd45.flowD.ind <- cd45.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the CD45+ population are being saved with respect to f
              # plotDens(lymph, channels = c(14,15), main = "Lymphocytes", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(lymph.flowD@flow.frame@exprs[cd45.flowD.ind,c(14,15)], col=2, cex = 0.1)
              
              all.events[q2,5] <- cd45.flowD@cell.count 
              all.props[q2,5] <- cd45.flowD@proportion
              
              # # Indices in the CD45 Population. Will be needed for the MFI calculation
              # cd45.ind <- lymph.ind[cd45.flowD.temp.ind[cd45.flowD.ind]]
              
              ###############################################################################################
              # Gating CD45+. Highlighting the Autofluorescence and then getting rid of the Autofluorescence and grabbing the NOT(P2) population. Plotting TCRd_CD4
              theta = atan(1)
              cd45s.flowD <- cd45.flowD
              cd45s.flowD@flow.frame<-rotate.data(cd45.flowD@flow.frame,c(18,16),theta = -pi/4)$data
              R <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)),2, 2)

              returnTry <- tryCatch({
                cd45sTemp  <- flowDensity(cd45s.flowD, channels=c(18,16),position=c(T, T), gates=c(gthres[10], gthres[9]))
                #cd45sTemp.ind <- cd45sTemp@index
                cd45sTempb <- flowDensity(cd45sTemp, channels=c(18,16),position=c(F, T), gates=c(gthres[11], gthres[9]) )
                #cd45sTempb.ind <- cd45sTempb@index
                cd45autoflour <- getflowFrame(cd45sTempb)
                cd45autoflour@exprs[,c(18,16)] <- t(t(R) %*% t(cd45autoflour@exprs[,c(18,16)]))
                TCRDCD4points <- t(t(R) %*% t(cd45sTempb@filter))
                
                # Grabbing the NOT(P2) population using notSubFrame() function
                cd45.NOTP2.flowD  <- notSubFrame(cd45s.flowD@flow.frame, channels=c(18,16), position= "logical", gates = "missing", cd45sTempb@filter)
                
                # Rotating back the Autofluorescence population
                cd45sTempb@filter <- rotate.data(cd45sTempb@filter,c(18,16),theta = pi/4)$data
                cd45sTempb@flow.frame <- rotate.data(cd45sTempb@flow.frame,c(18,16),theta = pi/4)$data
                cd45sTempb@proportion <- (cd45sTempb@cell.count/cd45.flowD@cell.count)*100
                cd45autoflour <- getflowFrame(cd45sTempb)
                
                # Rotating back the NOT(P2) population
                cd45.NOTP2.flowD@filter <- rotate.data(cd45.NOTP2.flowD@filter,c(18,16),theta = pi/4)$data
                cd45.NOTP2.flowD@flow.frame <- rotate.data(cd45.NOTP2.flowD@flow.frame,c(18,16),theta = pi/4)$data
                cd45.NOTP2.flowD@proportion <- (cd45.NOTP2.flowD@cell.count/cd45.flowD@cell.count)*100
                cd45.NOTP2 <- getflowFrame(cd45.NOTP2.flowD)
                
                # Indices in the Autofluorescence Population with respect to f. Will be needed for the MFI calculation
                cd45sTempb.ind <- cd45sTempb@index
                
                # Indices in the CD45+ (without the Autofluorescence) Population with respect to f. Will be needed for the MFI calculation
                cd45.NOTP2.flowD.ind <- cd45.NOTP2.flowD@index
                cd45.NOTP2 <- getflowFrame(cd45.NOTP2.flowD)
                
                },error = function(err) {
                print("Error with Normal"); return(0)
              })
              
              # # Optional Plot: Done to observe if the correct cells of the Autofluorescence population are being saved with respect to f
              # plotDens(cd45, channels = c(18,16), main="CD45+_Autofluorecence", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd45.flowD@flow.frame@exprs[cd45sTempb.ind,c(18,16)], col=2, cex = 0.1)
              # 
              # # Optional Plot: Done to observe if the correct cells of the CD45+ (without the Autofluorescence) population are being saved with respect to f
              # plotDens(cd45, channels = c(18,16), main="CD45+_Autofluorecence", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd45.flowD@flow.frame@exprs[cd45.NOTP2.flowD.ind,c(18,16)], col=2, cex = 0.1)
  
              all.events[q2,6] <- cd45sTempb@cell.count
              all.props[q2,6] <- cd45sTempb@proportion
              
              all.events[q2,7] <- cd45.NOTP2.flowD@cell.count
              all.props[q2,7] <- cd45.NOTP2.flowD@proportion
              
              
              ###############################################################################################
              ## Gating CD45+ (without the autofluorescence) to obtain gd T-cells. Plotting TCRd_CD4 to obtain gd T-cells
              gd.Tcells.flowD.temp <- flowDensity(cd45.NOTP2.flowD, channels = c(18,16), position = c(T,F), gates = c(gthres[12], gthres[14]))
              #gd.Tcells.flowD.temp.ind <- gd.Tcells.flowD.temp@index
              gd.Tcells.flowD <- flowDensity(gd.Tcells.flowD.temp, channels = c(18,16), position = c(F,F), gates = c(gthres[13], gthres[14]))
              gd.Tcells.flowD@proportion <- (gd.Tcells.flowD@cell.count/cd45.NOTP2.flowD@cell.count)*100
              gd.Tcells <- getflowFrame(gd.Tcells.flowD)
              
              # Indices in the CD45+ (NOT(P2)) to get the gd T-cells Population with respect to f. Will be needed for the MFI calculation
              gd.Tcells.flowD.ind <- gd.Tcells.flowD@index
            
              all.events[q2,8] <- gd.Tcells.flowD@cell.count 
              all.props[q2,8] <- gd.Tcells.flowD@proportion
              
              # # Optional Plot: Done to observe if the correct cells of the gd T-cell population are being saved with respect to f
              # plotDens(cd45.NOTP2, channels = c(18,16), main="CD45+_Autofluorecence", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd45.NOTP2.flowD@flow.frame@exprs[gd.Tcells.flowD.ind,c(18,16)], col=2, cex = 0.1)
              
              ###############################################################################################
              ## Gating CD45+ (without the autofluorescence) to obtain NOT gd T-cells.
              not.gd.Tcells.flowD <- notSubFrame(cd45.NOTP2.flowD@flow.frame, channels = c(18,16), position= "logical", gates = "missing", gd.Tcells.flowD@filter)
              not.gd.Tcells.flowD@proportion <- (not.gd.Tcells.flowD@cell.count/cd45.NOTP2.flowD@cell.count)*100
              not.gd.Tcells <- getflowFrame(not.gd.Tcells.flowD)
              
              # Indices in the CD45+ (NOT(P2)) to get the NOT gd T-cells Population. Will be needed for the MFI calculation
              not.gd.Tcells.ind <- not.gd.Tcells.flowD@index

              all.events[q2,9] <- not.gd.Tcells.flowD@cell.count
              all.props[q2,9] <- not.gd.Tcells.flowD@proportion
              
              # # Optional Plot: Done to observe if the correct cells of the NOT(gd T-cell) population are being saved with respect to f
              # plotDens(cd45.NOTP2, channels = c(18,16), main="CD45+_Autofluorecence", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd45.NOTP2.flowD@flow.frame@exprs[not.gd.Tcells.ind,c(18,16)], col=2, cex = 0.1)
              
              
              
              ###############################################################################################
              ## Gating gd T-cells.  Plotting CD62L_CD44 to obtain gd Resting, gd Effector, and gd Naive
              gd.Resting.flowD <- flowDensity(gd.Tcells.flowD, channels = c(8,7), position = c(T,T), gates = c(gthres[15], gthres[16]))
              gd.Resting.flowD@proportion <- (gd.Resting.flowD@cell.count/gd.Tcells.flowD@cell.count)*100
              # Indices in the gd T-cells to get the gd Resting Population with respect to f. Will be needed for the MFI calculation
              gd.Resting.flowD.ind <- gd.Resting.flowD@index
              
              all.events[q2,10] <- gd.Resting.flowD@cell.count
              all.props[q2,10] <- gd.Resting.flowD@proportion
              
              # # Optional Plot: Done to observe if the correct cells of the gd Resting population are being saved with respect to f
              # plotDens(gd.Tcells, channels = c(8,7), main="gd Resting", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(gd.Tcells.flowD@flow.frame@exprs[gd.Resting.flowD.ind,c(8,7)], col=2, cex = 0.1)
              
              
              
              gd.Effector.flowD <- flowDensity(gd.Tcells.flowD, channels = c(8,7), position = c(F,T), gates = c(gthres[15], gthres[16]))
              gd.Effector.flowD@proportion <- (gd.Effector.flowD@cell.count/gd.Tcells.flowD@cell.count)*100
              # Indices in the gd T-cells to get the gd Effector Population with respect to f. Will be needed for the MFI calculation
              gd.Effector.flowD.ind <- gd.Effector.flowD@index
              
              all.events[q2,11] <- gd.Effector.flowD@cell.count
              all.props[q2,11] <- gd.Effector.flowD@proportion
              # # Optional Plot: Done to observe if the correct cells of the gd Effector population are being saved with respect to f
              # plotDens(gd.Tcells, channels = c(8,7), main="gd Effector", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(gd.Tcells.flowD@flow.frame@exprs[gd.Effector.flowD.ind,c(8,7)], col=2, cex = 0.1)
              
                       
              
              gd.Naive.flowD <- flowDensity(gd.Tcells.flowD, channels = c(8,7), position = c(T,F), gates = c(gthres[15], gthres[16]))
              gd.Naive.flowD@proportion <- (gd.Naive.flowD@cell.count/gd.Tcells.flowD@cell.count)*100
              # Indices in the gd T-cells to get the gd Naive Population with respect to f. Will be needed for the MFI calculation
              gd.Naive.flowD.ind <- gd.Naive.flowD@index
              
              all.events[q2,12] <- gd.Naive.flowD@cell.count
              all.props[q2,12] <- gd.Naive.flowD@proportion
              # # Optional Plot: Done to observe if the correct cells of the gd Naive population are being saved with respect to f
              # plotDens(gd.Tcells, channels = c(8,7), main="gd Naive", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(gd.Tcells.flowD@flow.frame@exprs[gd.Naive.flowD.ind,c(8,7)], col=2, cex = 0.1)
              
                      
              
              #############################################################################################
            
              ## Gating gd T-cells. Plotting KLRG1_GITR to obtain GD KLRG1+ and GITR gd T-cells 
              gd.klrg1.flowD <- flowDensity(gd.Tcells.flowD, channels = c(12,17), position = c(T,F), gates = c(gthres[17], gthres[18]))
              gd.klrg1.flowD@proportion <- (gd.klrg1.flowD@cell.count/gd.Tcells.flowD@cell.count)*100
              # Indices in the gd T-cells to get the gd KLRG1+ Population with respect to f. Will be needed for the MFI calculation
              gd.klrg1.flowD.ind <- gd.klrg1.flowD@index
              
              all.events[q2,13] <- gd.klrg1.flowD@cell.count
              all.props[q2,13] <- gd.klrg1.flowD@proportion
              # # Optional Plot: Done to observe if the correct cells of the gd KLRG1+ population are being saved with respect to f
              # plotDens(gd.Tcells, channels = c(12,17), main="gd KLRG1+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(gd.Tcells.flowD@flow.frame@exprs[gd.klrg1.flowD.ind,c(12,17)], col=2, cex = 0.1)
              
         
              
              gd.gitr.flowD <- flowDensity(gd.Tcells.flowD, channels = c(12,17), position = c(F,T), gates = c(gthres[17], gthres[18]))
              gd.gitr.flowD@proportion <- (gd.gitr.flowD@cell.count/gd.Tcells.flowD@cell.count)*100
              # Indices in the gd T-cells to get the gd GITR Population with respect to f. Will be needed for the MFI calculation
              gd.gitr.flowD.ind <- gd.gitr.flowD@index
              
              all.events[q2,14] <- gd.gitr.flowD@cell.count
              all.props[q2,14] <- gd.gitr.flowD@proportion
              # # Optional Plot: Done to observe if the correct cells of the gd GITR population are being saved with respect to f
              # plotDens(gd.Tcells, channels = c(12,17), main="gd GITR", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(gd.Tcells.flowD@flow.frame@exprs[gd.gitr.flowD.ind,c(12,17)], col=2, cex = 0.1)
              
     
              
              #############################################################################################
              ## Gating gd T-cells. Plotting CD5_CD44 to obtain GD CD5+
              gd.cd5.flowD <- flowDensity(gd.Tcells.flowD, channels = c(13,7), position = c(T,NA), gates = c(gthres[19], NA))
              gd.cd5.flowD@proportion <- (gd.cd5.flowD@cell.count/gd.Tcells.flowD@cell.count)*100
              # Indices in the gd T-cells to get the gd GITR Population with respect to f. Will be needed for the MFI calculation
              gd.cd5.flowD.ind <- gd.cd5.flowD@index
              
              all.events[q2,15] <- gd.cd5.flowD@cell.count
              all.props[q2,15] <- gd.cd5.flowD@proportion
              # # Optional Plot: Done to observe if the correct cells of the gd CD5+ population are being saved with respect to f
              # plotDens(gd.Tcells, channels = c(13,7), main="gd CD5+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(gd.Tcells.flowD@flow.frame@exprs[gd.cd5.flowD.ind,c(13,7)], col=2, cex = 0.1)
              
              
              ##############################################################################################
              ## Gating NOT gd T-cells. Plotting CD161_CD5 to obtain CD5+ and NK-cells
              # CD5+
              cd5.flowD <- flowDensity(not.gd.Tcells.flowD, channels = c(15,13), position = c(NA,T), gates = c(NA,gthres[19]))
              cd5.flowD@proportion <- (cd5.flowD@cell.count/not.gd.Tcells.flowD@cell.count)*100
              cd5 <- getflowFrame(cd5.flowD)
              # Indices in the NOT gd T-cells to get the CD5+ Population with respect to f. Will be needed for the MFI calculation
              cd5.flowD.ind <- cd5.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the gd CD5+ population are being saved with respect to f
              # plotDens(not.gd.Tcells, channels = c(15,13), main="NOT(gd T-cells) CD5+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(not.gd.Tcells.flowD@flow.frame@exprs[cd5.flowD.ind,c(15,13)], col=2, cex = 0.1)
              
              all.events[q2,16] <- cd5.flowD@cell.count
              all.props[q2,16] <- cd5.flowD@proportion
             
              
              # NK-cells
              NKcells.flowD.temp <- flowDensity(not.gd.Tcells.flowD, channels = c(15,13), position = c(T,F), gates = c(gthres[20],gthres[19]))
              #NKcells.flowD.temp.ind <- NKcells.flowD.temp@index
              NKcells.flowD <- flowDensity(NKcells.flowD.temp, channels = c(15,13), position = c(F,F), gates = c(gthres[21], gthres[19]))
              NKcells.flowD@proportion <- (NKcells.flowD@cell.count/not.gd.Tcells.flowD@cell.count)*100
              NKcells <- getflowFrame(NKcells.flowD)
              # Indices in the NOT gd T-cells to get the NK-cells Population with respect to f. Will be needed for the MFI calculation
              NKcells.flowD.ind <- NKcells.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the NK-cells population are being saved with respect to f
              # plotDens(not.gd.Tcells, channels = c(15,13), main="NOT(gd T-cells) NK-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(not.gd.Tcells.flowD@flow.frame@exprs[NKcells.flowD.ind,c(15,13)], col=2, cex = 0.1)
              
              
              all.events[q2,17] <- NKcells.flowD@cell.count
              all.props[q2,17] <- NKcells.flowD@proportion
              
              ###############################################################################################
              ## Gating NK-cells. Plotting CD62L_CD44 to obtain NK Resting/Naive and NK Effector
              NK.Resting.Naive.flowD <- flowDensity(NKcells.flowD, channels = c(8,7), position = c(T,NA), gates = c(gthres[15],NA))
              NK.Resting.Naive.flowD@proportion <- (NK.Resting.Naive.flowD@cell.count/NKcells.flowD@cell.count)*100
              # Indices in the NK-cells to get the NK Resting/Naive Population with respect to f. Will be needed for the MFI calculation
              NK.Resting.Naive.flowD.ind <- NK.Resting.Naive.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the NK-cells:Resting/Naive population are being saved with respect to f
              # plotDens(NKcells, channels = c(8,7), main="NK-cells: Resting/Naive", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(NKcells.flowD@flow.frame@exprs[NK.Resting.Naive.flowD.ind,c(8,7)], col=2, cex = 0.1)
              
              all.events[q2,18] <- NK.Resting.Naive.flowD@cell.count
              all.props[q2,18] <- NK.Resting.Naive.flowD@proportion
              
              
              
              NK.Effector.flowD <- flowDensity(NKcells.flowD, channels = c(8,7), position = c(F,T), gates = c(gthres[15], gthres[16]))
              NK.Effector.flowD@proportion <- (NK.Effector.flowD@cell.count/NKcells.flowD@cell.count)*100
              # Indices in the NK-cells to get the NK Effector Population with respect to f. Will be needed for the MFI calculation
              NK.Effector.flowD.ind <- NK.Effector.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the NK-cells:Effector population are being saved with respect to f
              # plotDens(NKcells, channels = c(8,7), main="NK-cells: Effector", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(NKcells.flowD@flow.frame@exprs[NK.Effector.flowD.ind,c(8,7)], col=2, cex = 0.1)
              
              all.events[q2,19] <- NK.Effector.flowD@cell.count
              all.props[q2,19] <- NK.Effector.flowD@proportion
              
              
              ###############################################################################################
              ## Gating NK-cells. Plotting KLRG1_CD44 to obtain NK KLRG1+
              NK.klrg1.flowD <- flowDensity(NKcells.flowD, channels = c(12,7), position = c(T,T), gates = c(gthres[17], gthres[16]))
              NK.klrg1.flowD@proportion <- (NK.klrg1.flowD@cell.count/NKcells.flowD@cell.count)*100
              # Indices in the NK-cells to get the NK KLRG1+ Population with respect to f. Will be needed for the MFI calculation
              NK.klrg1.flowD.ind <- NK.klrg1.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the NK-cells:KLRG1+ population are being saved with respect to f
              # plotDens(NKcells, channels = c(12,7), main="NK-cells: KLRG1+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(NKcells.flowD@flow.frame@exprs[NK.klrg1.flowD.ind,c(12,7)], col=2, cex = 0.1)
              
              all.events[q2,20] <- NK.klrg1.flowD@cell.count
              all.props[q2,20] <- NK.klrg1.flowD@proportion
              
              ###########################################################################################
              # Gating CD5+. Plotting CD161_CD4 to obtain P2a, CD4- NKT-cells, and CD4+ NKT-cells
              # P2a
              P2a.flowD <- flowDensity(cd5.flowD, channels = c(15,16), position = c(F, NA), gates = c(gthres[20],NA))
              P2a.flowD@proportion <- (P2a.flowD@cell.count/cd5.flowD@cell.count)*100
              P2a <- getflowFrame(P2a.flowD)
              # Indices in the CD5+ to get the P2a Population with respect to f. Will be needed for the MFI calculation
              P2a.flowD.ind <- P2a.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the P2a population are being saved with respect to f
              # plotDens(cd5, channels = c(15,16), main="CD5+: P2a", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd5.flowD@flow.frame@exprs[P2a.flowD.ind,c(15,16)], col=2, cex = 0.1)
              
              all.events[q2,21] <- P2a.flowD@cell.count
              all.props[q2,21] <- P2a.flowD@proportion
              
              
              # CD4- NKT-cells
              cd4neg.NKTcells.flowD <- flowDensity(cd5.flowD, channels = c(15,16), position = c(T,F), gates = c(gthres[20], gthres[14]))
              cd4neg.NKTcells.flowD@proportion <- (cd4neg.NKTcells.flowD@cell.count/cd5.flowD@cell.count)*100
              cd4neg.NKTcells <- getflowFrame(cd4neg.NKTcells.flowD)
              # Indices in the CD5+ to get the CD4- NKT Population with respect to f. Will be needed for the MFI calculation
              cd4neg.NKTcells.flowD.ind <- cd4neg.NKTcells.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the CD4- population are being saved with respect to f
              # plotDens(cd5, channels = c(15,16), main="CD5+: CD4- NKT", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd5.flowD@flow.frame@exprs[cd4neg.NKTcells.flowD.ind,c(15,16)], col=2, cex = 0.1)
              
              all.events[q2,22] <- cd4neg.NKTcells.flowD@cell.count
              all.props[q2,22] <- cd4neg.NKTcells.flowD@proportion
              
              
              
              # CD4+ NKT-cells
              cd4pos.NKTcells.flowD <- flowDensity(cd5.flowD, channels = c(15,16), position = c(T,T), gates = c(gthres[20], gthres[14]))
              cd4pos.NKTcells.flowD@proportion <- (cd4pos.NKTcells.flowD@cell.count/cd5.flowD@cell.count)*100
              cd4pos.NKTcells <- getflowFrame(cd4pos.NKTcells.flowD)
              # Indices in the CD5+ to get the CD4+ NKT Population with respect to f. Will be needed for the MFI calculation
              cd4pos.NKTcells.flowD.ind <- cd4pos.NKTcells.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the CD4+ population are being saved with respect to f
              # plotDens(cd5, channels = c(15,16), main="CD5+: CD4+ NKT", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd5.flowD@flow.frame@exprs[cd4pos.NKTcells.flowD.ind,c(15,16)], col=2, cex = 0.1)
              
              all.events[q2,23] <- cd4pos.NKTcells.flowD@cell.count
              all.props[q2,23] <- cd4pos.NKTcells.flowD@proportion
              
              ###########################################################################################
              # Gating CD4- NKT-cells. Plotting CD62L_CD44 to obtain CD4- NKT Effector and CD4- NKT Resting/Naive
              cd4neg.NKT.Effector.flowD <- flowDensity(cd4neg.NKTcells.flowD, channels = c(8,7), position = c(F,T), gates = c(gthres[15], gthres[16]))
              cd4neg.NKT.Effector.flowD@proportion <- (cd4neg.NKT.Effector.flowD@cell.count/cd4neg.NKTcells.flowD@cell.count)*100
              # Indices in the CD4- NKT-cells to get the CD4- NKT Effector Population with respect to f. Will be needed for the MFI calculation
              cd4neg.NKT.Effector.flowD.ind <- cd4neg.NKT.Effector.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the CD4- NKT Effector population are being saved with respect to f
              # plotDens(cd4neg.NKTcells, channels = c(8,7), main="CD4- NKT Effector", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd4neg.NKTcells.flowD@flow.frame@exprs[cd4neg.NKT.Effector.flowD.ind,c(8,7)], col=2, cex = 0.1)
              
              all.events[q2,24] <- cd4neg.NKT.Effector.flowD@cell.count
              all.props[q2,24] <- cd4neg.NKT.Effector.flowD@proportion
            
              
              cd4neg.NKT.Resting.flowD <- flowDensity(cd4neg.NKTcells.flowD, channels = c(8,7), position = c(T,NA), gates = c(gthres[15], NA))
              cd4neg.NKT.Resting.flowD@proportion <- (cd4neg.NKT.Resting.flowD@cell.count/cd4neg.NKTcells.flowD@cell.count)*100
              # Indices in the CD4- NKT-cells to get the CD4- NKT Resting/Naive Population with respect to f. Will be needed for the MFI calculation
              cd4neg.NKT.Resting.flowD.ind <- cd4neg.NKT.Resting.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the CD4- NKT Resting/Naive population are being saved with respect to f
              # plotDens(cd4neg.NKTcells, channels = c(8,7), main="CD4- NKT Resting/Naive", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd4neg.NKTcells.flowD@flow.frame@exprs[cd4neg.NKT.Resting.flowD.ind,c(8,7)], col=2, cex = 0.1)
              
              all.events[q2,25] <- cd4neg.NKT.Resting.flowD@cell.count
              all.props[q2,25] <- cd4neg.NKT.Resting.flowD@proportion
                      
              ###########################################################################################
              # Gating CD4- NKT-cells. Plotting KLRG1_CD44 to obtain CD4- NKT KLRG1+
              cd4neg.NKT.klrg1.flowD <- flowDensity(cd4neg.NKTcells.flowD, channels = c(12,7), position = c(T,T), gates = c(gthres[17], gthres[16]))
              cd4neg.NKT.klrg1.flowD@proportion <- (cd4neg.NKT.klrg1.flowD@cell.count/cd4neg.NKTcells.flowD@cell.count)*100
              # Indices in the CD4- NKT-cells to get the CD4- NKT KLRG1+ Population with respect to f. Will be needed for the MFI calculation
              cd4neg.NKT.klrg1.flowD.ind <- cd4neg.NKT.klrg1.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the CD4- NKT KLRG1+ population are being saved with respect to f
              # plotDens(cd4neg.NKTcells, channels = c(12,7), main="CD4- NKT KLRG1+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd4neg.NKTcells.flowD@flow.frame@exprs[cd4neg.NKT.klrg1.flowD.ind,c(12,7)], col=2, cex = 0.1)
              
              all.events[q2,26] <- cd4neg.NKT.klrg1.flowD@cell.count
              all.props[q2,26] <- cd4neg.NKT.klrg1.flowD@proportion
              
              
              #############################################################################################
              # Gating CD4+ NKT-cells. Plotting CD62L_CD44 to obtain CD4+ NKT Effector and CD4+ NKT Resting/Naive
              cd4pos.NKT.Effector.flowD <- flowDensity(cd4pos.NKTcells.flowD, channels = c(8,7), position = c(F,T), gates = c(gthres[15], gthres[16]))
              cd4pos.NKT.Effector.flowD@proportion <- (cd4pos.NKT.Effector.flowD@cell.count/cd4pos.NKTcells.flowD@cell.count)*100
              # Indices in the CD4+ NKT-cells to get the CD4+ NKT Effector Population with respect to f. Will be needed for the MFI calculation
              cd4pos.NKT.Effector.flowD.ind <- cd4pos.NKT.Effector.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the CD4+ NKT Effector population are being saved with respect to f
              # plotDens(cd4pos.NKTcells, channels = c(8,7), main="CD4+ NKT Effector", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd4pos.NKTcells.flowD@flow.frame@exprs[cd4pos.NKT.Effector.flowD.ind,c(8,7)], col=2, cex = 0.1)
              
              all.events[q2,27] <- cd4pos.NKT.Effector.flowD@cell.count
              all.props[q2,27] <- cd4pos.NKT.Effector.flowD@proportion
              
                            
              
              cd4pos.NKT.Resting.flowD <- flowDensity(cd4pos.NKTcells.flowD, channels = c(8,7), position = c(T,NA), gates = c(gthres[15], NA))
              cd4pos.NKT.Resting.flowD@proportion <- (cd4pos.NKT.Resting.flowD@cell.count/cd4pos.NKTcells.flowD@cell.count)*100
              # Indices in the CD4+ NKT-cells to get the CD4+ NKT Resting/Naive Population with respect to f. Will be needed for the MFI calculation
              cd4pos.NKT.Resting.flowD.ind <- cd4pos.NKT.Resting.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the CD4+ NKT Resting/Naive population are being saved with respect to f
              # plotDens(cd4pos.NKTcells, channels = c(8,7), main="CD4+ NKT Resting/Naive", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd4pos.NKTcells.flowD@flow.frame@exprs[cd4pos.NKT.Resting.flowD.ind,c(8,7)], col=2, cex = 0.1)
              
              all.events[q2,28] <- cd4pos.NKT.Resting.flowD@cell.count
              all.props[q2,28] <- cd4pos.NKT.Resting.flowD@proportion
           
              
              ###########################################################################################
              # Gating CD4+ NKT-cells. Plotting KLRG1_CD44 to obtain CD4+ NKT KLRG1+
              cd4pos.NKT.klrg1.flowD <- flowDensity(cd4pos.NKTcells.flowD, channels = c(12,7), position = c(T,T), gates = c(gthres[17], gthres[16]))
              cd4pos.NKT.klrg1.flowD@proportion <- (cd4pos.NKT.klrg1.flowD@cell.count/cd4pos.NKTcells.flowD@cell.count)*100
              # Indices in the CD4+ NKT-cells to get the CD4+ NKT KLRG1+ Population with respect to f. Will be needed for the MFI calculation
              cd4pos.NKT.klrg1.flowD.ind <- cd4pos.NKT.klrg1.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the CD4+ NKT KLRG1+ population are being saved with respect to f
              # plotDens(cd4pos.NKTcells, channels = c(12,7), main="CD4+ NKT KLRG1+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd4pos.NKTcells.flowD@flow.frame@exprs[cd4pos.NKT.klrg1.flowD.ind,c(12,7)], col=2, cex = 0.1)
              
              all.events[q2,29] <- cd4pos.NKT.klrg1.flowD@cell.count
              all.props[q2,29] <- cd4pos.NKT.klrg1.flowD@proportion
              
     
              ###########################################################################################
              # Gating P2a. Plotting Cd8a_CD4 to obtain CD5+ CD4/CD8 and abT-cells.
              # CD5+ CD4/CD8
              cd5.cd4cd8.flowD <- flowDensity(P2a.flowD, channels = c(10,16), position = c(F,F), gates = c(gthres[22], gthres[14]))
              cd5.cd4cd8.flowD@proportion <- (cd5.cd4cd8.flowD@cell.count/P2a.flowD@cell.count)*100
              # Indices in the P2a cells to get the CD5+ CD4/CD8 Population. Will be needed for the MFI calculation
              cd5.cd4cd8.flowD.ind <- cd5.cd4cd8.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the P2a CD5+ CD4/CD8 population are being saved with respect to f
              # plotDens(P2a, channels = c(10,16), main="P2a: CD5+ CD4/CD8", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(P2a.flowD@flow.frame@exprs[cd5.cd4cd8.flowD.ind,c(10,16)], col=2, cex = 0.1)

              all.events[q2,30] <- cd5.cd4cd8.flowD@cell.count
              all.props[q2,30] <- cd5.cd4cd8.flowD@proportion
              

              
              # abT-cells
              abTcells.flowD <- notSubFrame(P2a.flowD@flow.frame, channels = c(10,16), position= "logical", gates = "missing", cd5.cd4cd8.flowD@filter)
              abTcells.flowD@proportion <- (abTcells.flowD@cell.count/P2a.flowD@cell.count)*100
              abTcells <- getflowFrame(abTcells.flowD)
              # Indices in the P2a cells to get the abT-cells Population. Will be needed for the MFI calculation
              abTcells.flowD.ind <- abTcells.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the P2a abT-cell population are being saved with respect to f
              # plotDens(P2a, channels = c(10,16), main="P2a: abT-cell", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(P2a.flowD@flow.frame@exprs[abTcells.flowD.ind,c(10,16)], col=2, cex = 0.1)
              
              all.events[q2,31] <- abTcells.flowD@cell.count
              all.props[q2,31] <- abTcells.flowD@proportion
              

              ###########################################################################################
              # Gating abT-cells. Plotting CD8a_CD4 to obtain CD8a+ T-cells and CD4+ T-cells
              # CD8a+ T-cells
              cd8a.Tcells.flowD <- flowDensity(abTcells.flowD, channels = c(10,16), position = c(T,F), gates = c(gthres[22], gthres[14]) )
              cd8a.Tcells.flowD@proportion <- (cd8a.Tcells.flowD@cell.count/abTcells.flowD@cell.count)*100
              cd8a.Tcells <- getflowFrame(cd8a.Tcells.flowD)
              # Indices in the abT-cells to get the CD8a+ T-cells Population. Will be needed for the MFI calculation
              cd8a.Tcells.flowD.ind <- cd8a.Tcells.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the abT-cell: CD8a+ population are being saved with respect to f
              # plotDens(abTcells, channels = c(10,16), main="abT-cell: CD8a+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(abTcells.flowD@flow.frame@exprs[cd8a.Tcells.flowD.ind,c(10,16)], col=2, cex = 0.1)
              
              all.events[q2,32] <- cd8a.Tcells.flowD@cell.count
              all.props[q2,32] <- cd8a.Tcells.flowD@proportion
              
             
              
              # CD4+ T-cells
              cd4.Tcells.flowD <- flowDensity(abTcells.flowD, channels = c(10,16), position = c(F,T), gates = c(gthres[22], gthres[14]) )
              cd4.Tcells.flowD@proportion <- (cd4.Tcells.flowD@cell.count/abTcells.flowD@cell.count)*100
              cd4.Tcells <- getflowFrame(cd4.Tcells.flowD)
              # Indices in the abT-cells to get the CD4+ T-cells Population with respect to f. Will be needed for the MFI calculation
              cd4.Tcells.flowD.ind <- cd4.Tcells.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the abT-cell: CD4+ population are being saved with respect to f
              # plotDens(abTcells, channels = c(10,16), main="abT-cell: CD4+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(abTcells.flowD@flow.frame@exprs[cd4.Tcells.flowD.ind,c(10,16)], col=2, cex = 0.1)
              
              all.events[q2,33] <- cd4.Tcells.flowD@cell.count
              all.props[q2,33] <- cd4.Tcells.flowD@proportion
              
              
              ##############################################################################################
              # Gating CD8a+ T-cells. Plotting CD62L_CD44 to obtain CD8 Naive, CD8 Effector, and CD8 Resting
              # CD8 Naive
              cd8.Naive.flowD <- flowDensity(cd8a.Tcells.flowD, channels = c(8,7), position = c(T,F), gates = c(gthres[15], gthres[23]))
              cd8.Naive.flowD@proportion <- (cd8.Naive.flowD@cell.count/cd8a.Tcells.flowD@cell.count)*100
              # Indices in the CD8a+ T-cells to get the CD8 Naive Population with respect to f. Will be needed for the MFI calculation
              cd8.Naive.flowD.ind <- cd8.Naive.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the CD8 Naive Population population are being saved with respect to f
              # plotDens(cd8a.Tcells, channels = c(8,7), main="CD8 Naive Population", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd8a.Tcells.flowD@flow.frame@exprs[cd8.Naive.flowD.ind,c(8,7)], col=2, cex = 0.1)
              
              all.events[q2,34] <- cd8.Naive.flowD@cell.count
              all.props[q2,34] <- cd8.Naive.flowD@proportion
              

              
              # CD8 Effector
              cd8.Effector.flowD <- flowDensity(cd8a.Tcells.flowD, channels = c(8,7), position = c(F,T), gates = c(gthres[15], gthres[16]))
              cd8.Effector.flowD@proportion <- (cd8.Effector.flowD@cell.count/cd8a.Tcells.flowD@cell.count)*100
              # Indices in the CD8a+ T-cells to get the CD8 Effector Population with respect to f. Will be needed for the MFI calculation
              cd8.Effector.flowD.ind <- cd8.Effector.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the CD8 Effector Population population are being saved with respect to f
              # plotDens(cd8a.Tcells, channels = c(8,7), main="CD8 Effector Population", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd8a.Tcells.flowD@flow.frame@exprs[cd8.Effector.flowD.ind,c(8,7)], col=2, cex = 0.1)
              
              all.events[q2,35] <- cd8.Effector.flowD@cell.count
              all.props[q2,35] <- cd8.Effector.flowD@proportion
       
                           
              
              # CD8 Resting
              cd8.Resting.flowD <- flowDensity(cd8a.Tcells.flowD, channels = c(8,7), position = c(T,T), gates = c(gthres[15], gthres[23]))
              cd8.Resting.flowD@proportion <- (cd8.Resting.flowD@cell.count/cd8a.Tcells.flowD@cell.count)*100
              # Indices in the CD8a+ T-cells to get the CD8 Resting Population with respect to f. Will be needed for the MFI calculation
              cd8.Resting.flowD.ind <- cd8.Resting.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the CD8 Resting Population population are being saved with respect to f
              # plotDens(cd8a.Tcells, channels = c(8,7), main="CD8 Resting Population", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd8a.Tcells.flowD@flow.frame@exprs[cd8.Resting.flowD.ind,c(8,7)], col=2, cex = 0.1)
              
              all.events[q2,36] <- cd8.Resting.flowD@cell.count
              all.props[q2,36] <- cd8.Resting.flowD@proportion
             

              ################################################################################################
              # Gating CD8a+ T-cells. Plotting KLRG1_CD44 to obtain CD8 KLRG1+
              cd8.KLRG1.flowD <- flowDensity(cd8a.Tcells.flowD, channels = c(12,7), position = c(T,T), gates = c(gthres[17],gthres[16]))
              cd8.KLRG1.flowD@proportion <- (cd8.KLRG1.flowD@cell.count/cd8a.Tcells.flowD@cell.count)*100
              # Indices in the CD8a+ T-cells to get the CD8 KLRG1+ Population with respect to f. Will be needed for the MFI calculation
              cd8.KLRG1.flowD.ind <- cd8.KLRG1.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the CD8 KLRG1+ Population population are being saved with respect to f
              # plotDens(cd8a.Tcells, channels = c(12,7), main="CD8 KLRG1+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd8a.Tcells.flowD@flow.frame@exprs[cd8.KLRG1.flowD.ind,c(12,7)], col=2, cex = 0.1)
              
              all.events[q2,37] <- cd8.KLRG1.flowD@cell.count
              all.props[q2,37] <- cd8.KLRG1.flowD@proportion
             
    
              ##################################################################################################
              # Gating CD4+ T-cells. Plotting CD25_GITR to obtain T-helper cells and Tregs
              # T-helper cells
              T.helper.flowD.temp <- flowDensity(cd4.Tcells.flowD, channels = c(9,17), position = c(F,NA), gates = c(gthres[25], gthres[18]))
              #T.helper.flowD.temp.ind <- T.helper.flowD.temp@index
              T.helper.flowD <- flowDensity(T.helper.flowD.temp, channels = c(9,17), position = c(T,NA), gates = c(gthres[24], gthres[18]))
              T.helper.flowD@proportion <- (T.helper.flowD@cell.count/cd4.Tcells.flowD@cell.count)*100
              T.helper <- getflowFrame(T.helper.flowD)
              # Indices in the CD4+ T-cells to get the T-helper Population with respect to f. Will be needed for the MFI calculation
              T.helper.flowD.ind <- T.helper.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the T-helper Population population are being saved with respect to f
              # plotDens(cd4.Tcells, channels = c(9,17), main="T-helper", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd4.Tcells.flowD@flow.frame@exprs[T.helper.flowD.ind,c(9,17)], col=2, cex = 0.1)
              
              all.events[q2,38] <- T.helper.flowD@cell.count
              all.props[q2,38] <- T.helper.flowD@proportion
              
              
                
                
              # Tregs
              Tregs.flowD <- flowDensity(cd4.Tcells.flowD, channels = c(9,17), position = c(T,T), gates = c(gthres[25], gthres[18]))
              Tregs.flowD@proportion <- (Tregs.flowD@cell.count/cd4.Tcells.flowD@cell.count)*100
              Tregs <- getflowFrame(Tregs.flowD)
              # Indices in the CD4+ T-cells to get the Tregs Population with respect to f. Will be needed for the MFI calculation
              Tregs.flowD.ind <- Tregs.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the Tregs Population population are being saved with respect to f
              # plotDens(cd4.Tcells, channels = c(9,17), main="Tregs", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(cd4.Tcells.flowD@flow.frame@exprs[Tregs.flowD.ind,c(9,17)], col=2, cex = 0.1)
              
              all.events[q2,39] <- Tregs.flowD@cell.count
              all.props[q2,39] <- Tregs.flowD@proportion
              
              
              
               ####################################################################################################
              # Gating T-helper cells. Plotting CD62L_CD44 to obtain CD4 Resting/Naive and CD4 Effector
              # CD4 Resting/Naive
              cd4.Resting.flowD <- flowDensity(T.helper.flowD, channels = c(8,7), position = c(T,NA), gates = c(gthres[15],NA))
              cd4.Resting.flowD@proportion <- (cd4.Resting.flowD@cell.count/T.helper.flowD@cell.count)*100
              # Indices in the T-helper cells to get the CD4 Resting/Naive Population with respect to f. Will be needed for the MFI calculation
              cd4.Resting.flowD.ind <- cd4.Resting.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the CD4 Resting/Naive Population population are being saved with respect to f
              # plotDens(T.helper, channels = c(8,7), main="CD4 Resting/Naive", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(T.helper.flowD@flow.frame@exprs[cd4.Resting.flowD.ind,c(8,7)], col=2, cex = 0.1)
              
              all.events[q2,40] <- cd4.Resting.flowD@cell.count
              all.props[q2,40] <- cd4.Resting.flowD@proportion
              
              
              
              # CD4 Effector
              cd4.Effector.flowD <- flowDensity(T.helper.flowD, channels = c(8,7), position = c(F,T), gates = c(gthres[15], gthres[16]))
              cd4.Effector.flowD@proportion <- (cd4.Effector.flowD@cell.count/T.helper.flowD@cell.count)*100
              # Indices in the T-helper cells to get the CD4 Effector Population with respect to f. Will be needed for the MFI calculation
              cd4.Effector.flowD.ind <- cd4.Effector.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the CD4 Effector Population population are being saved with respect to f
              # plotDens(T.helper, channels = c(8,7), main="CD4 Effector", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(T.helper.flowD@flow.frame@exprs[cd4.Effector.flowD.ind,c(8,7)], col=2, cex = 0.1)
              
              
              all.events[q2,41] <- cd4.Effector.flowD@cell.count
              all.props[q2,41] <- cd4.Effector.flowD@proportion
              
             

              ####################################################################################################
              # Gating T-helper cells. Plotting KLRG1_CD44 to obtain CD4+ KLRG1+ T-cells
              cd4.KLRG1.flowD <- flowDensity(T.helper.flowD, channels = c(12,7), position = c(T,T), gates = c(gthres[17], gthres[23]))
              cd4.KLRG1.flowD@proportion <- (cd4.KLRG1.flowD@cell.count/T.helper.flowD@cell.count)*100
              # Indices in the T-helper cells to get the CD4+ KLRG1+ T-cells Population with respect to f. Will be needed for the MFI calculation
              cd4.KLRG1.flowD.ind <- cd4.KLRG1.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the CD4+ KLRG1+ T-cells Population population are being saved with respect to f
              # plotDens(T.helper, channels = c(12,7), main="CD4+ KLRG1+ T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(T.helper.flowD@flow.frame@exprs[cd4.KLRG1.flowD.ind,c(12,7)], col=2, cex = 0.1)
              
              
              all.events[q2,42] <- cd4.KLRG1.flowD@cell.count
              all.props[q2,42] <- cd4.KLRG1.flowD@proportion
              
              ####################################################################################################
              # Gating Tregs. Plotting CD62L_CD44 to obtain Tregs Resting/Naive and Tregs Effector
              # Tregs Resting/Naive
              Tregs.Resting.flowD <- flowDensity(Tregs.flowD, channels = c(8,7), position = c(T,NA), gates = c(gthres[15],NA))
              Tregs.Resting.flowD@proportion <- (Tregs.Resting.flowD@cell.count/Tregs.flowD@cell.count)*100
              # Indices in the Tregs cells to get the Tregs Resting/Naive Population with respect to f. Will be needed for the MFI calculation
              Tregs.Resting.flowD.ind <- Tregs.Resting.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the Tregs Resting/Naive Population population are being saved with respect to f
              # plotDens(Tregs, channels = c(8,7), main="Tregs Resting/Naive", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(Tregs.flowD@flow.frame@exprs[Tregs.Resting.flowD.ind,c(8,7)], col=2, cex = 0.1)
              
              all.events[q2,43] <- Tregs.Resting.flowD@cell.count
              all.props[q2,43] <- Tregs.Resting.flowD@proportion
              
              
              
              
              # Tregs Effector
              Tregs.Effector.flowD <- flowDensity(Tregs.flowD, channels = c(8,7), position = c(F,T), gates = c(gthres[15], gthres[16]))
              Tregs.Effector.flowD@proportion <- (Tregs.Effector.flowD@cell.count/Tregs.flowD@cell.count)*100
              # Indices in the Tregs cells to get the Tregs Effector Population with respect to f. Will be needed for the MFI calculation
              Tregs.Effector.flowD.ind <- Tregs.Effector.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the Tregs Effector Population population are being saved with respect to f
              # plotDens(Tregs, channels = c(8,7), main="Tregs Effector", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(Tregs.flowD@flow.frame@exprs[Tregs.Effector.flowD.ind,c(8,7)], col=2, cex = 0.1)
              
              all.events[q2,44] <- Tregs.Effector.flowD@cell.count
              all.props[q2,44] <- Tregs.Effector.flowD@proportion
              
                           
              
              ####################################################################################################
              # Gating Tregs. Plotting KLRG1_CD44 to obtain Treg KLRG1+
              Tregs.KLRG1.flowD <- flowDensity(Tregs.flowD, channels = c(12,7), position = c(T,T), gates = c(gthres[17], gthres[23]))
              Tregs.KLRG1.flowD@proportion <- (Tregs.KLRG1.flowD@cell.count/Tregs.flowD@cell.count)*100
              # Indices in the Tregs cells to get the Tregs KLRG1+ Population with respect to f. Will be needed for the MFI calculation
              Tregs.KLRG1.flowD.ind <- Tregs.KLRG1.flowD@index
              # # Optional Plot: Done to observe if the correct cells of the Tregs KLRG1+ Population population are being saved with respect to f
              # plotDens(Tregs, channels = c(12,7), main="Tregs KLRG1+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              # points(Tregs.flowD@flow.frame@exprs[Tregs.KLRG1.flowD.ind,c(12,7)], col=2, cex = 0.1)
              
              all.events[q2,45] <- Tregs.KLRG1.flowD@cell.count
              all.props[q2,45] <- Tregs.KLRG1.flowD@proportion
              
             
              
              
              
              ###############################################################################################################################
              ###############################################################################################################################
              
              ## Saving the MFIs (mean) of Original Data (Data before Transformed)
              singlets.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[singlets.flowD.l.ind, x], na.rm = T)})
              live.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[live.flowD.ind, x], na.rm = T)})
              lymph.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[lymph.flowD.ind, x], na.rm = T)})
              cd45.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd45.flowD.ind, x], na.rm = T)})
              autoFluor.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd45sTempb.ind, x], na.rm = T)})
              cd45.NOTP2.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd45.NOTP2.flowD.ind, x], na.rm = T)})
              gd.Tcells.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[gd.Tcells.flowD.ind, x], na.rm = T)})
              not.gd.Tcells.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[not.gd.Tcells.ind, x], na.rm = T)})
              gd.Resting.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[gd.Resting.flowD.ind, x], na.rm = T)})
              gd.Effector.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[gd.Effector.flowD.ind, x], na.rm = T)})
              gd.Naive.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[gd.Naive.flowD.ind, x], na.rm = T)})
              gd.klrg1.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[gd.klrg1.flowD.ind, x], na.rm = T)})
              gd.gitr.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[gd.gitr.flowD.ind, x], na.rm = T)})
              gd.cd5.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[gd.cd5.flowD.ind, x], na.rm = T)})
              cd5.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd5.flowD.ind, x], na.rm = T)})
              NKcells.MFI <-  sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[NKcells.flowD.ind, x], na.rm = T)})
              NK.Resting.Naive.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[NK.Resting.Naive.flowD.ind, x], na.rm = T)})
              NK.Effector.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[NK.Effector.flowD.ind, x], na.rm = T)}) 
              NK.klrg1.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[NK.klrg1.flowD.ind, x], na.rm = T)}) 
              P2a.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[P2a.flowD.ind, x], na.rm = T)}) 
              cd4neg.NKTcells.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd4neg.NKTcells.flowD.ind, x], na.rm = T)}) 
              cd4pos.NKTcells.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd4pos.NKTcells.flowD.ind, x], na.rm = T)}) 
              cd4neg.NKT.Effector.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd4neg.NKT.Effector.flowD.ind, x], na.rm = T)}) 
              cd4neg.NKT.Resting.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd4neg.NKT.Resting.flowD.ind, x], na.rm = T)}) 
              cd4neg.NKT.klrg1.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd4neg.NKT.klrg1.flowD.ind, x], na.rm = T)}) 
              cd4pos.NKT.Effector.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd4pos.NKT.Effector.flowD.ind, x], na.rm = T)}) 
              cd4pos.NKT.Resting.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd4pos.NKT.Resting.flowD.ind, x], na.rm = T)}) 
              cd4pos.NKT.klrg1.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd4pos.NKT.klrg1.flowD.ind, x], na.rm = T)}) 
              cd5.cd4cd8.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd5.cd4cd8.flowD.ind, x], na.rm = T)}) 
              abTcells.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[abTcells.flowD.ind, x], na.rm = T)}) 
              cd8a.Tcells.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd8a.Tcells.flowD.ind, x], na.rm = T)}) 
              cd4.Tcells.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd4.Tcells.flowD.ind, x], na.rm = T)}) 
              cd8.Naive.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd8.Naive.flowD.ind, x], na.rm = T)}) 
              cd8.Effector.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd8.Effector.flowD.ind, x], na.rm = T)}) 
              cd8.Resting.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd8.Resting.flowD.ind, x], na.rm = T)}) 
              cd8.KLRG1.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd8.KLRG1.flowD.ind, x], na.rm = T)}) 
              T.helper.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[T.helper.flowD.ind, x], na.rm = T)}) 
              Tregs.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[Tregs.flowD.ind, x], na.rm = T)})
              cd4.Resting.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd4.Resting.flowD.ind, x], na.rm = T)}) 
              cd4.Effector.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd4.Effector.flowD.ind, x], na.rm = T)}) 
              cd4.KLRG1.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd4.KLRG1.flowD.ind, x], na.rm = T)}) 
              Tregs.Resting.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[Tregs.Resting.flowD.ind, x], na.rm = T)}) 
              Tregs.Effector.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[Tregs.Effector.flowD.ind, x], na.rm = T)}) 
              Tregs.KLRG1.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[Tregs.KLRG1.flowD.ind, x], na.rm = T)}) 
              
              MFI.mean.OriginalData[q2, ] <- c(singlets.MFI, live.MFI, lymph.MFI, cd45.MFI, autoFluor.MFI, cd45.NOTP2.MFI, gd.Tcells.MFI, not.gd.Tcells.MFI,
                                              gd.Resting.MFI, gd.Effector.MFI, gd.Naive.MFI, gd.klrg1.MFI, gd.gitr.MFI, gd.cd5.MFI, cd5.MFI, NKcells.MFI, NK.Resting.Naive.MFI,
                                              NK.Effector.MFI, NK.klrg1.MFI, P2a.MFI, cd4neg.NKTcells.MFI, cd4pos.NKTcells.MFI, cd4neg.NKT.Effector.MFI,
                                              cd4neg.NKT.Resting.MFI, cd4neg.NKT.klrg1.MFI, cd4pos.NKT.Effector.MFI, cd4pos.NKT.Resting.MFI, cd4pos.NKT.klrg1.MFI,
                                              cd5.cd4cd8.MFI, abTcells.MFI, cd8a.Tcells.MFI, cd4.Tcells.MFI, cd8.Naive.MFI, cd8.Effector.MFI, cd8.Resting.MFI, cd8.KLRG1.MFI,
                                              T.helper.MFI, Tregs.MFI, cd4.Resting.MFI, cd4.Effector.MFI, cd4.KLRG1.MFI, Tregs.Resting.MFI, Tregs.Effector.MFI, Tregs.KLRG1.MFI)
              
              
              ## Saving the MFIs (median) of Original Data (Data before Transformed)
              singlets.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[singlets.flowD.l.ind, x], na.rm = T)})
              live.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[live.flowD.ind, x], na.rm = T)})
              lymph.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[lymph.flowD.ind, x], na.rm = T)})
              cd45.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd45.flowD.ind, x], na.rm = T)})
              autoFluor.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd45sTempb.ind, x], na.rm = T)})
              cd45.NOTP2.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd45.NOTP2.flowD.ind, x], na.rm = T)})
              gd.Tcells.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[gd.Tcells.flowD.ind, x], na.rm = T)})
              not.gd.Tcells.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[not.gd.Tcells.ind, x], na.rm = T)})
              gd.Resting.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[gd.Resting.flowD.ind, x], na.rm = T)})
              gd.Effector.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[gd.Effector.flowD.ind, x], na.rm = T)})
              gd.Naive.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[gd.Naive.flowD.ind, x], na.rm = T)})
              gd.klrg1.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[gd.klrg1.flowD.ind, x], na.rm = T)})
              gd.gitr.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[gd.gitr.flowD.ind, x], na.rm = T)})
              gd.cd5.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[gd.cd5.flowD.ind, x], na.rm = T)})
              cd5.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd5.flowD.ind, x], na.rm = T)})
              NKcells.MFI <-  sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[NKcells.flowD.ind, x], na.rm = T)})
              NK.Resting.Naive.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[NK.Resting.Naive.flowD.ind, x], na.rm = T)})
              NK.Effector.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[NK.Effector.flowD.ind, x], na.rm = T)}) 
              NK.klrg1.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[NK.klrg1.flowD.ind, x], na.rm = T)}) 
              P2a.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[P2a.flowD.ind, x], na.rm = T)}) 
              cd4neg.NKTcells.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd4neg.NKTcells.flowD.ind, x], na.rm = T)}) 
              cd4pos.NKTcells.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd4pos.NKTcells.flowD.ind, x], na.rm = T)}) 
              cd4neg.NKT.Effector.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd4neg.NKT.Effector.flowD.ind, x], na.rm = T)}) 
              cd4neg.NKT.Resting.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd4neg.NKT.Resting.flowD.ind, x], na.rm = T)}) 
              cd4neg.NKT.klrg1.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd4neg.NKT.klrg1.flowD.ind, x], na.rm = T)}) 
              cd4pos.NKT.Effector.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd4pos.NKT.Effector.flowD.ind, x], na.rm = T)}) 
              cd4pos.NKT.Resting.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd4pos.NKT.Resting.flowD.ind, x], na.rm = T)}) 
              cd4pos.NKT.klrg1.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd4pos.NKT.klrg1.flowD.ind, x], na.rm = T)}) 
              cd5.cd4cd8.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd5.cd4cd8.flowD.ind, x], na.rm = T)}) 
              abTcells.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[abTcells.flowD.ind, x], na.rm = T)}) 
              cd8a.Tcells.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd8a.Tcells.flowD.ind, x], na.rm = T)}) 
              cd4.Tcells.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd4.Tcells.flowD.ind, x], na.rm = T)}) 
              cd8.Naive.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd8.Naive.flowD.ind, x], na.rm = T)}) 
              cd8.Effector.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd8.Effector.flowD.ind, x], na.rm = T)}) 
              cd8.Resting.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd8.Resting.flowD.ind, x], na.rm = T)}) 
              cd8.KLRG1.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd8.KLRG1.flowD.ind, x], na.rm = T)}) 
              T.helper.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[T.helper.flowD.ind, x], na.rm = T)}) 
              Tregs.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[Tregs.flowD.ind, x], na.rm = T)})
              cd4.Resting.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd4.Resting.flowD.ind, x], na.rm = T)}) 
              cd4.Effector.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd4.Effector.flowD.ind, x], na.rm = T)}) 
              cd4.KLRG1.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd4.KLRG1.flowD.ind, x], na.rm = T)}) 
              Tregs.Resting.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[Tregs.Resting.flowD.ind, x], na.rm = T)}) 
              Tregs.Effector.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[Tregs.Effector.flowD.ind, x], na.rm = T)}) 
              Tregs.KLRG1.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[Tregs.KLRG1.flowD.ind, x], na.rm = T)}) 
              
              MFI.median.OriginalData[q2, ] <- c(singlets.MFI, live.MFI, lymph.MFI, cd45.MFI, autoFluor.MFI, cd45.NOTP2.MFI, gd.Tcells.MFI, not.gd.Tcells.MFI,
                                               gd.Resting.MFI, gd.Effector.MFI, gd.Naive.MFI, gd.klrg1.MFI, gd.gitr.MFI, gd.cd5.MFI, cd5.MFI, NKcells.MFI, NK.Resting.Naive.MFI,
                                               NK.Effector.MFI, NK.klrg1.MFI, P2a.MFI, cd4neg.NKTcells.MFI, cd4pos.NKTcells.MFI, cd4neg.NKT.Effector.MFI,
                                               cd4neg.NKT.Resting.MFI, cd4neg.NKT.klrg1.MFI, cd4pos.NKT.Effector.MFI, cd4pos.NKT.Resting.MFI, cd4pos.NKT.klrg1.MFI,
                                               cd5.cd4cd8.MFI, abTcells.MFI, cd8a.Tcells.MFI, cd4.Tcells.MFI, cd8.Naive.MFI, cd8.Effector.MFI, cd8.Resting.MFI, cd8.KLRG1.MFI,
                                               T.helper.MFI, Tregs.MFI, cd4.Resting.MFI, cd4.Effector.MFI, cd4.KLRG1.MFI, Tregs.Resting.MFI, Tregs.Effector.MFI, Tregs.KLRG1.MFI)


              
              ## Saving the MFIs (mean) of Transformed (Data after Transformation)
              singlets.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[singlets.flowD.l.ind, x], na.rm = T)})
              live.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[live.flowD.ind, x], na.rm = T)})
              lymph.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[lymph.flowD.ind, x], na.rm = T)})
              cd45.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd45.flowD.ind, x], na.rm = T)})
              autoFluor.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd45sTempb.ind, x], na.rm = T)})
              cd45.NOTP2.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd45.NOTP2.flowD.ind, x], na.rm = T)})
              gd.Tcells.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[gd.Tcells.flowD.ind, x], na.rm = T)})
              not.gd.Tcells.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[not.gd.Tcells.ind, x], na.rm = T)})
              gd.Resting.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[gd.Resting.flowD.ind, x], na.rm = T)})
              gd.Effector.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[gd.Effector.flowD.ind, x], na.rm = T)})
              gd.Naive.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[gd.Naive.flowD.ind, x], na.rm = T)})
              gd.klrg1.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[gd.klrg1.flowD.ind, x], na.rm = T)})
              gd.gitr.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[gd.gitr.flowD.ind, x], na.rm = T)})
              gd.cd5.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[gd.cd5.flowD.ind, x], na.rm = T)})
              cd5.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd5.flowD.ind, x], na.rm = T)})
              NKcells.MFI <-  sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[NKcells.flowD.ind, x], na.rm = T)})
              NK.Resting.Naive.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[NK.Resting.Naive.flowD.ind, x], na.rm = T)})
              NK.Effector.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[NK.Effector.flowD.ind, x], na.rm = T)}) 
              NK.klrg1.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[NK.klrg1.flowD.ind, x], na.rm = T)}) 
              P2a.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[P2a.flowD.ind, x], na.rm = T)}) 
              cd4neg.NKTcells.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd4neg.NKTcells.flowD.ind, x], na.rm = T)}) 
              cd4pos.NKTcells.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd4pos.NKTcells.flowD.ind, x], na.rm = T)}) 
              cd4neg.NKT.Effector.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd4neg.NKT.Effector.flowD.ind, x], na.rm = T)}) 
              cd4neg.NKT.Resting.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd4neg.NKT.Resting.flowD.ind, x], na.rm = T)}) 
              cd4neg.NKT.klrg1.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd4neg.NKT.klrg1.flowD.ind, x], na.rm = T)}) 
              cd4pos.NKT.Effector.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd4pos.NKT.Effector.flowD.ind, x], na.rm = T)}) 
              cd4pos.NKT.Resting.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd4pos.NKT.Resting.flowD.ind, x], na.rm = T)}) 
              cd4pos.NKT.klrg1.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd4pos.NKT.klrg1.flowD.ind, x], na.rm = T)}) 
              cd5.cd4cd8.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd5.cd4cd8.flowD.ind, x], na.rm = T)}) 
              abTcells.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[abTcells.flowD.ind, x], na.rm = T)}) 
              cd8a.Tcells.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd8a.Tcells.flowD.ind, x], na.rm = T)}) 
              cd4.Tcells.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd4.Tcells.flowD.ind, x], na.rm = T)}) 
              cd8.Naive.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd8.Naive.flowD.ind, x], na.rm = T)}) 
              cd8.Effector.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd8.Effector.flowD.ind, x], na.rm = T)}) 
              cd8.Resting.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd8.Resting.flowD.ind, x], na.rm = T)}) 
              cd8.KLRG1.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd8.KLRG1.flowD.ind, x], na.rm = T)}) 
              T.helper.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[T.helper.flowD.ind, x], na.rm = T)}) 
              Tregs.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[Tregs.flowD.ind, x], na.rm = T)})
              cd4.Resting.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd4.Resting.flowD.ind, x], na.rm = T)}) 
              cd4.Effector.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd4.Effector.flowD.ind, x], na.rm = T)}) 
              cd4.KLRG1.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd4.KLRG1.flowD.ind, x], na.rm = T)}) 
              Tregs.Resting.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[Tregs.Resting.flowD.ind, x], na.rm = T)}) 
              Tregs.Effector.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[Tregs.Effector.flowD.ind, x], na.rm = T)}) 
              Tregs.KLRG1.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[Tregs.KLRG1.flowD.ind, x], na.rm = T)}) 
              
              MFI.mean.TransformedData[q2, ] <- c(singlets.MFI, live.MFI, lymph.MFI, cd45.MFI, autoFluor.MFI, cd45.NOTP2.MFI, gd.Tcells.MFI, not.gd.Tcells.MFI,
                                               gd.Resting.MFI, gd.Effector.MFI, gd.Naive.MFI, gd.klrg1.MFI, gd.gitr.MFI, gd.cd5.MFI, cd5.MFI, NKcells.MFI, NK.Resting.Naive.MFI,
                                               NK.Effector.MFI, NK.klrg1.MFI, P2a.MFI, cd4neg.NKTcells.MFI, cd4pos.NKTcells.MFI, cd4neg.NKT.Effector.MFI,
                                               cd4neg.NKT.Resting.MFI, cd4neg.NKT.klrg1.MFI, cd4pos.NKT.Effector.MFI, cd4pos.NKT.Resting.MFI, cd4pos.NKT.klrg1.MFI,
                                               cd5.cd4cd8.MFI, abTcells.MFI, cd8a.Tcells.MFI, cd4.Tcells.MFI, cd8.Naive.MFI, cd8.Effector.MFI, cd8.Resting.MFI, cd8.KLRG1.MFI,
                                               T.helper.MFI, Tregs.MFI, cd4.Resting.MFI, cd4.Effector.MFI, cd4.KLRG1.MFI, Tregs.Resting.MFI, Tregs.Effector.MFI, Tregs.KLRG1.MFI)
              
              
              ## Saving the MFIs (median) of Transformed (Data after Transformation)
              singlets.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[singlets.flowD.l.ind, x], na.rm = T)})
              live.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[live.flowD.ind, x], na.rm = T)})
              lymph.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[lymph.flowD.ind, x], na.rm = T)})
              cd45.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd45.flowD.ind, x], na.rm = T)})
              autoFluor.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd45sTempb.ind, x], na.rm = T)})
              cd45.NOTP2.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd45.NOTP2.flowD.ind, x], na.rm = T)})
              gd.Tcells.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[gd.Tcells.flowD.ind, x], na.rm = T)})
              not.gd.Tcells.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[not.gd.Tcells.ind, x], na.rm = T)})
              gd.Resting.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[gd.Resting.flowD.ind, x], na.rm = T)})
              gd.Effector.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[gd.Effector.flowD.ind, x], na.rm = T)})
              gd.Naive.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[gd.Naive.flowD.ind, x], na.rm = T)})
              gd.klrg1.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[gd.klrg1.flowD.ind, x], na.rm = T)})
              gd.gitr.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[gd.gitr.flowD.ind, x], na.rm = T)})
              gd.cd5.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[gd.cd5.flowD.ind, x], na.rm = T)})
              cd5.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd5.flowD.ind, x], na.rm = T)})
              NKcells.MFI <-  sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[NKcells.flowD.ind, x], na.rm = T)})
              NK.Resting.Naive.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[NK.Resting.Naive.flowD.ind, x], na.rm = T)})
              NK.Effector.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[NK.Effector.flowD.ind, x], na.rm = T)}) 
              NK.klrg1.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[NK.klrg1.flowD.ind, x], na.rm = T)}) 
              P2a.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[P2a.flowD.ind, x], na.rm = T)}) 
              cd4neg.NKTcells.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd4neg.NKTcells.flowD.ind, x], na.rm = T)}) 
              cd4pos.NKTcells.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd4pos.NKTcells.flowD.ind, x], na.rm = T)}) 
              cd4neg.NKT.Effector.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd4neg.NKT.Effector.flowD.ind, x], na.rm = T)}) 
              cd4neg.NKT.Resting.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd4neg.NKT.Resting.flowD.ind, x], na.rm = T)}) 
              cd4neg.NKT.klrg1.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd4neg.NKT.klrg1.flowD.ind, x], na.rm = T)}) 
              cd4pos.NKT.Effector.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd4pos.NKT.Effector.flowD.ind, x], na.rm = T)}) 
              cd4pos.NKT.Resting.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd4pos.NKT.Resting.flowD.ind, x], na.rm = T)}) 
              cd4pos.NKT.klrg1.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd4pos.NKT.klrg1.flowD.ind, x], na.rm = T)}) 
              cd5.cd4cd8.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd5.cd4cd8.flowD.ind, x], na.rm = T)}) 
              abTcells.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[abTcells.flowD.ind, x], na.rm = T)}) 
              cd8a.Tcells.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd8a.Tcells.flowD.ind, x], na.rm = T)}) 
              cd4.Tcells.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd4.Tcells.flowD.ind, x], na.rm = T)}) 
              cd8.Naive.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd8.Naive.flowD.ind, x], na.rm = T)}) 
              cd8.Effector.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd8.Effector.flowD.ind, x], na.rm = T)}) 
              cd8.Resting.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd8.Resting.flowD.ind, x], na.rm = T)}) 
              cd8.KLRG1.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd8.KLRG1.flowD.ind, x], na.rm = T)}) 
              T.helper.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[T.helper.flowD.ind, x], na.rm = T)}) 
              Tregs.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[Tregs.flowD.ind, x], na.rm = T)})
              cd4.Resting.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd4.Resting.flowD.ind, x], na.rm = T)}) 
              cd4.Effector.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd4.Effector.flowD.ind, x], na.rm = T)}) 
              cd4.KLRG1.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd4.KLRG1.flowD.ind, x], na.rm = T)}) 
              Tregs.Resting.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[Tregs.Resting.flowD.ind, x], na.rm = T)}) 
              Tregs.Effector.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[Tregs.Effector.flowD.ind, x], na.rm = T)}) 
              Tregs.KLRG1.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[Tregs.KLRG1.flowD.ind, x], na.rm = T)}) 
              
              MFI.median.TransformedData[q2, ] <- c(singlets.MFI, live.MFI, lymph.MFI, cd45.MFI, autoFluor.MFI, cd45.NOTP2.MFI, gd.Tcells.MFI, not.gd.Tcells.MFI,
                                                  gd.Resting.MFI, gd.Effector.MFI, gd.Naive.MFI, gd.klrg1.MFI, gd.gitr.MFI, gd.cd5.MFI, cd5.MFI, NKcells.MFI, NK.Resting.Naive.MFI,
                                                  NK.Effector.MFI, NK.klrg1.MFI, P2a.MFI, cd4neg.NKTcells.MFI, cd4pos.NKTcells.MFI, cd4neg.NKT.Effector.MFI,
                                                  cd4neg.NKT.Resting.MFI, cd4neg.NKT.klrg1.MFI, cd4pos.NKT.Effector.MFI, cd4pos.NKT.Resting.MFI, cd4pos.NKT.klrg1.MFI,
                                                  cd5.cd4cd8.MFI, abTcells.MFI, cd8a.Tcells.MFI, cd4.Tcells.MFI, cd8.Naive.MFI, cd8.Effector.MFI, cd8.Resting.MFI, cd8.KLRG1.MFI,
                                                  T.helper.MFI, Tregs.MFI, cd4.Resting.MFI, cd4.Effector.MFI, cd4.KLRG1.MFI, Tregs.Resting.MFI, Tregs.Effector.MFI, Tregs.KLRG1.MFI)
              
              
              ####################################################################################################
              ####################################################################################################
              
             
              ## Saving the Plots
              #--------Start Big Png------------
              png ( file = paste("Results/Figures/ScatterPlots_Updated/", store.allFCS[q2,1], "/", "Total_", store.allFCS[q2,2], ".png", sep = "" ), width=2100, height=2100*7/6)
              par(mfrow=c(7,6),mar=(c(5, 5, 4, 2) + 0.1))
              
              # 1st method of plotting FSC-A_SSC-W
              plotDens(f,c("FSC-A","SSC-W"),main="Ungated", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              lines(singlets.flowD.l@filter,lwd=2)
              
              # Plotting Live/Dead_SSC-A
              plotDens(singlets,  c(11,4), main= "Singlets", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=gthres[3], lwd=2);
              
              # Plotting FSC-A_SSC-A
              plotDens(live, c(1,4), main = "Live", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=gthres[5], lwd=2); abline(v=gthres[4], lwd=2); abline(h=gthres[6], lwd=2)
              
              # Plotting CD45_CD161
              plotDens(lymph, channels = c(14,15), main = "Lymphocytes", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=gthres[7], lwd=2); abline(v=gthres[8], lwd=2); abline(h=gthres[21], lwd=2)
              
              # Plotting TCRd_CD4. Highlighting the Autofluoresence part
              plotDens(cd45, channels = c(18,16), main="CD45+_Autofluorecence", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              if ( nrow(cd45autoflour@exprs) != 0) {
                points(TCRDCD4points, type="l", col="black", lty=2, lwd=2)
              }
              
              # Plotting TCRd_CD4. Getting rid of the Autofluoresence part
              plotDens(cd45.NOTP2,channels = c(18,16), main="CD45+_Without Autofluorecence", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              
              
              ## Plotting TCRd_CD4 to obtain gd T-cells
              plotDens(cd45.NOTP2, channels = c(18,16), main = "CD45+_gd T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2);
              lines(gd.Tcells.flowD@filter, lwd=2)
              
              ## Plotting CD62L_CD44 to obtain gd Resting, gd Effector, and gd Naive
              min.x <- min(exprs(gd.Tcells)[,8])-1
              max.x <- max(exprs(gd.Tcells)[,8])+1
              min.y <- min(exprs(gd.Tcells)[,7])-1
              max.y <- max(exprs(gd.Tcells)[,7])+1
              plotDens(gd.Tcells, c(8,7), main = "gd T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x,max.x), ylim = c(min.y, max.y)); abline(v=gthres[15], lwd=2); abline(h=gthres[16], lwd=2)
              
              ## Plotting KLRG1_GITR to obtain GITR gd T-cells and GD KLRG1+
              min.x <- min(exprs(gd.Tcells)[,12])
              max.x <- max(exprs(gd.Tcells)[,12])+1
              plotDens(gd.Tcells, c(12,17), main = "gd T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim =c(min.x, max.x)); abline(v=gthres[17], lwd=2); abline(h=gthres[18], lwd=2)
              
              ## Plotting CD5_CD44 to obtain GD CD5+
              plotDens(gd.Tcells, c(13,7), main = "gd T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=gthres[19], lwd=2)
              
              plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
              plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
              
              ## Plotting TCRd_CD4 to obtain NOT(P2) (without the gd T-cells part)
              plotDens(not.gd.Tcells, channels = c(18,16), main = "NOT(P2)", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              
              ## Plotting CD161_CD5 to obtain CD5+ and NK-cells
              plotDens(not.gd.Tcells, c(15,13), main = "NOT(gd T-cells)", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(h=gthres[19], lwd=2)
              #lines(x=c(cd161.gate.low, cd5.gate), y=c(cd161.gate.high, cd5.gate), lwd=2)
              lines(NKcells.flowD@filter, lwd=2)
              
              #Plotting CD62L_CD44 to obtain NK Effector and NK Resting/Naive
              min.x <- min(exprs(NKcells)[,8])-1
              min.y <- min(exprs(NKcells)[,7])-1
              max.y <-  max(exprs(NKcells)[,7])
              plotDens(NKcells, c(8,7), main = "NK-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, ylim =c(min.y, max.y)); abline(v=gthres[15], lwd=2)
              lines(x=c(min.x, gthres[15]), y=c(gthres[16],gthres[16]), lwd=2)
              
              # Plotting KLRG1_CD44 to obtain NK +KLRG1
              min.x <- min(exprs(NKcells)[,12])
              max.x <- max(exprs(NKcells)[,12])+1
              min.y <- min(exprs(NKcells)[,7])-1
              max.y <-  max(exprs(NKcells)[,7])+0.5
              plotDens(NKcells, c(12,7), main = "NK-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x,max.x), ylim = c(min.y,max.y))
              lines(x=c(gthres[17], gthres[17]), y=c(gthres[16], max.y+2), lwd=2)
              lines(x=c(gthres[17], max.x+1), y=c(gthres[16], gthres[16]), lwd=2)
              
              plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
              plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
              plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
              
              # Plotting CD161_CD4 to obtain P2a, CD4- NKT-cells, and CD4+ NKT-cells
              min.x <- min(exprs(cd5)[,15])
              max.x <- max(exprs(cd5)[,15])+1
              min.y <- min(exprs(cd5)[,16])-0.5
              max.y <-  max(exprs(cd5)[,16])+0.5
              plotDens(cd5, c(15,16), main = "CD5+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x,max.x), ylim = c(min.y,max.y)); abline(v=gthres[20], lwd=2)
              lines(x=c(gthres[20], max.x+1), y=c(gthres[14],gthres[14]), lwd=2)
              
              # Plotting CD62L_CD44 to obtain CD4- NKT Effector and CD4- NKT Resting/Naive
              min.x <- min(exprs(cd4neg.NKTcells)[,8])-1
              max.x <- max(exprs(cd4neg.NKTcells)[,8])+1
              min.y <- min(exprs(cd4neg.NKTcells)[,7])-2
              max.y <-  max(exprs(cd4neg.NKTcells)[,7])
              plotDens(cd4neg.NKTcells, c(8,7), main = "CD4- NKT-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x,max.x), ylim = c(min.y,max.y)); abline(v=gthres[15], lwd=2)
              lines(x=c(min.x-1,gthres[15]), y=c(gthres[16],gthres[16]), lwd=2)
              
              # Plotting KLRG1_CD44 to obtain CD4- NKT KLRG1+
              min.x <- min(exprs(cd4neg.NKTcells)[,12])
              max.x <- max(exprs(cd4neg.NKTcells)[,12])+1
              min.y <- min(exprs(cd4neg.NKTcells)[,7])-2
              max.y <-  max(exprs(cd4neg.NKTcells)[,7])
              plotDens(cd4neg.NKTcells, c(12,7), main = "CD4- NKT-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x,max.x), ylim = c(min.y,max.y))
              lines(x=c(gthres[17],max.x+1), y=c(gthres[16],gthres[16]),lwd=2)
              lines(x=c(gthres[17], gthres[17]), y=c(gthres[16], max.y+1),lwd=2)
              
              # Plotting CD62L_CD44 to obtain CD4+ NKT Effector and CD4+ NKT Resting/Naive
              min.x <- min(exprs(cd4pos.NKTcells)[,8])-1
              max.x <- max(exprs(cd4pos.NKTcells)[,8])+1
              min.y <- min(exprs(cd4pos.NKTcells)[,7])-2
              max.y <-  max(exprs(cd4pos.NKTcells)[,7])
              plotDens(cd4pos.NKTcells, c(8,7), main = "CD4+ NKT-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x,max.x), ylim = c(min.y,max.y)); abline(v=gthres[15], lwd=2)
              lines(x=c(min.x-1,gthres[15]), y=c(gthres[16],gthres[16]), lwd=2)
              
              # Plotting KLRG1_CD44 to obtain CD4+ NKT KLRG1+
              min.x <- min(exprs(cd4pos.NKTcells)[,12])
              max.x <- max(exprs(cd4pos.NKTcells)[,12])+1
              min.y <- min(exprs(cd4pos.NKTcells)[,7])-2
              max.y <-  max(exprs(cd4pos.NKTcells)[,7])
              plotDens(cd4pos.NKTcells, c(12,7), main = "CD4+ NKT-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x,max.x), ylim = c(min.y,max.y))
              lines(x=c(gthres[17],max.x+1), y=c(gthres[16],gthres[16]),lwd=2)
              lines(x=c(gthres[17], gthres[17]), y=c(gthres[16], max.y+1),lwd=2)
              
              plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
              
              # Plotting CD8a_CD4 to obtain CD5+ CD4/CD8 and abT-cell.
              min.x <- min(exprs(P2a)[,10])-1
              max.x <- max(exprs(P2a)[,10])
              min.y <- min(exprs(P2a)[,16])-1
              max.y <-  max(exprs(P2a)[,16])
              plotDens(P2a, c(10,16), main = "P2a", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              lines(x=c(gthres[22],gthres[22]), y=c(min.y,gthres[14]), lwd=2)
              lines(x=c(min.x,gthres[22]), y=c(gthres[14], gthres[14]), lwd=2)
              
              # Plotting CD8a_CD4 to obtain CD8a+ T-cells and CD4+ T-cells
              plotDens(abTcells, c(10,16), main = "abT-CELL", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(h=gthres[14], lwd=2); abline(v=gthres[22],lwd=2)
              
              # Plotting CD62L_CD44 to obtain CD8 Naive, CD8 Effector, and CD8 Resting
              min.x <- min(exprs(cd8a.Tcells)[,8])-1
              max.x <- max(exprs(cd8a.Tcells)[,8])+1
              min.y <- min(exprs(cd8a.Tcells)[,7])
              max.y <-  max(exprs(cd8a.Tcells)[,7])
              plotDens(cd8a.Tcells, c(8,7), main = "CD8a+ T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=gthres[15], lwd=2)
              lines(x=c(min.x,gthres[15]), y=c(gthres[16],gthres[16]), lwd=2)
              lines(x=c(gthres[15],max.x), y=c(gthres[23],gthres[23]), lwd=2)
              
              # Plotting KLRG1_CD44 to obtain CD8 KLRG1+
              min.x <- min(exprs(cd8a.Tcells)[,12])
              max.x <- max(exprs(cd8a.Tcells)[,12])+1
              min.y <- min(exprs(cd8a.Tcells)[,7])-1
              max.y <-  max(exprs(cd8a.Tcells)[,7])
              plotDens(cd8a.Tcells, c(12,7), main = "CD8a+ T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim=c(min.x,max.x), ylim =c(min.y,max.y)); abline(v=gthres[17], lwd=2)
              lines(x=c(gthres[17],max.x+1), y=c(gthres[16],gthres[16]), lwd=2)
              
              plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
              plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
              plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
              
              # Plotting CD25_GITR to obtain T-helper cells and Tregs
              min.x <- min(exprs(cd4.Tcells)[,9])
              max.x <- max(exprs(cd4.Tcells)[,9])+1
              min.y <- min(exprs(cd4.Tcells)[,17])
              max.y <-  max(exprs(cd4.Tcells)[,17])+1
              plotDens(cd4.Tcells, c(9,17), main = "CD4+ T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim=c(min.x,max.x), ylim =c(min.y,max.y)); abline(v=gthres[24], lwd=2); abline(v=gthres[25], lwd=2)
              lines(x=c(gthres[25],max.x+1), y=c(gthres[18],gthres[18]), lwd=2)
              
              # Plotting CD62L_CD44 to obtain CD4 Effector and CD4 Resting/Naive
              min.x <- min(exprs(T.helper)[,8])
              max.x <- max(exprs(T.helper)[,8])+0.5
              min.y <- min(exprs(T.helper)[,7])
              max.y <-  max(exprs(T.helper)[,7])
              plotDens(T.helper, c(8,7), main = "T-helper cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim=c(min.x,max.x), ylim =c(min.y,max.y)); abline(v=gthres[15], lwd=2)
              lines(x=c(min.x-1, gthres[15]), y=c(gthres[16],gthres[16]), lwd=2)
              
              # Plotting KLRG1_CD44 to obtain CD4+ KLRG1+ T-cells
              min.x <- min(exprs(T.helper)[,12])
              max.x <- max(exprs(T.helper)[,12])
              min.y <- min(exprs(T.helper)[,7])
              max.y <-  max(exprs(T.helper)[,7])
              plotDens(T.helper, c(12,7), main = "T-helper cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              lines(x=c(gthres[17], max.x+1), y=c(gthres[23], gthres[23]), lwd=2)
              lines(x=c(gthres[17], gthres[17]), y=c(max.y+1, gthres[23]),lwd=2)
              
              plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
              plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
              plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
              plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
              
              #Plotting CD62L_CD44 to obtain Tregs Resting/Naive and Tregs Effector
              min.x <- min(exprs(Tregs)[,8])-1
              max.x <- max(exprs(Tregs)[,8])+1
              min.y <- min(exprs(Tregs)[,7])-1
              max.y <-  max(exprs(Tregs)[,7])
              plotDens(Tregs, c(8,7), main = "Tregs", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim=c(min.x,max.x), ylim =c(min.y,max.y)); abline(v= gthres[15], lwd=2)
              lines(x=c(min.x-1, gthres[15]), y=c(gthres[16],gthres[16]), lwd=2)
              
              # Plotting KLRG1_CD44 to obtain Tregs KLRG1+
              min.x <- min(exprs(Tregs)[,12])
              max.x <- max(exprs(Tregs)[,12])
              min.y <- min(exprs(Tregs)[,7])
              max.y <-  max(exprs(Tregs)[,7])
              plotDens(Tregs, c(12,7), main = "Tregs", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
              lines(x=c(gthres[17], max.x+1), y=c(gthres[23], gthres[23]), lwd=2)
              lines(x=c(gthres[17], gthres[17]), y=c(max.y+1, gthres[23]),lwd=2)
              
              dev.off()
              par(mfrow=c(1,1))
  
  
    },error = function(err) {
            print(paste("ERROR in Gating/Plotting Message:  ",err, "Index: ", q2));
            #next
            errorFileIndex <- c(errorFileIndex,q2)
            return(errorFileIndex)
  }
) # end of tryCatch


    if (Do_flowType == T) {
          print("Start FlowType")
          # ## flowType on CD45 population. Live marker is excluded but CD45 marker is included. Also the maximum number of markers allowed is 6.
          # flowType.res <- flowType(Frame = cd45.NOTP2 , PropMarkers= as.vector(channels.ind.NoLive), MaxMarkersPerPop = 6, PartitionsPerMarker=2,
          #                         Methods='Thresholds', Thresholds= listGthres.NoLive, verbose=F, MemLimit=400)
          # save ( flowType.res, file =  paste("Results/FlowType/", store.allFCS[q2,1],"/FT_", store.allFCS[q2,2],".Rdata",sep="") )
          
          ## flowType on CD45 population. Live marker and CD45 markers are included. Also the maximum number of markers allowed is 8.
          ## For FlowType we start with the CD45 population without the Autofluorecence part
          flowType.res <- flowType(Frame = cd45.NOTP2 , PropMarkers= as.vector(channels.ind.NoLiveNoCD45), MaxMarkersPerPop = 8, PartitionsPerMarker=2,
                                   Methods='Thresholds', Thresholds= listGthres.NoLiveNoCD45, verbose=F, MemLimit=400)
          
          save ( flowType.res, file =  paste("Results/FlowType/", store.allFCS[q2,1],"/FT_", store.allFCS[q2,2],".Rdata",sep="") )
          
      }
        cat("Time is: ",TimeOutput(start2),"\n",sep="")
}# end of for-loop

colnames(all.props) <- c("All Events-%Parent", "Singlets-%Parent", "Live-%Parent", "Lymphocytes-%Parent", "CD45-%Parent", "Autofluorecence-%Parent", "NOT(P2)-%Parent", "gd T-cells-%Parent", "NOT(gd T-cells)-%Parent",
                             "GD Resting-%Parent","GD Effector-%Parent", "GD-Naive-%Parent", "GD KLRG1+-%Parent", "GITR GD T-cells-%Parent", "GD CD5+-%Parent", "CD5+-%Parent", "NK-cells-%Parent", "NK Resting/Naive-%Parent",
                             "NK Effector-%Parent", "NK KLRG1-%Parent", "P2a-%Parent", "CD4- NKT-cells-%Parent", "CD4+ NKT-cells-%Parent", "CD4- NKT Effector-%Parent", "CD4- NKT Resting-%Parent", "CD4- NKT KLRG1+-%Parent",
                             "CD4+ NKT Effector-%Parent", "CD4+ NKT Resting-%Parent", "CD4+ KLRG1+-%Parent", "CD5+ CD4/CD8-%Parent", "ab T-cell-%Parent", "CD8a+ T-cells-%Parent", "CD4+ T-cells-%Parent",
                             "CD8 Naive-%Parent", "CD8 Effector-%Parent", "CD8 Resting-%Parent", "Cd8 KLRG1-%Parent", "T-helper cells-%Parent", "Tregs-%Parent", "CD4 Resting-%Parent", "CD4 Effector-%Parent",
                             "CD4 KLRG1-%Parent", "Tregs Resting-%Parent", "Tregs Effector-%Parent", "Tregs KLRG1-%Parent")
rownames(all.props) <- 1:length(all.props[,1])

colnames(all.events) <- c("All Events-#Events", "Singlets-#Events", "Live-#Events", "Lymphocytes-#Events", "CD45-#Events", "Autofluorecence-#Events", "NOT(P2)-#Events", "gd T-cells-#Events", "NOT(gd T-cells)-#Events",
                              "GD Resting-#Events","GD Effector-#Events", "GD-Naive-#Events", "GD KLRG1+-#Events", "GITR GD T-cells-#Events", "GD CD5+-#Events", "CD5+-#Events", "NK-cells-#Events", "NK Resting/Naive-#Events",
                              "NK Effector-#Events", "NK KLRG1-#Events", "P2a-#Events", "CD4- NKT-cells-#Events", "CD4+ NKT-cells-#Events", "CD4- NKT Effector-#Events", "CD4- NKT Resting-#Events", "CD4- NKT KLRG1+-#Events",
                              "CD4+ NKT Effector-#Events", "CD4+ NKT Resting-#Events", "CD4+ KLRG1+-#Events", "CD5+ CD4/CD8-#Events", "ab T-cell-#Events", "CD8a+ T-cells-#Events", "CD4+ T-cells-#Events",
                              "CD8 Naive-#Events", "CD8 Effector-#Events", "CD8 Resting-#Events", "Cd8 KLRG1-#Events", "T-helper cells-#Events", "Tregs-#Events", "CD4 Resting-#Events", "CD4 Effector-#Events",
                              "CD4 KLRG1-#Events", "Tregs Resting-#Events", "Tregs Effector-#Events", "Tregs KLRG1-#Events")
rownames(all.events) <- 1:length(all.events[,1])

Events_Proportions_Table <- NULL
for(i in 1:ncol(all.props)){
  Events_Proportions_Table <- cbind(Events_Proportions_Table,all.events[,i],all.props[,i])
}
colnames(Events_Proportions_Table) <- c("All Events-#Events", "All Events-%Parent", "Singlets-#Events", "Singlets-%Parent", "Live-#Events", "Live-%Parent", "Lymphocytes-#Events", "Lymphocytes-%Parent", "CD45-#Events", "CD45-%Parent",
                                        "Autofluorecence-#Events", "Autofluorecence-%Parent", "NOT(P2)-#Events", "NOT(P2)-%Parent", "gd T-cells-#Events", "gd T-cells-%Parent", "NOT(gd T-cells)-#Events", "NOT(gd T-cells)-%Parent",
                                        "GD Resting-#Events", "GD Resting-%Parent", "GD Effector-#Events", "GD Effector-%Parent", "GD-Naive-#Events", "GD-Naive-%Parent", "GD KLRG1+-#Events", "GD KLRG1+-%Parent", "GITR GD T-cells-#Events", "GITR GD T-cells-%Parent",
                                        "GD CD5+-#Events", "GD CD5+-%Parent", "CD5+-#Events", "CD5+-%Parent", "NK-cells-#Events", "NK-cells-%Parent", "NK Resting/Naive-#Events", "NK Resting/Naive-%Parent",
                                        "NK Effector-#Events", "NK Effector-%Parent", "NK KLRG1-#Events", "NK KLRG1-%Parent", "P2a-#Events", "P2a-%Parent", "CD4- NKT-cells-#Events", "CD4- NKT-cells-%Parent", "CD4+ NKT-cells-#Events", "CD4+ NKT-cells-%Parent",
                                        "CD4- NKT Effector-#Events", "CD4- NKT Effector-%Parent", "CD4- NKT Resting-#Events", "CD4- NKT Resting-%Parent", "CD4- NKT KLRG1+-#Events", "CD4- NKT KLRG1+-%Parent",
                                        "CD4+ NKT Effector-#Events", "CD4+ NKT Effector-%Parent", "CD4+ NKT Resting-#Events", "CD4+ NKT Resting-%Parent", "CD4+ KLRG1+-#Events", "CD4+ KLRG1+-%Parent", "CD5+ CD4/CD8-#Events", "CD5+ CD4/CD8-%Parent",
                                        "ab T-cell-#Events", "ab T-cell-%Parent", "CD8a+ T-cells-#Events", "CD8a+ T-cells-%Parent", "CD4+ T-cells-#Events", "CD4+ T-cells-%Parent",
                                        "CD8 Naive-#Events", "CD8 Naive-%Parent", "CD8 Effector-#Events", "CD8 Effector-%Parent", "CD8 Resting-#Events", "CD8 Resting-%Parent", "Cd8 KLRG1-#Events", "Cd8 KLRG1-%Parent", "T-helper cells-#Events", "T-helper cells-%Parent",
                                        "Tregs-#Events", "Tregs-%Parent", "CD4 Resting-#Events", "CD4 Resting-%Parent", "CD4 Effector-#Events", "CD4 Effector-%Parent",
                                        "CD4 KLRG1-#Events", "CD4 KLRG1-%Parent", "Tregs Resting-#Events", "Tregs Resting-%Parent", "Tregs Effector-#Events", "Tregs Effector-%Parent", "Tregs KLRG1-#Events", "Tregs KLRG1-%Parent")
rownames(Events_Proportions_Table) <- 1:length(Events_Proportions_Table[,1])

# Combining the FCS files names, Genotype, The Number of cells before cleaning, number of channels, and the Assay Date with the Event_Proportion Table
Events_Proportions_Table <- cbind(store.allFCS, Events_Proportions_Table)


## Column names for all the Matrices storing the MFIs
colnames(MFI.mean.OriginalData) <- c(sapply(1:(ncol(f@exprs)-1), function(x){paste("Singlets:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("Live:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("Lymphocytes:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("CD45+:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("Autofluorecence:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("CD45.NOT(P2):", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("gd T-cells:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("NOT(gd T-cells):", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("gd Resting:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("gd Effector:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("gd Naive:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("gd KLRG1+:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("gd GITR:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("gd CD5+:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("CD5:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("NK-cells:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("NK Resting/Naive:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("NK Effector:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("NK KLRG1:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("CD5+ P2a:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("CD5+ CD4- NKT-cells:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("CD5+ CD4+ NKT-cells:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("CD4- NKT-cells Effector:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("CD4- NKT-cells Resting/Naive:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("CD4- NKT-cells KLRG1+: KLRG1:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("CD4+ NKT-cells Effector:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("CD4+ NKT-cells Resting/Naive:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("CD4+ NKT-cells KLRG1+:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("P2a CD5+ CD4CD8:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("P2a abT-cell:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("abT-cell CD8a+:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("abT-cell CD4+:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("CD8a+ Naive:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("CD8a+ Effector:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("CD8a+ Resting:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("CD8a+ KLRG1+:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("CD4+ T-helper:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("CD4+ Tregs:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("T-helper CD4 Resting/Naive:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("T-helper CD4 Effector:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("T-helper CD4 KLRG1+:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("Tregs Resting/Naive:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("Tregs Effector:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                  sapply(1:(ncol(f@exprs)-1), function(x){paste("Tregs KLRG1+:", sub(" .*", "", f@parameters@data$desc[[x]]))}))

colnames(MFI.mean.TransformedData) <- colnames(MFI.median.TransformedData) <- colnames(MFI.median.OriginalData) <- colnames(MFI.mean.OriginalData)


# Combining the FCS files names, Genotype, The Number of cells before cleaning, number of channels, and the Assay Date with the MFI.mean.OriginalData
MFI.mean.OriginalData <- cbind(store.allFCS, MFI.mean.OriginalData)

# Combining the FCS files names, Genotype, The Number of cells before cleaning, number of channels, and the Assay Date with the MFI.median.OriginalData
MFI.median.OriginalData <- cbind(store.allFCS, MFI.median.OriginalData)

# Combining the FCS files names, Genotype, The Number of cells before cleaning, number of channels, and the Assay Date with the MFI.mean.TransformedData
MFI.mean.TransformedData <- cbind(store.allFCS, MFI.mean.TransformedData)

# Combining the FCS files names, Genotype, The Number of cells before cleaning, number of channels, and the Assay Date with the MFI.median.TransformedData
MFI.median.TransformedData <- cbind(store.allFCS, MFI.median.TransformedData)


write.csv(Events_Proportions_Table, file =  paste("Results/Events_Proportions_Table.csv",sep=""))
write.csv(all.props, file =  paste("Results/allProportions_Table.csv",sep=""))
write.csv(all.events, file =  paste("Results/allEvents_Table.csv",sep=""))

write.csv(MFI.mean.OriginalData, file =  paste("Results/MFI_mean_OriginalData.csv",sep=""))
write.csv(MFI.median.OriginalData, file =  paste("Results/MFI_median_OriginalData.csv",sep=""))
write.csv(MFI.mean.TransformedData, file =  paste("Results/MFI_mean_TransformedData.csv",sep=""))
write.csv(MFI.median.TransformedData, file =  paste("Results/MFI_median_TransformedData.csv",sep=""))

save(channels.ind.NoLive, file = paste("Results/channels.ind.NoLive.Rdata", sep = ""))
save(channels.ind.NoLiveNoCD45, file = paste("Results/channels.ind.NoLiveNoCD45.Rdata", sep = ""))

cat("Total time is: ",TimeOutput(start),"\n",sep="")

