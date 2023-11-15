# Developed by Albina Rahim
# Date: June 09, 2016

## Re-gating and FlowType of Bone Marrow Panel cells
remove(list=ls())

Do_flowType <- T #For doing flowType at the end of the code

setwd("/code/Projects/3i/Panel_BM-cell/")

library("flowCore")
library("flowBin")
library("flowDensity")
library("flowType")
library("stringr")


source("Codes/3iTcellfunctions-BM.R")
source("Codes/rotate.data-BM.R")
source("Codes/getPeaks.R")

suppressWarnings ( dir.create ( "Results/FlowType") )
suppressWarnings ( dir.create ( "Results/Figures/ScatterPlots_Updated") )

load("Results/uniqueGT.Rdata")
load("Results/channels.ind.Rdata")
load("Results/store.allFCS.Rdata")
#load("Results/channels.ind.Rdata")

verbose_debris <- T
verbose_margin <- T

start <- Sys.time()

all.props <- matrix(nrow = nrow(store.allFCS), ncol = 21)
all.events <- matrix(nrow = nrow(store.allFCS), ncol = 21)
MFI.mean.OriginalData <- matrix(nrow =nrow(store.allFCS), ncol = 360)
MFI.median.OriginalData <- matrix(nrow = nrow(store.allFCS), ncol = 360)
MFI.mean.TransformedData <- matrix(nrow = nrow(store.allFCS), ncol = 360)
MFI.median.TransformedData <- matrix(nrow = nrow(store.allFCS), ncol = 360)


errorFileIndex <- NULL

# path to the FCS_Groups files
pathFCS <- paste("/code/Projects/3i/Panel_BM-cell/FCS_Groups")

# Reads all folders and files in current path folder and makes a list of all of their paths
allFCS <- dir(pathFCS, full.names=T, recursive=T) 

# Create directories
invisible(sapply(1:length(uniqueGT), function(x){
  suppressWarnings(dir.create ( paste("Results/FlowType/", uniqueGT[x], sep="") ))
  suppressWarnings(dir.create ( paste("Results/Figures/ScatterPlots_Updated/", uniqueGT[x], sep="") ))}))

load(file =  paste("Results/After_Clean/", store.allFCS[1,1],"/AftClean_", store.allFCS[1,2],".Rdata",sep="") )

## We will not include the Live marker while doing the FlowType
markers<- c("IgD","CD43","CD24","GR1","IgM","CD11b","CD138","CD3","BP1","B220", "CD45")
channels.ind.NoLive <- Find.markers(frame=f,marker.list=markers)
channels.ind.NoLive <- sort(channels.ind.NoLive)

## We will not include the Live marker and CD45 marker while doing the FlowType
markers<- c("IgD","CD43","CD24","GR1","IgM","CD11b","CD138","CD3","BP1","B220")
channels.ind.NoLiveNoCD45 <- Find.markers(frame=f,marker.list=markers)
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
        
        load (file =  paste("Results/After_Clean/", store.allFCS[q2,1],"/AftClean_", store.allFCS[q2,2],".Rdata",sep="") )
        load (file =  paste("Results/Gating_Thresholds_Updated/", store.allFCS[q2,1],"/Gthres_", store.allFCS[q2,2],".Rdata",sep="") ) 
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
  
  ## Gating All Events -----
  ## Gating All Events to obtain the Singlets. Plotting SSC-A_FSC-W
  singlets.flowD.h <- flowDensity(f, channels = c("SSC-A", "FSC-W"), position = c(NA, F), gates = c(NA, gthres[1]))
  singlets.flowD.l <- flowDensity(singlets.flowD.h, channels = c("SSC-A", "FSC-W"), position = c(NA, T), gates = c(NA, gthres[2]))
  singlets.flowD.l@proportion <- (singlets.flowD.l@cell.count/nrow(f))*100
  singlets <- getflowFrame(singlets.flowD.l)
  
  # Indices in the Singlets Population with respect to f. Will be needed for the MFI calculation
  singlets.flowD.l.ind <- singlets.flowD.l@index
  # # Optional Plot: Done to observe if the correct cells of the singlets population are being saved with respect to f
  # plotDens(f,c("SSC-A","FSC-W"),main="Ungated", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(f@exprs[singlets.flowD.l.ind,c(4,3)], col=2, cex = 0.1)

  all.events[q2,2] <- singlets.flowD.l@cell.count
  all.props[q2,2] <- singlets.flowD.l@proportion
  
  
  ###############################################################################################
  ## Gating Singlets to obtain the Live population. Plotting Live/Dead_SSC-A -----
  live.flowD <- flowDensity(singlets.flowD.l, channels = c(11,4), position = c(F,NA), gates = c(gthres[3], NA))
  live.flowD@proportion <- (live.flowD@cell.count/singlets.flowD.l@cell.count)*100
  live <- getflowFrame(live.flowD)
  
  # Indices in the Live Population with respect to f. Will be needed for the MFI calculation
  live.flowD.ind <- live.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the Live population are being saved with respect to f
  # plotDens(singlets,  c(11,4), main= "Singlets", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(singlets.flowD.l@flow.frame@exprs[live.flowD.ind,c(11,4)], col=2, cex = 0.1)
  
  all.events[q2,3] <- live.flowD@cell.count
  all.props[q2,3] <- live.flowD@proportion
  
  ###############################################################################################
  ## Gating Live to obtain the Lymphocyte population. Plotting FSC-A_SSC-A
  lymph.flowD.temp <- flowDensity(live.flowD, channels = c(1,4), position = c(T, F), gates = c(gthres[5], gthres[6]))
  lymph.flowD <- flowDensity(lymph.flowD.temp, channels = c(1,4), position = c(F,F), gates = c(gthres[4], gthres[6]))
  lymph.flowD@proportion <- (lymph.flowD@cell.count/live.flowD@cell.count)*100
  lymph <- getflowFrame(lymph.flowD)
  
  # Indices in the Lymphocytes Population with respect to f. Will be needed for the MFI calculation
  lymph.flowD.ind <- lymph.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the Lymph population are being saved with respect to f
  # plotDens(live, c(1,4), main = "Live", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(live.flowD@flow.frame@exprs[lymph.flowD.ind,c(1,4)], col=2, cex = 0.1)
  
  all.events[q2,4] <- lymph.flowD@cell.count
  all.props[q2,4] <- lymph.flowD@proportion
  ###############################################################################################
  
  ## Gating Lymphocytes to obtain the CD45+ population. Plotting CD45_CD43 (CD45+ population)
  theta = atan(1)
  lymph.temp.flowD <- lymph.flowD
  lymph.temp.flowD@flow.frame <-rotate.data(lymph.flowD@flow.frame,c(14,8),theta = pi/6.50)$data
 
  cd45.flowD.temp <- flowDensity(lymph.temp.flowD, channels = c(14,8), position = c(T,F), gates = c(gthres[7], gthres[9]))
  cd45.flowD <- flowDensity(cd45.flowD.temp, channels = c(14,8), position = c(F,F), gates = c(gthres[8], gthres[9]))

  cd45.flowD@filter <-  rotate.data(cd45.flowD@filter,c(14,8),theta = -pi/6.50)$data
  cd45.flowD@flow.frame <- rotate.data(cd45.flowD@flow.frame,c(14,8),theta = -pi/6.50)$data
  cd45.flowD@proportion <- (cd45.flowD@cell.count/lymph.flowD@cell.count)*100
  cd45 <- getflowFrame(cd45.flowD)
  
  # Indices in the CD45+ Population with respect to f. Will be needed for the MFI calculation
  cd45.flowD.ind <- cd45.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the CD45+ population are being saved with respect to f
  # plotDens(lymph, channels = c(14,8), main = "Lymphocytes", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # lines(cd45.flowD@filter)
  # points(lymph.flowD@flow.frame@exprs[cd45.flowD.ind,c(14,8)], col=2, cex = 0.1)

  all.events[q2,5] <- cd45.flowD@cell.count
  all.props[q2,5] <- cd45.flowD@proportion
  
  ###############################################################################################
  ## Gating CD45+ to obtain Granulocyte Pre and NOT(Granulocyte Pre) populations. Plotting GR1_CD43
  # NOT(Granulocyte Pre)
  NOT.granulocyte.flowD <- flowDensity(cd45.flowD, channels = c(10, 8), position = c(F,NA), gates = c(gthres[10], gthres[9]))
  NOT.granulocyte.flowD@proportion <- (NOT.granulocyte.flowD@cell.count/cd45.flowD@cell.count)*100
  NOT.granulocyte <- getflowFrame(NOT.granulocyte.flowD)
  
  # Indices in the NOT(Granulocyte) population with respect to f. Will be needed for the MFI calculation
  NOT.granulocyte.flowD.ind <- NOT.granulocyte.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the NOT(Granulocyte) population are being saved with respect to f
  # plotDens(cd45, channels = c(10,8), main = "CD45+: NOT(Granulocyte)", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(cd45.flowD@flow.frame@exprs[NOT.granulocyte.flowD.ind,c(10,8)], col=2, cex = 0.1)
  
  all.events[q2,6] <- NOT.granulocyte.flowD@cell.count
  all.props[q2,6] <- NOT.granulocyte.flowD@proportion
  
  
  
  # Granulocyte Pre
  #granulocyte.flowD <- notSubFrame(cd45.flowD@flow.frame, channels = c(10, 8), position = "logical", gates = "missing", NOT.granulocyte.flowD@filter)
  granulocyte.flowD <- flowDensity(cd45.flowD, channels = c(10, 8), position = c(T,NA), gates = c(gthres[10], gthres[9]))
  granulocyte.flowD@proportion <- (granulocyte.flowD@cell.count/cd45.flowD@cell.count)*100
  granulocyte <- getflowFrame(granulocyte.flowD)
  
  # Indices in the Granulocyte population with respect to f. Will be needed for the MFI calculation
  granulocyte.flowD.ind <- granulocyte.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the Granulocyte population are being saved with respect to f
  # plotDens(cd45, channels = c(10,8), main = "CD45+: Granulocyte", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # lines(granulocyte.flowD@filter)
  # points(cd45.flowD@flow.frame@exprs[granulocyte.flowD.ind,c(10,8)], col=2, cex = 0.1)
  
  all.events[q2,7] <- granulocyte.flowD@cell.count
  all.props[q2,7] <- granulocyte.flowD@proportion

  ###############################################################################################
  
  ## Gating NOT(Granulocyte Pre) population to obtain CD3 T-cells and NOT(CD3 T-cells). Plotting B220_CD3
  # CD3 T-cells
  cd3.Tcell.flowD <- flowDensity(NOT.granulocyte.flowD, channels = c(18,16), position = c(F, T), gates = c(gthres[11], gthres[12]))
  cd3.Tcell.flowD@proportion <- (cd3.Tcell.flowD@cell.count/NOT.granulocyte.flowD@cell.count)*100
  cd3.Tcell <- getflowFrame(cd3.Tcell.flowD)
  
  # Indices in the NOT(Granulocyte): CD3 T-cell population with respect to f. Will be needed for the MFI calculation
  cd3.Tcell.flowD.ind <- cd3.Tcell.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the NOT(Granulocyte): CD3 T-cell population are being saved with respect to f
  # plotDens(NOT.granulocyte, channels = c(18,16), main = "NOT(Granulocyte): CD3 T-cell", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(NOT.granulocyte.flowD@flow.frame@exprs[cd3.Tcell.flowD.ind,c(18,16)], col=2, cex = 0.1)
  
  all.events[q2,8] <- cd3.Tcell.flowD@cell.count
  all.props[q2,8] <- cd3.Tcell.flowD@proportion
  
  
  # NOT(CD3 T-cells)
  NOT.cd3.Tcell.flowD <- notSubFrame(NOT.granulocyte.flowD@flow.frame, channels = c(18,16), position = "logical", gates = "missing", cd3.Tcell.flowD@filter)
  NOT.cd3.Tcell.flowD@proportion <- (NOT.cd3.Tcell.flowD@cell.count/NOT.granulocyte.flowD@cell.count)*100
  NOT.cd3.Tcell <- getflowFrame(NOT.cd3.Tcell.flowD)
  
  # Indices in the NOT(Granulocyte): NOT(CD3 T-cells) population with respect to f. Will be needed for the MFI calculation
  NOT.cd3.Tcell.flowD.ind <- NOT.cd3.Tcell.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the NOT(Granulocyte): NOT(CD3 T-cell) population are being saved with respect to f
  # plotDens(NOT.granulocyte, channels = c(18,16), main = "NOT(Granulocyte): NOT(CD3 T-cell)", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(NOT.granulocyte.flowD@flow.frame@exprs[NOT.cd3.Tcell.flowD.ind,c(18,16)], col=2, cex = 0.1)
  
  all.events[q2,9] <- NOT.cd3.Tcell.flowD@cell.count
  all.props[q2,9] <- NOT.cd3.Tcell.flowD@proportion

  ###############################################################################################
 
  ## Gating NOT(CD3 T-cell) population to obtain Plasma and NOT Plasma. Plotting CD138_B220
  #Plasma
  plasma.flowD <- flowDensity(NOT.cd3.Tcell.flowD, channels = c(15,18), position = c(T, F), gates = c(gthres[13], gthres[11]))
  plasma.flowD@proportion <- (plasma.flowD@cell.count/NOT.cd3.Tcell.flowD@cell.count)*100
  
  # Indices in the NOT(CD3 T-cells): Plasma population with respect to f. Will be needed for the MFI calculation
  plasma.flowD.ind <- plasma.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the NOT(CD3 T-cell): Plasma population are being saved with respect to f
  # plotDens(NOT.cd3.Tcell, channels = c(15,18), main = "NOT(CD3 T-cell): Plasma", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(NOT.cd3.Tcell.flowD@flow.frame@exprs[plasma.flowD.ind,c(15,18)], col=2, cex = 0.1)
  
  all.events[q2,10] <- plasma.flowD@cell.count
  all.props[q2,10] <- plasma.flowD@proportion
  
  
  
  # NOT Plasma
  NOT.plasma.flowD <- notSubFrame(NOT.cd3.Tcell.flowD@flow.frame, channels = c(15,18), position = "logical", gates = "missing", plasma.flowD@filter)
  NOT.plasma.flowD@proportion <- (NOT.plasma.flowD@cell.count/NOT.cd3.Tcell.flowD@cell.count)*100
  NOT.plasma <- getflowFrame(NOT.plasma.flowD)
  # Indices in the NOT(CD3 T-cells): NOT(Plasma) population with respect to f. Will be needed for the MFI calculation
  NOT.plasma.flowD.ind <- NOT.plasma.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the NOT(CD3 T-cell): NOT(Plasma) population are being saved with respect to f
  # plotDens(NOT.cd3.Tcell, channels = c(15,18), main = "NOT(CD3 T-cell): NOT(Plasma)", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(NOT.cd3.Tcell.flowD@flow.frame@exprs[NOT.plasma.flowD.ind,c(15,18)], col=2, cex = 0.1)
  
  all.events[q2,11] <- NOT.plasma.flowD@cell.count
  all.props[q2,11] <- NOT.plasma.flowD@proportion
  
  ###############################################################################################
  
  ## Gating NOT Plasma population to obtain Myeloid Pre and B-cells. Plotting B220_CD11b
  # B-cells
  Bcell.flowD <- flowDensity(NOT.plasma.flowD, channels = c(18,13), position = c(T, NA), gates = c(gthres[11], gthres[14]))
  Bcell.flowD@proportion <- (Bcell.flowD@cell.count/NOT.plasma.flowD@cell.count)*100
  Bcell <- getflowFrame(Bcell.flowD)
  
  # Indices in the NOT(Plasma): B-cells population with respect to f. Will be needed for the MFI calculation
  Bcell.flowD.ind <- Bcell.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the NOT(Plasma): B-cells population are being saved with respect to f
  # plotDens(NOT.plasma, channels = c(18,13), main = "NOT(Plasma): B-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(NOT.plasma.flowD@flow.frame@exprs[Bcell.flowD.ind,c(18,13)], col=2, cex = 0.1)
  
  all.events[q2,12] <- Bcell.flowD@cell.count
  all.props[q2,12] <- Bcell.flowD@proportion
  
  
  # Myeloid Pre
  myeloid.flowD <- flowDensity(NOT.plasma.flowD, channels = c(18,13), position = c(F, NA), gates = c(gthres[11], gthres[14]))
  myeloid.flowD@proportion <- (myeloid.flowD@cell.count/NOT.plasma.flowD@cell.count)*100
  # Indices in the NOT(Plasma): Myeloid Pre population with respect to f. Will be needed for the MFI calculation
  myeloid.flowD.ind <- myeloid.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the NOT(Plasma): Myeloid Pre population are being saved with respect to f
  # plotDens(NOT.plasma, channels = c(18,13), main = "NOT(Plasma): Myeloid Pre", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(NOT.plasma.flowD@flow.frame@exprs[myeloid.flowD.ind,c(18,13)], col=2, cex = 0.1)
  
  all.events[q2,13] <- myeloid.flowD@cell.count
  all.props[q2,13] <- myeloid.flowD@proportion
  
  ###############################################################################################
 
  ## Gating B-cells population to obtain CD43+ and CD43-. Plotting B220_CD43
  # CD43+
  cd43.pos.flowD <- flowDensity(Bcell.flowD, channels = c(18,8), position = c(NA, T), gates = c(gthres[11], gthres[15]))
  cd43.pos.flowD@proportion <- (cd43.pos.flowD@cell.count/Bcell.flowD@cell.count)*100
  cd43.pos <- getflowFrame(cd43.pos.flowD)
  
  # Indices in the B-cells: CD43+ population with respect to f. Will be needed for the MFI calculation
  cd43.pos.flowD.ind <- cd43.pos.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the B-cells: CD43+ population are being saved with respect to f
  # plotDens(Bcell, channels = c(18,8), main = "B-cells: CD43+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(Bcell.flowD@flow.frame@exprs[cd43.pos.flowD.ind,c(18,8)], col=2, cex = 0.1)
  
  all.events[q2,14] <- cd43.pos.flowD@cell.count
  all.props[q2,14] <- cd43.pos.flowD@proportion
  
  
  #CD43-
  cd43.neg.flowD <- notSubFrame(Bcell.flowD@flow.frame, channels = c(18,8), position = "logical", gates = "missing", cd43.pos.flowD@filter)
  cd43.neg.flowD@proportion <- (cd43.neg.flowD@cell.count/Bcell.flowD@cell.count)*100
  cd43.neg <- getflowFrame(cd43.neg.flowD)
  
  # Indices in the B-cells: CD43- population with respect to f. Will be needed for the MFI calculation
  cd43.neg.flowD.ind <- cd43.neg.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the B-cells: CD43- population are being saved with respect to f
  # plotDens(Bcell, channels = c(18,8), main = "B-cells: CD43-", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(Bcell.flowD@flow.frame@exprs[cd43.neg.flowD.ind,c(18,8)], col=2, cex = 0.1)
  
  all.events[q2,15] <- cd43.neg.flowD@cell.count
  all.props[q2,15] <- cd43.neg.flowD@proportion
  
  #############################################################################################

  ## Gating CD43+ population to obtain HFA, HFB, and HFC. Plotting CD24_BP1
  # HFA
  HFA.flowD <- flowDensity(cd43.pos.flowD, channels = c(9, 17), position = c(F,F), gates = c(gthres[16], gthres[19]))
  HFA.flowD@proportion <- (HFA.flowD@cell.count/cd43.pos.flowD@cell.count)*100
  
  # Indices in the CD43+: HFA population with respect to f. Will be needed for the MFI calculation
  HFA.flowD.ind <- HFA.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the CD43+: HFA population are being saved with respect to f
  # plotDens(cd43.pos, channels = c(9,17), main = "CD43+: HFA", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(cd43.pos.flowD@flow.frame@exprs[HFA.flowD.ind,c(9,17)], col=2, cex = 0.1)
  
  all.events[q2,16] <- HFA.flowD@cell.count
  all.props[q2,16] <- HFA.flowD@proportion
  
  
  # HFB
  HFB.flowD <- flowDensity(cd43.pos.flowD, channels = c(9, 17), position = c(T,F), gates = c(gthres[16], gthres[19]))
  HFB.flowD@proportion <- (HFB.flowD@cell.count/cd43.pos.flowD@cell.count)*100
  
  # Indices in the CD43+: HFB population with respect to f. Will be needed for the MFI calculation
  HFB.flowD.ind <- HFB.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the CD43+: HFB population are being saved with respect to f
  # plotDens(cd43.pos, channels = c(9,17), main = "CD43+: HFB", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(cd43.pos.flowD@flow.frame@exprs[HFB.flowD.ind,c(9,17)], col=2, cex = 0.1)
  
  all.events[q2,17] <- HFB.flowD@cell.count
  all.props[q2,17] <- HFB.flowD@proportion
  
  
  # HFC
  HFC.flowD <- flowDensity(cd43.pos.flowD, channels = c(9, 17), position = c(T,T), gates = c(gthres[16], gthres[19]))
  HFC.flowD@proportion <- (HFC.flowD@cell.count/cd43.pos.flowD@cell.count)*100
  
  # Indices in the CD43+: HFC population with respect to f. Will be needed for the MFI calculation
  HFC.flowD.ind <- HFC.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the CD43+: HFC population are being saved with respect to f
  # plotDens(cd43.pos, channels = c(9,17), main = "CD43+: HFC", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(cd43.pos.flowD@flow.frame@exprs[HFC.flowD.ind,c(9,17)], col=2, cex = 0.1)
  
  all.events[q2,18] <- HFC.flowD@cell.count
  all.props[q2,18] <- HFC.flowD@proportion
  
  #############################################################################################
 
  ## Gating CD43- population to obtain HFD, HFE, and HFF. Plotting IgM_IgD
  # HFD
  HFD.flowD <- flowDensity(cd43.neg.flowD, channels = c(12,7), position = c(F,F), gates = c(gthres[17], gthres[18]))
  HFD.flowD@proportion <- (HFD.flowD@cell.count/cd43.neg.flowD@cell.count)*100
  
  # Indices in the CD43-: HFD population with respect to f. Will be needed for the MFI calculation
  HFD.flowD.ind <- HFD.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the CD43-: HFD population are being saved with respect to f
  # plotDens(cd43.neg, channels = c(12,7), main = "CD43-: HFD", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(cd43.neg.flowD@flow.frame@exprs[HFD.flowD.ind,c(12,7)], col=2, cex = 0.1)
  
  all.events[q2,19] <- HFD.flowD@cell.count
  all.props[q2,19] <- HFD.flowD@proportion
  
  
  # HFE
  HFE.flowD <- flowDensity(cd43.neg.flowD, channels = c(12,7), position = c(T,F), gates = c(gthres[17], gthres[18]))
  HFE.flowD@proportion <- (HFE.flowD@cell.count/cd43.neg.flowD@cell.count)*100
  
  # Indices in the CD43-: HFE population with respect to f. Will be needed for the MFI calculation
  HFE.flowD.ind <- HFE.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the CD43-: HFE population are being saved with respect to f
  # plotDens(cd43.neg, channels = c(12,7), main = "CD43-: HFE", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(cd43.neg.flowD@flow.frame@exprs[HFE.flowD.ind,c(12,7)], col=2, cex = 0.1)
  
  all.events[q2,20] <- HFE.flowD@cell.count
  all.props[q2,20] <- HFE.flowD@proportion
  
  
  # HFF
  HFF.flowD <- flowDensity(cd43.neg.flowD, channels = c(12,7), position = c(NA,T), gates = c(gthres[17], gthres[18]))
  HFF.flowD@proportion <- (HFF.flowD@cell.count/cd43.neg.flowD@cell.count)*100
  
  # Indices in the CD43-: HFF population with respect to f. Will be needed for the MFI calculation
  HFF.flowD.ind <- HFF.flowD@index
  # # Optional Plot: Done to observe if the correct cells of the CD43-: HFF population are being saved with respect to f
  # plotDens(cd43.neg, channels = c(12,7), main = "CD43-: HFF", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  # points(cd43.neg.flowD@flow.frame@exprs[HFF.flowD.ind,c(12,7)], col=2, cex = 0.1)
  
  all.events[q2,21] <- HFF.flowD@cell.count
  all.props[q2,21] <- HFF.flowD@proportion
  

  ############################################################################################################
  ############################################################################################################
  
  
  ## Saving the MFIs (mean) of Original Data (Data before Transformed)
  singlets.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[singlets.flowD.l.ind, x], na.rm = T)})
  live.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[live.flowD.ind, x], na.rm = T)})
  lymph.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[lymph.flowD.ind, x], na.rm = T)})
  cd45.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd45.flowD.ind, x], na.rm = T)})
  NOT.granulocyte.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[NOT.granulocyte.flowD.ind, x], na.rm = T)})
  granulocyte.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[granulocyte.flowD.ind, x], na.rm = T)})
  cd3.Tcell.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd3.Tcell.flowD.ind, x], na.rm = T)})
  NOT.cd3.Tcell.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[NOT.cd3.Tcell.flowD.ind, x], na.rm = T)})
  plasma.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[plasma.flowD.ind, x], na.rm = T)})
  NOT.plasma.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[NOT.plasma.flowD.ind, x], na.rm = T)})
  Bcell.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[Bcell.flowD.ind, x], na.rm = T)})
  myeloid.MFI <-sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[myeloid.flowD.ind, x], na.rm = T)})
  cd43.pos.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd43.pos.flowD.ind, x], na.rm = T)}) 
  cd43.neg.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[cd43.neg.flowD.ind, x], na.rm = T)})
  HFA.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[HFA.flowD.ind, x], na.rm = T)})
  HFB.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[HFB.flowD.ind, x], na.rm = T)})
  HFC.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[HFC.flowD.ind, x], na.rm = T)})
  HFD.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[HFD.flowD.ind, x], na.rm = T)})
  HFE.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[HFE.flowD.ind, x], na.rm = T)})
  HFF.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){mean(f.beforeTransform@exprs[HFF.flowD.ind, x], na.rm = T)})
  
  MFI.mean.OriginalData[q2, ] <- c(singlets.MFI, live.MFI, lymph.MFI, cd45.MFI, NOT.granulocyte.MFI, granulocyte.MFI, cd3.Tcell.MFI, NOT.cd3.Tcell.MFI,
                                   plasma.MFI, NOT.plasma.MFI, Bcell.MFI, myeloid.MFI, cd43.pos.MFI, cd43.neg.MFI, HFA.MFI, HFB.MFI, HFC.MFI, HFD.MFI, HFE.MFI, HFF.MFI)
  
  
  ## Saving the MFIs (median) of Original Data (Data before Transformed)
  singlets.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[singlets.flowD.l.ind, x], na.rm = T)})
  live.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[live.flowD.ind, x], na.rm = T)})
  lymph.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[lymph.flowD.ind, x], na.rm = T)})
  cd45.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd45.flowD.ind, x], na.rm = T)})
  NOT.granulocyte.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[NOT.granulocyte.flowD.ind, x], na.rm = T)})
  granulocyte.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[granulocyte.flowD.ind, x], na.rm = T)})
  cd3.Tcell.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd3.Tcell.flowD.ind, x], na.rm = T)})
  NOT.cd3.Tcell.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[NOT.cd3.Tcell.flowD.ind, x], na.rm = T)})
  plasma.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[plasma.flowD.ind, x], na.rm = T)})
  NOT.plasma.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[NOT.plasma.flowD.ind, x], na.rm = T)})
  Bcell.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[Bcell.flowD.ind, x], na.rm = T)})
  myeloid.MFI <-sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[myeloid.flowD.ind, x], na.rm = T)})
  cd43.pos.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd43.pos.flowD.ind, x], na.rm = T)}) 
  cd43.neg.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[cd43.neg.flowD.ind, x], na.rm = T)})
  HFA.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[HFA.flowD.ind, x], na.rm = T)})
  HFB.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[HFB.flowD.ind, x], na.rm = T)})
  HFC.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[HFC.flowD.ind, x], na.rm = T)})
  HFD.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[HFD.flowD.ind, x], na.rm = T)})
  HFE.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[HFE.flowD.ind, x], na.rm = T)})
  HFF.MFI <- sapply(1:(ncol(f.beforeTransform@exprs)-1), function(x){median(f.beforeTransform@exprs[HFF.flowD.ind, x], na.rm = T)})
  
  MFI.median.OriginalData[q2, ] <- c(singlets.MFI, live.MFI, lymph.MFI, cd45.MFI, NOT.granulocyte.MFI, granulocyte.MFI, cd3.Tcell.MFI, NOT.cd3.Tcell.MFI,
                                   plasma.MFI, NOT.plasma.MFI, Bcell.MFI, myeloid.MFI, cd43.pos.MFI, cd43.neg.MFI, HFA.MFI, HFB.MFI, HFC.MFI, HFD.MFI, HFE.MFI, HFF.MFI)
  
  
  ## Saving the MFIs (mean) of Transformed Data (Data after Transformation)
  singlets.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[singlets.flowD.l.ind, x], na.rm = T)})
  live.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[live.flowD.ind, x], na.rm = T)})
  lymph.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[lymph.flowD.ind, x], na.rm = T)})
  cd45.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd45.flowD.ind, x], na.rm = T)})
  NOT.granulocyte.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[NOT.granulocyte.flowD.ind, x], na.rm = T)})
  granulocyte.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[granulocyte.flowD.ind, x], na.rm = T)})
  cd3.Tcell.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd3.Tcell.flowD.ind, x], na.rm = T)})
  NOT.cd3.Tcell.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[NOT.cd3.Tcell.flowD.ind, x], na.rm = T)})
  plasma.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[plasma.flowD.ind, x], na.rm = T)})
  NOT.plasma.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[NOT.plasma.flowD.ind, x], na.rm = T)})
  Bcell.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[Bcell.flowD.ind, x], na.rm = T)})
  myeloid.MFI <-sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[myeloid.flowD.ind, x], na.rm = T)})
  cd43.pos.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd43.pos.flowD.ind, x], na.rm = T)}) 
  cd43.neg.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[cd43.neg.flowD.ind, x], na.rm = T)})
  HFA.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[HFA.flowD.ind, x], na.rm = T)})
  HFB.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[HFB.flowD.ind, x], na.rm = T)})
  HFC.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[HFC.flowD.ind, x], na.rm = T)})
  HFD.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[HFD.flowD.ind, x], na.rm = T)})
  HFE.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[HFE.flowD.ind, x], na.rm = T)})
  HFF.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){mean(f@exprs[HFF.flowD.ind, x], na.rm = T)})
  
  MFI.mean.TransformedData[q2, ] <- c(singlets.MFI, live.MFI, lymph.MFI, cd45.MFI, NOT.granulocyte.MFI, granulocyte.MFI, cd3.Tcell.MFI, NOT.cd3.Tcell.MFI,
                                   plasma.MFI, NOT.plasma.MFI, Bcell.MFI, myeloid.MFI, cd43.pos.MFI, cd43.neg.MFI, HFA.MFI, HFB.MFI, HFC.MFI, HFD.MFI, HFE.MFI, HFF.MFI)
  
  
  ## Saving the MFIs (median) of Transformed Data (Data after Transformation)
  singlets.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[singlets.flowD.l.ind, x], na.rm = T)})
  live.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[live.flowD.ind, x], na.rm = T)})
  lymph.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[lymph.flowD.ind, x], na.rm = T)})
  cd45.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd45.flowD.ind, x], na.rm = T)})
  NOT.granulocyte.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[NOT.granulocyte.flowD.ind, x], na.rm = T)})
  granulocyte.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[granulocyte.flowD.ind, x], na.rm = T)})
  cd3.Tcell.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd3.Tcell.flowD.ind, x], na.rm = T)})
  NOT.cd3.Tcell.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[NOT.cd3.Tcell.flowD.ind, x], na.rm = T)})
  plasma.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[plasma.flowD.ind, x], na.rm = T)})
  NOT.plasma.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[NOT.plasma.flowD.ind, x], na.rm = T)})
  Bcell.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[Bcell.flowD.ind, x], na.rm = T)})
  myeloid.MFI <-sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[myeloid.flowD.ind, x], na.rm = T)})
  cd43.pos.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd43.pos.flowD.ind, x], na.rm = T)}) 
  cd43.neg.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[cd43.neg.flowD.ind, x], na.rm = T)})
  HFA.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[HFA.flowD.ind, x], na.rm = T)})
  HFB.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[HFB.flowD.ind, x], na.rm = T)})
  HFC.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[HFC.flowD.ind, x], na.rm = T)})
  HFD.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[HFD.flowD.ind, x], na.rm = T)})
  HFE.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[HFE.flowD.ind, x], na.rm = T)})
  HFF.MFI <- sapply(1:(ncol(f@exprs)-1), function(x){median(f@exprs[HFF.flowD.ind, x], na.rm = T)})
  
  MFI.median.TransformedData[q2, ] <- c(singlets.MFI, live.MFI, lymph.MFI, cd45.MFI, NOT.granulocyte.MFI, granulocyte.MFI, cd3.Tcell.MFI, NOT.cd3.Tcell.MFI,
                                      plasma.MFI, NOT.plasma.MFI, Bcell.MFI, myeloid.MFI, cd43.pos.MFI, cd43.neg.MFI, HFA.MFI, HFB.MFI, HFC.MFI, HFD.MFI, HFE.MFI, HFF.MFI)
  

  #############################################################################################################
  #############################################################################################################
  
  ## Saving the Plots
  #--------Start Big Png------------
  png ( file = paste("Results/Figures/ScatterPlots_Updated/", store.allFCS[q2,1], "/", "Total_", store.allFCS[q2,2], ".png", sep = "" ), width=2100, height=2100*4/4)
  #png ( file = paste("Results/tempPlots/", "Total_", store.allFCS[q2,2], ".png", sep = "" ), width=2100, height=2100*4/4)
  
  
  par(mfrow=c(4,4),mar=(c(5, 5, 4, 2) + 0.1))
  
  # 1st method of plotting FSC-A_SSC-W
  plotDens(f,c("SSC-A","FSC-W"),main="Ungated", devn = F, cex.lab = 2, cex.axis = 2, cex.main=2)
  lines(singlets.flowD.l@filter,lwd=2)
  
  # Plotting Live/Dead_SSC-A 
  plotDens(singlets,  c(11,4), main= "Singlets", devn = F, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=gthres[3], lwd=2); 
  
  # Plotting FSC-A_SSC-A
  plotDens(live, c(1,4), main = "Live", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=gthres[5], lwd=2); abline(v=gthres[4], lwd=2); abline(h=gthres[6], lwd=2)
  
  # Plotting CD45_CD43 to illustrate the Lymphocyte population
  plotDens(lymph, channels = c(14,8), main = "Lymphocytes", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
  
  # Plotting CD45_CD43 to obtain the CD45+ population from the rotated Lymphocyte population
  plotDens(lymph.temp.flowD@flow.frame, channels = c(14,8), main = "Lymphocytes (Rotated)", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=gthres[7], lwd=2); abline(v=gthres[8], lwd=2); abline(h=gthres[9], lwd=2)
  
  # Plotting GR1_CD43
  min.y <- min(exprs(cd45)[,8])
  max.y <- max(exprs(cd45)[,8])+0.5
  plotDens(cd45, channels = c(10,8),  main = "CD45+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y))
  
  # Plotting GR1_CD43 to obtain the Granulocyte Pre and NOT(Granulocyte Pre) populations
  min.y <- min(exprs(cd45)[,8])
  max.y <- max(exprs(cd45)[,8])+0.5
  plotDens(cd45, channels = c(10,8),  main = "CD45+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y)); abline(v=gthres[10], lwd=2)
  
  # Plotting B220_CD3 to obtain CD3 T-cells and NOT(CD3 T-cells)
  min.x <- min(exprs(NOT.granulocyte)[,18])
  max.x <- max(exprs(NOT.granulocyte)[,18])
  min.y <- min(exprs(NOT.granulocyte)[,16])
  max.y <-  max(exprs(NOT.granulocyte)[,16])
  plotDens(NOT.granulocyte, channels = c(18,16),  main = "NOT(Granulocyte Pre)", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); 
  lines(cd3.Tcell.flowD@filter, lwd=2)
  
  # Plotting CD138_B220 to obtain Plasma and NOT Plasma
  plotDens(NOT.cd3.Tcell, channels = c(15,18),  main = "NOT(CD3 T-cell)", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); 
  lines(plasma.flowD@filter, lwd=2)
  
  #  Plotting B220_CD11b to obtain Myeloid Pre and B-cells
  min.x <- min(exprs(NOT.plasma)[,18])
  max.x <- max(exprs(NOT.plasma)[,18])
  min.y <- min(exprs(NOT.plasma)[,13])
  max.y <-  max(exprs(NOT.plasma)[,13])
  plotDens(NOT.plasma, channels = c(18, 13),  main = "NOT Plasma", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=gthres[11], lwd=2)
  
  
  # Plotting B220_CD43 to obtain CD43+ and CD43-
  min.x <- min(exprs(Bcell)[,18])
  max.x <- max(exprs(Bcell)[,18])
  min.y <- min(exprs(Bcell)[,8])
  max.y <-  max(exprs(Bcell)[,8])
  plotDens(Bcell, channels = c(18, 8),  main = "B-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x-1.5, max.x)); abline(h=gthres[15], lwd=2); abline(v=min.x, lwd=2)
  
  # Plotting CD24_BP1 to obtain HFA, HFB, and HFC. 
  min.x <- min(exprs(cd43.pos)[,9])
  max.x <- max(exprs(cd43.pos)[,9])
  min.y <- min(exprs(cd43.pos)[,17])
  max.y <-  max(exprs(cd43.pos)[,17])
  plotDens(cd43.pos, channels = c(9, 17),  main = "CD43+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=gthres[16], lwd=2); abline(h=gthres[19], lwd=2)
  
  
  # Plotting IgM_IgD to obtain HFD, HFE, and HFF. 
  min.y <- min(exprs(cd43.neg)[,7])
  plotDens(cd43.neg, channels = c(12, 7),  main = "CD43-", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(h=gthres[18], lwd=2)
  lines(x=c(gthres[17], gthres[17]), y=c(gthres[18], min.y-1), lwd=2)
  
  
  dev.off()
  par(mfrow=c(1,1))
  #--------End Big Png------------

},error = function(err) {
  print(paste("ERROR in Gating/Plotting Message:  ",err, "Index: ", q2));
  errorFileIndex <- c(errorFileIndex,q2)
  return(errorFileIndex)
}) # end of tryCatch

        ## Start of flowType
          if (Do_flowType == T) {
          print("Start FlowType")
            # ## flowType on CD45 population. Live marker is excluded but CD45 marker is included. Also the maximum number of markers allowed is 6.
            # flowType.res <- flowType(Frame = cd45, PropMarkers= as.vector(channels.ind.NoLive), MaxMarkersPerPop = 6, PartitionsPerMarker=2,
            #                          Methods='Thresholds', Thresholds= listGthres.NoLive, verbose=F, MemLimit=400)
            # 
            # save ( flowType.res, file =  paste("Results/FlowType/", store.allFCS[q2,1],"/FT_", store.allFCS[q2,2],".Rdata",sep="") )
                    
            ## flowType on CD45 population. Live marker and CD45 marker is excluded. Also the maximum number of markers allowed is 8.
            flowType.res <- flowType(Frame = cd45, PropMarkers= as.vector(channels.ind.NoLiveNoCD45), MaxMarkersPerPop = 8, PartitionsPerMarker=2,
                                     Methods='Thresholds', Thresholds= listGthres.NoLiveNoCD45, verbose=F, MemLimit=400)
            
            save ( flowType.res, file =  paste("Results/FlowType/", store.allFCS[q2,1],"/FT_", store.allFCS[q2,2],".Rdata",sep="") )
            
        }
        
        cat("Time is: ",TimeOutput(start2),"\n",sep="")
} # end of for-loop

colnames(all.events) <- c("All Events-#Events", "Singlets-#Events", "Live-#Events", "Lymphocytes-#Events", "CD45+-#Events", "NOT(Granulocyte Pre)-#Events", "Granulocyte Pre-#Events",
                          "CD3 T-cell-#Events", "NOT(CD3 T-cell)-#Events", "Plasma-#Events", "NOT Plasma-#Events", "B-cells-#Events", "Myeloid Pre-#Events", "CD43+ #Events", "CD43- #Events",
                          "HFA-#Events", "HFB-#Events", "HFC-#Events", "HFD-#Events", "HFE-#Events", "HFF-#Events")
rownames(all.events) <- 1:length(all.events[,1])

colnames(all.props) <- c("All Events-%Parent", "Singlets-%Parent", "Live-%Parent", "Lymphocytes-%Parent", "CD45+-%Parent", "NOT(Granulocyte Pre)-%Parent", "Granulocyte Pre-%Parent",
                         "CD3 T-cell-%Parent", "NOT(CD3 T-cell)-%Parent", "Plasma-%Parent", "NOT Plasma-%Parent", "B-cells-%Parent", "Myeloid Pre-%Parent", "CD43+ %Parent", "CD43- %Parent",
                         "HFA-%Parent", "HFB-%Parent", "HFC-%Parent", "HFD-%Parent", "HFE-%Parent", "HFF-%Parent")
rownames(all.props) <- 1:length(all.props[,1])

Events_Proportions_Table <- NULL
for(i in 1:ncol(all.props)){
  Events_Proportions_Table <- cbind(Events_Proportions_Table,all.events[ ,i],all.props[ ,i])
}

colnames(Events_Proportions_Table) <- c("All Events-#Events", "All Events-%Parent", "Singlets-#Events", "Singlets-%Parent", "Live-#Events", "Live-%Parent", "Lymphocytes-#Events", "Lymphocytes-%Parent", 
                                        "CD45+-#Events", "CD45+-%Parent", "NOT(Granulocyte Pre)-#Events", "NOT(Granulocyte Pre)-%Parent", "Granulocyte Pre-#Events", "Granulocyte Pre-%Parent", 
                                        "CD3 T-cell-#Events", "CD3 T-cell-%Parent", "NOT(CD3 T-cell)-#Events", "NOT(CD3 T-cell)-%Parent", "Plasma-#Events", "Plasma-%Parent", "NOT Plasma-#Events", "NOT Plasma-%Parent", 
                                        "B-cells-#Events", "B-cells-%Parent", "Myeloid Pre-#Events", "Myeloid Pre-%Parent", "CD43+ #Events",  "CD43+ %Parent", "CD43- #Events", "CD43- %Parent", 
                                        "HFA-#Events", "HFA-%Parent", "HFB-#Events", "HFB-%Parent", "HFC-#Events", "HFC-%Parent", "HFD-#Events", "HFD-%Parent", "HFE-#Events", "HFE-%Parent", "HFF-#Events", "HFF-%Parent")


## Code for removing duplicate FCS files based on their barcodes
Barcodes <- str_extract(store.allFCS[,2],"L[0-9]+")
store.allFCS.unique <- cbind(store.allFCS, Barcodes)
duplicate.index <- which(duplicated(Barcodes)==TRUE)
duplicate.FCS <- store.allFCS[duplicate.index, 1:2]
store.allFCS.unique <- store.allFCS.unique[!duplicated(store.allFCS.unique[,c('Barcodes')]),]
rownames(store.allFCS.unique) <- 1:length(store.allFCS.unique[,1])
# Combining the FCS files names, Genotype, The Number of cells before cleaning, number of channels, and the Assay Date with the Event_Proportion Table
Events_Proportions_Table <- cbind(store.allFCS.unique, Events_Proportions_Table)


## Column names for all the Matrices storing the MFIs
colnames(MFI.mean.OriginalData) <- c(sapply(1:(ncol(f@exprs)-1), function(x){paste("Singlets:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("Live:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("Lymphocytes:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("CD45+:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("NOT(Granulocyte Pre):", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("Granulocyte Pre:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("CD3 T-cells:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("NOT(CD3 T-cells):", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("Plasma:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("NOT(Plasma):", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("B-cells:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("Myeloid Pre:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("CD43+:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("CD43-:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("HFA:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("HFB:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("HFC:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("HFD:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("HFE:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                     sapply(1:(ncol(f@exprs)-1), function(x){paste("HFF:", sub(" .*", "", f@parameters@data$desc[[x]]))}))

# Combining the FCS files names, Genotype, The Number of cells before cleaning, number of channels, and the Assay Date with the MFI.mean.OriginalData
MFI.mean.OriginalData <- cbind(store.allFCS, MFI.mean.OriginalData)

colnames(MFI.median.OriginalData) <-c(sapply(1:(ncol(f@exprs)-1), function(x){paste("Singlets:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("Live:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("Lymphocytes:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("CD45+:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("NOT(Granulocyte Pre):", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("Granulocyte Pre:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("CD3 T-cells:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("NOT(CD3 T-cells):", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("Plasma:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("NOT(Plasma):", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("B-cells:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("Myeloid Pre:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("CD43+:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("CD43-:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("HFA:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("HFB:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("HFC:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("HFD:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("HFE:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                      sapply(1:(ncol(f@exprs)-1), function(x){paste("HFF:", sub(" .*", "", f@parameters@data$desc[[x]]))}))

# Combining the FCS files names, Genotype, The Number of cells before cleaning, number of channels, and the Assay Date with the MFI.median.OriginalData
MFI.median.OriginalData <- cbind(store.allFCS, MFI.median.OriginalData)

colnames(MFI.mean.TransformedData) <- c(sapply(1:(ncol(f@exprs)-1), function(x){paste("Singlets:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("Live:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("Lymphocytes:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("CD45+:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("NOT(Granulocyte Pre):", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("Granulocyte Pre:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("CD3 T-cells:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("NOT(CD3 T-cells):", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("Plasma:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("NOT(Plasma):", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("B-cells:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("Myeloid Pre:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("CD43+:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("CD43-:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("HFA:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("HFB:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("HFC:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("HFD:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("HFE:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                        sapply(1:(ncol(f@exprs)-1), function(x){paste("HFF:", sub(" .*", "", f@parameters@data$desc[[x]]))}))

# Combining the FCS files names, Genotype, The Number of cells before cleaning, number of channels, and the Assay Date with the MFI.mean.TransformedData
MFI.mean.TransformedData <- cbind(store.allFCS, MFI.mean.TransformedData)

colnames(MFI.median.TransformedData) <- c(sapply(1:(ncol(f@exprs)-1), function(x){paste("Singlets:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("Live:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("Lymphocytes:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("CD45+:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("NOT(Granulocyte Pre):", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("Granulocyte Pre:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("CD3 T-cells:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("NOT(CD3 T-cells):", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("Plasma:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("NOT(Plasma):", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("B-cells:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("Myeloid Pre:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("CD43+:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("CD43-:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("HFA:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("HFB:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("HFC:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("HFD:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("HFE:", sub(" .*", "", f@parameters@data$desc[[x]]))}),
                                          sapply(1:(ncol(f@exprs)-1), function(x){paste("HFF:", sub(" .*", "", f@parameters@data$desc[[x]]))}))

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
save(store.allFCS.unique, file = paste("Results/store.allFCS.unique.Rdata", sep = ""))

cat("Total time is: ",TimeOutput(start),"\n",sep="")

