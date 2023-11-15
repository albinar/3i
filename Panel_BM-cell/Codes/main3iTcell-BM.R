# Developed by Albina Rahim
# Date: May 24, 2016

# Gating of Bone Marrow files (August 2015+March 2016+October 2016)
remove(list=ls())
setwd("/code/Projects/3i/Panel_BM-cell/")

library("flowCore")
library("flowDensity")
library("stringr")
library("plyr")
library("doMC")

source("Codes/3iTcellfunctions-BM.R")
source("Codes/rotate.data-BM.R")
source("Codes/getPeaks.R")


results.dir <- "/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results"

load(paste0(results.dir,"/lgl.Rdata"))
load(paste0(results.dir,"/Genotype.Rdata"))
load(paste0(results.dir,"/uniqueGT.Rdata"))
load(paste0(results.dir,"/store.allFCS.Rdata"))
load(paste0(results.dir,"/ind.marg.neg.clean.all.Rdata"))
load(paste0(results.dir, "/res.clean.Rdata"))

# Create directories
suppressWarnings(dir.create(paste0(results.dir,"/Figures/ScatterPlots/")))
invisible(sapply(1:length(uniqueGT), function(x){
  suppressWarnings(dir.create(paste0(results.dir,"/Figures/ScatterPlots/", uniqueGT[x])))
}))

start <- Sys.time()

## Reading the TvsF flagged csv file which was sent to Adam for feedback
flagged.FCS <- as.matrix(read.csv(paste0(results.dir, "/flagged.FCS.BoneMarrow.Feedback.csv"), sep = ","))
## List the indices of the flagged FCS files based on TvsF algorithm
flagged.FCS.index <- which(res.clean[c('Has the file passed'),] == "F")

rownames(flagged.FCS) <- flagged.FCS.index

## Manually including some flagged files for further analysis
manualExclude.flagged.FCS <- flagged.FCS[-which(flagged.FCS[,c('Final.Action')] == "Keep"),]
manualExclude.flagged.FCS.index <- as.integer(rownames(manualExclude.flagged.FCS))

## Removing flagged files based on Adam's feedback from further analysis after recieving confirmation from Adam (Jan 04, 2017)
if ( length(manualExclude.flagged.FCS.index) != 0){
  store.allFCS <- store.allFCS[-manualExclude.flagged.FCS.index,]
  res.clean <- res.clean[,-manualExclude.flagged.FCS.index]
  ind.marg.neg.clean.all[manualExclude.flagged.FCS.index] <- NULL
}

## This part of the script was taken from Sibyl for the purpose of parallelizing the execution of this script
no_cores <- detectCores() - 2
registerDoMC(no_cores)

file.names <- data.frame(store.allFCS, stringsAsFactors = F)

print("Starting Gating & Plotting")
#props.events.gates <- ldply(1:25, function(i){

props.events.gates <- ldply(1:nrow(file.names), function(i){ 
  
  x <- file.names[i,]
  
  all.events <- matrix(nrow = 1, ncol = 21, data = NA) # matrix for saving the Event counts
  all.props <- matrix(nrow = 1, ncol = 21, data = NA) # matrix for saving the Proportions
  gthres <- matrix(nrow = 1, ncol = 19, data = NA) # matrix for saving the thresholds
  
  tryCatch({
    
    # Load FCS files
    f <- read.FCS(filename = paste0(x$Path,"/", x$FCS.file))
    
    # Compensate, remove indices, and then transform
    cat("Starting Compensation.\n")
    if( det(f@description$SPILL)==1 ){
      message("Check the spillover matrix, it's probably an identity matrix!")
    } else {
      f <- compensate(f, f@description$SPILL)
    }
    if(length(ind.marg.neg.clean.all[[i]]$ind.marg.neg) != 0){
      f <- f[-(ind.marg.neg.clean.all[[i]]$ind.marg.neg),]
    }
    #f.marg.neg <- f
    if(length(ind.marg.neg.clean.all[[i]]$ind.clean) != 0){
      f <- f[-(ind.marg.neg.clean.all[[i]]$ind.clean),]
    }
    f <- transform(f, lgl)
    #f.marg.neg <- f
    
    ###############################################################################################
    ## Cell counts/proportions after cleaning and margin event removal
    
    all.events[1] <- nrow(f)
    all.props[1] <- (nrow(f)/nrow(f))*100 
    
    ###############################################################################################
    
    ## Gating All Events -----
    ## Gating All Events to obtain the Singlets. Plotting SSC-A_FSC-W
    singlets.gate.h <- deGate(f, channel = c("FSC-W"), use.percentile = T, percentile = 0.85)
    singlets.gate.l <- deGate(f, channel = c("FSC-W"), use.percentile = T, percentile = 0.0001)
    #deGate(f, channel = c("FSC-W"), upper = F)
    names(singlets.gate.h) <- names(singlets.gate.l) <- identifier(f)
    
    singlets.flowD.h <- flowDensity(f, channels = c("SSC-A", "FSC-W"), position = c(NA, F), gates = c(NA, singlets.gate.h[identifier(f)]))
    singlets.flowD.l <- flowDensity(singlets.flowD.h, channels = c("SSC-A", "FSC-W"), position = c(NA, T), gates = c(NA, singlets.gate.l[identifier(f)]))
    singlets.flowD.l@proportion <- (singlets.flowD.l@cell.count/nrow(f))*100
    singlets.flowD.l.ind <- singlets.flowD.l@index
    singlets <- getflowFrame(singlets.flowD.l)
    
    gthres[1] <- singlets.gate.h
    gthres[2] <- singlets.gate.l
    all.events[2] <- singlets.flowD.l@cell.count
    all.props[2] <- singlets.flowD.l@proportion
    ######################################################################################################
    
    ## Gating Singlets to obtain the Live population. Plotting Live/Dead_SSC-A -----
    live.gate <- deGate(singlets, channel = c(11), use.upper = T, upper = T, tinypeak.removal = 0.90)
    names(live.gate) <- identifier(singlets)
    live.flowD <- flowDensity(singlets, channels = c(11,4), position = c(F,NA), gates = c(live.gate,NA))
    live <- getflowFrame(live.flowD)
    
    gthres[3] <- live.gate
    all.events[3] <- live.flowD@cell.count
    all.props[3] <- live.flowD@proportion
    
    ##############################################################################################
    
    ## Gating Live to obtain the Lymphocyte population. Plotting FSC-A_SSC-A
    peaks <- getPeaks(live, c(1), tinypeak.removal=0.01)
    numPeaks <- length(peaks$Peaks)
    if(numPeaks > 2){
      fsc.a.gate.low <- deGate(live, channel = c(1), tinypeak.removal = 0.01)
    } else{
      fsc.a.gate.low <- deGate(live, channel = c(1), use.upper = T, upper = FALSE)
    }

    fsc.a.gate.high <- deGate(live, channel = c(1), use.percentile = T, percentile = 0.99999)
    ssc.a.gate <- deGate(live, channel = c(4), use.percentile = T, percentile = 0.9999999)
    names(fsc.a.gate.low) <- names(fsc.a.gate.high) <- names(ssc.a.gate) <- identifier(live)
    
    lymph.flowD.temp <- flowDensity(live, channels = c(1,4), position = c(T, F), gates = c(fsc.a.gate.low, ssc.a.gate))
    lymph.flowD <- flowDensity(lymph.flowD.temp, channels = c(1,4), position = c(F,F), gates = c(fsc.a.gate.high, ssc.a.gate))
    lymph.flowD@proportion <- (lymph.flowD@cell.count/live.flowD@cell.count)*100
    lymph <- getflowFrame(lymph.flowD)
    
    gthres[4] <- fsc.a.gate.high
    gthres[5] <- fsc.a.gate.low
    gthres[6] <- ssc.a.gate
    all.events[4] <- lymph.flowD@cell.count
    all.props[4] <- lymph.flowD@proportion
    
    #############################################################################################
    
    ## Gating Lymphocytes to obtain the CD45+ population. Plotting CD45_CD43 (CD45+ population)
     
    theta = atan(1)
    lymph.temp <- lymph
    #lymph.temp <-rotate.data(lymph,c(14,8),theta = pi/6.50)$data
    lymph.temp <-rotate.data(lymph,c(14,8),theta = pi)$data
    cd45.gate.low <- deGate(lymph.temp, channel = c(14), use.upper = T, upper = F)-0.5
    #cd45.gate.low <- deGate(lymph.temp, channel=c(14),upper = T, tinypeak.removal = 0.9, graph = T)
   
    # cd45.gate.low.temp1 <- min(deGate(lymph.temp, channel = c(14), upper = F, tinypeak.removal = 0.1, all.cut = T, graph = T))
    # cd45.gate.low <- cd45.gate.low.temp1
    # print(paste("List of CD45gate.low:", deGate(lymph.temp, channel = c(14), upper = F, tinypeak.removal = 0.1, all.cut = T)))
    # 
    # if(cd45.gate.low.temp1 > 2.6){
    #     print("Recalculating CD45.gate.low")
    #     cd45.gate.low.temp2 <- min(deGate(lymph.temp, channel = c(14), upper = F, tinypeak.removal = 0.1))
    #     cd45.gate.low.temp3 <- deGate(lymph.temp, channel = c(14))
    #     print(paste("List of CD45gate.low:", cd45.gate.low.temp2 ,"and", cd45.gate.low.temp3))
    #     if((cd45.gate.low.temp2 < 1) | (cd45.gate.low.temp3 < 1)){
    #         cd45.gate.low <- mean(c(cd45.gate.low.temp2, cd45.gate.low.temp3))
    #     }else{
    #         cd45.gate.low <- min(cd45.gate.low.temp2, cd45.gate.low.temp3)
    #     }
    # }
    # 
    # # if(cd45.gate.low.temp1 > 2.6){
    # #       print("Recalculating CD45.gate.low")
    # #       cd45.gate.low.temp2 <- min(deGate(lymph.temp, channel = c(14), upper = F, tinypeak.removal = 0.1))
    # #       cd45.gate.low <- cd45.gate.low.temp2
    # #       print(paste("List of CD45gate.low:", deGate(lymph.temp, channel = c(14), upper = F, tinypeak.removal = 0.1)))
    # #       if(cd45.gate.low.temp2 > 2.6){
    # #         print("CD45.low.gate is still > 2.6. Recalculating again.")
    # #         cd45.gate.low.temp3 <- deGate(lymph.temp, channel = c(14))
    # #         cd45.gate.low <- cd45.gate.low.temp3
    # #         print(paste("List of CD45gate.low:", deGate(lymph.temp, channel = c(14))))
    # #         if(cd45.gate.low.temp3 > 2.6){
    # #           cd45.gate.low <- min(cd45.gate.low.temp1, cd45.gate.low.temp2, cd45.gate.low.temp3)
    # #         }
    # #       }
    # # }
    # print(paste("Final CD45gate.low:", cd45.gate.low))
    # remove(cd45.gate.low.temp1, cd45.gate.low.temp2, cd45.gate.low.temp3)
      #median(c(deGate(lymph.temp, channel = c(14), upper = F, tinypeak.removal = 0.1, all.cut = T), deGate(lymph.temp, channel = c(14), all.cut = T)))
    
      #median(c(deGate(lymph.temp, channel = c(14), upper = F, tinypeak.removal = 0.1, all.cut = T)))
     
    cd45.gate.high <- deGate(lymph.temp, channel=c(14),all.cut=T)
    # cd45.gate.high <- deGate(lymph.temp, channel = c(14), use.percentile = TRUE, percentile = 0.9999, graph = T)
    cd43.gate.high <- deGate(lymph.temp, channel = c(8), use.percentile = T, percentile = 0.99999)
    
      
    names(cd45.gate.low) <- names(cd45.gate.high) <-  names(cd43.gate.high) <- identifier(lymph.temp)
    
    cd45.flowD.temp <- flowDensity(lymph.temp, channels = c(14,8), position = c(T,F), gates = c(cd45.gate.low, cd43.gate.high))
    cd45.flowD <- flowDensity(cd45.flowD.temp, channels = c(14,8), position = c(F,F), gates = c(cd45.gate.high, cd43.gate.high))
    
    # cd45.flowD@filter <-  rotate.data(cd45.flowD@filter,c(14,8),theta = -pi/6.50)$data
    # cd45.flowD@flow.frame <- rotate.data(getflowFrame(cd45.flowD),c(14,8),theta = -pi/6.50)$data
    cd45.flowD@filter <-  rotate.data(cd45.flowD@filter,c(14,8),theta = -pi)$data
    cd45.flowD@flow.frame <- rotate.data(getflowFrame(cd45.flowD),c(14,8),theta = -pi)$data
    cd45.flowD@proportion <- (cd45.flowD@cell.count/lymph.flowD@cell.count)*100
    
    print(paste("CD45 proportion:",  cd45.flowD@proportion))
    #plotDens(lymph.temp, channels = c(14,8)); abline(v=cd45.gate.low, lwd=2); abline(v=cd45.gate.high, lwd=2); abline(h=cd43.gate.high, lwd=2)
    
    # ## NOTE: I had to add this condition because for the Lymphocyte population, it was seen that there was a small population sometimes
    # ## with the CD45+ population, which I wanted to grab. So in order to grab the population and not cut off too many cells
    # ## I added the following if-else condition.
    # #if(ceiling(cd45.flowD@proportion) < 85.00){
    #   if(cd45.flowD@proportion < 80.00){
    #       print("CD45 proportion is less than 80%")
    #   
    #       cd45.gate.low.temp1 <- deGate(lymph.temp, channel = c(14))
    #       
    #       if(cd45.gate.low <= cd45.gate.low.temp1) {
    #         print("Recalculate CD45.low.gate")
    #         cd45.gate.low <- deGate(lymph.temp, channel = c(14), upper = F, tinypeak.removal = 0.1)
    #       } else{
    #         cd45.gate.low <- cd45.gate.low.temp1
    #       }
    #       
    #       print(paste("New Final CD45gate.low:", cd45.gate.low))
    #       
    #       names(cd45.gate.low) <- identifier(lymph.temp)
    #       
    #       cd45.flowD.temp <- flowDensity(lymph.temp, channels = c(14,8), position = c(T,F), gates = c(cd45.gate.low, cd43.gate.high))
    #       cd45.flowD <- flowDensity(cd45.flowD.temp, channels = c(14,8), position = c(F,F), gates = c(cd45.gate.high, cd43.gate.high))
    #       cd45.flowD@filter <-  rotate.data(cd45.flowD@filter,c(14,8),theta = -pi/6.50)$data
    #       cd45.flowD@flow.frame <- rotate.data(getflowFrame(cd45.flowD),c(14,8),theta = -pi/6.50)$data
    #       cd45.flowD@proportion <- (cd45.flowD@cell.count/lymph.flowD@cell.count)*100
    #       print(paste("New CD45 proportion:",  cd45.flowD@proportion))
    # }
    
    cd45 <- getflowFrame(cd45.flowD)
    
    gthres[7] <- cd45.gate.low
    gthres[8] <- cd45.gate.high
    gthres[9] <- cd43.gate.high
    all.events[5] <- cd45.flowD@cell.count
    all.props[5] <- cd45.flowD@proportion
    
    #############################################################################################
    
    ## Gating CD45+ to obtain Granulocyte Pre and NOT(Granulocyte Pre) populations. Plotting GR1_CD43
    #gr1.gate <- deGate(cd45, channel = c(10))+0.1
    # peaks <- getPeaks(cd45, c(1), tinypeak.removal=0.1)
    # numPeaks <- length(peaks$Peaks)
    gr1.gate <- deGate(cd45, channel = c(10))
      #deGate(cd45, channel=c(10), use.percentile = T, percentile = 0.51)
    names(gr1.gate) <- identifier(cd45)
    gthres[10] <- gr1.gate
    
    NOT.granulocyte.flowD <- flowDensity(cd45, channels = c(10, 8), position = c(F,NA), gates = c(gr1.gate, cd43.gate.high))
    NOT.granulocyte <- getflowFrame(NOT.granulocyte.flowD)
    all.events[6] <- NOT.granulocyte.flowD@cell.count
    all.props[6] <- NOT.granulocyte.flowD@proportion
    
    granulocyte <- cd45
    granulocyte@exprs <- granulocyte@exprs[-NOT.granulocyte.flowD@index,]
    all.events[7] <- nrow(granulocyte)
    all.props[7] <- (nrow(granulocyte)/nrow(cd45))*100
    
    #############################################################################################
    
    ## Gating NOT(Granulocyte Pre) population to obtain CD3 T-cells and NOT(CD3 T-cells). Plotting B220_CD3
    b220.gate <- min(deGate(NOT.granulocyte, channel = c(18), all.cut = T))
    #cd3.gate <- min(deGate(NOT.granulocyte, channel = c(16), all.cut = T))-0.15
    cd3.gate <- min(deGate(NOT.granulocyte, channel = c(16), all.cut = T))
    
    names(b220.gate) <- names(cd3.gate) <- identifier(NOT.granulocyte)
    gthres[11] <- b220.gate
    gthres[12] <- cd3.gate
    
    cd3.Tcell.flowD <- flowDensity(NOT.granulocyte, channels = c(18,16), position = c(F, T), gates = c(b220.gate, cd3.gate))
    all.events[8] <- cd3.Tcell.flowD@cell.count
    all.props[8] <- cd3.Tcell.flowD@proportion
    
    NOT.cd3.Tcell <- NOT.granulocyte
    NOT.cd3.Tcell@exprs <- NOT.cd3.Tcell@exprs[-cd3.Tcell.flowD@index,]
    all.events[9] <- nrow(NOT.cd3.Tcell)
    all.props[9] <- (nrow(NOT.cd3.Tcell)/nrow(NOT.granulocyte))*100
    
    ##############################################################################################
    
    ## Gating NOT(CD3 T-cell) population to obtain Plasma and NOT Plasma. Plotting CD138_B220
    #cd138.gate <- deGate(NOT.cd3.Tcell, channel = c(15), use.percentile = T, percentile = 0.94)
    cd138.gate <- deGate(NOT.cd3.Tcell, channel = c(15), use.upper = T, upper = T)-0.5
    names(cd138.gate) <- identifier(NOT.cd3.Tcell)
    gthres[13] <- cd138.gate
    
    plasma.flowD <- flowDensity(NOT.cd3.Tcell, channels = c(15,18), position = c(T, F), gates = c(cd138.gate, b220.gate))
    all.events[10] <- plasma.flowD@cell.count
    all.props[10] <- plasma.flowD@proportion
    
    NOT.plasma <- NOT.cd3.Tcell
    NOT.plasma@exprs <-  NOT.plasma@exprs[-plasma.flowD@index,]
    all.events[11] <- nrow(NOT.plasma)
    all.props[11] <- (nrow(NOT.plasma)/nrow(NOT.cd3.Tcell))*100
    
    ##############################################################################################
    
    ## Gating NOT Plasma population to obtain Myeloid Pre and B-cells. Plotting B220_CD11b
    cd11b.gate <- deGate(NOT.plasma, channel = c(13), use.percentile = T, percentile = 0.4)
    names(cd11b.gate) <- identifier(NOT.plasma)
    gthres[14] <- cd11b.gate
    
    #Bcell.flowD <- flowDensity(NOT.plasma, channels = c(18,13), position = c(T, NA), gates = c(b220.gate-0.5, cd11b.gate))
    Bcell.flowD <- flowDensity(NOT.plasma, channels = c(18,13), position = c(T, NA), gates = c(b220.gate, cd11b.gate))
    Bcell <- getflowFrame(Bcell.flowD)
    all.events[12] <- Bcell.flowD@cell.count
    all.props[12] <- Bcell.flowD@proportion
    
    #myeloid.flowD <- flowDensity(NOT.plasma, channels = c(18,13), position = c(F, NA), gates = c(b220.gate-0.5, cd11b.gate))
    myeloid.flowD <- flowDensity(NOT.plasma, channels = c(18,13), position = c(F, NA), gates = c(b220.gate, cd11b.gate))
    all.events[13] <- myeloid.flowD@cell.count
    all.props[13] <- myeloid.flowD@proportion
    
    ##############################################################################################
    
    ## Gating B-cells population to obtain CD43+ and CD43-. Plotting B220_CD43
    cd43.gate <-  deGate(Bcell, channel = c(8), use.percentile = T, percentile = 0.64)
    names(cd43.gate) <- identifier(Bcell)
    gthres[15] <- cd43.gate
    
    #deGate(Bcell, channel = c(8), upper = T, tinypeak.removal = 0.9)
    cd43.pos.flowD <- flowDensity(Bcell, channels = c(18,8), position = c(NA, T), gates = c(b220.gate, cd43.gate))
    cd43.pos <- getflowFrame(cd43.pos.flowD)
    all.events[14] <- cd43.pos.flowD@cell.count
    all.props[14] <- cd43.pos.flowD@proportion
    
    cd43.neg <- Bcell
    cd43.neg@exprs <- cd43.neg@exprs[-cd43.pos.flowD@index,]
    all.events[15] <- nrow(cd43.neg)
    all.props[15] <- (nrow(cd43.neg)/Bcell.flowD@cell.count)*100
    
    ##############################################################################################
    
    ## Gating CD43+ population to obtain HFA, HFB, and HFC. Plotting CD24_BP1
    
    cd24.gate <- deGate(cd43.pos, channel = c(9), tinypeak.removal = 0.001)
    if(cd24.gate < 1){
      print(paste("CD24.gate: ", cd24.gate))
      print("Recalculate CD24.gate")
      cd24.gate <- deGate(cd43.pos, channel = c(9))
    }
    print(paste("Final CD24.gate: ", cd24.gate))
    
    bp1.gate <- deGate(cd43.pos, channel = c(17), use.percentile = T, percentile = 0.65)
    names(cd24.gate) <- names(bp1.gate) <- identifier(cd43.pos)
    gthres[16] <- cd24.gate
    gthres[17] <- bp1.gate
    
    
    #deGate(cd43.pos, channel = c(9), use.percentile = T, percentile = 0.3)
    #deGate(cd43.pos, channel = c(9), use.percentile = T, percentile = 0.3)
    #deGate(cd43.pos, channel = c(17), upper = T, tinypeak.removal = 0.9)
   
    HFA.flowD <- flowDensity(cd43.pos, channels = c(9, 17), position = c(F,F), gates = c(cd24.gate, bp1.gate))
    all.events[16] <- HFA.flowD@cell.count
    all.props[16] <- HFA.flowD@proportion
    
    HFB.flowD <- flowDensity(cd43.pos, channels = c(9, 17), position = c(T,F), gates = c(cd24.gate, bp1.gate))
    all.events[17] <- HFB.flowD@cell.count
    all.props[17] <- HFB.flowD@proportion
    
    HFC.flowD <- flowDensity(cd43.pos, channels = c(9, 17), position = c(T,T), gates = c(cd24.gate, bp1.gate))
    all.events[18] <- HFC.flowD@cell.count
    all.props[18] <- HFC.flowD@proportion
    
    ##############################################################################################
    
    ## Gating CD43- population to obtain HFD, HFE, and HFF. Plotting IgM_IgD
    
    igm.gate <- deGate(cd43.neg, channel = c(12))+0.2
    igd.gate <- deGate(cd43.neg, channel = c(7), use.upper = T, upper = T, tinypeak.removal = 0.9)+0.2
    names(igm.gate) <- names(igd.gate) <- identifier(cd43.neg)
    gthres[18] <- igm.gate
    gthres[19] <- igd.gate
    
    HFD.flowD <- flowDensity(cd43.neg, channels = c(12,7), position = c(F,F), gates = c(igm.gate, igd.gate))
    all.events[19] <- HFD.flowD@cell.count
    all.props[19] <- HFD.flowD@proportion
    
    HFE.flowD <- flowDensity(cd43.neg, channels = c(12,7), position = c(T,F), gates = c(igm.gate, igd.gate))
    all.events[20] <- HFE.flowD@cell.count
    all.props[20] <- HFE.flowD@proportion
    
    HFF.flowD <- flowDensity(cd43.neg, channels = c(12,7), position = c(NA,T), gates = c(igm.gate, igd.gate))
    all.events[21] <- HFF.flowD@cell.count
    all.props[21] <- HFF.flowD@proportion
    
    ##############################################################################################
    ##############################################################################################
    
    # ## Saving all the Gates
    # gthres <- c(singlets.gate.h, singlets.gate.l, live.gate, fsc.a.gate.high, fsc.a.gate.low, ssc.a.gate, cd45.gate.low, cd45.gate.high, cd43.gate.high,
    #             gr1.gate, b220.gate, cd3.gate, cd138.gate, cd11b.gate, cd43.gate, cd24.gate, igm.gate, igd.gate, bp1.gate)
    # save ( gthres, file =  paste("Results/Gating_Thresholds/", store.allFCS[i,1],"/Gthres_", store.allFCS[i,2],".Rdata",sep="") )
    # 
    # ## Saving the CD45 gates separately
    # gates.cd45 <- cbind(igd.gate, cd43.gate, cd24.gate, gr1.gate, live.gate, igm.gate, cd11b.gate, cd45.gate.low, cd138.gate, cd3.gate, bp1.gate, b220.gate)
    # colnames(gates.cd45) <- c("IgD", "CD43", "CD24", "GR1", "Live", "IgM", "CD11b", "CD45", "CD138", "CD3", "BP1", "B220")
    #   
    # na.ind <-which(apply(gates.cd45,2,function(x) all(is.na(x))))
    # if(length(na.ind)!=0)
    # {
    #   gates.cd45 <- gates.cd45[,-na.ind]
    # }
    # save ( gates.cd45, file =  paste("Results/Gating_Thresholds_CD45/", store.allFCS[i,1],"/Gthres_", store.allFCS[i,2],".Rdata",sep="") )
    
    
    ##############################################################################################
    ##############################################################################################
    
    ## Saving the Plots
    #--------Start Big Png------------
    png ( file = paste0(results.dir,"/Figures/ScatterPlots/", x$Genotype, "/", "Total_", x$FCS.files, ".png"), width=2100, height=2100*4/4)
    #png ( file = paste("Results/tempPlots/", "Total_", store.allFCS[i,2], ".png", sep = "" ), width=2100, height=2100*4/4)
    
    
    par(mfrow=c(4,4),mar=(c(5, 5, 4, 2) + 0.1))
    
    # 1st method of plotting FSC-A_SSC-W
    plotDens(f,c("SSC-A","FSC-W"),main="Ungated", devn = F, cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(singlets.flowD.l@filter,lwd=2)
    
    # Plotting Live/Dead_SSC-A 
    plotDens(singlets,  c(11,4), main= "Singlets", devn = F, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=live.gate, lwd=2); 
    
    # Plotting FSC-A_SSC-A
    plotDens(live, c(1,4), main = "Live", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=fsc.a.gate.low, lwd=2); abline(v=fsc.a.gate.high, lwd=2); abline(h=ssc.a.gate, lwd=2)
    
    # Plotting CD45_CD43 to illustrate the Lymphocyte population
    plotDens(lymph, channels = c(14,8), main = "Lymphocytes", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    
    # Plotting CD45_CD43 to obtain the CD45+ population from the rotated Lymphocyte population
    plotDens(lymph.temp, channels = c(14,8), main = "Lymphocytes (Rotated)", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=cd45.gate.low, lwd=2); abline(v=cd45.gate.high, lwd=2); abline(h=cd43.gate.high, lwd=2)
    
    # plotDens(lymph, channels = c(14,8), main = "Lymphocytes_Rotated", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    # points(cd45points, type="l", col="black", lty=2, lwd=2)  
    
    # Plotting GR1_CD43
    min.y <- min(exprs(cd45)[,8])
    max.y <- max(exprs(cd45)[,8])+0.5
    plotDens(cd45, channels = c(10,8),  main = "CD45+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y))
    
    # Plotting GR1_CD43 to obtain the Granulocyte Pre and NOT(Granulocyte Pre) populations
    min.y <- min(exprs(cd45)[,8])
    max.y <- max(exprs(cd45)[,8])+0.5
    plotDens(cd45, channels = c(10,8),  main = "CD45+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, ylim = c(min.y, max.y)); abline(v=gr1.gate, lwd=2)
    
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
    plotDens(NOT.plasma, channels = c(18, 13),  main = "NOT Plasma", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=b220.gate, lwd=2); #abline(v=b220.gate-0.5, lwd=2)
    
    
    # Plotting B220_CD43 to obtain CD43+ and CD43-
    min.x <- min(exprs(Bcell)[,18])
    max.x <- max(exprs(Bcell)[,18])
    min.y <- min(exprs(Bcell)[,8])
    max.y <-  max(exprs(Bcell)[,8])
    plotDens(Bcell, channels = c(18, 8),  main = "B-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x-1.5, max.x)); abline(h=cd43.gate, lwd=2); abline(v=min.x, lwd=2)
    
    # Plotting CD24_BP1 to obtain HFA, HFB, and HFC. 
    min.x <- min(exprs(cd43.pos)[,9])
    max.x <- max(exprs(cd43.pos)[,9])
    min.y <- min(exprs(cd43.pos)[,17])
    max.y <-  max(exprs(cd43.pos)[,17])
    plotDens(cd43.pos, channels = c(9, 17),  main = "CD43+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=cd24.gate, lwd=2); abline(h=bp1.gate, lwd=2)
  
    
    # Plotting IgM_IgD to obtain HFD, HFE, and HFF. 
    min.y <- min(exprs(cd43.neg)[,7])
    plotDens(cd43.neg, channels = c(12, 7),  main = "CD43-", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(h=igd.gate, lwd=2)
    lines(x=c(igm.gate, igm.gate), y=c(igd.gate, min.y-1), lwd=2)
    
    # # Plotting the TvsF Worst Fluorescence Channel with the parts removed by TvsF
    # plotDens(f, channels = c(19, as.numeric(res.clean["Worst Channel",i])), cex.main=2, cex.lab=2, cex.axis=2, main=paste0("TvsF Removal"))
    # points(f@exprs[ind.marg.neg.clean.all[[i]]$ind.clean, c(19,as.numeric(res.clean["Worst Channel",i]))],col=1,pch=".")
    
    
    dev.off()
    par(mfrow=c(1,1))
    #--------End Big Png------------
    
  },error = function(err) {
    print(paste("ERROR in Gating/Plotting Message:  ",err, "Index: ", i));
    # errorFileIndex <- c(errorFileIndex,i)
    # return(errorFileIndex)
    return(0)
  }
  ) # end of tryCatch
  
  data.frame(x$FCS.files, all.props, all.events, gthres)
  #data.frame(all.props, all.events, gthres)

#   cat("Time is: ",TimeOutput(start2),"\n",sep="")
# } # end of for-loop

}, .parallel = TRUE) # end ldply

#cat("Time is: ",TimeOutput(start),"\n",sep="")
print("Finished Gating & Plotting")

# Saving the big dataframe of Proportions, Events, and Gating thresholds 
save(props.events.gates, file = paste0(results.dir,"/PropsEventsGates0.Rdata"))  

# Checking if any files failed the gating by finding if there is any NA
failedGating.files.index <- which(is.na(props.events.gates[,ncol(props.events.gates)]))

if(length(failedGating.files.index) != 0){
  print("Files which failed Gating:")
  print(props.events.gates[failedGating.files.index,'x.FCS.files'])
  failedGating.files <- store.allFCS[failedGating.files.index,]
  rownames(failedGating.files) <- failedGating.files.index
  save(failedGating.files, file = paste0(results.dir,"/failedGating.files.Rdata"))  
  ## For sending a spreadsheet containing list of FCS files which Failed the Gating to center
  write.csv(failedGating.files, file = paste0(results.dir, "/failedGating.BoneMarrow.csv"), row.names = FALSE)
}else{
  print("No files failed the Gating")
}


# Extracting the Proportions from the large dataframe that was returned drom ldply 
all.props <- props.events.gates[,2:22] 
all.props <- as.matrix(all.props)
colnames(all.props) <- c("All Events-%Parent", "Singlets-%Parent", "Live-%Parent", "Lymphocytes-%Parent", "CD45+-%Parent", "NOT(Granulocyte Pre)-%Parent", "Granulocyte Pre-%Parent",
                         "CD3 T-cell-%Parent", "NOT(CD3 T-cell)-%Parent", "Plasma-%Parent", "NOT Plasma-%Parent", "B-cells-%Parent", "Myeloid Pre-%Parent", "CD43+ %Parent", "CD43- %Parent",
                         "HFA-%Parent", "HFB-%Parent", "HFC-%Parent", "HFD-%Parent", "HFE-%Parent", "HFF-%Parent")

# Extracting the #Events from the large dataframe that was returned drom ldply 
all.events <- props.events.gates[,23:43]
all.events <- as.matrix(all.events)
colnames(all.events) <- c("All Events-#Events", "Singlets-#Events", "Live-#Events", "Lymphocytes-#Events", "CD45+-#Events", "NOT(Granulocyte Pre)-#Events", "Granulocyte Pre-#Events",
                          "CD3 T-cell-#Events", "NOT(CD3 T-cell)-#Events", "Plasma-#Events", "NOT Plasma-#Events", "B-cells-#Events", "Myeloid Pre-#Events", "CD43+ #Events", "CD43- #Events",
                          "HFA-#Events", "HFB-#Events", "HFC-#Events", "HFD-#Events", "HFE-#Events", "HFF-#Events")

# Extracting the Gating Thresholds from the large dataframe that was returned drom ldply
all.gthres <- props.events.gates[,44:ncol(props.events.gates)]
all.gthres <- as.matrix(all.gthres)
colnames(all.gthres) <- c("singlets.gate.h", "singlets.gate.l", "live.gate", "fsc.a.gate.high", "fsc.a.gate.low", "ssc.a.gate", "cd45.gate.low", "cd45.gate.high", "cd43.gate.high",
             "gr1.gate", "b220.gate", "cd3.gate", "cd138.gate", "cd11b.gate", "cd43.gate", "cd24.gate", "bp1.gate", "igm.gate", "igd.gate")
## Saving the Gating Thresholds as numerical data, since we will need this for detecting gating threshold outliers
save(all.gthres, file = paste0(results.dir, "/all.gthres.Rdata"))
all.gthres <- cbind(store.allFCS[,c('Panel/Organ', 'Genotype', 'FCS files', 'Barcodes', 'Assay Date', 'Gender', 'Number of Channels', 'Number of Cells')], all.gthres)


## Combining Events & Proportions into one large matrix
Events_Proportions_Table <- NULL
for(i in 1:ncol(all.props)){
  Events_Proportions_Table <- cbind(Events_Proportions_Table,all.events[,i],all.props[,i])
}

colnames(Events_Proportions_Table) <- c("All Events-#Events", "All Events-%Parent", "Singlets-#Events", "Singlets-%Parent", "Live-#Events", "Live-%Parent", "Lymphocytes-#Events", "Lymphocytes-%Parent", 
                                        "CD45+-#Events", "CD45+-%Parent", "NOT(Granulocyte Pre)-#Events", "NOT(Granulocyte Pre)-%Parent", "Granulocyte Pre-#Events", "Granulocyte Pre-%Parent", 
                                        "CD3 T-cell-#Events", "CD3 T-cell-%Parent", "NOT(CD3 T-cell)-#Events", "NOT(CD3 T-cell)-%Parent", "Plasma-#Events", "Plasma-%Parent", "NOT Plasma-#Events", "NOT Plasma-%Parent", 
                                        "B-cells-#Events", "B-cells-%Parent", "Myeloid Pre-#Events", "Myeloid Pre-%Parent", "CD43+ #Events",  "CD43+ %Parent", "CD43- #Events", "CD43- %Parent", 
                                        "HFA-#Events", "HFA-%Parent", "HFB-#Events", "HFB-%Parent", "HFC-#Events", "HFC-%Parent", "HFD-#Events", "HFD-%Parent", "HFE-#Events", "HFE-%Parent", "HFF-#Events", "HFF-%Parent")

# Combining the store.allFCS information with the Events_Proportions Table, all.props, and all.events 
Events_Proportions_Table <- cbind(store.allFCS[,c('Panel/Organ', 'Genotype', 'FCS files', 'Barcodes', 'Assay Date', 'Gender', 'Number of Channels', 'Number of Cells')], Events_Proportions_Table)
all.props <- cbind(store.allFCS[,c('Panel/Organ', 'Genotype', 'FCS files', 'Barcodes', 'Assay Date', 'Gender', 'Number of Channels', 'Number of Cells')], all.props)
all.events <- cbind(store.allFCS[,c('Panel/Organ', 'Genotype', 'FCS files', 'Barcodes', 'Assay Date', 'Gender', 'Number of Channels', 'Number of Cells')], all.events)

write.csv(all.props, file =  paste0(results.dir, "/allProportions_Table_Initial.csv"))
write.csv(all.events, file =  paste0(results.dir, "/allEvents_Table_Initial.csv"))
write.csv(all.gthres, file =  paste0(results.dir, "/allGthres_Table_Initial.csv"))
write.csv(Events_Proportions_Table, file =  paste0(results.dir, "/Events_Proportions_Table_Initial.csv"))

cat("Total time is: ",TimeOutput(start),"\n",sep="")

