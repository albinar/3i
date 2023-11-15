# Developed by Albina Rahim
# Date: March 21, 2016

# Gating of T-cell Panel (MLN) organ (August 2015 + March 2016 + October 2016) using the updated Gating Thresholds
remove(list=ls())
setwd("/code/Projects/3i/Panel_T-cell_MLN/")

library("flowCore")
library("flowBin")
library("flowDensity")
library("stringr")
library("plyr")
library("doMC")


source("Codes/3iTcellfunctions.R")
source("Codes/rotate.data.R")

results.dir <- "/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results"

load(paste0(results.dir,"/lgl.Rdata"))
load(paste0(results.dir,"/Genotype.Rdata"))
load(paste0(results.dir,"/uniqueGT.Rdata"))
load(paste0(results.dir,"/store.allFCS.Rdata"))
load(paste0(results.dir,"/ind.marg.neg.clean.all.Rdata"))
load(paste0(results.dir, "/res.clean.Rdata"))
load(paste0(results.dir, "/failedGating.files.Rdata"))
load(paste0(results.dir, "/all.gthres.Updated.Rdata"))

# Create directories
suppressWarnings(dir.create(paste0(results.dir,"/Figures/ScatterPlotsUpdated/")))
invisible(sapply(1:length(uniqueGT), function(x){
  suppressWarnings(dir.create(paste0(results.dir,"/Figures/ScatterPlotsUpdated/", uniqueGT[x])))
}))


## Reading the TvsF flagged csv file which was sent to Adam for feedback
flagged.FCS <- as.matrix(read.csv(paste0(results.dir, "/flagged.FCS.TcellMLN.Feedback.csv"), sep = ","))
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
## There are files which failed the Gating in our previous step. So we exclude those files
if (nrow(failedGating.files) != 0){
  failedGating.files.index <- as.integer(rownames(failedGating.files))
  store.allFCS <- store.allFCS[-failedGating.files.index,]
  res.clean <- res.clean[,-failedGating.files.index]
  ind.marg.neg.clean.all[failedGating.files.index] <- NULL
}

rownames(store.allFCS) <- 1:nrow(store.allFCS)

start <- Sys.time()


# # Optional: Added for checking the gating thresholds for the gd.Tcells population 
# gd.Tcells.gthres <- matrix(nrow = nrow(store.allFCS.unique), ncol = 3)
# colnames(gd.Tcells.gthres) <- c("tcrd.gate.low", "tcrd.gate.high", "cd4.gate")

## This part of the script was taken from Sibyl for the purpose of parallelizing the execution of this script
no_cores <- detectCores() - 6
registerDoMC(no_cores)

file.names <- data.frame(store.allFCS, stringsAsFactors = F)

print("Starting Gating & Plotting")

#props.events <- ldply(1:10, function(i){ 
props.events <- ldply(1:nrow(file.names), function(i){ 
  
  x<- file.names[i,]
  
  all.events <- matrix(nrow = 1, ncol = 48, data = NA)# Previously it was 45
  all.props <- matrix(nrow = 1, ncol = 53, data = NA)# Previously it was 45
  
  
  tryCatch({
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
    f.marg.neg <- f
    if(length(ind.marg.neg.clean.all[[i]]$ind.clean) != 0){
      f <- f[-(ind.marg.neg.clean.all[[i]]$ind.clean),]
    }
    f <- transform(f, lgl)
    
    ###############################################################################################
    ## Cell counts/proportions after cleaning and margin event removal

    all.events[1] <- nrow(f)
    all.props[1] <- (nrow(f)/nrow(f))*100
    
    ###############################################################################################
    
    ## Gating All Events -----
    ## 1st method of plotting FSC-A_SSC-W
    singlets.flowD.h <- flowDensity(f, channels = c("FSC-A", "SSC-W"), position = c(NA, F), gates = c(NA, all.gthres.Updated[i,c('singlets.gate.h')]))
    singlets.flowD.l <- flowDensity(singlets.flowD.h, channels = c("FSC-A", "SSC-W"), position = c(NA, T), gates = c(NA, all.gthres.Updated[i,c('singlets.gate.l')]))
    singlets.flowD.l@proportion <- (singlets.flowD.l@cell.count/nrow(f))*100
    singlets <- getflowFrame(singlets.flowD.l)
    
    all.events[2] <- singlets.flowD.l@cell.count
    all.props[2] <- singlets.flowD.l@proportion 
    
    ######################################################################################################
    
    ## Gating Singlets. Plotting Live/Dead_SSC-A -----
    live.flowD <- flowDensity(singlets, channels = c(11,4), position = c(F,NA), gates = c(all.gthres.Updated[i,c('live.gate')],NA))
    live <- getflowFrame(live.flowD)
    
    all.events[3] <- live.flowD@cell.count
    all.props[3] <- live.flowD@proportion
    ##############################################################################################
    
    ## Gating Live. Plotting FSC-A_SSC-A
    lymph.flowD.temp <- flowDensity(live, channels = c(1,4), position = c(T, F), gates = c(all.gthres.Updated[i,c('fsc.a.gate.low')], all.gthres.Updated[i,c('ssc.a.gate')]))
    lymph.flowD <- flowDensity(lymph.flowD.temp, channels = c(1,4), position = c(F,F), gates = c(all.gthres.Updated[i,c('fsc.a.gate.high')], all.gthres.Updated[i,c('ssc.a.gate')]))
    lymph.flowD@proportion <- (lymph.flowD@cell.count/live.flowD@cell.count)*100
    lymph <- getflowFrame(lymph.flowD)
    
    all.events[4] <- lymph.flowD@cell.count
    all.props[4] <- lymph.flowD@proportion   
    
    #############################################################################################
    
    ## Gating Lymphocytes. Plotting CD45_CD161
    cd45.flowD.temp <- flowDensity(lymph, channels = c(14,15), position = c(T,F), gates = c(all.gthres.Updated[i,c('cd45.gate.low')], all.gthres.Updated[i,c('cd161.gate.high')]))
    cd45.flowD <- flowDensity(cd45.flowD.temp, channels = c(14,15), position = c(F,F), gates = c(all.gthres.Updated[i,c('cd45.gate.high')], all.gthres.Updated[i, c('cd161.gate.high')]))
    cd45.flowD@proportion <- (cd45.flowD@cell.count/lymph.flowD@cell.count)*100
    cd45 <- getflowFrame(cd45.flowD)
  
    all.events[5] <- cd45.flowD@cell.count 
    all.props[5] <- cd45.flowD@proportion
    
    
    #############################################################################################
    
    ## Gating CD45+. Highlighting the Autofluorescence. Plotting TCRd_CD4
    ## This part of the code was taken from Justin
    theta = atan(1)
    cd45s <- cd45
    cd45s<-rotate.data(cd45,c(18,16),theta = -pi/4)$data
    R <- matrix( c(cos(theta), sin(theta), -sin(theta), cos(theta)),2 ,2)
    
    returnTry <- tryCatch({
      cd45sTemp  <- flowDensity(cd45s, channels=c(18,16),position=c(T, T), gates=c(all.gthres.Updated[i,c('tcrd.gate.slanta')], all.gthres.Updated[i, c('cd4.gate.slant')]))
      cd45sTempb <- flowDensity(cd45sTemp, channels=c(18,16),position=c(F, T), gates=c(all.gthres.Updated[i, c('tcrd.gate.slantb')], all.gthres.Updated[i,c('cd4.gate.slant')]))
      cd45autoflour <- getflowFrame(cd45sTempb)
      cd45autoflour@exprs[,c(18,16)] <- t(t(R) %*% t(cd45autoflour@exprs[,c(18,16)]))
      TCRDCD4points <- t(t(R) %*% t(cd45sTempb@filter))
      
      cd45sTempb@filter <- rotate.data(cd45sTempb@filter,c(18,16),theta = pi/4)$data
      cd45sTempb@flow.frame <- rotate.data(getflowFrame(cd45sTempb),c(18,16),theta = pi/4)$data
      cd45sTempb@proportion <- (cd45sTempb@cell.count/nrow(cd45))*100
      cd45autoflour <- getflowFrame(cd45sTempb)
      
      # ## Optional:
      # autoFluor.Prop[i] <- cd45sTempb@proportion 
      # autoFluor.Events[i] <- cd45sTempb@cell.count
      
    },error = function(err) { 
      #print("Error with Normal"); 
      return(0)
    })

   
    all.events[6] <- cd45sTempb@cell.count
    all.props[6] <- cd45sTempb@proportion
    
    # plotDens(cd45, channels = c(18,16), main="CD45+_Autofluorecence", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    # if ( nrow(cd45autoflour@exprs) != 0) {
    #   points(TCRDCD4points, type="l", col="black", lty=2, lwd=2)
    # }
    ############################################################################################
    
    ## Gating CD45+. Getting rid of the Autofluorescence. Plotting TCRd_CD4
    cd45modified <- cd45
    cd45modified@exprs <- cd45modified@exprs[-cd45sTempb@index,]
    
    all.events[7] <- nrow(cd45modified)
    all.props[7] <- (nrow(cd45modified)/nrow(cd45))*100
    
    
    ############################################################################################
    
    ## Gating CD45+ (without the autofluorescence) to obtain gd T-cells. Plotting TCRd_CD4 to obtain gd T-cells
    
    gd.Tcells.flowD.temp <- flowDensity(cd45modified, channels = c(18,16), position = c(T,F), gates = c(all.gthres.Updated[i, c('tcrd.gate.low')], (all.gthres.Updated[i,c('cd4.gate')]-0.5)))
    gd.Tcells.flowD <- flowDensity(gd.Tcells.flowD.temp, channels = c(18,16), position = c(F,F), gates = c(all.gthres.Updated[i,c('tcrd.gate.high')], (all.gthres.Updated[i, c('cd4.gate')]-0.5)))
    gd.Tcells.flowD@proportion <- (gd.Tcells.flowD@cell.count/nrow(cd45modified))*100
    gd.Tcells <- getflowFrame(gd.Tcells.flowD)
    

    all.events[8] <- gd.Tcells.flowD@cell.count 
    all.props[8] <- gd.Tcells.flowD@proportion
    # ## Optional
    # gd.Tcells.Prop[i] <- nrow(gd.Tcells)/nrow(cd45modified)*100
    # gd.Tcells.Events[i] <- nrow(gd.Tcells)
    
    ############################################################################################
    
    ## Gating CD45+ (without the autofluorescence) to obtain NOT gd T-cells.
    not.gd.Tcells <- cd45modified
    not.gd.Tcells@exprs <- not.gd.Tcells@exprs[-gd.Tcells.flowD@index,]

    all.events[9] <- nrow(not.gd.Tcells)
    all.props[9] <- (nrow(not.gd.Tcells)/nrow(cd45modified))*100
    
    
    ############################################################################################

    ## Gating gd T-cells.  Plotting CD62L_CD44 to obtain gd Resting, gd Effector, and gd Naive

    gd.Resting.flowD <- flowDensity(gd.Tcells, channels = c(8,7), position = c(T,T), gates = c(all.gthres.Updated[i,c('cd62.gate')], all.gthres.Updated[i,c('cd44.gate')]))
    ##gd.Resting <- getflowFrame(gd.Resting.flowD)

    all.events[10] <- gd.Resting.flowD@cell.count
    all.props[10] <- gd.Resting.flowD@proportion

    gd.Effector.flowD <- flowDensity(gd.Tcells, channels = c(8,7), position = c(F,T), gates = c(all.gthres.Updated[i,c('cd62.gate')], all.gthres.Updated[i,c('cd44.gate')]))
    #gd.Effector <- getflowFrame(gd.Effector.flowD)

    all.events[11] <- gd.Effector.flowD@cell.count
    all.props[11] <- gd.Effector.flowD@proportion
    
    gd.Naive.flowD <- flowDensity(gd.Tcells, channels = c(8,7), position = c(T,F), gates = c(all.gthres.Updated[i,c('cd62.gate')], all.gthres.Updated[i,c('cd44.gate')]))
    #gd.Naive <- getflowFrame(gd.Naive.flowD)

    all.events[12] <- gd.Naive.flowD@cell.count
    all.props[12] <- gd.Naive.flowD@proportion
    
    #############################################################################################

    ## Gating gd T-cells. Plotting KLRG1_GITR to obtain GITR gd T-cells and GD KLRG1+

    gd.klrg1.flowD <- flowDensity(gd.Tcells, channels = c(12,17), position = c(T,F), gates = c(all.gthres.Updated[i,c('klrg1.gate')], all.gthres.Updated[i,c('gitr.gate')]))
    #gd.klrg1 <- getflowFrame(gd.klrg1.flowD)

    all.events[13] <- gd.klrg1.flowD@cell.count
    all.props[13] <- gd.klrg1.flowD@proportion

    gd.gitr.flowD <- flowDensity(gd.Tcells, channels = c(12,17), position = c(F,T), gates = c(all.gthres.Updated[i,c('klrg1.gate')], all.gthres.Updated[i,c('gitr.gate')]))
    #gd.gitr <- getflowFrame(gd.gitr.flowD)

    all.events[14] <- gd.gitr.flowD@cell.count
    all.props[14] <- gd.gitr.flowD@proportion
    
    #############################################################################################

    ## Gating gd T-cells. Plotting CD5_CD44 to obtain GD CD5+

    gd.cd5.flowD <- flowDensity(gd.Tcells, channels = c(13,7), position = c(T,NA), gates = c(all.gthres.Updated[i, c('cd5.gate')], NA))
    #gd.cd5 <- getflowFrame(gd.cd5.flowD)

    all.events[15] <- gd.cd5.flowD@cell.count
    all.props[15] <- gd.cd5.flowD@proportion
    
    #############################################################################################
    ## Gating NOT gd T-cells. Plotting CD161_CD5 to obtain CD5+ and NK-cells
    # CD5+
    cd5.flowD <- flowDensity(not.gd.Tcells, channels = c(15,13), position = c(NA,T), gates = c(NA,all.gthres.Updated[i,c('cd5.gate.new')]))
    cd5 <- getflowFrame(cd5.flowD)

    all.events[16] <- cd5.flowD@cell.count
    all.props[16] <- cd5.flowD@proportion

    # NK-cells
    NKcells.flowD.temp <- flowDensity(not.gd.Tcells, channels = c(15,13), position = c(T,F), gates = c(all.gthres.Updated[i,c('cd161.gate.low')],all.gthres.Updated[i,c('cd5.gate.new')]))
    NKcells.flowD <- flowDensity(NKcells.flowD.temp, channels = c(15,13), position = c(F,F), gates = c(all.gthres.Updated[i,c('cd161.gate.high')],all.gthres.Updated[i,c('cd5.gate.new')]))
    NKcells.flowD@proportion <- (NKcells.flowD@cell.count/nrow(not.gd.Tcells))*100
    NKcells <- getflowFrame(NKcells.flowD)

    all.events[17] <- NKcells.flowD@cell.count
    all.props[17] <- NKcells.flowD@proportion

    ###############################################################################################

    ## Gating NK-cells. Plotting CD62L_CD44 to obtain NK Effector and NK Resting/Naive
    NK.Resting.Naive.flowD <- flowDensity(NKcells, channels = c(8,7), position = c(T,NA), gates = c(all.gthres.Updated[i,c('cd62.gate')],NA))
    #NK.Resting.Naive <- getflowFrame(NK.Resting.Naive.flowD)

    all.events[18] <- NK.Resting.Naive.flowD@cell.count
    all.props[18] <- NK.Resting.Naive.flowD@proportion
    
    NK.Effector.flowD <- flowDensity(NKcells, channels = c(8,7), position = c(F,T), gates = c(all.gthres.Updated[i,c('cd62.gate')], all.gthres.Updated[i,c('cd44.gate')]))
    #NK.Effector <- getflowFrame(NK.Effector.flowD)

    all.events[19] <- NK.Effector.flowD@cell.count
    all.props[19] <- NK.Effector.flowD@proportion
    
    ###############################################################################################

    ## Gating NK-cells. Plotting KLRG1_CD44 to obtain NK KLRG1+
    NK.klrg1.flowD <- flowDensity(NKcells, channels = c(12,7), position = c(T,T), gates = c(all.gthres.Updated[i,c('klrg1.gate')], all.gthres.Updated[i,c('cd44.gate')]))
    #NK.klrg1 <- getflowFrame(NK.klrg1.flowD)

    all.events[20] <- NK.klrg1.flowD@cell.count
    all.props[20] <- NK.klrg1.flowD@proportion
    
    ###########################################################################################

    # Gating CD5+. Plotting CD161_CD4 to obtain P2a, CD4- NKT-cells, and CD4+ NKT-cells
    # P2a
    P2a.flowD <- flowDensity(cd5, channels = c(15,16), position = c(F, NA), gates = c(all.gthres.Updated[i,c('cd161.gate')],NA))
    P2a <- getflowFrame(P2a.flowD)

    all.events[21] <- P2a.flowD@cell.count
    all.props[21] <- P2a.flowD@proportion

    # CD4- NKT-cells
    cd4neg.NKTcells.flowD <- flowDensity(cd5, channels = c(15,16), position = c(T,F), gates = c(all.gthres.Updated[i,c('cd161.gate')], all.gthres.Updated[i,c('cd4.gate')]))
    cd4neg.NKTcells <- getflowFrame(cd4neg.NKTcells.flowD)

    all.events[22] <- cd4neg.NKTcells.flowD@cell.count
    all.props[22] <- cd4neg.NKTcells.flowD@proportion

    # CD4+ NKT-cells
    cd4pos.NKTcells.flowD <- flowDensity(cd5, channels = c(15,16), position = c(T,T), gates = c(all.gthres.Updated[i,c('cd161.gate')], all.gthres.Updated[i,c('cd4.gate')]))
    cd4pos.NKTcells <- getflowFrame(cd4pos.NKTcells.flowD)

    all.events[23] <- cd4pos.NKTcells.flowD@cell.count
    all.props[23] <- cd4pos.NKTcells.flowD@proportion
    
    ###########################################################################################

    # Gating CD4- NKT-cells. Plotting CD62L_CD44 to obtain CD4- NKT Effector and CD4- NKT Resting/Naive
    cd4neg.NKT.Effector.flowD <- flowDensity(cd4neg.NKTcells, channels = c(8,7), position = c(F,T), gates = c(all.gthres.Updated[i,c('cd62.gate')], all.gthres.Updated[i,c('cd44.gate')]))
    #cd4neg.NKT.Effector <- getflowFrame(cd4neg.NKT.Effector.flowD)

    all.events[24] <- cd4neg.NKT.Effector.flowD@cell.count
    all.props[24] <- cd4neg.NKT.Effector.flowD@proportion

    cd4neg.NKT.Resting.flowD <- flowDensity(cd4neg.NKTcells, channels = c(8,7), position = c(T,NA), gates = c(all.gthres.Updated[i,c('cd62.gate')], NA))
    #cd4neg.NKT.Resting <- getflowFrame(cd4neg.NKT.Resting.flowD)

    all.events[25] <- cd4neg.NKT.Resting.flowD@cell.count
    all.props[25] <- cd4neg.NKT.Resting.flowD@proportion
    
    ###########################################################################################

    # Gating CD4- NKT-cells. Plotting KLRG1_CD44 to obtain CD4- NKT KLRG1+
    cd4neg.NKT.klrg1.flowD <- flowDensity(cd4neg.NKTcells, channels = c(12,7), position = c(T,T), gates = c(all.gthres.Updated[i,c('klrg1.gate')], all.gthres.Updated[i,c('cd44.gate')]))
    #cd4neg.NKT.klrg1 <- getflowFrame(cd4neg.NKT.klrg1.flowD)

    all.events[26] <- cd4neg.NKT.klrg1.flowD@cell.count
    all.props[26] <- cd4neg.NKT.klrg1.flowD@proportion
    
    #############################################################################################

    # Gating CD4+ NKT-cells. Plotting CD62L_CD44 to obtain CD4+ NKT Effector and CD4+ NKT Resting/Naive
    cd4pos.NKT.Effector.flowD <- flowDensity(cd4pos.NKTcells, channels = c(8,7), position = c(F,T), gates = c(all.gthres.Updated[i,c('cd62.gate')], all.gthres.Updated[i,c('cd44.gate')]))
    #cd4pos.NKT.Effector <- getflowFrame(cd4pos.NKT.Effector.flowD)

    all.events[27] <- cd4pos.NKT.Effector.flowD@cell.count
    all.props[27] <- cd4pos.NKT.Effector.flowD@proportion

    cd4pos.NKT.Resting.flowD <- flowDensity(cd4pos.NKTcells, channels = c(8,7), position = c(T,NA), gates = c(all.gthres.Updated[i,c('cd62.gate')], NA))
    #cd4pos.NKT.Resting <- getflowFrame(cd4pos.NKT.Resting.flowD)

    all.events[28] <- cd4pos.NKT.Resting.flowD@cell.count
    all.props[28] <- cd4pos.NKT.Resting.flowD@proportion
    
    ###########################################################################################

    # Gating CD4+ NKT-cells. Plotting KLRG1_CD44 to obtain CD4+ NKT KLRG1+
    cd4pos.NKT.klrg1.flowD <- flowDensity(cd4pos.NKTcells, channels = c(12,7), position = c(T,T), gates = c(all.gthres.Updated[i,c('klrg1.gate.cd4posNKT')], all.gthres.Updated[i,c('cd44.gate')]))
    #cd4pos.NKT.klrg1 <- getflowFrame(cd4pos.NKT.klrg1.flowD)

    all.events[29] <- cd4pos.NKT.klrg1.flowD@cell.count
    all.props[29] <- cd4pos.NKT.klrg1.flowD@proportion

    ###########################################################################################

    # Gating P2a. Plotting Cd8a_CD4 to obtain CD5+ CD4/CD8 and abT-cells.
    # CD5+ CD4/CD8
    
    cd5.cd4cd8.flowD <- flowDensity(P2a, channels = c(10,16), position = c(F,F), gates = c(all.gthres.Updated[i,c('cd8.gate')], all.gthres.Updated[i,c('cd4.gate')]))
    #cd5.cd4cd8 <- getflowFrame(cd5.cd4cd8.flowD)

    all.events[30] <- cd5.cd4cd8.flowD@cell.count
    all.props[30] <- cd5.cd4cd8.flowD@proportion

    # abT-cells
    abTcells <- P2a
    abTcells@exprs <- abTcells@exprs[-cd5.cd4cd8.flowD@index,]

    all.events[31] <- nrow(abTcells)
    all.props[31] <- (nrow(abTcells)/nrow(P2a))*100
    
    ###########################################################################################
    # Gating abT-cells. Plotting CD8a_CD4 to obtain CD8a+ T-cells and CD4+ T-cells
    # CD8a+ T-cells
    cd8a.Tcells.flowD <- flowDensity(abTcells, channels = c(10,16), position = c(T,F), gates = c(all.gthres.Updated[i,c('cd8.gate')], all.gthres.Updated[i,c('cd4.gate')]))
    cd8a.Tcells <- getflowFrame(cd8a.Tcells.flowD)

    all.events[32] <- cd8a.Tcells.flowD@cell.count
    all.props[32] <- cd8a.Tcells.flowD@proportion

    # CD4+ T-cells
    cd4.Tcells.flowD <- flowDensity(abTcells, channels = c(10,16), position = c(F,T), gates = c(all.gthres.Updated[i,c('cd8.gate')], all.gthres.Updated[i,c('cd4.gate')]))
    cd4.Tcells <- getflowFrame(cd4.Tcells.flowD)

    all.events[33] <- cd4.Tcells.flowD@cell.count
    all.props[33] <- cd4.Tcells.flowD@proportion
    
    ##############################################################################################

    # Gating CD8a+ T-cells. Plotting CD62L_CD44 to obtain CD8 Naive, CD8 Effector, and CD8 Resting
    ## NOTE: cd44 gate seems have to a different value for gating the CD8 Resting/Naive population. Therefore, I calculated cd44.gate.high just for gating the CD8a+ T-cells.
  
    # CD8 Naive
    cd8.Naive.flowD <- flowDensity(cd8a.Tcells, channels = c(8,7), position = c(T,F), gates = c(all.gthres.Updated[i,c('cd62.gate')], all.gthres.Updated[i,c('cd44.gate.high')]))
    #cd8.Naive <- getflowFrame(cd8.Naive.flowD)

    all.events[34] <- cd8.Naive.flowD@cell.count
    all.props[34] <- cd8.Naive.flowD@proportion

    # CD8 Effector
    cd8.Effector.flowD <- flowDensity(cd8a.Tcells, channels = c(8,7), position = c(F,T), gates = c(all.gthres.Updated[i,c('cd62.gate')], all.gthres.Updated[i,c('cd44.gate')]))
    #cd8.Effector <- getflowFrame(cd8.Effector.flowD)

    all.events[35] <- cd8.Effector.flowD@cell.count
    all.props[35] <- cd8.Effector.flowD@proportion

    # CD8 Resting
    cd8.Resting.flowD <- flowDensity(cd8a.Tcells, channels = c(8,7), position = c(T,T), gates = c(all.gthres.Updated[i,c('cd62.gate')], all.gthres.Updated[i,c('cd44.gate.high')]))
    #cd8.Resting <- getflowFrame(cd8.Resting.flowD)

    all.events[36] <- cd8.Resting.flowD@cell.count
    all.props[36] <- cd8.Resting.flowD@proportion
    
    ################################################################################################

    # Gating CD8a+ T-cells. Plotting KLRG1_CD44 to obtain CD8 KLRG1+
    cd8.KLRG1.flowD <- flowDensity(cd8a.Tcells, channels = c(12,7), position = c(T,T), gates = c(all.gthres.Updated[i,c('klrg1.gate')], all.gthres.Updated[i,c('cd44.gate')]))
    #cd8.KLRG1 <- getflowFrame(cd8.KLRG1.flowD)

    all.events[37] <- cd8.KLRG1.flowD@cell.count
    all.props[37] <- cd8.KLRG1.flowD@proportion
    
    ##################################################################################################

    # Gating CD4+ T-cells. Plotting CD25_GITR to obtain T-helper cells and Tregs

    # T-helper cells
    T.helper.flowD.temp <- flowDensity(cd4.Tcells, channels = c(9,17), position = c(F,NA), gates = c(all.gthres.Updated[i,c('cd25.gate.high')], all.gthres.Updated[i,c('gitr.gate')]))
    T.helper.flowD <- flowDensity(T.helper.flowD.temp, channels = c(9,17), position = c(T,NA), gates = c(all.gthres.Updated[i,c('cd25.gate.low')], all.gthres.Updated[i,c('gitr.gate')]))
    T.helper.flowD@proportion <- (T.helper.flowD@cell.count/cd4.Tcells.flowD@cell.count)*100
    T.helper <- getflowFrame(T.helper.flowD)

    all.events[38] <- T.helper.flowD@cell.count
    all.props[38] <- T.helper.flowD@proportion

    # Tregs
    Tregs.flowD <- flowDensity(cd4.Tcells, channels = c(9,17), position = c(T,T), gates = c(all.gthres.Updated[i,c('cd25.gate.high')], all.gthres.Updated[i,c('gitr.gate')]))
    Tregs <- getflowFrame(Tregs.flowD)

    all.events[39] <- Tregs.flowD@cell.count
    all.props[39] <- Tregs.flowD@proportion
    
    ####################################################################################################

    # Gating T-helper cells. Plotting CD62L_CD44 to obtain CD4 Effector and CD4 Resting/Naive
    # CD4 Resting/Naive
    cd4.Resting.flowD <- flowDensity(T.helper, channels = c(8,7), position = c(T,NA), gates = c(all.gthres.Updated[i,c('cd62.gate')],NA))
    #cd4.Resting <- getflowFrame(cd4.Resting.flowD)
    all.events[40] <- cd4.Resting.flowD@cell.count
    all.props[40] <- cd4.Resting.flowD@proportion

    # CD4 Effector
    cd4.Effector.flowD <- flowDensity(T.helper, channels = c(8,7), position = c(F,T), gates = c(all.gthres.Updated[i,c('cd62.gate')], all.gthres.Updated[i,c('cd44.gate')]))
    #cd4.Effector <- getflowFrame(cd4.Effector.flowD)
    all.events[41] <- cd4.Effector.flowD@cell.count
    all.props[41] <- cd4.Effector.flowD@proportion

    ####################################################################################################

    # Gating T-helper cells. Plotting KLRG1_CD44 to obtain CD4+ KLRG1+ T-cells
    cd4.KLRG1.flowD <- flowDensity(T.helper, channels = c(12,7), position = c(T,T), gates = c(all.gthres.Updated[i,c('klrg1.gate')], all.gthres.Updated[i,c('cd44.gate.high')]))
    #cd4.KLRG1 <- getflowFrame(cd4.KLRG1.flowD)
    all.events[42] <- cd4.KLRG1.flowD@cell.count
    all.props[42] <- cd4.KLRG1.flowD@proportion

    ####################################################################################################

    # Gating Tregs. Plotting CD62L_CD44 to obtain Tregs Resting/Naive and Tregs Effector
    # Tregs Resting/Naive
    Tregs.Resting.flowD <- flowDensity(Tregs, channels = c(8,7), position = c(T,NA), gates = c(all.gthres.Updated[i,c('cd62.gate')],NA))
    #Tregs.Resting <- getflowFrame(Tregs.Resting.flowD)
    all.events[43] <- Tregs.Resting.flowD@cell.count
    all.props[43] <- Tregs.Resting.flowD@proportion
    

    # Tregs Effector
    Tregs.Effector.flowD <- flowDensity(Tregs, channels = c(8,7), position = c(F,T), gates = c(all.gthres.Updated[i,c('cd62.gate')], all.gthres.Updated[i,c('cd44.gate')]))
    #Tregs.Effector <- getflowFrame(Tregs.Effector.flowD)
    all.events[44] <- Tregs.Effector.flowD@cell.count
    all.props[44] <- Tregs.Effector.flowD@proportion
    
    ####################################################################################################

    # Gating Tregs. Plotting KLRG1_CD44 to obtain Treg KLRG1+
    Tregs.KLRG1.flowD <- flowDensity(Tregs, channels = c(12,7), position = c(T,T), gates = c(all.gthres.Updated[i,c('klrg1.gate')], all.gthres.Updated[i,c('cd44.gate.high')]))
    #Tregs.KLRG1 <- getflowFrame(Tregs.KLRG1.flowD)
    all.events[45] <- Tregs.KLRG1.flowD@cell.count
    all.props[45] <- Tregs.KLRG1.flowD@proportion

    ##############################################################################################
    
    ## Adding the new subset of populations CD161+ gd T cells (November 29, 2016)
    ## Gating gd T-cells to obtain the CD161+ CD5+ gd T cells and CD161+ KLRG1+ gd T cells populations.
    
    gd.cd161cd5pos.flowD <- flowDensity(gd.Tcells, channels = c(15,13), position = c(T,NA), gates = c(all.gthres.Updated[i,c('cd161.gate.low')], NA))
    all.events[46] <- gd.cd161cd5pos.flowD@cell.count # CD161+ gd T cells - Events
    all.props[46] <- gd.cd161cd5pos.flowD@proportion # CD161+ gd T cells - as % of gd T cells
    all.props[47] <- (gd.cd161cd5pos.flowD@cell.count/nrow(cd45modified))*100 # cd161+ gd T-cells as % of CD45 population (without the Autofluoresence)
    
    
    gd.cd161cd5.temp.flowD <- flowDensity(gd.Tcells, channels = c(15,13), position = c(NA,T), gates = c(NA, all.gthres.Updated[i,c('cd5.gate')]))
    gd.cd161cd5.flowD <- flowDensity(gd.Tcells, channels = c(15,13), position = c(T,T), gates = c(all.gthres.Updated[i,c('cd161.gate.low')], all.gthres.Updated[i,c('cd5.gate')]))
    #gd.cd5 <- getflowFrame(gd.cd5.flowD)
    all.events[47] <- gd.cd161cd5.flowD@cell.count # CD161+ CD5+ gd T cells - Events
    all.props[48] <- gd.cd161cd5.flowD@proportion # CD161+ CD5+ gd T cells - as % of gd T cells
    all.props[49] <- (gd.cd161cd5.flowD@cell.count/gd.cd161cd5.temp.flowD@cell.count)*100 # CD161+ CD5+ gd T cells - as % of CD5+ gd T cells
    all.props[50] <- (gd.cd161cd5.flowD@cell.count/nrow(cd45modified))*100 # CD161+ CD5+ gd T cells - as % of CD45 population (without the Autofluoresence)
    
    # plotDens(gd.Tcells, c(15,13), main = "gd T-cells (CD161+ gd T cells)", devn = T, cex.lab = 1, cex.axis = 1, cex.main=1); abline(v=cd161.gate.low, lwd=2, col = "blue"); abline(h=cd5.gate, lwd=2, col = "red")
    # lines(gd.cd161cd5.flowD@filter, lwd=2)
    
    # gd.cd161klrg1pos.flowD <- flowDensity(gd.Tcells, channels = c(15,12), position = c(T,NA), gates = c(cd161.gate.low, NA))
    # all.events[i,46] <- gd.cd161cd5pos.flowD@cell.count
    # all.props[i,46] <- gd.cd161cd5pos.flowD@proportion
    # all.props[i,47] <- (gd.cd161cd5pos.flowD@cell.count/nrow(cd45modified))*100 # cd161+ gd T-cells as % of CD45 population (without the Autofluoresence)
    
    gd.cd161klrg1.temp.flowD <- flowDensity(gd.Tcells, channels = c(15,12), position = c(NA,T), gates = c(NA, all.gthres.Updated[i,c('klrg1.gate.new')]))
    #gd.cd161klrg1.flowD <- flowDensity(gd.Tcells, channels = c(15,12), position = c(T,T), gates = c(all.gthres.Updated[i,c('cd161.gate.low')], all.gthres.Updated[i,c('klrg1.gate')]))
    gd.cd161klrg1.flowD <- flowDensity(gd.Tcells, channels = c(15,12), position = c(T,T), gates = c(all.gthres.Updated[i,c('cd161.gate.low')], all.gthres.Updated[i,c('klrg1.gate.new')]))
    all.events[48] <- gd.cd161klrg1.flowD@cell.count # CD161+ KLRG1+ gd T cells - Events
    all.props[51] <- gd.cd161klrg1.flowD@proportion # CD161+ KLRG1+ gd T cells - as % of gd T cells
    all.props[52] <- (gd.cd161klrg1.flowD@cell.count/gd.cd161klrg1.temp.flowD@cell.count)*100 # CD161+ KLRG1+ gd T cells - as % of KLRG1+ gd T cells
    all.props[53] <- (gd.cd161klrg1.flowD@cell.count/nrow(cd45modified))*100 # CD161+ KLRG1+ gd T cells - as % of CD45 population (without the Autofluoresence)
    

    #################################################################################################
    #################################################################################################
    
    ## Saving the Plots
    #--------Start Big Png------------
     #png ( file = paste("Results/Figures/ScatterPlots/", x$Genotype, "/", "Total_", x$FCS.file, ".png", sep = "" ), width=2100, height=2100*7/6)
    png ( file = paste0(results.dir,"/Figures/ScatterPlotsUpdated/", x$Genotype, "/", "Total_", x$FCS.files, ".png"), width=2100, height=2100*7/6)
     
    par(mfrow=c(7,6),mar=(c(5, 5, 4, 2) + 0.1))
    
    # 1st method of plotting FSC-A_SSC-W
    plotDens(f, c("FSC-A","SSC-W"), main="Ungated", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(singlets.flowD.l@filter,lwd=2)
    
    # Plotting Live/Dead_SSC-A 
    plotDens(singlets,  c(11,4), main= "Singlets", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=all.gthres.Updated[i,c('live.gate')], lwd=2); 
    
    # Plotting FSC-A_SSC-A
    plotDens(live, c(1,4), main = "Live", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=all.gthres.Updated[i,c('fsc.a.gate.low')], lwd=2); abline(v=all.gthres.Updated[i,c('fsc.a.gate.high')], lwd=2); abline(h=all.gthres.Updated[i,c('ssc.a.gate')], lwd=2)
    
    # Plotting CD45_CD161
    plotDens(lymph, channels = c(14,15), main = "Lymphocytes", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=all.gthres.Updated[i,c('cd45.gate.low')], lwd=2); abline(v=all.gthres.Updated[i,c('cd45.gate.high')], lwd=2); abline(h=all.gthres.Updated[i,c('cd161.gate.high')], lwd=2)
    
    # Plotting TCRd_CD4. Highlighting the Autofluoresence part
    plotDens(cd45, channels = c(18,16), main="CD45+_Autofluorecence", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    if ( nrow(cd45autoflour@exprs) != 0) {
      points(TCRDCD4points, type="l", col="black", lty=2, lwd=2)    
    }
    
    # Plotting TCRd_CD4. Getting rid of the Autofluoresence part
    plotDens(cd45modified,channels = c(18,16), main="CD45+_Without Autofluorecence", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    
    
    ## Plotting TCRd_CD4 to obtain gd T-cells
    plotDens(cd45modified, channels = c(18,16), main = "CD45+_gd T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2);
    lines(gd.Tcells.flowD@filter, lwd=2)
    
    ## Plotting CD62L_CD44 to obtain gd Resting, gd Effector, and gd Naive 
    min.x <- min(exprs(gd.Tcells)[,8])-1
    max.x <- max(exprs(gd.Tcells)[,8])+1
    min.y <- min(exprs(gd.Tcells)[,7])-1
    max.y <- max(exprs(gd.Tcells)[,7])+1
    plotDens(gd.Tcells, c(8,7), main = "gd T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x,max.x), ylim = c(min.y, max.y)); abline(v=all.gthres.Updated[i,c('cd62.gate')], lwd=2); abline(h=all.gthres.Updated[i,c('cd44.gate')], lwd=2)
    
    ## Plotting KLRG1_GITR to obtain GITR gd T-cells and GD KLRG1+
    min.x <- min(exprs(gd.Tcells)[,12])
    max.x <- max(exprs(gd.Tcells)[,12])+1
    plotDens(gd.Tcells, c(12,17), main = "gd T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim =c(min.x, max.x)); abline(v=all.gthres.Updated[i,c('klrg1.gate')], lwd=2); abline(h=all.gthres.Updated[i,c('gitr.gate')], lwd=2)
    
    ## Plotting CD5_CD44 to obtain GD CD5+
    plotDens(gd.Tcells, c(13,7), main = "gd T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=all.gthres.Updated[i,c('cd5.gate')], lwd=2)
    
    ## New population - added in Nov 29, 2016
    ## Plotting CD161_CD5 to obtain CD161+ gd T cell population and CD161+CD5+ gd T cell
    plotDens(gd.Tcells, c(15,13), main = "CD161+ gd T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=all.gthres.Updated[i,c('cd161.gate.low')], lwd=2, col = "blue"); abline(h=all.gthres.Updated[i,c('cd5.gate')], lwd=2, col = "red")
    lines(gd.cd161cd5.flowD@filter, lwd=2)

    ## New population - added in Nov 29, 2016
    ## Plotting CD161_KLRG1 to obtain CD161+ gd T cell population and CD161+KLRG1+ gd T cell
    #plotDens(gd.Tcells, c(15,12), main = "CD161+ gd T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=all.gthres.Updated[i,c('cd161.gate.low')], lwd=2, col = "blue"); abline(h=all.gthres.Updated[i,c('klrg1.gate')], lwd=2, col = "red")
    plotDens(gd.Tcells, c(15,12), main = "CD161+ gd T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=all.gthres.Updated[i,c('cd161.gate.low')], lwd=2, col = "blue"); abline(h = all.gthres.Updated[i,c('klrg1.gate.new')], lwd=2, col = "red")
    # plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
    # plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
    
    ## Plotting TCRd_CD4 to obtain NOT(P2) (without the gd T-cells part)
    plotDens(not.gd.Tcells, channels = c(18,16), main = "NOT(P2)", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2);
    #lines(gd.Tcells.flowD@filter, lwd=2)

    ## Plotting CD161_CD5 to obtain CD5+ and NK-cells
    plotDens(not.gd.Tcells, c(15,13), main = "NOT(gd T-cells)", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(h=all.gthres.Updated[i,c('cd5.gate.new')], lwd=2)
    #lines(x=c(cd161.gate.low, cd5.gate), y=c(cd161.gate.high, cd5.gate), lwd=2)
    lines(NKcells.flowD@filter, lwd=2)

    #Plotting CD62L_CD44 to obtain NK Effector and NK Resting/Naive
    min.x <- min(exprs(NKcells)[,8])-1
    min.y <- min(exprs(NKcells)[,7])-1
    max.y <-  max(exprs(NKcells)[,7])
    plotDens(NKcells, c(8,7), main = "NK-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, ylim =c(min.y, max.y)); abline(v=all.gthres.Updated[i,c('cd62.gate')], lwd=2)
    lines(x=c(min.x, all.gthres.Updated[i,c('cd62.gaklrg1.gate.cd4posNKTte')]), y=c(all.gthres.Updated[i,c('cd44.gate')],all.gthres.Updated[i,c('cd44.gate')]), lwd=2)

    # Plotting KLRG1_CD44 to obtain NK +KLRG1
    min.x <- min(exprs(NKcells)[,12])
    max.x <- max(exprs(NKcells)[,12])+1
    min.y <- min(exprs(NKcells)[,7])-1
    max.y <-  max(exprs(NKcells)[,7])+0.5
    plotDens(NKcells, c(12,7), main = "NK-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x,max.x), ylim = c(min.y,max.y)); #abline(h=cd44.gate, lwd=2)
    lines(x=c(all.gthres.Updated[i,c('klrg1.gate')], all.gthres.Updated[i,c('klrg1.gate')]), y=c(all.gthres.Updated[i,c('cd44.gate')], max.y+2), lwd=2)
    lines(x=c(all.gthres.Updated[i,c('klrg1.gate')], max.x+1), y=c(all.gthres.Updated[i,c('cd44.gate')], all.gthres.Updated[i,c('cd44.gate')]), lwd=2)

    plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
    plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
    plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank

    # Plotting CD161_CD4 to obtain P2a, CD4- NKT-cells, and CD4+ NKT-cells
    min.x <- min(exprs(cd5)[,15])
    max.x <- max(exprs(cd5)[,15])+1
    min.y <- min(exprs(cd5)[,16])-0.5
    max.y <-  max(exprs(cd5)[,16])+0.5
    plotDens(cd5, c(15,16), main = "CD5+", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x,max.x), ylim = c(min.y,max.y)); abline(v=all.gthres.Updated[i,c('cd161.gate')], lwd=2)
    lines(x=c(all.gthres.Updated[i,c('cd161.gate')], max.x+1), y=c(all.gthres.Updated[i,c('cd4.gate')],all.gthres.Updated[i,c('cd4.gate')]), lwd=2)

    # Plotting CD62L_CD44 to obtain CD4- NKT Effector and CD4- NKT Resting/Naive
    min.x <- min(exprs(cd4neg.NKTcells)[,8])-1
    max.x <- max(exprs(cd4neg.NKTcells)[,8])+1
    min.y <- min(exprs(cd4neg.NKTcells)[,7])-2
    max.y <-  max(exprs(cd4neg.NKTcells)[,7])
    plotDens(cd4neg.NKTcells, c(8,7), main = "CD4- NKT-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x,max.x), ylim = c(min.y,max.y)); abline(v=all.gthres.Updated[i,c('cd62.gate')], lwd=2)
    lines(x=c(min.x-1,all.gthres.Updated[i,c('cd62.gate')]), y=c(all.gthres.Updated[i,c('cd44.gate')],all.gthres.Updated[i,c('cd44.gate')]), lwd=2)

    # Plotting KLRG1_CD44 to obtain CD4- NKT KLRG1+
    min.x <- min(exprs(cd4neg.NKTcells)[,12])
    max.x <- max(exprs(cd4neg.NKTcells)[,12])+1
    min.y <- min(exprs(cd4neg.NKTcells)[,7])-2
    max.y <-  max(exprs(cd4neg.NKTcells)[,7])
    plotDens(cd4neg.NKTcells, c(12,7), main = "CD4- NKT-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x,max.x), ylim = c(min.y,max.y))
    lines(x=c(all.gthres.Updated[i,c('klrg1.gate')],max.x+1), y=c(all.gthres.Updated[i,c('cd44.gate')],all.gthres.Updated[i,c('cd44.gate')]),lwd=2)
    lines(x=c(all.gthres.Updated[i,c('klrg1.gate')], all.gthres.Updated[i,c('klrg1.gate')]), y=c(all.gthres.Updated[i,c('cd44.gate')], max.y+1),lwd=2)

    # Plotting CD62L_CD44 to obtain CD4+ NKT Effector and CD4+ NKT Resting/Naive
    min.x <- min(exprs(cd4pos.NKTcells)[,8])-1
    max.x <- max(exprs(cd4pos.NKTcells)[,8])+1
    min.y <- min(exprs(cd4pos.NKTcells)[,7])-2
    max.y <-  max(exprs(cd4pos.NKTcells)[,7])
    plotDens(cd4pos.NKTcells, c(8,7), main = "CD4+ NKT-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x,max.x), ylim = c(min.y,max.y)); abline(v=all.gthres.Updated[i,c('cd62.gate')], lwd=2)
    lines(x=c(min.x-1,all.gthres.Updated[i,c('cd62.gate')]), y=c(all.gthres.Updated[i,c('cd44.gate')],all.gthres.Updated[i,c('cd44.gate')]), lwd=2)

    # Plotting KLRG1_CD44 to obtain CD4+ NKT KLRG1+
    min.x <- min(exprs(cd4pos.NKTcells)[,12])
    max.x <- max(exprs(cd4pos.NKTcells)[,12])+1
    min.y <- min(exprs(cd4pos.NKTcells)[,7])-2
    max.y <-  max(exprs(cd4pos.NKTcells)[,7])
    plotDens(cd4pos.NKTcells, c(12,7), main = "CD4+ NKT-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(min.x,max.x), ylim = c(min.y,max.y))
    lines(x=c(all.gthres.Updated[i,c('klrg1.gate.cd4posNKT')],max.x+1), y=c(all.gthres.Updated[i,c('cd44.gate')],all.gthres.Updated[i,c('cd44.gate')]),lwd=2)
    lines(x=c(all.gthres.Updated[i,c('klrg1.gate.cd4posNKT')], all.gthres.Updated[i,c('klrg1.gate.cd4posNKT')]), y=c(all.gthres.Updated[i,c('cd44.gate')], max.y+1),lwd=2)

    plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank

    # Plotting CD8a_CD4 to obtain CD5+ CD4/CD8 and abT-cell.
    min.x <- min(exprs(P2a)[,10])-1
    max.x <- max(exprs(P2a)[,10])
    min.y <- min(exprs(P2a)[,16])-1
    max.y <-  max(exprs(P2a)[,16])
    plotDens(P2a, c(10,16), main = "P2a", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(x=c(all.gthres.Updated[i,c('cd8.gate')],all.gthres.Updated[i,c('cd8.gate')]), y=c(min.y,all.gthres.Updated[i,c('cd4.gate')]), lwd=2)
    lines(x=c(min.x,all.gthres.Updated[i,c('cd8.gate')]), y=c(all.gthres.Updated[i,c('cd4.gate')], all.gthres.Updated[i,c('cd4.gate')]), lwd=2)

    # Plotting CD8a_CD4 to obtain CD8a+ T-cells and CD4+ T-cells
    plotDens(abTcells, c(10,16), main = "abT-CELL", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(h=all.gthres.Updated[i,c('cd4.gate')], lwd=2); abline(v=all.gthres.Updated[i,c('cd8.gate')],lwd=2)

    # Plotting CD62L_CD44 to obtain CD8 Naive, CD8 Effector, and CD8 Resting
    min.x <- min(exprs(cd8a.Tcells)[,8])-1
    max.x <- max(exprs(cd8a.Tcells)[,8])+1
    min.y <- min(exprs(cd8a.Tcells)[,7])
    max.y <-  max(exprs(cd8a.Tcells)[,7])
    plotDens(cd8a.Tcells, c(8,7), main = "CD8a+ T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2); abline(v=all.gthres.Updated[i,c('cd62.gate')], lwd=2)
    lines(x=c(min.x,all.gthres.Updated[i,c('cd62.gate')]), y=c(all.gthres.Updated[i,c('cd44.gate')],all.gthres.Updated[i,c('cd44.gate')]), lwd=2)
    lines(x=c(all.gthres.Updated[i,c('cd62.gate')],max.x), y=c(all.gthres.Updated[i,c('cd44.gate.high')],all.gthres.Updated[i,c('cd44.gate.high')]), lwd=2)

    # Plotting KLRG1_CD44 to obtain CD8 KLRG1+
    min.x <- min(exprs(cd8a.Tcells)[,12])
    max.x <- max(exprs(cd8a.Tcells)[,12])+1
    min.y <- min(exprs(cd8a.Tcells)[,7])-1
    max.y <-  max(exprs(cd8a.Tcells)[,7])
    plotDens(cd8a.Tcells, c(12,7), main = "CD8a+ T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim=c(min.x,max.x), ylim =c(min.y,max.y)); abline(v=all.gthres.Updated[i,c('klrg1.gate')], lwd=2)
    lines(x=c(all.gthres.Updated[i,c('klrg1.gate')],max.x+1), y=c(all.gthres.Updated[i,c('cd44.gate')],all.gthres.Updated[i,c('cd44.gate')]), lwd=2)

    plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
    plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
    plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank

    # Plotting CD25_GITR to obtain T-helper cells and Tregs
    min.x <- min(exprs(cd4.Tcells)[,9])
    max.x <- max(exprs(cd4.Tcells)[,9])+1
    min.y <- min(exprs(cd4.Tcells)[,17])
    max.y <-  max(exprs(cd4.Tcells)[,17])+1
    plotDens(cd4.Tcells, c(9,17), main = "CD4+ T-cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim=c(min.x,max.x), ylim =c(min.y,max.y)); abline(v=all.gthres.Updated[i,c('cd25.gate.low')], lwd=2); abline(v=all.gthres.Updated[i,c('cd25.gate.high')], lwd=2)
    lines(x=c(all.gthres.Updated[i,c('cd25.gate.high')],max.x+1), y=c(all.gthres.Updated[i,c('gitr.gate')],all.gthres.Updated[i,c('gitr.gate')]), lwd=2)

    # Plotting CD62L_CD44 to obtain CD4 Effector and CD4 Resting/Naive
    min.x <- min(exprs(T.helper)[,8])
    max.x <- max(exprs(T.helper)[,8])+0.5
    min.y <- min(exprs(T.helper)[,7])
    max.y <-  max(exprs(T.helper)[,7])
    plotDens(T.helper, c(8,7), main = "T-helper cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim=c(min.x,max.x), ylim =c(min.y,max.y)); abline(v=all.gthres.Updated[i,c('cd62.gate')], lwd=2)
    lines(x=c(min.x-1, all.gthres.Updated[i,c('cd62.gate')]), y=c(all.gthres.Updated[i,c('cd44.gate')],all.gthres.Updated[i,c('cd44.gate')]), lwd=2)

    # Plotting KLRG1_CD44 to obtain CD4+ KLRG1+ T-cells
    min.x <- min(exprs(T.helper)[,12])
    max.x <- max(exprs(T.helper)[,12])
    min.y <- min(exprs(T.helper)[,7])
    max.y <-  max(exprs(T.helper)[,7])
    plotDens(T.helper, c(12,7), main = "T-helper cells", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(x=c(all.gthres.Updated[i,c('klrg1.gate')], max.x+1), y=c(all.gthres.Updated[i,c('cd44.gate.high')], all.gthres.Updated[i,c('cd44.gate.high')]), lwd=2)
    lines(x=c(all.gthres.Updated[i,c('klrg1.gate')], all.gthres.Updated[i,c('klrg1.gate')]), y=c(max.y+1, all.gthres.Updated[i,c('cd44.gate.high')]),lwd=2)

    plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
    plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
    plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
    plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank

    #Plotting CD62L_CD44 to obtain Tregs Resting/Naive and Tregs Effector
    min.x <- min(exprs(Tregs)[,8])-1
    max.x <- max(exprs(Tregs)[,8])+1
    min.y <- min(exprs(Tregs)[,7])-1
    max.y <-  max(exprs(Tregs)[,7])
    plotDens(Tregs, c(8,7), main = "Tregs", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2, xlim=c(min.x,max.x), ylim =c(min.y,max.y)); abline(v= all.gthres.Updated[i,c('cd62.gate')], lwd=2)
    lines(x=c(min.x-1, all.gthres.Updated[i,c('cd62.gate')]), y=c(all.gthres.Updated[i,c('cd44.gate')],all.gthres.Updated[i,c('cd44.gate')]), lwd=2)

    # Plotting KLRG1_CD44 to obtain Tregs KLRG1+
    min.x <- min(exprs(Tregs)[,12])
    max.x <- max(exprs(Tregs)[,12])
    min.y <- min(exprs(Tregs)[,7])
    max.y <-  max(exprs(Tregs)[,7])
    plotDens(Tregs, c(12,7), main = "Tregs", devn = T, cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(x=c(all.gthres.Updated[i,c('klrg1.gate')], max.x+1), y=c(all.gthres.Updated[i,c('cd44.gate.high')], all.gthres.Updated[i,c('cd44.gate.high')]), lwd=2)
    lines(x=c(all.gthres.Updated[i,c('klrg1.gate')], all.gthres.Updated[i,c('klrg1.gate')]), y=c(max.y+1, all.gthres.Updated[i,c('cd44.gate.high')]),lwd=2)
    
    dev.off()
    par(mfrow=c(1,1))
    
    #--------End Big Png------------
    
  },error = function(err) {
    print(paste("ERROR in Gating/Plotting Message:  ",err, "Index: ", i));
    #errorFileIndex <- c(errorFileIndex,i)
    #return(errorFileIndex)
    return(0)
  }
  ) # end of tryCatch
  
  data.frame(x$FCS.files, all.props, all.events)
  # cat("Time is: ",TimeOutput(start2),"\n",sep="")
#} # end of for-loop
  
}, .parallel = TRUE) # end ldply

print("Finished Gating & Plotting")

# Saving the big dataframe of Proportions, Events, and Gating thresholds 
save(props.events, file = paste0(results.dir,"/PropsEvents.Updated.Rdata"))  

# Checking if any files failed the gating by finding if there is any NA
failedGating.files.index <- which(is.na(props.events[,ncol(props.events)]))

if(length(failedGating.files.index) != 0){
  print("Files which failed Gating:")
  print(props.events[failedGating.files.index,'x.FCS.files'])
  failedGating.files <- store.allFCS[failedGating.files.index,]
  rownames(failedGating.files) <- failedGating.files.index
  save(failedGating.files, file = paste0(results.dir,"/failedGating.files.Updated.Rdata"))  
  ## For sending a spreadsheet containing list of FCS files which Failed the Gating to center
  write.csv(failedGating.files, file = paste0(results.dir, "/failedGatingUpdated.TcellMLN.csv"), row.names = FALSE)
}else{
  print("No files failed the Gating")
}

# ## Arranging the output from ddply in order of the file.names or store.allFCS.unique matrix
# props.events.gates <- join(file.names, props.events.gates)

# Extracting the Proportions from the large dataframe that was returned drom ldply
all.props <- props.events[,2:54] 
all.props <- as.matrix(all.props)
colnames(all.props) <- c("All Events-%Parent", "Singlets-%Parent", "Live-%Parent", "Lymphocytes-%Parent", "CD45-%Parent", "Autofluorecence-%Parent", "NOT(P2)-%Parent", "gd T-cells-%Parent", "NOT(gd T-cells)-%Parent",
                               "GD Resting-%Parent","GD Effector-%Parent", "GD-Naive-%Parent", "GD KLRG1+-%Parent", "GITR GD T-cells-%Parent", "GD CD5+-%Parent", "CD5+-%Parent", "NK-cells-%Parent", "NK Resting/Naive-%Parent",
                               "NK Effector-%Parent", "NK KLRG1-%Parent", "P2a-%Parent", "CD4- NKT-cells-%Parent", "CD4+ NKT-cells-%Parent", "CD4- NKT Effector-%Parent", "CD4- NKT Resting-%Parent", "CD4- NKT KLRG1+-%Parent",
                               "CD4+ NKT Effector-%Parent", "CD4+ NKT Resting-%Parent", "CD4+ KLRG1+-%Parent", "CD5+ CD4/CD8-%Parent", "ab T-cell-%Parent", "CD8a+ T-cells-%Parent", "CD4+ T-cells-%Parent",
                               "CD8 Naive-%Parent", "CD8 Effector-%Parent", "CD8 Resting-%Parent", "Cd8 KLRG1-%Parent", "T-helper cells-%Parent", "Tregs-%Parent", "CD4 Resting-%Parent", "CD4 Effector-%Parent",
                               "CD4 KLRG1-%Parent", "Tregs Resting-%Parent", "Tregs Effector-%Parent", "Tregs KLRG1-%Parent", "CD161+ gd T cells-as % of gd T cells", "CD161+ gd T cells-as % of CD45+ cells", 
                               "CD161+ CD5+ gd T cells-as % of gd T cells", "CD161+ CD5+ gd T cells-as % of CD5+ gd T cells", "CD161+ CD5+ gd T cells - as % of CD45+ cells", "CD161+ KLRG1+ gd T cells - as % of gd T cells",
                               "CD161+ KLRG1+ gd T cells-as % of KLRG1+ gd T cells", "CD161+ KLRG1+ gd T cells-as % of CD45+ cells")
all.props <- cbind(store.allFCS[,c('Panel/Organ', 'Genotype', 'FCS files', 'Barcodes', 'Assay Date', 'Gender', 'Number of Channels', 'Number of Cells')], all.props)


# Extracting the #Events from the large dataframe that was returned drom ddply
all.events <- props.events[,55:ncol(props.events)]
all.events <- as.matrix(all.events)
colnames(all.events) <- c("All Events-#Events", "Singlets-#Events", "Live-#Events", "Lymphocytes-#Events", "CD45-#Events", "Autofluorecence-#Events", "NOT(P2)-#Events", "gd T-cells-#Events", "NOT(gd T-cells)-#Events",
                                "GD Resting-#Events","GD Effector-#Events", "GD-Naive-#Events", "GD KLRG1+-#Events", "GITR GD T-cells-#Events", "GD CD5+-#Events", "CD5+-#Events", "NK-cells-#Events", "NK Resting/Naive-#Events",
                                "NK Effector-#Events", "NK KLRG1-#Events", "P2a-#Events", "CD4- NKT-cells-#Events", "CD4+ NKT-cells-#Events", "CD4- NKT Effector-#Events", "CD4- NKT Resting-#Events", "CD4- NKT KLRG1+-#Events",
                                "CD4+ NKT Effector-#Events", "CD4+ NKT Resting-#Events", "CD4+ KLRG1+-#Events", "CD5+ CD4/CD8-#Events", "ab T-cell-#Events", "CD8a+ T-cells-#Events", "CD4+ T-cells-#Events",
                                "CD8 Naive-#Events", "CD8 Effector-#Events", "CD8 Resting-#Events", "Cd8 KLRG1-#Events", "T-helper cells-#Events", "Tregs-#Events", "CD4 Resting-#Events", "CD4 Effector-#Events",
                                "CD4 KLRG1-#Events", "Tregs Resting-#Events", "Tregs Effector-#Events", "Tregs KLRG1-#Events", "CD161+ gd Tcells-#Events", "CD161+ CD5+ gd T cells-#Events", "CD161+ KLRG1+ gd T cells-#Events")
all.events <- cbind(store.allFCS[,c('Panel/Organ', 'Genotype', 'FCS files', 'Barcodes', 'Assay Date', 'Gender', 'Number of Channels', 'Number of Cells')], all.events)


## There are files which failed the Gating in second run after the threshold outliers have been fixed. So we exclude those files
if (nrow(failedGating.files) != 0){
  failedGating.files.index <- which(is.na(props.events[,ncol(props.events)]))
  all.props <- all.props[-failedGating.files.index,]
  all.events <- all.events[-failedGating.files.index,]
}

rownames(all.props) <- 1:nrow(all.props)
rownames(all.events) <- 1:nrow(all.events)

write.csv(all.props, file =  paste0(results.dir, "/allProportions_Table.csv"), row.names = FALSE)
write.csv(all.events, file =  paste0(results.dir, "/allEvents_Table.csv"), row.names = FALSE)


cat("Total time is: ",TimeOutput(start),"\n",sep="")

