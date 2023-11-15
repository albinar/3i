# Originally written by Justin Meskas 
# Modified & Updated by Albina Rahim
# Date: March 21, 2016

remove(list=ls())

#setwd("/code/Projects/3i/Panel_T-cell/")
setwd("/data/Panel_M-cell")

library("flowCore")
library("flowBin")
library("flowDensity")
library("flowType")
# library("flowClean")
# library("snowfall")
# library("foreach")
library("stringr")
library('pracma') 
library('MASS')

source("Codes/3iTcellfunctions-M.R")
source("Codes/rotate.data-M.R")

load("Results_MLN/store.allFCS.Rdata")

start <- Sys.time()
errorFileIndex <- NULL


# Ly6Ghi.SSCAlo.files.barcodes <- c(92105, 109527,92123, 96523,96531, 96537, 55606, 57930, 70281, 84494, 88737, 97854, 88740, 88741, 89595, 93968,
#                                   97849, 93985, 97840, 64391, 97838, 97848, 108033, 108038, 109526, 109527, 90429)
# loop.idx <- unlist(lapply(Ly6Ghi.SSCAlo.files.barcodes, function(x){ grep(x, store.allFCS[,2])}))
write.data <- F # write cell count & proportion data to cvs file
do.flowType <- T

# Some of the files are pre-gated CD45+, Live, single cells. 
pregate.mln.barcodes <- c(53519, 53521:53524, 53527:53531) 
pregate.mln.barcodes <- paste(pregate.mln.barcodes, collapse = "|")
#pregate.idx <- lapply(pregate.mln.barcodes, function(x){ grep(x, store.allFCS[,2])})

# all.props <- matrix(nrow = nrow(store.allFCS), ncol = 20)
# all.events <- matrix(nrow = nrow(store.allFCS), ncol = 20)
# gthres <- matrix(nrow = nrow(store.allFCS), ncol = 22)

library(plyr)
library(doMC)
no_cores <- detectCores() - 1
registerDoMC(no_cores)

file.names <- data.frame(store.allFCS, stringsAsFactors = F)

props.events.gates <- ddply(file.names[1:2,], "FCS.file", function(x){
  
#for(i in loop.idx){
  
  all.props <- matrix(nrow = 1, ncol = 6)
  all.events <- matrix(nrow = 1, ncol = 6)
  gthres <- matrix(nrow = 1, ncol = 10)

  possibleError <- tryCatch({
    #start2 <- Sys.time()
    
    #print(paste("Starting ", x$Genotype, "/", x$FCS.file, sep = ""))
    load(file = paste("Results_MLN/After_Clean/", x$Genotype, "/AftClean_", x$FCS.file, ".Rdata", sep=""))
    
    all.props[1] <- nrow(f)/nrow(f)*100
    all.events[1] <- nrow(f)
    
    ###############################################################################################
    ## Gating All Events. Plotting SSC-A_SSC-W to obtain singlets
    
    # Find location of minima (valley) between the largest SSC-W peak and the next highest peak in the +SSC-W direction
    maxDens <- density(f@exprs[, 6])
    temp <- findpeaks(maxDens$y)
    n = length(temp[, 2])
    first.peak <- sort(temp[, 1], partial = n-1)[n]
    first.peak.idx <- temp[which(temp[, 1] == first.peak),2]
    # look for the next peak only for higher SSC-W values
    temp <- temp[which(temp[,2] > (first.peak.idx + 12)),]    
    n = length(temp[, 2])
    second.peak <- sort(temp[, 1], partial = n-1)[n]
    valley = min(maxDens$y[which(maxDens$y == first.peak):which(maxDens$y == second.peak)])
    singlets.gate.h = maxDens$x[which(maxDens$y == valley)]
    # if the valley and second peak are essentially the same
    if ((which(maxDens$y == second.peak) - which(maxDens$y == valley)) < 5){ 
      temp <- deGate(f, "SSC-W", use.upper = T, upper =T)
      if(abs(singlets.gate.h - 80842) > abs(temp - 80842)){ #only use temp if it is closer to the mean
        singlets.gate.h <- temp
      }
    }
    
    singlets.gate.l <- deGate(f, channel = c("SSC-W"), use.upper = T, upper = F)
    
    names(singlets.gate.h) <- names(singlets.gate.l) <- identifier(f)
    
    singlets.flowD.h <- flowDensity(f, channels = c("SSC-A", "SSC-W"), position = c(NA, F), gates = c(NA, singlets.gate.h[identifier(f)]))
    singlets.flowD.l <- flowDensity(singlets.flowD.h, channels = c("SSC-A", "SSC-W"), position = c(NA, T), gates = c(NA, singlets.gate.l[identifier(f)]))
    singlets.flowD.l@proportion <- (singlets.flowD.l@cell.count/nrow(f))*100
    singlets.flowD.l.ind <- singlets.flowD.l@index
    singlets <- getflowFrame(singlets.flowD.l)
    
    gthres[1] <- singlets.gate.h
    gthres[2] <- singlets.gate.l
    all.props[2] <- singlets.flowD.l@proportion
    all.events[2] <- singlets.flowD.l@cell.count
    
    ######################################################################################################
    ## Gating Singlets. Plotting Live/Dead_SSC-A 
    
    live.gate <- deGate(singlets, channel = c(11))
    if(live.gate > 3.5){ # rel. small live population
      temp <- deGate(singlets, channel = c(11), use.upper = T, upper = F) - 0.2
      temp <- c(temp, deGate(singlets, channel = c(11), all.cut = T))
      live.gate <- temp[which.min(abs(temp - 2.71))]        # find gate closest to expected location
    }
    names(live.gate) <- identifier(singlets)
    
    live.flowD <- flowDensity(singlets, channels = c(11, 4), position = c(F, NA), gates = c(live.gate, NA))
    live.flowD.ind <- live.flowD@index
    live <- getflowFrame(live.flowD)
    
    gthres[3] <- live.gate
    all.props[3] <- live.flowD@proportion
    all.events[3] <- live.flowD@cell.count
    
    ##############################################################################################
    ## Gating Live. Plotting FSC-A_SSC-A
    
    fsc.a.gate.low <- max(deGate(live, channel = c(1), all.cut = T, use.upper = T, upper = F))
    # There are places where max(deGate()) is above the lymphocyte peak (~60,000)
    # Check if this is the case pick second highest gate if it is
    if(fsc.a.gate.low > 60000){
      temp <- deGate(live, channel = c(1), all.cut = T, use.upper = T, upper = F)
      if(length(temp) > 1){
        fsc.a.gate.low <-temp[length(temp) - 1]
      }
    }
    fsc.a.gate.high <- deGate(live, channel = c(1), use.percentile = T, percentile = 0.99999999)
    ssc.a.gate <- deGate(live, channel = c(4), use.percentile = T, percentile = 0.9999999)
    names(fsc.a.gate.low) <- names(fsc.a.gate.high) <- names(ssc.a.gate) <- identifier(live)
    
    temp <- flowDensity(live, channels = c(1, 4), position = c(T, F), gates = c(fsc.a.gate.low, ssc.a.gate))
    lymph.flowD <- flowDensity(temp, channels = c(1, 4), position = c(F, F), gates = c(fsc.a.gate.high, ssc.a.gate))
    lymph.flowD.ind <- lymph.flowD@index
    lymph.flowD@proportion <- (lymph.flowD@cell.count/nrow(live))*100
    lymph <- getflowFrame(lymph.flowD)
    
    gthres[4] <- fsc.a.gate.high
    gthres[5] <- fsc.a.gate.low
    gthres[6] <- ssc.a.gate
    all.props[4] <- lymph.flowD@proportion
    all.events[4] <- lymph.flowD@cell.count
    
    #############################################################################################
    ## Gating Lymphocytes. Plotting CD45_Lin(CD3, cd19, CD161) to get CD45 cells
    
    # CD45
    cd45.gate.low <- c(deGate(lymph, channel = c(14), all.cut = T, use.upper = T, upper =F, alpha = 0.05),
                       deGate(lymph, channel = c(14), use.upper = T, upper =F, alpha = 0.05)) # in case there is a local min in the CD45 peak
    cd45.gate.low <- cd45.gate.low[which.min(abs(cd45.gate.low - 2))]
    if(cd45.gate.low < 1){ 
      # A handful of the MLN files contain only a cd45 lo distn
      # it is not always caught by the outlier check, so I am setting it equal to the mean manually
      cd45.gate.low <- 2.0 
    }
    # cd45.gate.high <- deGate(lymph, channel = c(14), use.percentile = T, percentile = 1)
    cd45.gate.high <- deGate(lymph, channel = c(14), use.percentile = T, percentile = 0.999) + 0.2
    if((cd45.gate.high - deGate(lymph, channel = c(14), use.upper = T, upper = T)) > 1){ # suggests outliers/edge effects
      cd45.gate.high <-  deGate(lymph, channel = c(14), use.upper = T, upper = T) + 0.4
    }
    # gating below suggested by Albina to mimic her T cell cd161 gating
    lincd3cd19cd161.gate.high <- deGate(lymph, channel = c(12), use.percentile = T, percentile = 0.9999) + 0.1 
    
    names(cd45.gate.low) <- names(cd45.gate.high) <-  names(lincd3cd19cd161.gate.high) <- identifier(lymph)
    
    cd45.flowD.temp <- flowDensity(lymph, channels = c(14, 12), position = c(T, F), gates = c(cd45.gate.low, lincd3cd19cd161.gate.high))
    cd45.flowD <- flowDensity(cd45.flowD.temp, channels = c(14, 12), position = c(F, F), gates = c(cd45.gate.high, lincd3cd19cd161.gate.high))
    cd45.flowD.ind <- cd45.flowD@index
    cd45.flowD@proportion <- (cd45.flowD@cell.count/nrow(lymph))*100
    cd45 <- getflowFrame(cd45.flowD)
    
    gthres[7] <- cd45.gate.low
    gthres[8] <- cd45.gate.high
    gthres[9] <- lincd3cd19cd161.gate.high
    all.props[5] <- cd45.flowD@proportion
    all.events[5] <- cd45.flowD@cell.count
    
    #############################################################################################
    # Some files are pregated. This just undoes the above for the pregated files
    # This is not the most efficient way to do this, could in fact skip code above in this case
    
    if(length(grep(pregate.mln.barcodes, x$FCS.file)) > 0){ 
    #if(i %in% pregate.idx){
      
      singlets.gate.h <- deGate(f, "SSC-W", use.percentile = T, percentile = 1)
      singlets.gate.l <- deGate(f, "SSC-W", use.percentile = T, percentile = 0)
      gthres[1] <- singlets.gate.h
      gthres[2] <- singlets.gate.l
      all.props[2] <- nrow(f)/nrow(f)*100
      all.events[2] <- nrow(f)
      singlets <- f
      
      live.gate <- deGate(singlets, channel = c(11), use.percentile = T, percentile = 1)
      gthres[3] <- live.gate
      all.props[3] <- nrow(f)/nrow(f)*100
      all.events[3] <- nrow(f)
      live <- f
      
      fsc.a.gate.low <- deGate(live, channel = c(1), use.percentile = T, percentile = 0)
      gthres[5] <- fsc.a.gate.low
      all.props[4] <- nrow(f)/nrow(f)*100
      all.events[4] <- nrow(f)
      lymph <- f
      
      cd45.gate.low <- deGate(lymph, channel = c(14), use.percentile = T, percentile = 0)
      gthres[7] <- cd45.gate.low
      all.props[5] <- nrow(f)/nrow(f)*100
      all.events[5] <- nrow(f)
      cd45 <- f
    
    }
    #############################################################################################
    ## Gating CD45 population. Plotting CD45_Lin(CD3, cd19, CD161) to get Lineage negative cells
    
    # Rotate CD45 data and threshold
    theta0 = -pi*25/180
    cd45s <- cd45
    cd45s <- rotate.data(cd45, c(14, 12), theta = -theta0)$data
    R <- matrix(c(cos(theta0), sin(theta0), -sin(theta0), cos(theta0)), 2, 2)
    
    # Find location the minimum (valley) between the Lin(CD3,CD19,cd161) hi and low peaks
    maxDens <- density(cd45s@exprs[, 12]) 
    temp <- findpeaks(maxDens$y[min(which(maxDens$x > 1.3)):length(maxDens$y)])
    n = length(temp[, 2])
    first.peak <- sort(temp[, 1], partial = n)[n]
    first.peak.idx <- which(maxDens$y == first.peak)
    #temp <- findpeaks(maxDens$y[1:min(which(maxDens$x > 1.19))])
    maxDensy.smooth <- as.numeric(filter(maxDens$y[1:min(which(maxDens$x > 1.19))], rep(1/11,11), sides = 2))
    temp <- findpeaks(maxDensy.smooth)
    n = length(temp[, 2])
    if(n > 0){
      # second.peak <- sort(temp[, 1], partial = n)[n]
      # second.peak.idx <- which(maxDens$y == second.peak)
      if(n > 1){
        second.peak.idx <- temp[which.max(temp[,1]),2]
      }else{
        second.peak.idx <- temp[2]
      }
    }else{ # no peaks found
      second.peak.idx <- which.min(abs(maxDens$x - 0.75))
    }
    # Check to see if there are are two peaks detected in the Lin(CD3,CD19,cd161) hi population 
    if(abs(maxDens$x[first.peak.idx] - maxDens$x[second.peak.idx]) < 0.4){ 
      second.peak <- sort(temp[, 1], partial = n-1)[n-1]
      second.peak.idx <- which(maxDens$y == second.peak)
    }
    if((maxDens$x[second.peak.idx] > 1.6) || (maxDens$x[second.peak.idx] < 0.1) ){
      second.peak.idx <- which.min(abs(maxDens$x-0.8))
    }
    valley <- which.min(maxDens$y[second.peak.idx:first.peak.idx])
    valley.yval <- maxDens$y[second.peak.idx + valley]
    
    lin.gate.slant <- max(maxDens$x[which(maxDens$y[1:first.peak.idx] < 1.1*valley.yval)])
    
    if((lin.gate.slant - maxDens$x[second.peak.idx]) < 0.1){
      lin.gate.slant <- c(lin.gate.slant,  deGate(cd45s, channel = c(12), use.upper = T, upper = F, alpha = 0.05, tinypeak.removal = 0.9))
      lin.gate.slant <- lin.gate.slant[which.min(abs(lin.gate.slant - maxDens$x[first.peak.idx] + 0.55))]
    }

    cd45sTempb <- flowDensity(cd45s, channels=c(14, 12),position=c(NA, F), gates = c(NA, lin.gate.slant))
    cd45sTempb.ind <- cd45sTempb@index
    lin.neg <- getflowFrame(cd45sTempb)
    lin.neg@exprs[, c(14, 12)] <- t(t(R) %*% t(lin.neg@exprs[, c(14, 12)]))
    
    cd45sTempb@filter <- rotate.data(cd45sTempb@filter, c(14, 12),theta = theta0)$data
    cd45sTempb@flow.frame <- rotate.data(getflowFrame(cd45sTempb), c(14, 12), theta = theta0)$data
    cd45sTempb@proportion <- (cd45sTempb@cell.count/nrow(cd45))*100
    lin.neg <- getflowFrame(cd45sTempb)
    
    # Get the Lineage positive population
    
    linpos.flowD <- flowDensity(cd45s, channels=c(14, 12),position=c(NA, T), gates = c(NA, lin.gate.slant))
    linpos.flowD.ind <- linpos.flowD@index
    lin.pos <- getflowFrame(linpos.flowD)
    lin.pos@exprs[, c(14, 12)] <- t(t(R) %*% t(lin.pos@exprs[, c(14, 12)]))
    
    linpos.flowD@filter <- rotate.data(linpos.flowD@filter, c(14, 12),theta = theta0)$data
    linpos.flowD@flow.frame <- rotate.data(getflowFrame(linpos.flowD), c(14, 12), theta = theta0)$data
    linpos.flowD@proportion <- (linpos.flowD@cell.count/nrow(cd45))*100
    lin.pos <- getflowFrame(linpos.flowD)
    
    gthres[10] <-  lin.gate.slant
    all.props[6] <- cd45sTempb@proportion
    all.events[6] <- cd45sTempb@cell.count
    
    #############################################################################################
    ## Gating CD45 population (i.e. the Singlet - Live - Lymphocyte - CD45+ population)
    # Linear Gates as Adam Liang suggested
    
    mhcII.gate <- deGate(cd45, channel = c(7))
    F4_80.gate <- deGate(cd45, channel = c(8), use.upper = T, upper = T, tinypeak.removal = 0.1)
    Ly6G.gate <- deGate (cd45, channel = c(9), use.upper = T, upper = T, tinypeak.removal = 0.1)
    CD11b.gate <- min(deGate (cd45, channel = c(13), all.cut = T, use.upper = T, upper = T, alpha = 0.025))
    CD317.gate <- min(deGate (cd45, channel = c(15), all.cut = T, use.upper = T, upper = T, alpha = 0.025))
    CD11c.gate <- min(deGate (cd45, channel = c(16), all.cut = T, use.upper = T, upper = T, alpha = 0.025))
    CD86.gate <- min(deGate(cd45, channel = c(18), all.cut = T, use.upper = T, upper = T))
    
    #--------Start Big Png------------
    png ( file = paste("Results_MLN/Figures/ScatterPlots_LinearGates/", x$Genotype, "/", "Total_", x$FCS.file, ".png", sep = "" ), width=2500, height=2500*3/5)
    par(mfrow=c(3,5),mar=(c(5, 5, 4, 2) + 0.1))
    #par(mfrow=c(5,4))
    
    plotDens(f,c("SSC-A","SSC-W"), main="Ungated", cex.lab = 2, cex.axis = 2, cex.main=2)
    #lines(singlets.flowD.l@filter,lwd =1)
    abline(h = singlets.gate.l, lwd =1)
    abline(h = singlets.gate.h, lwd =1)
    
    # Plotting Live/Dead_SSC-A 
    plotDens(singlets,  c(11,4), main= "Singlets", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v=live.gate, lwd =1); 
    
    # Plotting FSC-A_SSC-A
    plotDens(live, c(1,4), main = "Live", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v=fsc.a.gate.low, lwd =1); abline(v=fsc.a.gate.high, lwd =1)
    abline(h=ssc.a.gate, lwd =1)
    
    # Plotting CD45_Lin(CD3,CD19,CD161)
    plotDens(lymph, channels = c(14,12), main='Lymphocytes', cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v=cd45.gate.low, lwd =1); abline(v=cd45.gate.high, lwd =1)
    abline(h= lincd3cd19cd161.gate.high, lwd =1)
    #lines(cd45.flowD@filter,lwd =1)
    
    plotDens(cd45, channels = c(14,12), main="CD45", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(cd45sTempb@filter, lwd =1)
    
    plotDens(lin.pos, channels = c(7,13), main = "lin.pos", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = mhcII.gate)
    
    plotDens(lin.pos, channels = c(14,12), main = "lin.pos", cex.lab = 2, cex.axis = 2, cex.main=2)
    
    plotDens(lin.pos, channels = c(8,13), main = "lin.pos", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = F4_80.gate)
    
    plotDens(lin.pos, channels = c(9,13), main = "lin.pos", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = Ly6G.gate)
    
    plotDens(lin.pos, channels = c(13,10), main = "lin.pos", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = CD11b.gate)
    
    plotDens(lin.pos, channels = c(15,10), main = "lin.pos", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = CD317.gate)
    
    plotDens(lin.pos, channels = c(16,10), main = "lin.pos", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = CD11c.gate)
    
    plotDens(lin.pos, channels = c(18,17), main = "lin.pos", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = CD86.gate)
    
    dev.off()
    
    
    #--------End Big Png------------
    #cat("Time is: ",TimeOutput(start2),"\n",sep="")
    
    if (do.flowType == T) {
      #print("Start FlowType")
      listGthres <- list(mhcII.gate, F4_80.gate, Ly6G.gate, CD11b.gate, CD317.gate, CD11c.gate,  CD86.gate)
      
      # Run flowType on the Lineage Negative population
      flowType.res <- flowType(Frame = lin.pos, PropMarkers= c(7, 8, 9, 13, 15, 16, 18), MaxMarkersPerPop = NULL, PartitionsPerMarker=2,
                               Methods='Thresholds', Thresholds= listGthres, verbose=F, MemLimit=400)
      
      save ( flowType.res, file =  paste("Results_MLN/FlowType_LineagePositive/", x$Genotype,"/FT_", x$FCS.file,".Rdata",sep="") )
      
      # Run flowType on the Lineage Positive Population
    }
    
    data.frame(all.props, all.events, gthres)
    
    
  },error = function(err) {
    err
  }
  ) # end of tryCatch
  
}, .parallel = TRUE) # end ddply
# } # end of for-loop

################################################
# Write cell count and propoortion data to files
################################################
# if(write.data){
#   
#   colnames(gthres) <- c('singlets.gate.h', 'singlets.gate.l', 'live.gate', 'fsc.a.gate.high', 'fsc.a.gate.low', 'ssc.a.gate', 'cd45.gate.low', 
#                         'cd45.gate.high', 'lincd3cd19cd161.gate.high', 'lin.gate.slant', 'linneg.gate.slant', 'linneg.macneg.gate.slanta', 
#                         'linneg.macneg.gate.slantb', 'Ly6c.gate', 'cd11b.gate', 'Ly6c.gate2', 'cd317.gate', 'mhcII.gate', 'cd11c.gate', 
#                         'cDC.gate.slant', 'cd103.gate', 'linneg.macneg.gate.slanta2')
#   rownames(gthres) <- 1:length(gthres[,1])
#   save(gthres, file = paste("Results_MLN/ThreshTable.Rdata",sep="") )  
#   save(errorFileIndex, file = paste("Results_MLN/store.Gate.Error.Rdata",sep=""))
#   
#   colnames(all.props) <- c("All Events-%Parent", "Singlets-%Parent", "Live-%Parent", "Lymphocytes-%Parent", "CD45-%Parent", "Lineage neg-%Parent", 
#                            "Lin-Mac--%Parent", "RP macrophages", "Ly6G--%Parent", "Neutrophils 2-%Parent", "Eosinophils 3-%Parent", 
#                            "Monocytes Ly6c hi-%Parent","NOT(Monocytes Ly6c hi)-%Parent", "pDC-%Parent", "NOT pDC-%Parent", "cDC-%Parent", 
#                            "Misc-%Parent", "CD8A Type DC-%Parent","CD11B+ CD86Lo-%Parent", "CD103+ DC-%Parent")
#   rownames(all.props) <- 1:length(all.props[,1])
#   
#   colnames(all.events) <- c("All Events-#Events", "Singlets-#Events", "Live-#Events", "Lymphocytes-#Events", "CD45-#Events", "Lineage neg-#Events", 
#                             "Lin-Mac--#Events", "RP macrophages", "Ly6G--#Events", "Neutrophils 2-#Events", "Eosinophils 3-#Events", 
#                             "Monocytes Ly6c hi-#Events","NOT(Monocytes Ly6c hi)-#Events", "pDC-#Events", "NOT pDC-#Events", "cDC-%#Events", 
#                             "Misc-#Events", "CD8A Type DC-#Events","CD11B+ CD86Lo-#Events", "CD103+ DC-#Events")
#   rownames(all.events) <- 1:length(all.props[,1])
#   
#   Events_Proportions_Table <- NULL
#   for(i in 1:ncol(all.props)){
#     Events_Proportions_Table <- cbind(Events_Proportions_Table,all.events[,i],all.props[,i])
#   }
#   colnames(Events_Proportions_Table) <- c(rbind(colnames(all.events), colnames(all.props)))
#   rownames(Events_Proportions_Table) <- 1:length(Events_Proportions_Table[,1])
# 
#   all.props <- cbind(store.allFCS[,1:2], all.props)
#   all.events <- cbind(store.allFCS[,1:3], all.events)
#   
#   date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv")
#   write.csv(Events_Proportions_Table, file =  paste("Results_MLN/Events_Proportions_Table", date.time, sep=""))
#   write.csv(all.props, file =  paste("Results_MLN/allProportions_Table", date.time, sep=""))
#   write.csv(all.events, file =  paste("Results_MLN/allEvents_Table", date.time, sep=""))
# }

cat("Total time is: ",TimeOutput(start),"\n",sep="")

