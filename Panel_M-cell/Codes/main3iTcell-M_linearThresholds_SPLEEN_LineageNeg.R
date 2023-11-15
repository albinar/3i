# Originally written by Justin Meskas 
# Modified & Updated for the T-cell panel by Albina Rahim
# Date: Aug 17, 2016

remove(list=ls())


setwd("/data/Panel_M-cell")

library("flowCore")
library("flowBin")
library("flowDensity")
library("flowType")
# library("flowClean")
#library("snowfall")
#library("foreach")
library("stringr")
library('pracma')  
library('MASS')

source("Codes/3iTcellfunctions-M.R")
source("Codes/rotate.data-M.R")

# load("Results/lgl657.Rdata")
# load("Results/Genotype.Rdata")
# load("Results/uniqueGT.Rdata")
#load("Results/channels.ind.Rdata")
load("Results_SPLEEN/store.allFCS.Rdata")

start <- Sys.time()
write.data <- F
do.flowType <- T


# According to Adam, these files are pre-gated CD45+, Live, single cells. 
pregate.spln.barcodes <- c(53506, 53508:53518) 
pregate.spln.barcodes <- paste(pregate.spln.barcodes, collapse = "|")
# pregate.idx <- lapply(pregate.spln.barcodes, function(x){ grep(x, store.allFCS[,2])})


library(plyr)
library(doMC)
no_cores <- detectCores() - 1
registerDoMC(no_cores)

file.names <- data.frame(store.allFCS, stringsAsFactors = F)

props.events.gates <- ddply(file.names[1783:1785,], "FCS.file", function(x){
  
  # for(i in loop.idx){
  
  all.props <- matrix(nrow = 1, ncol = 6)
  all.events <- matrix(nrow = 1, ncol = 6)
  gthres <- matrix(nrow = 1, ncol = 26)
  
  possibleError <- tryCatch({
    # start2 <- Sys.time()
    
    print(paste("Starting ",  x$Genotype, "/", x$FCS.file, sep = ""))
    load(file = paste("Results_SPLEEN/After_Clean/",  x$Genotype, "/AftClean_", x$FCS.file, ".Rdata", sep=""))
    
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
      temp <- deGate(f, "SSC-W", use.upper = T, upper=T)
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
    names(live.gate) <- identifier(singlets)
    
    live.flowD <- flowDensity(singlets, channels = c(11, 4), position = c(F, NA), gates = c(live.gate, NA))
    live.flowD.ind <- live.flowD@index
    live <- getflowFrame(live.flowD)
    
    gthres[3] <- live.gate
    all.props[3] <- live.flowD@proportion
    all.events[3] <- live.flowD@cell.count
    
    ##############################################################################################
    ## Gating Live. Plotting FSC-A_SSC-A
    
    fsc.a.gate.low <- max(deGate(live, channel = c(1), all.cuts = T, use.upper = T, upper = F))
    # There are places where max(deGate()) is above the lymphocyte peak (~60,000)
    # Check if this is the case pick second highest gate if it is
    if(fsc.a.gate.low > 60000){
      temp <- deGate(live, channel = c(1), all.cuts = T, use.upper = T, upper = F)
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
    ## Gating Lymphocytes. Plotting CD45_Lin(CD3, CD109, CD161) to get CD45 cells
    
    # CD45
    cd45.gate.low <- max(deGate(lymph, channel = c(14), all.cuts = T, use.upper = T, upper=F))
    cd45.gate.high <- deGate(lymph, channel = c(14), use.percentile = T, percentile = 1)
    # gating below suggested by Albina to mimic her T cell cd161 gating
    lincd3cd109cd161.gate.high <- deGate(lymph, channel = c(12), use.percentile = T, percentile = 0.9999) + 0.1 
    
    names(cd45.gate.low) <- names(cd45.gate.high) <-  names(lincd3cd109cd161.gate.high) <- identifier(lymph)
    
    cd45.flowD.temp <- flowDensity(lymph, channels = c(14, 12), position = c(T, F), gates = c(cd45.gate.low, lincd3cd109cd161.gate.high))
    cd45.flowD <- flowDensity(cd45.flowD.temp, channels = c(14, 12), position = c(F, F), gates = c(cd45.gate.high, lincd3cd109cd161.gate.high))
    cd45.flowD.ind <- cd45.flowD@index
    cd45.flowD@proportion <- (cd45.flowD@cell.count/nrow(lymph))*100
    cd45 <- getflowFrame(cd45.flowD)
    
    gthres[7] <- cd45.gate.low
    gthres[8] <- cd45.gate.high
    gthres[9] <- lincd3cd109cd161.gate.high
    all.props[5] <- cd45.flowD@proportion
    all.events[5] <- cd45.flowD@cell.count
    
    #############################################################################################
    # Some files are pregated. This just undoes the above for the pregated files
    # This is not the most efficient way to do this, but is the simplest to delete if/when we get
    # the ungated files
    
    if(length(grep(pregate.spln.barcodes, x$FCS.file)) > 0){ 
      # if(i %in% pregate.idx){
      
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
    ## Gating CD45 population. Plotting CD45_Lin(CD3, CD109, CD161) to get Lineage negative cells
    
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
    temp <- findpeaks(maxDens$y[1:min(which(maxDens$x > 1.19))])
    n = length(temp[, 2])
    if(n > 0){
      second.peak <- sort(temp[, 1], partial = n)[n]
      second.peak.idx <- which(maxDens$y == second.peak)
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
    valley = which.min(maxDens$y[second.peak.idx:first.peak.idx])
    lin.gate.slant = maxDens$x[second.peak.idx + valley]
    
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
    # Split up into Lineage Positive and Lineage Negative
    
    # Linear Gates for the Lineage Positive population as Adam Liang suggested

    mhcII.gate <- deGate(lin.pos, channel = c(7))

    CD11b.gate <- deGate (lin.pos, channel = c(13), all.cuts = T, use.upper = T, upper = T, tinypeak.removal = 0.2)
    maxDens <- density(lin.pos@exprs[,c(13)])
    CD11b.gate <- CD11b.gate[which.min(abs(CD11b.gate - maxDens$x[which.max(maxDens$y) - 0.5]))]

    cd45.temp <- flowDensity(lin.pos, channels = c(8,10), position = c(NA, F), gates = c(NA, CD11b.gate))

    F4_80.gate <- deGate(cd45.temp@flow.frame, channel = c(8), all.cuts = T, use.upper = T, upper = T, alpha = 0.025, tinypeak.removal = 0.1)
    maxDens <- density(lin.pos@exprs[,c(8)])
    F4_80.gate <- F4_80.gate[which.min(abs(F4_80.gate - maxDens$x[which.max(maxDens$y) - 0.5]))]

    Ly6G.gate <- deGate (cd45.temp@flow.frame, channel = c(9), all.cuts = T, use.upper = T, upper = T, alpha = 0.025, tinypeak.removal = 0.1)
    maxDens <- density(lin.pos@exprs[,c(9)])
    Ly6G.gate <- Ly6G.gate[which.min(abs(Ly6G.gate - maxDens$x[which.max(maxDens$y) - 0.5]))]

    maxDens <- density(lin.pos@exprs[,c(15)])
    CD317.gate <- 2*maxDens$x[which.max(maxDens$y)] - deGate (lin.pos, channel = c(15), use.upper = T, upper = F)

    maxDens <- density(na.omit(lin.pos@exprs[,c(16)]), adjust = 0.1)
    peak.idx <- which.max(maxDens$y)
    CD11c.gate <- maxDens$x[min(which(maxDens$y[peak.idx:length(maxDens$y)] < 0.1*max(maxDens$y))) + peak.idx - 1] + 0.1

    maxDens <- density(lin.pos@exprs[,c(18)])
    CD86.gate <- deGate(lin.pos, channel = c(18), all.cuts = T, use.upper = T, upper = T)
    if(CD86.gate < maxDens$x[which.max(maxDens$y)]){
      CD86.gate <- deGate(lin.pos, channel = c(18), use.upper = T, upper = T)
    }

    
    # Linear Gates for the Lineage Negative population as Adam Liang suggested
    
    mhcII.ln.gate <-  deGate(lin.neg, all.cuts = T, channel = c(7), use.upper = T, upper = T)
    maxDens <- density(lin.neg@exprs[,c(7)])
    stop.idx <- min(which(maxDens$x > 2))
    peak.lcn <- maxDens$x[which.max(maxDens$y[1:stop.idx])]
    mhcII.ln.gate <- mhcII.ln.gate[which.min(abs(mhcII.ln.gate - peak.lcn - 0.5))]
    
    F4_80.ln.gate <- min(deGate(lin.neg, channel = c(8), use.upper = T, upper = T, tinypeak.removal = 0.9))
    
    maxDens <- density(lin.neg@exprs[,c(9)])
    stop.idx <- which.min(abs(maxDens$x - 2))
    Ly6G.ln.gate <- deGate(lin.neg, channel = c(9), all.cuts = T)
    Ly6G.ln.gate <- Ly6G.ln.gate[which.min(abs(Ly6G.ln.gate - maxDens$x[which.max(maxDens$y[1:stop.idx])] + 0.5))]
    if(Ly6G.ln.gate > 2.2){
      temp.flowD <- flowDensity(lin.neg, channels = c(14,9), position = c(NA, F), gates = c(NA, Ly6G.ln.gate))
      Ly6G.ln.gate <- c(Ly6G.ln.gate, 
                        deGate(temp.flowD@flow.frame, channel = c(9), all.cuts = T),
                        deGate(temp.flowD@flow.frame, channel = c(9), use.upper = T, upper = T, tinypeak.removal = 0.9))
      Ly6G.ln.gate <- Ly6G.ln.gate[which.min(abs(Ly6G.ln.gate - maxDens$x[which.max(maxDens$y[1:stop.idx])] + 0.5))]
    }
    
    CD11b.ln.gate <- deGate(lin.neg, channel = c(13), all.cuts = T)
    maxDens <- density(lin.neg@exprs[,c(13)])

    peak.lcns <- findpeaks(maxDens$y)
    peak.lcns <- peak.lcns[which(peak.lcns[,1] > 0.1*max(peak.lcns[,1])),]
    peak.lcns <- peak.lcns[which.min(abs(maxDens$x[peak.lcns[,2]] - 2)), 2]
    
    CD11b.ln.gate <-  CD11b.ln.gate[which.min(abs(maxDens$x[peak.lcns] + 0.25 - CD11b.ln.gate))]
                                                   
    temp <- lin.neg@exprs[, c(14,15)]
    temp <- temp[which(temp[,1] < 2.7),]
    maxDens2 <- density(temp[,2])
    
    start.idx <- max(which(maxDens2$x < 2.1))
    
    peak.lcn <- findpeaks(maxDens2$y[start.idx:length(maxDens2$y)])
    if(length(peak.lcn) > 4){
      peak.lcn <- peak.lcn[which.max(peak.lcn[,1]), 2]
    }else{
      peak.lcn <- peak.lcn[2]
    }

    temp.flowD <- flowDensity(lin.neg, channels = c(14,15), position = c(F, NA), gates = c(2.6, NA))
    CD317.ln.gate <- c(deGate(temp.flowD@flow.frame, channel = c(15), all.cuts = T),
                       deGate(temp.flowD@flow.frame, channel = c(15), use.upper = T, upper = T, tinypeak.removal = 0.2),
                       2*maxDens2$x[start.idx - 1 + peak.lcn] - deGate(temp.flowD@flow.frame, channel = c(15), use.upper = T, upper = T))
    CD317.ln.gate <- CD317.ln.gate[which.min(abs(CD317.ln.gate - maxDens2$x[start.idx - 1 + peak.lcn] + 0.25))]
    
    maxDens <- density(lin.neg@exprs[,c(16)])
    CD11c.ln.gate <- deGate(lin.neg, channel = c(16), use.upper = T, upper = T, tinypeak.removal = 0.9)
    CD11c.ln.gate <- CD11c.ln.gate[which.min(abs(CD11c.ln.gate - maxDens$x[which.max(maxDens$y)] - 0.5))]

    maxDens <- density(lin.neg@exprs[,c(17)])
    CD103.ln.gate <- deGate(lin.neg, channel = c(17), use.upper = T, upper = T, tinypeak.removal = 0.9, alpha = 0.5)
    CD103.ln.gate <-  CD103.ln.gate[which.min(abs(CD103.ln.gate- maxDens$x[which.max(maxDens$y)] - 0.5))]
    
    maxDens <- density(lin.neg@exprs[,c(18)])
    CD86.ln.gate <- deGate(lin.neg, channel = c(18), use.upper = T, upper = T, tinypeak.removal = 0.9, alpha = 0.5)
    CD86.ln.gate <- CD86.ln.gate[which(CD86.ln.gate > maxDens$x[which.max(maxDens$y)])]
    CD86.ln.gate <-  CD86.ln.gate[which.min(abs(CD86.ln.gate- maxDens$x[which.max(maxDens$y)] - 0.7))]
    if(CD86.ln.gate > 2.4){
      CD86.ln.gate <- c(CD86.ln.gate, deGate(lin.neg, channel = c(18), use.percentile = F, percentile = NA))
      CD86.ln.gate <- CD86.ln.gate[which(CD86.ln.gate > maxDens$x[which.max(maxDens$y)])]
      CD86.ln.gate <-  CD86.ln.gate[which.min(abs(CD86.ln.gate- maxDens$x[which.max(maxDens$y)] - 0.7))]
    }
    
    maxDens <- density(lin.neg@exprs[,c(4)])
    ssca.ln.gate <- c(deGate(lin.neg, channel = c(4), all.cuts = T),
                      deGate(lin.neg, channel = c(4), use.upper = T, upper = T, tinypeak.removal = 0.9))
    ssca.ln.gate <-  ssca.ln.gate[which.min(abs(ssca.ln.gate- maxDens$x[which.max(maxDens$y)] - 25000))]
    if(ssca.ln.gate > 90000){
      temp.flowD <-flowDensity(lin.neg, channels = c(14,4), position = c(NA, F), gates = c(NA, ssca.ln.gate))
      ssca.ln.gate <- deGate(temp.flowD@flow.frame, channel = c(4), all.cuts = T, percentile = NA)
      ssca.ln.gate <-  ssca.ln.gate[which.min(abs(ssca.ln.gate- maxDens$x[which.max(maxDens$y)] - 25000))]
    }
    
    gthres[11:26] <- c(mhcII.gate, F4_80.gate, Ly6G.gate, CD11b.gate, CD317.gate, CD11c.gate,  CD86.gate,
               mhcII.ln.gate, F4_80.ln.gate, Ly6G.ln.gate, CD11b.ln.gate, CD317.ln.gate, 
              CD11c.ln.gate, CD103.ln.gate, CD86.ln.gate, ssca.ln.gate)
    
    
    #--------Start Big Png------------
    png ( file = paste("Results_SPLEEN/Figures/ScatterPlots_LinearGates/", x$Genotype, "/", "Total_", x$FCS.file, ".png", sep = "" ), width=2500, height=2500)
    par(mfrow=c(5,5),mar=(c(5, 5, 4, 2) + 0.1))
    
    plotDens(f,c("SSC-A","SSC-W"), main="Ungated", cex.lab = 2, cex.axis = 2, cex.main=2)
    #lines(singlets.flowD.l@filter,lwd=1)
    abline(h = singlets.gate.l, lwd=1)
    abline(h = singlets.gate.h, lwd=1)
    
    # Plotting Live/Dead_SSC-A 
    plotDens(singlets,  c(11,4), main= "Singlets", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v=live.gate, lwd=1); 
    
    # Plotting FSC-A_SSC-A
    plotDens(live, c(1,4), main = "Live", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v=fsc.a.gate.low, lwd=1); abline(v=fsc.a.gate.high, lwd=1)
    abline(h=ssc.a.gate, lwd=1)
    
    # Plotting CD45_Lin(CD3,CD19,CD161)
    plotDens(lymph, channels = c(14,12), main='Lymphocytes', cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v=cd45.gate.low, lwd=1); abline(v=cd45.gate.high, lwd=1)
    abline(h= lincd3cd109cd161.gate.high, lwd=1)
    #lines(cd45.flowD@filter,lwd=1)
    
    plotDens(cd45, channels = c(14,12), main="CD45", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(cd45sTempb@filter, lwd =1)
    
    # Plots for the Lineage Negative Population
    
    plotDens(lin.neg, channels = c(14,7), main = "Lineage Negative", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(h = mhcII.ln.gate)
    
    plotDens(lin.neg, channels = c(14,8), main = "Lineage Negative", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(h = F4_80.ln.gate)
    
    plotDens(lin.neg, channels = c(14,9), main = "Lineage Negative", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(h = Ly6G.ln.gate)
    
    plotDens(lin.neg, channels = c(14,13), main = "Lineage Negative", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(h = CD11b.ln.gate)
    
    plotDens(lin.neg, channels = c(14,15), main = "Lineage Negative", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(h=CD317.ln.gate)
    
    plotDens(lin.neg, channels = c(14,16), main = "Lineage Negative", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(h=CD11c.ln.gate)
    
    plotDens(lin.neg, channels = c(14,17), main = "Lineage Negative", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(h = CD103.ln.gate)
    
    plotDens(lin.neg, channels = c(14,18), main = "Lineage Negative", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(h = CD86.ln.gate)
    
    maxDens <- density(lin.neg@exprs[,c(18)])
    plot(maxDens, cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = CD86.ln.gate)
    
    plotDens(lin.neg, channels = c(14,4), main = "Lineage Negative", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(h = ssca.ln.gate)
    
    #plot(1, type="n", axes=F, xlab="", ylab="")
  
    
    # Plots for the Lineage Positive Population
    
    plotDens(lin.pos, channels = c(7,13), main = "Lineage Postive", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = mhcII.gate)
    
    plotDens(lin.pos, channels = c(14,12), main = "Lineage Postive", cex.lab = 2, cex.axis = 2, cex.main=2)
    
    plotDens(lin.pos, channels = c(8,13), main = "Lineage Postive", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = F4_80.gate)
    
    plotDens(lin.pos, channels = c(9,13), main = "Lineage Postive", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = Ly6G.gate)
    
    plotDens(lin.pos, channels = c(13,10), main = "Lineage Postive", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = CD11b.gate)
    
    plotDens(lin.pos, channels = c(15,10), main = "Lineage Postive", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = CD317.gate)
    
    plotDens(lin.pos, channels = c(16,10), main = "Lineage Postive", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = CD11c.gate)
    
    plotDens(lin.pos, channels = c(18,17), main = "Lineage Postive", cex.lab = 2, cex.axis = 2, cex.main=2)
    abline(v = CD86.gate)
    
    plot(1, type="n", axes=F, xlab="", ylab="")
    plot(1, type="n", axes=F, xlab="", ylab="")
  
    dev.off()
    
    if (do.flowType == T) {
      
      if(nrow(lin.pos) > 2000){
      # Run flowType on the Lineage Positive population
      #print("Start FlowType")
      listGthres <- list(mhcII.gate, F4_80.gate, Ly6G.gate, CD11b.gate, CD317.gate, CD11c.gate,  CD86.gate)
      flowType.res <- flowType(Frame = lin.pos, PropMarkers= c(7, 8, 9, 13, 15, 16, 18), MaxMarkersPerPop = NULL, PartitionsPerMarker=2,
                               Methods='Thresholds', Thresholds= listGthres, verbose=F, MemLimit=400)
      save ( flowType.res, file =  paste("Results_SPLEEN/FlowType_LineagePositive/", x$Genotype,"/FT_", x$FCS.file,".Rdata",sep="") )
      }
      
      if(nrow(lin.neg) > 2000){
        # Run flowType on Lineage Negative population
        #print("Start FlowType")
        listGthres.ln <- list(mhcII.ln.gate, F4_80.ln.gate, Ly6G.ln.gate, CD11b.ln.gate, CD317.ln.gate, 
                              CD11c.ln.gate, CD103.ln.gate, CD86.ln.gate, ssca.ln.gate)
        flowType.res <- flowType(Frame = lin.neg, PropMarkers= c(7, 8, 9, 13, 15, 16, 17, 18, 4), MaxMarkersPerPop = NULL, PartitionsPerMarker=2,
                                 Methods='Thresholds', Thresholds= listGthres.ln, verbose=F, MemLimit=400)
        save ( flowType.res, file =  paste("Results_SPLEEN/FlowType_LineageNegative/", x$Genotype,"/FT_", x$FCS.file,".Rdata",sep="") )
      
      }
    }
    
  },error = function(err) {
    err
  }) # end of tryCatch
  
  return(data.frame(all.props, all.events, gthres))
  
}, .parallel = TRUE) # end ddply

# } # end of for-loop

# Check for errors?!
error.files.indices <- which(is.na(props.events.gates[,2]))
if(length(error.files.indices) > 0){
  print('The following files failed gating:')
  print(file.names[error.files.indices, 2])
}

if(write.data){
  
  colnames(props.events.gates) <- c('FCS.file', 
                                    "All Events-%Parent", "Singlets-%Parent", "Live-%Parent", "Lymphocytes-%Parent", "CD45-%Parent", "Lineage neg-%Parent",
                                    "All Events-#Events", "Singlets-#Events", "Live-#Events", "Lymphocytes-#Events", "CD45-#Events", "Lineage neg-#Events", 
                                    'singlets.gate.h', 'singlets.gate.l', 'live.gate', 'fsc.a.gate.high', 'fsc.a.gate.low', 'ssc.a.gate', 'cd45.gate.low', 
                                    'cd45.gate.high', 'lincd3cd109cd161.gate.high', 'lin.gate.slant', 'mhcII.gate', 'F4_80.gate', 'Ly6G.gate', 'CD11b.gate',
                                    'CD317.gate', 'CD11c.gate', 'CD86.gate', 'mhcII.ln.gate', 'F4_80.ln.gate', 'Ly6G.ln.gate', 'CD11b.ln.gate', 'CD317.ln.gate', 
                                    'CD11c.ln.gate', 'CD103.ln.gate', 'CD86.ln.gate', 'ssca.ln.gate')
  save(props.events.gates, file = paste("Results_SPLEEN/PropsEventsGates.Rdata",sep="") ) 
  
}


# ################################################
# # Write cell count and propoortion data to files
# ################################################
# if(write.data){
#   
#   colnames(gthres) <- c('singlets.gate.h', 'singlets.gate.l', 'live.gate', 'fsc.a.gate.high', 'fsc.a.gate.low', 'ssc.a.gate', 'cd45.gate.low', 
#                         'cd45.gate.high', 'lincd3cd109cd161.gate.high', 'lin.gate.slant', 'linneg.gate.slant', 'linneg.macneg.gate.slanta', 
#                         'linneg.macneg.gate.slantb', 'Ly6c.gate', 'cd11b.gate', 'Ly6c.gate2', 'cd317.gate', 'mhcII.gate', 'cd11c.gate', 
#                         'cDC.gate.slant', 'cd103.gate')
#   rownames(gthres) <- 1:length(gthres[,1])
#   save(gthres, file = paste("Results/ThreshTable.Rdata",sep="") )  
#   save(errorFileIndex, file = paste("Results/store.Gate.Error.Rdata",sep=""))
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
#   # 
#   # # # Albina constructed CSVfile.Rdata by binding rows of the following .csv files in Separate_File_Folders-M.R
#   # # # /mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2015-08-27/attachments/3i_IMPC_Data_Genotypes.csv
#   # # # /mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2016-03-08/attachments/SPLN data RB.csv
#   # # # I diffed these files - most, but not all of the lines in te first file are duplicated in the second,
#   # # # Thus, CSVfile.Rdata has many duplicated rows. 
#   # # 
#   # load("Results/CSVfile.Rdata")
#   # Barcodes <- str_extract(store.allFCS[, 2],"L[0-9]+")
#   # Genotype <- str_replace(store.allFCS[, 1], '_', '/')
#   # Mouse_Label <-  CSVfile[, 'Label.Barcode']   
#   # # Thus when I get the indices for the barcodes I just arbitrarily pick the largest because as far as I can tell
#   # # the rows in CSVfile pertaining to one barcode are the same
#   # indices <- sapply(1:length(Barcodes), function(xdat){current.index <- max(grep(Barcodes[xdat], Mouse_Label))})
#   # AssayDate <- CSVfile[indices,'Assay.Date']
#   # Events_Proportions_Table <- cbind(store.allFCS[,1:2], Barcodes, AssayDate, store.allFCS[,3], Events_Proportions_Table)
#   # colnames(Events_Proportions_Table)[5] <- "Total Number of Cells"
#   all.props <- cbind(store.allFCS[,1:2], all.props)
#   all.events <- cbind(store.allFCS[,1:3], all.events)
#   
#   date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv")
#   write.csv(Events_Proportions_Table, file =  paste("Results/Events_Proportions_Table", date.time, sep=""))
#   write.csv(all.props, file =  paste("Results/allProportions_Table", date.time, sep=""))
#   write.csv(all.events, file =  paste("Results/allEvents_Table", date.time, sep=""))
# }

cat("Total time is: ",TimeOutput(start),"\n",sep="")

