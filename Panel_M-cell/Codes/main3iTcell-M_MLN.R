# Developed by Albina Rahim
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

load("Results/store.allFCS.Rdata")

start <- Sys.time()
errorFileIndex <- NULL

loop.idx <- 1:nrow(store.allFCS)
# Ly6Ghi.SSCAlo.files.barcodes <- c(92105, 109527,92123, 96523,96531, 96537, 55606, 57930, 70281, 84494, 88737, 97854, 88740, 88741, 89595, 93968,
#                                   97849, 93985, 97840, 64391, 97838, 97848, 108033, 108038, 109526, 109527, 90429)
# loop.idx <- unlist(lapply(Ly6Ghi.SSCAlo.files.barcodes, function(x){ grep(x, store.allFCS[,2])}))
write.data <- T # write cell count & proportion data to cvs file

# Some of the files are pre-gated CD45+, Live, single cells. 
pregate.mln.barcodes <- c(53519, 53521:53524, 53527:53531) 
pregate.idx <- lapply(pregate.mln.barcodes, function(x){ grep(x, store.allFCS[,2])})

all.props <- matrix(nrow = nrow(store.allFCS), ncol = 20)
all.events <- matrix(nrow = nrow(store.allFCS), ncol = 20)
gthres <- matrix(nrow = nrow(store.allFCS), ncol = 22)

for(i in loop.idx){

  possibleError <- tryCatch({
    start2 <- Sys.time()
    
    print(paste(i, ": Starting ", store.allFCS[i,1], "/", store.allFCS[i,2], sep = ""))
    load(file = paste("Results/After_Clean/", store.allFCS[i,1], "/AftClean_", store.allFCS[i,2], ".Rdata", sep=""))
    
    all.props[i,1] <- nrow(f)/nrow(f)*100
    all.events[i,1] <- nrow(f)
    
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
    
    gthres[i,1] <- singlets.gate.h
    gthres[i,2] <- singlets.gate.l
    all.props[i,2] <- singlets.flowD.l@proportion
    all.events[i,2] <- singlets.flowD.l@cell.count
    
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
    
    gthres[i,3] <- live.gate
    all.props[i,3] <- live.flowD@proportion
    all.events[i,3] <- live.flowD@cell.count
    
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
    
    gthres[i,4] <- fsc.a.gate.high
    gthres[i,5] <- fsc.a.gate.low
    gthres[i,6] <- ssc.a.gate
    all.props[i,4] <- lymph.flowD@proportion
    all.events[i,4] <- lymph.flowD@cell.count
    
    #############################################################################################
    ## Gating Lymphocytes. Plotting CD45_Lin(CD3, CD109, CD161) to get CD45 cells
    
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
    lincd3cd109cd161.gate.high <- deGate(lymph, channel = c(12), use.percentile = T, percentile = 0.9999) + 0.1 
    
    names(cd45.gate.low) <- names(cd45.gate.high) <-  names(lincd3cd109cd161.gate.high) <- identifier(lymph)
    
    cd45.flowD.temp <- flowDensity(lymph, channels = c(14, 12), position = c(T, F), gates = c(cd45.gate.low, lincd3cd109cd161.gate.high))
    cd45.flowD <- flowDensity(cd45.flowD.temp, channels = c(14, 12), position = c(F, F), gates = c(cd45.gate.high, lincd3cd109cd161.gate.high))
    cd45.flowD.ind <- cd45.flowD@index
    cd45.flowD@proportion <- (cd45.flowD@cell.count/nrow(lymph))*100
    cd45 <- getflowFrame(cd45.flowD)
    
    gthres[i,7] <- cd45.gate.low
    gthres[i,8] <- cd45.gate.high
    gthres[i,9] <- lincd3cd109cd161.gate.high
    all.props[i,5] <- cd45.flowD@proportion
    all.events[i,5] <- cd45.flowD@cell.count
    
    #############################################################################################
    # Some files are pregated. This just undoes the above for the pregated files
    # This is not the most efficient way to do this, could in fact skip code above in this case
    
    if(i %in% pregate.idx){
      
      singlets.gate.h <- deGate(f, "SSC-W", use.percentile = T, percentile = 1)
      singlets.gate.l <- deGate(f, "SSC-W", use.percentile = T, percentile = 0)
      gthres[i,1] <- singlets.gate.h
      gthres[i,2] <- singlets.gate.l
      all.props[i,2] <- nrow(f)/nrow(f)*100
      all.events[i,2] <- nrow(f)
      singlets <- f
      
      live.gate <- deGate(singlets, channel = c(11), use.percentile = T, percentile = 1)
      gthres[i,3] <- live.gate
      all.props[i,3] <- nrow(f)/nrow(f)*100
      all.events[i,3] <- nrow(f)
      live <- f
      
      fsc.a.gate.low <- deGate(live, channel = c(1), use.percentile = T, percentile = 0)
      gthres[i,5] <- fsc.a.gate.low
      all.props[i,4] <- nrow(f)/nrow(f)*100
      all.events[i,4] <- nrow(f)
      lymph <- f
      
      cd45.gate.low <- deGate(lymph, channel = c(14), use.percentile = T, percentile = 0)
      gthres[i,7] <- cd45.gate.low
      all.props[i,5] <- nrow(f)/nrow(f)*100
      all.events[i,5] <- nrow(f)
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
    
    gthres[i,10] <-  lin.gate.slant
    all.props[i,6] <- cd45sTempb@proportion
    all.events[i,6] <- cd45sTempb@cell.count
    
    ############################################################################################
    ## Gating Lin-. Plotting F4/80_CD11b to obtain Lin-Mac- and RP Macrophages
    
    # Rotate Lin- data and threshold
    theta0 = pi*25/180
    lin.neg.s <- lin.neg
    lin.neg.s <- rotate.data(lin.neg, c(8, 13), theta = -theta0)$data
    R <- matrix(c(cos(theta0), sin(theta0), -sin(theta0), cos(theta0)), 2 , 2)
    
    linneg.gate.slant <- deGate(lin.neg.s, channel=c(8))
    if(!is.null(names(linneg.gate.slant))){ # using 95% cutoff 
      linneg.gate.slant <- deGate(lin.neg.s, channel=c(8), use.upper = T, upper = T, alpha = 0.05, tinypeak.removal = 0.9) # Works for MLN because insignif. RP Mac population
    }else if(linneg.gate.slant < 0.5){ # Gate lower than expected, likely b/c no RP mac pop
      linneg.gate.slant <- deGate(lin.neg.s, channel=c(8), use.upper = T, upper = T)
    }
    
    linneg.Tempb <- flowDensity(lin.neg.s, channels=c(8, 13),position=c(F, NA), gates = c(linneg.gate.slant, NA))
    linneg.Tempb.ind <- linneg.Tempb@index
    linneg.macneg <- getflowFrame(linneg.Tempb)
    linneg.macneg@exprs[, c(8, 13)] <- t(t(R) %*% t(linneg.macneg@exprs[, c(8, 13)]))
    
    linneg.Tempb@filter <- rotate.data(linneg.Tempb@filter, c(8, 13), theta = theta0)$data
    linneg.Tempb@flow.frame <- rotate.data(getflowFrame(linneg.Tempb), c(8, 13),theta = theta0)$data
    linneg.Tempb@proportion <- (linneg.Tempb@cell.count/nrow(lin.neg.s))*100
    linneg.macneg <- getflowFrame(linneg.Tempb)
    
    gthres[i,11] <-  linneg.gate.slant
    all.props[i,7] <- linneg.Tempb@proportion
    all.events[i,7] <- linneg.Tempb@cell.count

    # Cell count of RP Macrophages = Cell count of Lin- minus cell count of Lin-Mac- 
    RP.mac.cell.count <- nrow(lin.neg) - nrow(linneg.macneg)
    all.props[i,8] <- RP.mac.cell.count/nrow(lin.neg)*100
    all.events[i,8] <- RP.mac.cell.count

    ############################################################################################
    ## Gating Lin-Mac-. Plotting Ly6G_SSC-A to obtain Neutrophils2, Ly6G-, Eosinohils3
    
    # There are some differences between MLN and SPLEEN gating becuase I found a SSC-A- Ly6G+ population in some of the files
    #
    # I asked Adam Liang about how to gate this population
    # His response: This is very strange indeed! I have had a look at the original FCS files - The Ly6G+ population is definitely not biological,
    # I suspect a technical issue but i can’t put my finger on just what it is. I think the solution is to just include the entire SSA-
    # population as the Ly6G- gate. By eye at least, this does not appear to have any major impact on the downstream populations although we do
    # need to keep an eye on these files for outliers post analysis.
    
    # Get all three cell populations using linear thresholds a rotated frame 
    theta0 = atan(-150000)
    linneg.macneg.s <- linneg.macneg
    linneg.macneg.s <- rotate.data(linneg.macneg,c(9,4),theta = -theta0)$data
    R <- matrix(c(cos(theta0), sin(theta0), -sin(theta0), cos(theta0)), 2 , 2)

    # Gate dividing Neu2 from Ly6G-
    linneg.macneg.gate.slanta <- max(deGate(linneg.macneg.s, channel=c(4), all.cut = T, use.upper = T, upper = F, tinypeak.removal = 0.1))
    if(linneg.macneg.gate.slanta > -1.7){
      linneg.macneg.gate.slanta <- linneg.macneg.gate.slanta - 0.25
    }
    temp <- flowDensity(linneg.macneg.s, channels =c(9, 4), position = c(NA, T), gates = c(NA, linneg.macneg.gate.slanta)) 
    
    # Gate dividing Eosin3 from Ly6G-
    linneg.macneg.gate.slantb <- c(deGate(temp@flow.frame, channel=c(9), all.cut = T, use.upper = T, upper = T, alpha = 0.05, tinypeak.removal = 0.9),
                                   deGate(temp@flow.frame, channel=c(9), use.upper = T, upper = T, alpha= 0.05, tinypeak.removal = 0.9))
    maxDens <- density(na.omit(temp@flow.frame@exprs[,c(9)]))
    peak.lcn <- maxDens$x[which.max(maxDens$y)]
    peak.idx <- which.max(maxDens$y)
    if(min(linneg.macneg.gate.slantb) > 100000){
      linneg.macneg.gate.slantb <- c(linneg.macneg.gate.slantb, 
                                     maxDens$x[min(which(maxDens$y[(peak.idx+1):length(maxDens$y)] < 0.1*max(maxDens$y))) + peak.idx],
                                     deGate(temp@flow.frame, channel=c(9), use.upper = T, upper = T, alpha = 0.25))
    }
    linneg.macneg.gate.slantb <- linneg.macneg.gate.slantb[which.min(abs(linneg.macneg.gate.slantb - (peak.lcn + 50000)))]
    
    # Check that I did not place the gate at to the left of the Ly6G- population, and if so move to the right
    Ly6GnegTempb <-flowDensity(linneg.macneg.s, channels =c(9, 4), position = c(F, NA), gates = c(linneg.macneg.gate.slantb, NA)) 
    temp <- flowDensity(Ly6GnegTempb, channels =c(9, 4), position = c(NA, F), gates = c(NA, linneg.macneg.gate.slanta)) 
    if((temp@cell.count/linneg.Tempb@cell.count) > 0.14){ 
      linneg.macneg.gate.slanta <- deGate(linneg.macneg.s, channel=c(4), use.upper = T, upper = F, tinypeak.removal = 0.1)
    }
    
    # The gate to the right of the Ly6G- population may not work for the Eosinophils, in particular if there is the
    # not biological SSC-A- Ly6G+ population. In this case the Eosinophils get a different gate.
    if(linneg.macneg.gate.slanta < -2){
      temp <- flowDensity(linneg.macneg.s, channels =c(9, 4), position = c(T, NA), gates = c(linneg.macneg.gate.slantb, NA)) 
      linneg.macneg.gate.slanta2 <- deGate(temp@flow.frame, channel = c(4), use.upper = T, upper = F, alpha = 0.02, tinypeak.removal = 0.9)
    }else{
      linneg.macneg.gate.slanta2 <- linneg.macneg.gate.slanta
    }
    
    eosinTempb <-flowDensity(linneg.macneg.s, channels =c(9, 4), position = c(T, T), gates = c(linneg.macneg.gate.slantb, linneg.macneg.gate.slanta2)) 
    eosinTempb.ind <- eosinTempb@index
    eosin3 <- getflowFrame(eosinTempb)
    eosin3@exprs[, c(9, 4)] <- t(t(R) %*% t(eosin3@exprs[, c(9, 4)]))
    eosinTempb@filter <- rotate.data(eosinTempb@filter, c(9, 4),theta = theta0)$data
    eosinTempb@flow.frame <- rotate.data(getflowFrame(eosinTempb), c(9, 4), theta = theta0)$data
    eosinTempb@proportion <- (eosinTempb@cell.count/nrow(linneg.macneg))*100
    eosin3 <- getflowFrame(eosinTempb)
    
    Ly6GnegTempb <-flowDensity(linneg.macneg.s, channels =c(9, 4), position = c(F, T), gates = c(linneg.macneg.gate.slantb, linneg.macneg.gate.slanta)) 
    Ly6GnegTempb.ind <- Ly6GnegTempb@index
    Ly6Gneg <- getflowFrame(Ly6GnegTempb)
    Ly6Gneg@exprs[, c(9, 4)] <- t(t(R) %*% t(Ly6Gneg@exprs[, c(9, 4)]))
    Ly6GnegTempb@filter <- rotate.data(Ly6GnegTempb@filter, c(9, 4),theta = theta0)$data
    Ly6GnegTempb@flow.frame <- rotate.data(getflowFrame(Ly6GnegTempb), c(9, 4), theta = theta0)$data
    Ly6GnegTempb@proportion <- (Ly6GnegTempb@cell.count/nrow(linneg.macneg))*100
    Ly6Gneg <- getflowFrame(Ly6GnegTempb)
      
    # Take the Neu2 population to be all cells that is not Ly6G- or EOSIN3
    all.events[i,10] <- nrow(linneg.macneg) - nrow(Ly6Gneg) - nrow(eosin3)
    all.props[i,10] <- all.events[i,10]/nrow(linneg.macneg)*100

    gthres[i,22] <- linneg.macneg.gate.slanta2
    gthres[i,12] <- linneg.macneg.gate.slanta
    gthres[i,13] <- linneg.macneg.gate.slantb
    all.props[i,9] <- nrow(Ly6Gneg)/nrow(linneg.macneg)*100
    all.props[i,11] <- nrow(eosin3)/nrow(linneg.macneg)*100
    all.events[i,9] <- nrow(Ly6Gneg)
    all.events[i,11] <- nrow(eosin3)
    
    ############################################################################################
    ## Gating Ly6G-. Plotting Ly6C_CD11b to obtain Monocytes Ly6c hi and NOT(Monocytes Ly6c hi)
    
    # Monocytes Ly6c hi
    # set the cd11b gate using only Ly6c hi data
    temp <- flowDensity(Ly6Gneg, channels = c(10, 13), position = c(T, NA), gates = c(2, NA))
    cd11b.gate <- deGate(getflowFrame(temp), channel = c(13))
    if(cd11b.gate < 2){
      temp.gate <- deGate(getflowFrame(temp), channel = c(13), all.cut = T)
      if(!is.null(names(temp.gate))){ # using 95% cutoff 
        maxDens <- density(getflowFrame(temp)@exprs[, 13])
        temp.gate <- maxDens$x[which.max(maxDens$y)] + 0.7
      }
      cd11b.gate <- max(temp.gate)
    }
    # Set the Ly6c gate using only the CD11b hi data
    temp <- flowDensity(Ly6Gneg, channels = c(10, 13), position = c(NA, T), gates = c(NA, cd11b.gate))
    Ly6c.gate <- deGate(getflowFrame(temp), channel = c(10), all.cut = T)
    if(length(Ly6c.gate)>1){
      # find gate closest to expected gate location
      Ly6c.gate <- Ly6c.gate[which.min(abs(Ly6c.gate - 2.35))]
    }
    # If the Ly6c gate is low try re-gating using the data above this gate only 
    if(Ly6c.gate < 2.1){
      temp <- flowDensity(Ly6Gneg, channels = c(10, 13), position = c(T, T), gates = c(Ly6c.gate, cd11b.gate))
      temp.gate <- deGate(getflowFrame(temp), channel = c(10), all.cut = T)
      if(length(temp.gate)>1){
        temp.gate <- temp.gate[which.min(abs(temp.gate - 2.35))]
      }
      if(abs(temp.gate -2.35) < abs(Ly6c.gate-2.35)){
        Ly6c.gate <- temp.gate
      }
    }
    names(Ly6c.gate) <- names(cd11b.gate) <- identifier(Ly6Gneg)
    
    mono.Ly6c.hi.flowD <- flowDensity(Ly6Gneg, channels = c(10, 13), position = c(T, T), gates = c(Ly6c.gate, cd11b.gate))
    mono.Ly6c.hi <- getflowFrame(mono.Ly6c.hi.flowD)

    gthres[i,14] <- Ly6c.gate
    gthres[i,15] <- cd11b.gate
    all.props[i,12] <- mono.Ly6c.hi.flowD@proportion
    all.events[i,12] <- mono.Ly6c.hi.flowD@cell.count
    
    # Take complement of the mono.Ly6c.hi data to obtain NOT(Monocytes Ly6c hi)
    not.mono.Ly6c.hi <- Ly6Gneg
    not.mono.Ly6c.hi@exprs <- not.mono.Ly6c.hi@exprs[-mono.Ly6c.hi.flowD@index, ]
    
    all.props[i,13] <- nrow(not.mono.Ly6c.hi)/nrow(Ly6Gneg)*100
    all.events[i,13] <- nrow(not.mono.Ly6c.hi)
    
    ############################################################################################
    ## Gating NOT(Monocytes Ly6c hi). Plotting Ly6C_CD317 to obtain pDC and NOT pDC
    
    # pDC
    Ly6c.gate2 <- deGate(not.mono.Ly6c.hi, channel = c(10))
    temp <- flowDensity(not.mono.Ly6c.hi, channels = c(10, 15), position = c(T, NA), gates = c(Ly6c.gate2,  NA))
    cd137.gate <- deGate(getflowFrame(temp), channel = c(15))
    if(cd137.gate < 1.7){ # gate is lower than expected
      if(!is.null(names(cd137.gate))){ # using 95% cutoff - suggests a large cd137 lo population
        maxDens <- density(getflowFrame(temp)@exprs[, 15]) 
        if(max(maxDens$x)>2){
          cd137lo.peak <- max(maxDens$y[1:(min(which(maxDens$x>2))-1)])
          cd137hi.peak <- max(maxDens$y[min(which(maxDens$x>2)):length(maxDens$y)])
          cd137lo.peak.idx <- which(maxDens$y == cd137lo.peak)
          cd137hi.peak.idx <- which(maxDens$y == cd137hi.peak)
          cd137.gate <- 0.5*(maxDens$x[cd137lo.peak.idx] + maxDens$x[cd137hi.peak.idx]) 
        }else{
          cd137.gate <- max(maxDens$x)
        }
      }
    }
    names(Ly6c.gate2) <- names(cd137.gate) <- identifier(not.mono.Ly6c.hi)
    
    pDC.flowD <- flowDensity(not.mono.Ly6c.hi, channels = c(10, 15), position = c(T, T), gates = c(Ly6c.gate2,  cd137.gate))
    pDC <- getflowFrame(pDC.flowD)
    
    gthres[i,16] <- Ly6c.gate2
    gthres[i,17] <- cd137.gate
    all.props[i,14] <- pDC.flowD@proportion
    all.events[i,14] <- pDC.flowD@cell.count
    
    # Take complement of the pDC data to obtain NOT(pDC)
    not.pDC <- not.mono.Ly6c.hi
    if(pDC.flowD@cell.count > 0){ # this ensures that the code does not fail if there is no pDC population
      not.pDC@exprs <- not.pDC@exprs[-pDC.flowD@index, ]
    }
    
    all.props[i,15] <-  nrow(not.pDC)/nrow(not.mono.Ly6c.hi)*100
    all.events[i,15] <- nrow(not.pDC)

    ############################################################################################    
    ## Gating NOT(pDC). Plotting Ly6C_CD317 to obtain cDC and Misc
    
    mhcII.gate <- deGate(not.pDC, channel = c(7))
    cDC.flowD.temp <- flowDensity(not.pDC, channels = c(16, 7), position = c(NA, T), gates = c(NA,  mhcII.gate))
    temp <- deGate(getflowFrame(cDC.flowD.temp), channel = c(16), all.cut = T)
    cd11c.gate <- temp[which.min(abs(temp-1.7))]   # pick gate closest to expected gate position
    names(mhcII.gate) <- names(cd11c.gate) <- identifier(not.pDC)
    
    cDC.flowD <- flowDensity(not.pDC, channels = c(16, 7), position = c(T, T), gates = c(cd11c.gate,  mhcII.gate))
    cDC <- getflowFrame(cDC.flowD)
    
    gthres[i,18] <- mhcII.gate
    gthres[i,19] <- cd11c.gate
    all.props[i,16] <- cDC.flowD@proportion
    all.events[i,16] <- cDC.flowD@cell.count
    
    # Take complement of the cDC data to obtain Misc
    misc <- not.pDC
    misc@exprs <- misc@exprs[-cDC.flowD@index, ]
    
    all.props[i,17] <- nrow(misc)/nrow(not.pDC)*100
    all.events[i,17] <- nrow(misc)
    
    ############################################################################################
    
    ## Gating cDC. Plotting CD11b_CD86 to obtain CD8A Type DC and CD11B+ CD86Lo
    # CD8A Type DC
    theta0 = -pi*45/180
    cDC.s <- cDC
    cDC.s <- rotate.data(cDC, c(13, 18), theta = -theta0)$data
    R <- matrix(c(cos(theta0), sin(theta0), -sin(theta0), cos(theta0)), 2 , 2)
    
    cDC.gate.slant <- deGate(cDC.s, channel = c(18), percentile = NA)
    
    # In case deGate cannot find the gate, find the peaks in the density distribution
    data.new <- exprs(cDC.s)[,c(13,18)]
    z <- kde2d(data.new[, 1], data.new[, 2])
    
    # Find the cd11b+cd86lo peak
    ylo <- which.min(abs(z$y + 1.0))  
    peak.lo <- which(z$z == max(z$z[,1:ylo]), arr.ind = T)
    y.peak.lo <- z$y[peak.lo[2]]
    # Find the CD8A peak
    ylo2 <- which.min(abs(z$y - (y.peak.lo + 0.7)))
    yhi <- which.min(abs(z$y - 0.25)) 
    peak.hi <- which(z$z == max(z$z[,ylo2:yhi]), arr.ind = T)
    y.peak.hi <- z$y[peak.hi[2]]
    if(y.peak.hi < -0.69){ # if rel. small or no cd8A peak
      y.peak.hi <- -0.4
    }
    if((y.peak.hi-y.peak.lo) < 0.15){
      if(z$z[peak.hi[1],peak.hi[2]] > z$z[peak.lo[1],peak.lo[2]]){ #if rel. small or no  cd11b+cd86lo peak
        y.peak.lo <- y.peak.hi - 0.8
      }
    }
    
    # pick the gating closest to the expected location (0.5 above the cd11b+cd86lo peak)
    cDC.gate.slant <- c(cDC.gate.slant,(y.peak.lo+y.peak.hi)/2)
    cDC.gate.slant <- cDC.gate.slant[which.min(abs(cDC.gate.slant - (y.peak.lo+0.56)))]
    
    cd8a.typedc.Tempb <- flowDensity(cDC.s, channels=c(13, 18),position=c(NA, T), gates = c(NA, cDC.gate.slant))
    cd8a.typedc.Tempb.ind <- cd8a.typedc.Tempb@index
    cd8a.typedc <- getflowFrame(cd8a.typedc.Tempb)
    cd8a.typedc@exprs[, c(13, 18)] <- t(t(R) %*% t(cd8a.typedc@exprs[, c(13, 18)]))
    
    cd8a.typedc.Tempb@filter <- rotate.data(cd8a.typedc.Tempb@filter, c(13, 18), theta = theta0)$data
    cd8a.typedc.Tempb@flow.frame <- rotate.data(getflowFrame(cd8a.typedc.Tempb), c(13, 18), theta = theta0)$data
    cd8a.typedc.Tempb@proportion <- (cd8a.typedc.Tempb@cell.count/nrow(cDC))*100
    cd8a.typedc <- getflowFrame(cd8a.typedc.Tempb)
    
    gthres[i,20] <- cDC.gate.slant
    all.props[i,18] <- cd8a.typedc.Tempb@proportion
    all.events[i,18] <- cd8a.typedc.Tempb@cell.count

    # Take complement of the  CD8A Type DC to obtain CD11B+ CD86Lo cells
    cd11b.plus.cd86.lo <- cDC
    cd11b.plus.cd86.lo@exprs <- cd11b.plus.cd86.lo@exprs[-cd8a.typedc.Tempb@index, ]
    
    all.props[i,19] <- nrow(cd11b.plus.cd86.lo)/nrow(cDC)*100
    all.events[i,19] <- nrow(cd11b.plus.cd86.lo)
    
    #     cd11b.plus.cd86.lo2.flowD <-flowDensity(cDC, channels=c(13, 18), position=c(T, F),  ellip.gate=T)
    #     cd11b.plus.cd86.lo2 <- getflowFrame(cd11b.plus.cd86.lo2.flowD)
    
    ############################################################################################
    ## Gating CD8A Type DC. Plotting CD103_CD11b to obtain CD103+ DC
    
    # Find the peak of the not(CD103+DC) population using the 2d density 
    maxDens <- density(cd8a.typedc@exprs[,17])
    peak.idx <- which.max(maxDens$y[1:max(which(maxDens$x < 1.5))]) 

    # If the not(CD103+DC) peak is part of a larger flat peak, shift to the right if necessary
    if(maxDens$x[peak.idx] < 0.25){
      if(min(maxDens$y[peak.idx:max(which(maxDens$x < 0.8))]) > (maxDens$y[peak.idx]*.96)){
        peak.idx <- max(which(maxDens$x < 0.8))
      }
    }
    temp <- deGate(cd8a.typedc, channel = c(17), all.cut = T)
    
    # now assume the not(cd103+ DC) peak is circular 
    # find the distance from the peak to the bottom in the cd11b direction
    # and use this as the distance of the gate from the peak in the cd103 direction
    cd11b.gate.low <- deGate(cd8a.typedc, channel = c(13), use.upper = T, upper = F) + 0.1
    cd11b.gate.low2 <- deGate(cd8a.typedc, channel = c(13), use.percentile = T, percentile = 0.025)
    maxDens.cd11b <- density(cd8a.typedc@exprs[,13])
    peak.cd11b <- maxDens.cd11b$x[which.max(maxDens.cd11b$y)]
    temp <- c(temp, maxDens$x[peak.idx] + (peak.cd11b - cd11b.gate.low), maxDens$x[peak.idx] + (peak.cd11b - cd11b.gate.low2))
    cd103.gate <- temp[which.min(abs(temp-(maxDens$x[peak.idx] + 0.75)))]

    # if there is a larger cd103 hi peak and the gate is not close to the minimum between peaks, shift it
    peak2.idx <- which.max(maxDens$y) 
    if(maxDens$x[peak2.idx] > 2.15){
      valley.idx = which.min(maxDens$y[peak.idx:peak2.idx]) + peak.idx - 1
      if(abs(valley.idx - peak.idx) < 5){
        valley <- (maxDens$x[peak.idx] + maxDens$x[peak2.idx])/2
      }else{
        valley <- maxDens$x[valley.idx]
      }
      if(abs(cd103.gate - valley) > 0.1){
        cd103.gate <- valley
      }
    }
    
    # last check for outliers
    if((cd103.gate < (maxDens$x[peak.idx] + 0.35)) || (cd103.gate > 2.08)){
      cd103.gate <- maxDens$x[peak.idx] + 0.75
    }

    names(cd103.gate) <- identifier(cd8a.typedc)
    
    cd103plus.dc.flowD <- flowDensity(cd8a.typedc, channels = c(17, 13), position = c(T, NA), gates = c(cd103.gate, NA))
    cd103plus.dc <- getflowFrame(cd103plus.dc.flowD)
    
    gthres[i,21] <- cd103.gate
    all.props[i,20] <- cd103plus.dc.flowD@proportion
    all.events[i,20] <- cd103plus.dc.flowD@cell.count
    
  },error = function(err) {
      err
    }
  ) # end of tryCatch
  
  # if there is an error in gating do not save gates/figures
  if(inherits(possibleError, "error")){
    errorFileIndex <- c(errorFileIndex, i)
    print(paste0("Error in gating: ", possibleError))
    #next
  }

  ## Saving the Plots
  tryCatch({
    #--------Start Big Png------------
    png ( file = paste("Results/Figures/ScatterPlots/", store.allFCS[i,1], "/", "Total_", store.allFCS[i,2], ".png", sep = "" ), width=2500, height=2500*3/5)
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
    abline(h= lincd3cd109cd161.gate.high, lwd =1)
    #lines(cd45.flowD@filter,lwd =1)
    
    plotDens(cd45, channels = c(14,12), main="CD45", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(cd45sTempb@filter, lwd =1)
    
    # Plotting F4/80_CD11b
    plotDens(lin.neg, channels = c(8,13), main="Lineage neg", cex.lab = 2, cex.axis = 2, cex.main=2)
    lines(linneg.Tempb@filter,lwd =1)
    
    # Plotting Ly6G_SSC-A
    #     plotDens(linneg.macneg, channels = c(9,4), main="Lin-Mac-")
    #     lines(eosinTempb@filter, lwd =1)
    #     segments(Ly6G.gate,par('usr')[3],Ly6G.gate,ssc.a.gate2, lwd =1); segments(par('usr')[1],ssc.a.gate2,Ly6G.gate,ssc.a.gate2, lwd =1)
    #     neu2.gate <- Ly6G.gate
    #     segments(neu2.gate,par('usr')[3],neu2.gate,ssc.a.gate3, lwd =1); segments(neu2.gate,ssc.a.gate3,par('usr')[2],ssc.a.gate3, lwd =1)
    
    plotDens(linneg.macneg, channels = c(9,4), main="Lin-Mac-", cex.lab = 2, cex.axis = 2, cex.main=2)
      # lines(neu2Tempb@filter)
    lines(Ly6GnegTempb@filter)
    lines(eosinTempb@filter)
    
    #plotDens(linneg.macneg, channels = c(9,4), main="Lin-Mac-")
    #lines(Ly6Gneg.flowD@filter,lwd =1)
    #lines(eosinTempb@filter, lwd =1)
    #lines(neu2.flowD@filter, lwd =1)
    
    # Plotting Ly6C_CD11b
    plotDens(Ly6Gneg, channels = c(10,13), main="Ly6G-", cex.lab = 2, cex.axis = 2, cex.main=2)
    #lines(mono.Ly6c.hi.flowD@filter,lwd =1)
    segments(Ly6c.gate, cd11b.gate, Ly6c.gate, par('usr')[4], lwd =1)
    segments(Ly6c.gate, cd11b.gate, par('usr')[2], cd11b.gate, lwd =1)
    
    # Plotting Ly6C_CD317 
    plotDens(not.mono.Ly6c.hi, channels = c(10,15), main="NOT(Monocytes Ly6c hi)", cex.lab = 2, cex.axis = 2, cex.main=2)
    #lines(pDC.flowD@filter,lwd =1)
    segments(Ly6c.gate2, cd137.gate, Ly6c.gate2, par('usr')[4], lwd =1)
    segments(Ly6c.gate2, cd137.gate, par('usr')[2], cd137.gate, lwd =1)
    
    # Plotting Ly6C_CD317
    plotDens(not.pDC, channels = c(16, 7), main="NOT pDC", cex.lab = 2, cex.axis = 2, cex.main=2)
    #lines(cDC.flowD@filter,lwd =1)
    segments(cd11c.gate, mhcII.gate, cd11c.gate, par('usr')[4], lwd =1)
    segments(cd11c.gate, mhcII.gate, par('usr')[2], mhcII.gate, lwd =1)
    
#     # Plotting CD11b_CD86
#     plotDens(cDC, channels = c(13, 18), main="cDC", cex.lab = 2, cex.axis = 2, cex.main=2)
#     lines(cd8a.typedc.Tempb@filter, lwd =1)
    
#     plotDens(cDC.s, channels = c(13,18))
#     abline(h=y.peak.lo, lty =2)
#     abline(h=y.peak.hi, lty = 2)
#     abline(h=cDC.gate.slant)
    
    # Plotting CD11b_CD86 with contour
    plotDens(cDC, channels = c(13, 18), main="cDC", cex.lab = 2, cex.axis = 2, cex.main=2)
    data.new <- exprs(cDC)[,c(13,18)]
    z <- kde2d(data.new[, 1], data.new[, 2])
    contour(z, drawlabels = FALSE, add = TRUE, lty = 2)
    lines(cd8a.typedc.Tempb@filter, lwd =1)
    
#     # Plotting CD103_CD11b
#     plotDens(cd8a.typedc, channels = c(17, 13), main="CD8A Type DC", cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(-1,4), ylim = c(-1,4))
#     abline(v=cd103.gate, lwd =1)
    
    # Plotting CD103_CD11b with contour
    plotDens(cd8a.typedc, channels = c(17, 13), main="CD8A Type DC", cex.lab = 2, cex.axis = 2, cex.main=2, xlim = c(-1,4), ylim = c(-1,4))
    data.new <- exprs(cd8a.typedc)[,c(17,13)]
    z <- kde2d(data.new[, 1], data.new[, 2])
    contour(z, drawlabels = FALSE, add = TRUE, lty = 2)
    abline(v=cd103.gate, lwd =1)
    
    #     plot(1, type="n", axes=F, xlab="", ylab="") # plot a blank
    
    #graphics.off()
    dev.off()
    #par(mfrow=c(1,1))
  },error = function(err) {print(paste0("Error in plotting ", err))})
  
  #--------End Big Png------------
  cat("Time is: ",TimeOutput(start2),"\n",sep="")

} # end of for-loop

################################################
# Write cell count and propoortion data to files
################################################
if(write.data){
  
  colnames(gthres) <- c('singlets.gate.h', 'singlets.gate.l', 'live.gate', 'fsc.a.gate.high', 'fsc.a.gate.low', 'ssc.a.gate', 'cd45.gate.low', 
                        'cd45.gate.high', 'lincd3cd109cd161.gate.high', 'lin.gate.slant', 'linneg.gate.slant', 'linneg.macneg.gate.slanta', 
                        'linneg.macneg.gate.slantb', 'Ly6c.gate', 'cd11b.gate', 'Ly6c.gate2', 'cd317.gate', 'mhcII.gate', 'cd11c.gate', 
                        'cDC.gate.slant', 'cd103.gate', 'linneg.macneg.gate.slanta2')
  rownames(gthres) <- 1:length(gthres[,1])
  save(gthres, file = paste("Results/ThreshTable.Rdata",sep="") )  
  save(errorFileIndex, file = paste("Results/store.Gate.Error.Rdata",sep=""))
  
  colnames(all.props) <- c("All Events-%Parent", "Singlets-%Parent", "Live-%Parent", "Lymphocytes-%Parent", "CD45-%Parent", "Lineage neg-%Parent", 
                           "Lin-Mac--%Parent", "RP macrophages", "Ly6G--%Parent", "Neutrophils 2-%Parent", "Eosinophils 3-%Parent", 
                           "Monocytes Ly6c hi-%Parent","NOT(Monocytes Ly6c hi)-%Parent", "pDC-%Parent", "NOT pDC-%Parent", "cDC-%Parent", 
                           "Misc-%Parent", "CD8A Type DC-%Parent","CD11B+ CD86Lo-%Parent", "CD103+ DC-%Parent")
  rownames(all.props) <- 1:length(all.props[,1])
  
  colnames(all.events) <- c("All Events-#Events", "Singlets-#Events", "Live-#Events", "Lymphocytes-#Events", "CD45-#Events", "Lineage neg-#Events", 
                            "Lin-Mac--#Events", "RP macrophages", "Ly6G--#Events", "Neutrophils 2-#Events", "Eosinophils 3-#Events", 
                            "Monocytes Ly6c hi-#Events","NOT(Monocytes Ly6c hi)-#Events", "pDC-#Events", "NOT pDC-#Events", "cDC-%#Events", 
                            "Misc-#Events", "CD8A Type DC-#Events","CD11B+ CD86Lo-#Events", "CD103+ DC-#Events")
  rownames(all.events) <- 1:length(all.props[,1])
  
  Events_Proportions_Table <- NULL
  for(i in 1:ncol(all.props)){
    Events_Proportions_Table <- cbind(Events_Proportions_Table,all.events[,i],all.props[,i])
  }
  colnames(Events_Proportions_Table) <- c(rbind(colnames(all.events), colnames(all.props)))
  rownames(Events_Proportions_Table) <- 1:length(Events_Proportions_Table[,1])

  all.props <- cbind(store.allFCS[,1:2], all.props)
  all.events <- cbind(store.allFCS[,1:3], all.events)
  
  date.time <- strftime(Sys.time(), "_20%y%m%d_%H%M.csv")
  write.csv(Events_Proportions_Table, file =  paste("Results/Events_Proportions_Table", date.time, sep=""))
  write.csv(all.props, file =  paste("Results/allProportions_Table", date.time, sep=""))
  write.csv(all.events, file =  paste("Results/allEvents_Table", date.time, sep=""))
}

cat("Total time is: ",TimeOutput(start),"\n",sep="")

