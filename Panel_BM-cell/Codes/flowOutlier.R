#flowOutlier function

flowOutlier1 <- function(thresholds, channel, rank_upper_T, rank_upper_F, rank_tinypeak, rank_adjust, rank_matrix, densities, wp=5, po = 0.95, SD=1, verbose=FALSE){

    # Rank_Total  <- rank(rank_upper_T) + rank(rank_upper_F) + rank(1-rank_tinypeak) + rank(1-rank_adjust) # Find the files that are not well behaved.
    Rank_Total  <- rank(1-rank_tinypeak) + rank(10-rank_adjust) # Find the files that are not well behaved.

    cbind(rank(rank_upper_T), rank(rank_upper_F), rank(1-rank_tinypeak), rank(1-rank_adjust), Rank_Total)

    # plot(density(Rank_Total, adjust = 1))
    # plot(density(Rank_Total, adjust = 0.5))
    # plot(density(Rank_Total, adjust = 0.3))

    BadIndices  <- which(Rank_Total >= sort(Rank_Total)[round(length(densities)*po)])
    GoodIndices <- which(Rank_Total <  sort(Rank_Total)[round(length(densities)*po)])

    # print(paste0("flowOutlier has chosen to update the following files: ", paste(unlist(BadIndices), collapse=", ") ))

        thresholdsShortened <- thresholds[GoodIndices]

    SDchange <- NULL
    if ( verbose == TRUE )
        print("Similarity Outlier Detection")
    for (j in 1:length(BadIndices)) {  # third version
        similarGateValues <- c(rank_matrix[BadIndices[j],GoodIndices]+ rank_matrix[GoodIndices,BadIndices[j]])
        NewValue <- sum(thresholdsShortened*(1/similarGateValues^wp)) / sum(1/similarGateValues^wp) # weighted mean
        weightedSD <- sqrt( sum((thresholdsShortened-NewValue)^2*(1/similarGateValues^wp)) / (sum(1/similarGateValues^wp))) # weighted SD
        # print(j)
        # print(thresholds[BadIndices[j]])
        # print(NewValue)
        # print(SD)
        # print(weightedSD)
        if ( (thresholds[BadIndices[j]] <= (NewValue - SD*weightedSD )) || (thresholds[BadIndices[j]] >= (NewValue + SD*weightedSD )) ) { # if outside range
            if ( verbose == TRUE )
                print(paste0("Changed flowFrame ", BadIndices[j], ", channel ", channel,": ",round(thresholds[BadIndices[j]], digits=2)," to ", round(NewValue, digits=2)))
            thresholds[BadIndices[j]] <- NewValue # update threshold
        }
    }
    return( thresholds )
}
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
flowOutlier3 <- function(thresholds, channel, rank_matrix, wp=5, SD=3.25, verbose=FALSE){

    SDchange <- NULL
    if ( verbose == TRUE )
        print("Standard Outlier Detection")
    for (j in 1:length(thresholds)) {
        similarGateValues <- c(rank_matrix[j,]+ rank_matrix[,j])
        NewValue <- sum(thresholds[-j]*(1/similarGateValues[-j]^wp))/sum(1/similarGateValues[-j]^wp) # weighted mean
        weightedSD <- sqrt( sum((thresholds[-j]-NewValue)^2*(1/similarGateValues[-j]^wp)) / (sum(1/similarGateValues[-j]^wp))) # weighted SD
        if ( (thresholds[j] <= (NewValue - SD*weightedSD )) || (thresholds[j] >= (NewValue + SD*weightedSD )) ) { # if outside range
            if ( verbose == TRUE )
                print(paste0("Changed flowFrame ", j, ", channel ", channel,": ",round(thresholds[j], digits=2)," to ", round(NewValue, digits=2),
                         "  [", round((NewValue - SD*weightedSD ), digits=2),", ", round((NewValue + SD*weightedSD ), digits=2),"]" ))
            thresholds[j] <- NewValue # update threshold
        }
    }
    return( thresholds )
}
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# Create Rank Matrix
flowOutlierRankMatrix <- function(densities){

    rank_matrix <- matrix (0,length(densities),length(densities))
    for ( v1 in 1:(length(densities)-1)){
      for ( v2 in (v1+1):length(densities)){
          y2mod <- approx(densities[[v2]]$x, densities[[v2]]$y, densities[[v1]]$x) # shift/change x and y values to have the same x values as other function

          rank_matrix[v1,v2] <- (sum(abs(y2mod$y - densities[[v1]]$y),na.rm = T)/    # ************** Rank of how similar each function is to each other ************
                                mean(c(sum(y2mod$y, na.rm = T), sum(densities[[v1]]$y, na.rm = T) )) ) # One sided Reimann sum. dx's cancel
      }
    }
    return( rank_matrix )
}

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# Rank tinypeak.removal
flowOutlierRankTinyPeak <- function(fSet, channel){

    rank_tinypeak <- NULL
    for ( j in 1:length(fSet)){
        highP <- 1
        lowP <- 0
        perc95 <- tryCatch(deGate(fSet[[j]], channel, use.percentile=T, percentile = 0.95, verbose = F), error=function(e) print("error_found"))
        if(perc95 == "error_found"){
            lowP <- 0
        } else{
            repeat{ # newton iteration
                newP <- round((highP+lowP)/2, digits=2)
                tempDeG <- deGate(fSet[[j]], channel, tinypeak.removal=newP, percentile = 0.95, verbose = F)
                if(tempDeG == perc95){
                    highP <- newP
                } else {
                    lowP <- newP
                }
                if( round((highP - lowP), digits=2) <= 0.01){
                    break;
                }
            }
        }
        rank_tinypeak[j] <- lowP
    }
    return( rank_tinypeak )
}
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# Rank adjust.dens
flowOutlierRankAdjust <- function(fSet, channel){

    rank_adjust <- NULL
    for ( j in 1:length(fSet)){
        highAdj <- 10
        lowAdj <- 1
        perc95 <- tryCatch(deGate(fSet[[j]], channel, use.percentile=T, percentile = 0.95, verbose = F), error=function(e) print("error_found"))
        if(perc95 == "error_found"){
            lowAdj <- 0
        } else{
            repeat{ # newton iteration
                newAdj <- round((highAdj+lowAdj)/2, digits=2)
                tempDeG <- deGate(fSet[[j]], channel, adjust.dens = newAdj, percentile = 0.95, verbose = F)
                if(tempDeG == perc95){
                    highAdj <- newAdj
                } else {
                    lowAdj <- newAdj
                }
                if( round((highAdj - lowAdj), digits=2) <= 0.01){
                    break;
                }
            }
        }
        rank_adjust[j] <- lowAdj
    }
    return( rank_adjust )
}
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# Rank alpha upper T
flowOutlierRankAlphaUpper <- function(fSet, channel, upper.para=T){

    rank_upper <- NULL
    for ( j in 1:length(fSet)){
        tempUpper <- NULL
        for ( k in 1:40){
            tempUpper[k] <- deGate(fSet[[j]], channel, use.upper = T, upper = upper.para, tinypeak.removal=0.99, alpha = (k)/100, verbose = F)
        }
        rank_upper[j] <- sd(tempUpper) * mean(tempUpper) #***********  Rank for the upper  ****************
    }
    return( rank_upper )
}
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# Rank_deGate
flowOutlierDensities <- function(fSet, channel){
    densities <- list()
    for ( j in 1:length(fSet)){
        densities[[j]] <- density(fSet[[j]]@exprs[,channel]) # x and y values of density, used by flowOutlier1
    }
    return( densities )
}
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# Prints out the time since start_time. Used for optimizing code.
TimeOutput <- function(start_time) {
    start_time <- as.POSIXct(start_time)
    dt <- difftime(Sys.time(), start_time, units="secs")
    # Since you only want the H:M:S, we can ignore the date...
    # but you have to be careful about time-zone issues
    format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
    # thresholds <- thresholds.basic
    # fSet <- fs
    # channel=7
    # wp1=5
    # wp3=5
    # SD1=1
    # SD3=3.25
    # po=0.95
    # verbose=FALSE

flowOutlier <- function(thresholds, fSet, channel, wp1=5, wp3=5, SD1=1, SD3=3.25, po=0.95, verbose=FALSE){

    start <- Sys.time()
    densities     <- flowOutlierDensities(fSet, channel)
    rank_upper_T  <- flowOutlierRankAlphaUpper(fSet, channel, upper.para=T)
    rank_upper_F  <- flowOutlierRankAlphaUpper(fSet, channel, upper.para=F)
    rank_tinypeak <- flowOutlierRankTinyPeak(fSet, channel)
    rank_adjust   <- flowOutlierRankAdjust(fSet, channel)
    rank_matrix   <- flowOutlierRankMatrix(densities)
    thresholds    <- flowOutlier1(thresholds, channel, rank_upper_T, rank_upper_F, rank_tinypeak, rank_adjust, rank_matrix, densities, wp=wp1, po=po, SD=SD1, verbose=verbose)
    thresholds    <- flowOutlier3(thresholds, channel, rank_matrix, wp=wp3, SD=SD3, verbose=verbose)
    if (verbose == T){ cat("flowOutlier time: ", TimeOutput(start),"\n",sep="") }

    return( thresholds )
}