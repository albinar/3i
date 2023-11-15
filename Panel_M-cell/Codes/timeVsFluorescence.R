##################################################################################################
# removes groups of cells that show some statisitcally significant differences from other groups #
##################################################################################################
timeVsFluorescence <- function(f, segment, alpha, directory, name1, name2, Plt){

    f.org <- f
    PlotPlots <- T

    FSC.loc  <- grep("fsc", tolower(f@parameters@data$name));  names( FSC.loc) <- NULL
    SSC.loc  <- grep("ssc", tolower(f@parameters@data$name)); names( SSC.loc) <- NULL
    Time.loc <- which(tolower(f@parameters@data$name) == "time"); names(Time.loc) <- NULL
    OtherChan.loc <- (1:ncol(f))[-c(FSC.loc, SSC.loc, Time.loc)]

    # save the time plots before any removal
    if ( PlotPlots == F && Plt == T){
        for (x in OtherChan.loc){
          png ( file = paste0(directory, "/", name1, "_", name2, "_Before_", x, ".png"))
          plotDens(f, c(Time.loc,x))
          dev.off ( )
        }
    }

    listDeleteTotal <- NULL
    totalNumSeg <- floor(nrow(f.org@exprs)/segment)

    for ( j in OtherChan.loc){
    # create matrix from summary for the arrayOutliers function
        segSummary <- NULL
        for ( k in 1:totalNumSeg){
            if ( k == totalNumSeg)
                temp <- f@exprs[(segment*(k-1)):nrow(f@exprs), c(j)]
            if ( k != totalNumSeg)
                temp <- f@exprs[segment*(k-1)+(1:segment), c(j)]
            segSummary[[k]] <- summary(temp)[2:5]
        #     segSummary[[k]] <- summary(temp)

        }
        allSegSum <- NULL
        oneSegSum <- NULL

        for (m in 1:4){ # 1st Qu / Median / Mean / 3rd Qu
          for (l in 1:length(segSummary)){
            oneSegSum[l] <- segSummary[[l]][m]
          }
        allSegSum <- cbind(allSegSum,oneSegSum)
        }

        cellDelete <- arrayMvout:::ArrayOutliers(data.frame(allSegSum), alpha=alpha )

        if ( is.null(cellDelete$outl)) {next}
        listDelete <- sort(as.numeric(rownames(cellDelete$outl)),decreasing=T)
        listDeleteTotal <- c(listDeleteTotal, listDelete)

    } # end of for-loop

    # noSegmentRem <- F
    if (is.null(listDeleteTotal)){ # if nothing is deleted
        print("None deleted from TvsF segment removal.")
        removed.ind <- NULL
        # noSegmentRem <- T
    } else {
        listDeleteTotal <- sort( unique( listDeleteTotal ) , decreasing = T )
        cat(paste0("Removing sections ", paste0(sort(listDeleteTotal), collapse = ", "), " out of ", totalNumSeg, " sections.\n"))
        ####################################################
        # delete segments that are statistically different #
        ####################################################
        removed.ind <- NULL
        for (n in 1:length(listDeleteTotal)){
            if (listDeleteTotal[n] == totalNumSeg){
                removed.ind <- c(removed.ind, (segment*(listDeleteTotal[n]-1)+1):nrow(f@exprs))
                f@exprs <- f@exprs[-((segment*(listDeleteTotal[n]-1)+1):nrow(f@exprs)),]
            }
            if (listDeleteTotal[n] != totalNumSeg){
                removed.ind <- c(removed.ind, segment*(listDeleteTotal[n]-1)+(1:segment))
                f@exprs <- f@exprs[-(segment*(listDeleteTotal[n]-1)+(1:segment)),]
            }
        }
    }
    ################################################
    # remove sections that have a very low density #
    ################################################
    extraPoints1 <- min(f@exprs[ ,Time.loc])-100
    extraPoints2 <- max(f@exprs[ ,Time.loc])+100
    # add extra points at the beginning and the end to deal with the fact that density curves always have a zero point at the ends
    # then when the ranges for min to close to the right and max to just below can be removed because they are fake anyways.
    # This method allows the algorithm to remove events that have low density at the endpoints.
    dens.f <- density(c(rep(extraPoints1, 1000), f@exprs[ ,Time.loc], rep(extraPoints2, 1000)), n= nrow(f), adjust=0.1)
    low.dens.rem.ind <- which(dens.f$y <= 0.1*max(dens.f$y))

    range.low.dens <- list()
    count <- 1
    low <- 1
    for(b1 in 1:(length(low.dens.rem.ind)-1) ) {
        if( low.dens.rem.ind[b1] != (low.dens.rem.ind[b1+1]-1) ){
            range.low.dens[[count]] <- low.dens.rem.ind[low]:low.dens.rem.ind[b1]
            low <- b1+1
            count <- count + 1
        }
        if( b1 == (length(low.dens.rem.ind)-1) ) # end point
            range.low.dens[[count]] <- low.dens.rem.ind[low]:low.dens.rem.ind[b1+1]
    }

    # does [[1]] or [[N]] contain the first or last point
    if ( length(which(range.low.dens[[length(range.low.dens)]] == nrow(f)) >= 1) )
        range.low.dens[[length(range.low.dens)]] <- NULL
    if ( length(which(range.low.dens[[1]] == 1) >= 1) )
        range.low.dens[[1]] <- NULL

    if (length(range.low.dens) != 0 ){
        range.low.dens <- lapply(1:length(range.low.dens), function(x){ range(range.low.dens[[x]]) }) # change to range
        range.low.dens <- lapply(1:length(range.low.dens), function(x){ c(dens.f$x[range.low.dens[[x]][1]], dens.f$x[range.low.dens[[x]][2]]) } ) # change to time coordinates

        rem.ind <- NULL
        for ( b2 in 1:length(range.low.dens) ){
            rem.ind <- c(rem.ind, intersect(which(f@exprs[,Time.loc] >= range.low.dens[[b2]][1]), which(f@exprs[,Time.loc] <= range.low.dens[[b2]][2])) )
        }
        f.iterm <- f

        if (length(rem.ind) == 0 ){
            print("None deleted from TvsF low dens removal.")
            rem.ind <- NULL
            # if (noSegmentRem == T)
                # return(f)
        } else {
            cat(paste0("Removing low density ranges ", paste0(round(unlist(range.low.dens), digits = 2), collapse = ", "), "\n"))
            f@exprs <- f@exprs[-rem.ind,]
        }
    } else {
        print("None deleted from TvsF low dens removal.")
         rem.ind <- NULL
        # if (noSegmentRem == T)
            # return(f)
    }
    #################################################################################
    # save the time plots with black points indicating which events will be deleted #
    #################################################################################
    alphaChar <-as.character(alpha)
    alphaChar <- sub("[.]", "-", alphaChar)
    if ( PlotPlots == T && Plt == T){
        z <- ceiling(length(OtherChan.loc)/4)
        png ( file = paste(save.dir_TvsF, "/", segment, "_", alphaChar, "_", name1, "_", name2, "_Black.png", sep = "" ), width = 4*400, height = z*400)
        par(mfrow=c(z,4),mar=c(5,5,4,2))
            for (x in OtherChan.loc){
                plotDens(f.org, c(Time.loc,x), cex.main=2, cex.lab=2, cex.axis=2)
                if ( length(removed.ind) != 0 )
                    points(exprs(f.org)[removed.ind,c(Time.loc,x)],col=1,pch=".")
                if ( (length(range.low.dens) != 0) && (length(rem.ind) != 0 ) )
                    points(exprs(f.iterm)[rem.ind,c(Time.loc,x)],col=1,pch=".", cex=1)
            }
        dev.off ( )
    }

    # save the time plots after the removal
    if ( PlotPlots == T && Plt == T){

        z <- ceiling(length(OtherChan.loc)/4)
        png ( file = paste0(save.dir_TvsF, "/", name1, "_", name2, "_After.png"), width = 4*400, height = z*400)
        par(mfrow=c(z,4),mar=c(5,5,4,2))
            for (x in OtherChan.loc){
                plotDens(f, c(Time.loc,x))
            }
        dev.off ( )
    }

    # organize the indices to be removed
    if(  is.null(removed.ind ) &&  is.null(rem.ind))
        to.be.removed <- NULL
    if(  is.null(removed.ind ) && !is.null(rem.ind))
        to.be.removed <- unique( (1:nrow(f.org))[rem.ind])
    if( !is.null(removed.ind ) &&  is.null(rem.ind))
        to.be.removed <- unique(removed.ind)
    if( !is.null(removed.ind ) && !is.null(rem.ind))
        to.be.removed <- unique(removed.ind, (1:nrow(f.org))[-removed.ind][rem.ind])
    if(is.null(to.be.removed)){
        to.be.removed <- NA
    }

return(list(frame=f,inds= to.be.removed))
}
