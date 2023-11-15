    remove(list=ls())
    setwd("/code/Projects/3i/Panel_BM-cell/")

    # load locally
    source("Codes/IMPC_Functions.R")
    source("Codes/timeVsFluorescence.R")
    source("Codes/SangerCleanFunc.R")
    
    results.dir <- "/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results"

    marker.names <- c("Live|I515*|Syto*|DAPI","IgD","CD43","CD24","GR1","IgM","CD11b","CD138","CD3","BP1","B220", "CD45")
    
    res.clean <- SangerCleanFunc(results.dir, marker.names, overwrite.channel.names = NULL)
    
    load(paste0(results.dir, "/store.allFCS.Rdata"))
    ## List the indices of the flagged FCS files based on TvsF algorithm
    flagged.FCS.index <- which(res.clean[c('Has the file passed'),] == "F")
    
    flagged.FCS <- store.allFCS[flagged.FCS.index,]
    rownames(flagged.FCS) <- flagged.FCS.index
    flagged.FCS.TvsF <- res.clean[,flagged.FCS.index]
    
    save (flagged.FCS, file =  paste0(results.dir,"/flagged.FCS.Rdata") )
    save (flagged.FCS.TvsF, file =  paste0(results.dir,"/flagged.FCS.TvsF.Rdata") )
    
    ## For sending a spreadsheet containing list of flagged FCS files to center
    write.csv(flagged.FCS, file = paste0(results.dir, "/flagged.FCS.BoneMarrow.csv"), row.names = FALSE)
    