# Developed by Albina Rahim
# Date: April 28, 2016

# This code creates lgl657.Rdata which is just the estimateLogicle transformation of GlobalFrame657.
remove(list=ls())
setwd("/code/Projects/3i/Panel_BM-cell/")


load (file =  paste("Results/GlobalFrame657.Rdata",sep="") )
verbose_debris <- T
verbose_margin <- T

library("flowCore")
library("flowDensity")
source("Codes/3iTcellfunctions-BM.R")

    ## Removing Margin events------------------------------------------------------------------------------------------
    scat.chans <- c(grep (colnames(g),pattern = "FSC*"),grep (colnames(g),pattern = "SSC*"))
    ## Removing margin events in Scatter channels
    g <- removeMargins(g, chans=scat.chans, verbose= verbose_margin)
    #Removing negative values in scatter channels
    g <- removeMargins(g, chans=scat.chans, debris=T, neg=T, verbose= verbose_debris) 
    

    # Finding markers for BM cells---------------------
     markers<- c("Live|I515*|Syto*|DAPI","IgD","CD43","CD24","GR1","IgM","CD11b","CD138","CD3","BP1","B220", "CD45")
    
    channels.ind <- Find.markers(frame=g,marker.list=markers)
    names(channels.ind)[grep(names(channels.ind),pattern = "Live*")] <- "Live"
    channels.ind <- sort(channels.ind)
    
    ## Tranformation
    lgl657<-estimateLogicle(g, channels= colnames(g)[channels.ind])
    
    
    
    save(channels.ind, file = paste("Results/channels.ind.Rdata", sep = ""))
    save( lgl657,    file =  paste("Results/lgl657.Rdata",sep=""))