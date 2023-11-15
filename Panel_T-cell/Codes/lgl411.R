# Developed by Albina Rahim
# Date: March 16, 2016

# This code creates lgl657.Rdata which is just the estimateLogicle transformation of GlobalFrame411.
remove(list=ls())
setwd("/code/Projects/3i/Panel_T-cell_MLN/")

load (file =  paste("Results/GlobalFrame657.Rdata",sep="") )
verbose_debris <- T
verbose_margin <- T

library("flowCore")
library("flowDensity")
source("Codes/3iTcellfunctions.R")

    ## Removing Margin events------------------------------------------------------------------------------------------
    scat.chans <- c(grep (colnames(g),pattern = "FSC*"),grep (colnames(g),pattern = "SSC*"))
    ## Removing margin events in Scatter channels
    g <- removeMargins(g, chans=scat.chans, verbose= verbose_margin)
    #Removing negative values in scatter channels
    g <- removeMargins(g, chans=scat.chans, debris=T, neg=T, verbose= verbose_debris) 
    

    # Finding markers (added the marker CD45 because dataset from 3i have this marker.)---------------------
    markers<- c("Live|I515*|Syto*|DAPI","CD5","CD161","CD62L","CD44","CD8","CD4","CD25","GITR","CD24","gdTCR|TCRd|Syto*","KLRG", "CD45")
    channels.ind <- Find.markers(frame=g,marker.list=markers)
    names(channels.ind)[grep(names(channels.ind),pattern = "Live*")] <- "Live"
    names(channels.ind)[grep(names(channels.ind),pattern = "TCRd*")] <- "TCRd"
    channels.ind <- sort(channels.ind)
    
    ## Tranformation
    lgl657<-estimateLogicle(g, channels= colnames(g)[channels.ind])
    #lgl657.temp<-estimateLogicle(g, channels= colnames(g)[7:18])
    
    
    save(channels.ind, file = paste("Results/channels.ind.Rdata", sep = ""))
    save( lgl657,    file =  paste("Results/lgl657.Rdata",sep=""))