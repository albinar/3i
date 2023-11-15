# Originally written by Nima Aghaeepour (UFO package) -------------------
# Re-edited by Albina Rahim (July 11, 2016)

remove(list=ls())

library("stringr")

# This code creates an .Rdata file for each .fcs file after running flowType

##Logicle transformation -> flowType on 10^4 random FCS files (increase this if you want) -> save the results

logicleTrans <- function(Frame, Markers){
      lgl <- estimateLogicle(Frame,channels= Frame@parameters@data$name[7:18])
      Frame <- transform(Frame, lgl)
  return(Frame);
}

setwd("/code/Projects/3i/Panel_BM-cell/")

load("Results/store.allFCS.Rdata")
## Code for removing duplicate FCS files based on their barcodes
Barcodes <- str_extract(store.allFCS[,2],"L[0-9]+")
store.allFCS.unique <- cbind(store.allFCS, Barcodes)
duplicate.index <- which(duplicated(Barcodes)==TRUE)
duplicate.FCS <- store.allFCS[duplicate.index, 1:2]
store.allFCS.unique <- store.allFCS.unique[!duplicated(store.allFCS.unique[,c('Barcodes')]),]
meta_data <- store.allFCS.unique

suppressWarnings ( dir.create ("Results/UFO/"))
suppressWarnings ( dir.create ("Results/UFO_WT/"))

# ## The following code can be used if we want to check batch effect on Wild types only:
# f.path1 <- dir("FCS_Groups/+_+", full.names = T, recursive = T)
# f.path2 <- dir("FCS_Groups/+_Y", full.names = T, recursive = T)
# f.path  <- c(f.path1, f.path2)

## The following code is used for determining batch effect on all the files (WTs + KOs)
# path to the FCS files
pathFCS <- paste("/code/Projects/3i/Panel_BM-cell/FCS_Groups")
# Read all folders and files in current path folder and makes a list of all of their paths
f.path <- dir(pathFCS, full.names = T, recursive = T)

# # this line is for finding batch effect for WTs only
# meta_data <- meta_data[match(File_Names, meta_data[,2]),]

File_Names <- meta_data[,"FCS file"]
WT_dates <- meta_data[,"Assay Date"]
Genotype <- meta_data[,"Genotype"]
gender <- meta_data[,"Gender"]
gender[which(gender == "Female")] <- 0
gender[which(gender == "Male")] <- 1
gender <- as.numeric((gender))

WT_dates2 <- strsplit(as.character(WT_dates), split="-")
WT_dates3 <- NULL
months <- c(31,28,31,30,31,30,31,31,30,31,30,31)
for ( q1 in 1:length(WT_dates2) ) {
        temp <- WT_dates2[[q1]]
        WT_dates3[q1] <- (as.numeric(temp[3])-14)*365 + sum(months[ 1:match(temp[2], month.abb)] ) - months[match(temp[2], month.abb)] + as.numeric(temp[1])
}

meta_data <- cbind(Genotype,File_Names, WT_dates, WT_dates3, rep(1,100), gender)
meta_data[which(meta_data[,1] == "+_+"), 5] <- 0
meta_data[which(meta_data[,1] == "+_Y"), 5] <- 0

colnames(meta_data) <- c("Genotype", "Mouse/File","Date", "Day", "Group", "Gender")


files = meta_data[,2]
datapath='/Results/UFO'
#CSN <- 1:5
#CSN <- 1:nrow(meta_dataWT)
tempfunc <- function(CSN){
    print(CSN)
    if (file.exists(paste0("Results/UFO/",File_Names[CSN],".Rdata")))
        return(1);
    library(flowCore)
    f = read.FCS(f.path[CSN])
    f <- compensate(f, f@description$SPILL)
    f <- logicleTrans(f, c(7:18))
    subsample = sample(seq(nrow(exprs(f))), 10^3) ## sub-sample of 10^3 cells
    x <- exprs(f)[subsample,]
    marker.names <- f@parameters@data$name
    #prop.markers=c(7:18)
    #MFI.markers=NULL
    f <- new('flowFrame',exprs=x)
    library(flowType)
    FT <- flowType(Frame=f, PropMarkers= c(7:18), MFIMarkers=NULL, Methods='kmeans',
                   MarkerNames=NULL,MemLimit=4,verbose=TRUE, MaxMarkersPerPop=8) ## The MaxMarkersPerPop has been changed from 3 to 8 (decided after Ania's email on August 06, 2016)
    save(list='FT', file=paste0("Results/UFO/",File_Names[CSN],".Rdata"))
}

library(snow)
library(Rmpi)
cl <- makeCluster(spec=8, type = "SOCK")
clusterExport(cl,list=c('datapath','files', 'logicleTrans', 'File_Names', 'f.path'))
# clusterApplyLB(cl, 1:10, tempfunc);
clusterApplyLB(cl, 1:dim(meta_data)[1], tempfunc);
stopCluster(cl)

save(WT_dates, file = paste("Results/WT_dates.Rdata", sep = ""))
save(WT_dates3, file = paste("Results/WT_dates3.Rdata", sep = ""))
save(File_Names, file = paste("Results/File_Names.Rdata", sep = ""))
save(meta_data, file = paste("Results/meta_data.Rdata", sep = ""))

########################################################################################

# UFO analysis on the WT samples only (This part of the script is optional)

## The following code can be used if we want to check batch effect on Wild types only:
f.path1 <- dir("FCS_Groups/+_+", full.names = T, recursive = T)
f.path2 <- dir("FCS_Groups/+_Y", full.names = T, recursive = T)
f.path  <- c(f.path1, f.path2)

meta_dataWT <- meta_data[which(meta_data[,5] == 0),]
File_NamesWT <- meta_dataWT[,"Mouse/File"]

files = meta_dataWT[,2]
datapath='/Results/UFO_WT'
#CSN <- 1
#CSN <- 1:5
#CSN <- 1:nrow(meta_dataWT)
tempfuncWT <- function(CSN){
  print(CSN)
  if (file.exists(paste0("Results/UFO_WT/",File_NamesWT[CSN],".Rdata")))
    return(1);
  library(flowCore)
  print(paste("Reading", f.path[CSN]))
  f = read.FCS(f.path[CSN])
  f <- compensate(f, f@description$SPILL)
  f <- logicleTrans(f, c(7:18))
  subsample = sample(seq(nrow(exprs(f))), 10^3) ## sub-sample of 10^3 cells
  x <- exprs(f)[subsample,]
  marker.names <- f@parameters@data$name
  #prop.markers=c(7:18)
  #MFI.markers=NULL
  f <- new('flowFrame',exprs=x)
  library(flowType)
  FT <- flowType(Frame=f, PropMarkers= c(7:18), MFIMarkers=NULL, Methods='kmeans',
                 MarkerNames=NULL,MemLimit=4,verbose=TRUE, MaxMarkersPerPop=8)
  save(list='FT', file=paste0("Results/UFO_WT/",File_NamesWT[CSN],".Rdata"))
}

library(snow)
library(Rmpi)
cl <- makeCluster(spec=8, type = "SOCK")
clusterExport(cl,list=c('datapath','files', 'logicleTrans', 'File_NamesWT', 'f.path'))
# clusterApplyLB(cl, 1:10, tempfunc);
clusterApplyLB(cl, 1:dim(meta_dataWT)[1], tempfuncWT);
stopCluster(cl)

save(File_NamesWT, file = paste("Results/File_NamesWT.Rdata", sep = ""))
save(meta_dataWT, file = paste("Results/meta_dataWT.Rdata", sep = ""))
save(duplicate.FCS, file = paste("Results/duplicate.FCS.Rdata", sep = ""))
save(store.allFCS.unique, file = paste("Results/store.allFCS.unique.Rdata", sep = ""))

