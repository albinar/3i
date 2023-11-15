# Developed by Albina Rahim
# Date: April 27, 2016
# This code is for BM Panel cells for datasets sent on August 27, 2015 and March 08, 2016


remove(list=ls())
setwd("/code/Projects/3i/Panel_BM-cell/")

library("flowCore")
library("flowBin")
library("flowDensity")
library("flowType")
library("stringr")

source("Codes/3iTcellfunctions-BM.R")

start <- Sys.time()

FCS.Groups <- "FCS_Groups"
dir.create (FCS.Groups)

## path to the files sent on October 11, 2016 for the Bone Marrow Panel
pathMice <- paste("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2016-10-11/BM Labelled")
# # path to the files sent on August 27, 2015
# pathMice <- paste("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2015-08-27/Bone Marrow Labelled")
# Reads all folders and files in current path folder and makes a list of all of their paths
allMice <- dir(pathMice, full.names=T, recursive=T, pattern = "*.fcs") 
store.allMice <- sapply(1:length(allMice), function(x){pathMice})
store.allMice <- cbind(store.allMice, sapply(1:length(allMice), function(x){unlist(strsplit(allMice[x], split = "/"))[length(unlist(strsplit(allMice[x], split = "/")))-1]}))
store.allMice <- cbind(store.allMice, sapply(1:length(allMice), function(x){unlist(strsplit(allMice[x], split = "/"))[length(unlist(strsplit(allMice[x], split = "/")))]}))
store.allMice <- cbind(store.allMice, str_extract(store.allMice[,3],"L[0-9]+"))

  
# path to the files sent on March 08, 2016
pathMice <- paste("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2016-03-08/BM 169+")
# Reads all folders and files in current path folder and makes a list of all of their paths
allMice <- dir(pathMice, full.names=T, recursive=T, pattern = "*.fcs")
store.allMice.temp <- sapply(1:length(allMice), function(x){pathMice})
store.allMice.temp <- cbind(store.allMice.temp, sapply(1:length(allMice), function(x){unlist(strsplit(allMice[x], split = "/"))[length(unlist(strsplit(allMice[x], split = "/")))-1]}))
store.allMice.temp <- cbind(store.allMice.temp, sapply(1:length(allMice), function(x){unlist(strsplit(allMice[x], split = "/"))[length(unlist(strsplit(allMice[x], split = "/")))]}))
store.allMice.temp <- cbind(store.allMice.temp, str_extract(store.allMice.temp[,3],"L[0-9]+"))
  
store.allMice <- rbind(store.allMice, store.allMice.temp)
remove(store.allMice.temp)

index.Remove <- which(is.na(store.allMice[,4]))
if(length(index.Remove != 0)){
    store.allMice<-store.allMice[-index.Remove,]
}

store.allMice <- cbind(store.allMice, NA, NA) # Column for Assay date & Gender
colnames(store.allMice) <- c("Path", "Folder", "FCS file", "Label Barcode", "Assay Date", "Gender")

## Reading the metadata spreadsheet for the October 2016 batch
CSVfile <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2016-10-11/attachments/download_phnBoneImmunophenotyping.csv")
CSVfile <- as.matrix(CSVfile)

# ## Reading the metadata spreadsheet for the August 2015 and March 2016 batches
# CSVfile <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2016-03-08/attachments/Bone Marrow data RB.csv")
# CSVfile <- as.matrix(CSVfile)

Genotype <- CSVfile[,c(11)]
Genotype <- sub("/","_", Genotype)
Mouse_Label <- CSVfile[,c(20)]
Assay_Date <- CSVfile[,c(15)]
Gender <- CSVfile[,c(12)]

uniqueGT <- unique(Genotype)

invisible(sapply(1:length(uniqueGT), function(x){
dir.create ( paste(FCS.Groups,"/", uniqueGT[x], sep="") )}))

count.FCSmoved <- 0
countX <-0
index.Remove <- 0
index.Skip <- 0

for(x in 1:nrow(store.allMice)){
    temp <- grep(store.allMice[x,4], Mouse_Label)
    if(length(temp) == 1){
      print(paste("Move:", store.allMice[x,3]))
      if (!file.exists(paste(FCS.Groups, "/", Genotype[temp], "/", store.allMice[x,3], sep="")) ) {
        print(paste("Copying", store.allMice[x,3]))
        count.FCSmoved <- count.FCSmoved+1
        store.allMice[x,5] <- Assay_Date[temp]
        store.allMice[x,6] <- Gender[temp]
        file.copy(paste(store.allMice[x,1], "/", store.allMice[x,3], sep=""), paste(FCS.Groups, "/", Genotype[temp], "/", store.allMice[x,3], sep=""))
      } else {
        print(paste("Skipping", store.allMice[x,3]))
        store.allMice[x,5] <- Assay_Date[temp]
        store.allMice[x,6] <- Gender[temp]
        index.Skip <- c(index.Skip,x)
      }
    } else {print(paste("Did not move", store.allMice[x,3]))
      index.Remove <- c(index.Remove,x)
      countX <- countX+1}
}
index.Remove <-index.Remove[index.Remove !=0]
if(length(index.Remove != 0)){
    store.allMice<-store.allMice[-index.Remove,]
}
#count.FCSmoved <- length(intersect(Mouse_Label, store.allMice[,3]))
print(paste("Copied and moved", count.FCSmoved, "files"))

suppressWarnings(dir.create ( "Results/"))
save ( store.allMice , file =  paste("Results/store.allMice.Rdata",sep="") )
save(Genotype, file = paste("Results/Genotype.Rdata", sep = ""))
save(uniqueGT, file = paste("Results/uniqueGT.Rdata", sep = ""))

cat("Total time is: ",TimeOutput(start),sep="")



