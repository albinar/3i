# Developed by Albina Rahim
# Date: May 03, 2016
# This code is for Myeloid Panel cells (for both SPLN and MLN organs) for datasets sent on August 27, 2015 and March 08, 2016


remove(list=ls())
#setwd("/code/Projects/3i/Panel_M-cell/")
setwd("/data/Panel_M-cell")

library("flowCore")
library("flowBin")
# library("flowDensity")
# library("flowType")
library("stringr")
library('plyr')

source("Codes/3iTcellfunctions-M.R")

start <- Sys.time()

FCS.Groups <- "FCS_Groups_MLN"
dir.create (FCS.Groups)

# path to the files sent on August 27, 2015
pathMice <- paste("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2015-08-27/Myeloid MLN Labelled")
# Reads all folders and files in current path folder and makes a list of all of their paths
allMice <- dir(pathMice, full.names=T, recursive=T, pattern = "*.fcs") 

store.allMice1 <- sapply(1:length(allMice), function(x){pathMice})
store.allMice1 <- cbind(store.allMice1, sapply(1:length(allMice), function(x){unlist(strsplit(allMice[x], split = "/"))[length(unlist(strsplit(allMice[x], split = "/")))-1]}))
store.allMice1 <- cbind(store.allMice1, sapply(1:length(allMice), function(x){unlist(strsplit(allMice[x], split = "/"))[length(unlist(strsplit(allMice[x], split = "/")))]}))
store.allMice1 <- cbind(store.allMice1, str_extract(store.allMice1[,3],"L[0-9]+"))
  
index.Remove <- which(is.na(store.allMice1[,4]))

if(length(index.Remove != 0)){
    store.allMice1 <- store.allMice1[-index.Remove,]
}

#CSVfile1 <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2015-08-27/attachments/3i_IMPC_Data_Genotypes.csv")
CSVfile1 <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2015-08-27/attachments/download_phnMLNImmuno.csv")
CSVfile1 <- as.matrix(CSVfile1)

Genotype1 <- CSVfile1[, "Genotype"]
Genotype1 <- sub("/","_", Genotype1)
Mouse_Label1 <- CSVfile1[, 'Label.Barcode']


# path to the files sent on March 08, 2016
pathMice <- paste("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2016-03-08/MLN M 169+")
# Reads all folders and files in current path folder and makes a list of all of their paths
allMice <- dir(pathMice, full.names=T, recursive=T, pattern = "*.fcs")
store.allMice2 <- sapply(1:length(allMice), function(x){pathMice})
store.allMice2 <- cbind(store.allMice2, sapply(1:length(allMice), function(x){unlist(strsplit(allMice[x], split = "/"))[length(unlist(strsplit(allMice[x], split = "/")))-1]}))
store.allMice2 <- cbind(store.allMice2, sapply(1:length(allMice), function(x){unlist(strsplit(allMice[x], split = "/"))[length(unlist(strsplit(allMice[x], split = "/")))]}))
store.allMice2 <- cbind(store.allMice2, str_extract(store.allMice2[,3],"L[0-9]+"))

index.Remove <- which(is.na(store.allMice2[,4]))
if(length(index.Remove != 0)){
    store.allMice2<-store.allMice2[-index.Remove,]
}

## Combining the spreadsheets from August and March. 
CSVfile2 <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2016-03-08/attachments/MLN data RB.csv")
CSVfile2 <- as.matrix(CSVfile2)
Genotype2 <- CSVfile2[,"Genotype"]
Genotype2 <- sub("/","_", Genotype2)
Mouse_Label2 <- CSVfile2[,'Label.Barcode']

# path to the files sent on Nov 11, 2016
pathMice <- paste("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2016-10-11/MLN M Labelled")
# Reads all folders and files in current path folder and makes a list of all of their paths
allMice <- dir(pathMice, full.names=T, recursive=T, pattern = "*.fcs")
store.allMice3 <- sapply(1:length(allMice), function(x){pathMice})
store.allMice3 <- cbind(store.allMice3, sapply(1:length(allMice), function(x){unlist(strsplit(allMice[x], split = "/"))[length(unlist(strsplit(allMice[x], split = "/")))-1]}))
store.allMice3 <- cbind(store.allMice3, sapply(1:length(allMice), function(x){unlist(strsplit(allMice[x], split = "/"))[length(unlist(strsplit(allMice[x], split = "/")))]}))
store.allMice3 <- cbind(store.allMice3, str_extract(store.allMice3[,3],"L[0-9]+"))

index.Remove <- which(is.na(store.allMice3[,4]))
if(length(index.Remove != 0)){
  store.allMice3<-store.allMice3[-index.Remove,]
}

## Combining the spreadsheets from August and March. 
CSVfile3 <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2016-10-11/attachments/download_phnMLNImmuno.csv")
CSVfile3 <- as.matrix(CSVfile3)
Genotype3 <- CSVfile3[,"Genotype"]
Genotype3 <- sub("/","_", Genotype3)
Mouse_Label3 <- CSVfile3[,'Label.Barcode']

Genotype <- c(Genotype1, Genotype2, Genotype3)
Mouse_Label <- c(Mouse_Label1, Mouse_Label2, Mouse_Label3)
uniqueGT <- unique(Genotype)
CSVfile <- rbind.fill.matrix(CSVfile1, CSVfile2, CSVfile3) 

invisible(sapply(1:length(uniqueGT), function(x){
dir.create ( paste(FCS.Groups,"/", uniqueGT[x], sep="") )}))

count.FCSmoved <- 0
countX <-0
index.Remove1 <- 0
index.Skip1 <- 0
for(x in 1:nrow(store.allMice1)){
    temp <- grep(store.allMice1[x,4], Mouse_Label1)
    if(length(temp) == 1){
      print(paste("Move:", store.allMice1[x,3]))
      if (!file.exists(paste(FCS.Groups, "/", Genotype1[temp], "/", store.allMice1[x,3], sep="")) ) {
        print(paste("Copying", store.allMice1[x,3]))
        count.FCSmoved <- count.FCSmoved+1
        file.copy(paste(store.allMice1[x,1], "/", store.allMice1[x,3], sep=""), paste(FCS.Groups, "/", Genotype1[temp], "/", store.allMice1[x,3], sep=""))
      } else {
        print(paste("Skipping", store.allMice1[x,3]))
        index.Skip1 <- c(index.Skip1,x)
      }
    } else{print(paste("Did not move", store.allMice1[x,3]))
      index.Remove1 <- c(index.Remove1,x)
      countX <- countX+1}
}
index.Remove1 <-index.Remove1[index.Remove1 !=0]
store.allMice1<-store.allMice1[-index.Remove1,]

index.Remove2 <- 0
index.Skip2 <- 0
for(x in 1:nrow(store.allMice2)){
  temp <- grep(store.allMice2[x,4], Mouse_Label2)
  if(length(temp) == 1){
    print(paste("Move:", store.allMice2[x,3]))
    if (!file.exists(paste(FCS.Groups, "/", Genotype2[temp], "/", store.allMice2[x,3], sep="")) ) {
      print(paste("Copying", store.allMice2[x,3]))
      count.FCSmoved <- count.FCSmoved+1
      file.copy(paste(store.allMice2[x,1], "/", store.allMice2[x,3], sep=""), paste(FCS.Groups, "/", Genotype2[temp], "/", store.allMice2[x,3], sep=""))
    } else {
      print(paste("Skipping", store.allMice2[x,3]))
      index.Skip2 <- c(index.Skip2,x)
    }
  } else{print(paste("Did not move", store.allMice2[x,3]))
    index.Remove2 <- c(index.Remove2,x)
    countX <- countX+1}
}
index.Remove2 <-index.Remove2[index.Remove2 !=0]
store.allMice2 <- store.allMice2[-index.Remove2,]

index.Remove3 <- 0
index.Skip3 <- 0
for(x in 1:nrow(store.allMice3)){
  temp <- grep(store.allMice3[x,4], Mouse_Label3)
  if(length(temp) == 1){
    print(paste("Move:", store.allMice3[x,3]))
    if (!file.exists(paste(FCS.Groups, "/", Genotype3[temp], "/", store.allMice3[x,3], sep="")) ) {
      print(paste("Copying", store.allMice3[x,3]))
      count.FCSmoved <- count.FCSmoved+1
      file.copy(paste(store.allMice3[x,1], "/", store.allMice3[x,3], sep=""), paste(FCS.Groups, "/", Genotype3[temp], "/", store.allMice3[x,3], sep=""))
    } else {
      print(paste("Skipping", store.allMice3[x,3]))
      index.Skip3 <- c(index.Skip3,x)
    }
  } else{print(paste("Did not move", store.allMice3[x,3]))
    index.Remove3 <- c(index.Remove3,x)
    countX <- countX+1}
}
index.Remove3 <-index.Remove3[index.Remove3 !=0]
store.allMice3 <- store.allMice3[-index.Remove3,]

## Combine information of Datasets from August & March in one matrix.
store.allMice <- rbind(store.allMice1, store.allMice2, store.allMice3)
#remove(store.allMice1, store.allMice2)
colnames(store.allMice) <- c("Path", "Folder", "FCS file", "Label Barcode")


print(paste("Copied and moved", count.FCSmoved, "files"))

suppressWarnings(dir.create ( "Results_MLN/"))
save ( store.allMice , file =  paste("Results_MLN/store.allMice.Rdata",sep="") )
save(Genotype, file = paste("Results_MLN/Genotype.Rdata", sep = ""))
save(uniqueGT, file = paste("Results_MLN/uniqueGT.Rdata", sep = ""))
save(Mouse_Label, file = paste("Results_MLN/Mouse_Label.Rdata", sep = ""))
## Saving the two combined spreadsheets
save(CSVfile, file = paste("Results_MLN/CSVfile.Rdata", sep = ""))

cat("Total time is: ",TimeOutput(start),sep="")



