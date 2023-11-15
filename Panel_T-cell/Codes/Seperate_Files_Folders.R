# Developed by Albina Rahim
# Date: March 10, 2016

remove(list=ls())

setwd("/code/Projects/3i/Panel_T-cell_MLN/")

library("flowCore")
library("flowBin")
library("flowDensity")
library("flowType")
library("stringr")

source("Codes/3iTcellfunctions.R")

start <- Sys.time()

FCS.Groups <- "FCS_Groups"
dir.create (FCS.Groups)

# # path to the files sent on August 27, 2015-MLN Organ
# pathMice <- paste("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2015-08-27/MLN T Labelled")

# path to the files sent on October 11, 2016-MLN Organ
pathMice <- paste("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2016-10-11/MLN T Labelled")

# Reads all folders and files in current path folder and makes a list of all of their paths
allMice <- dir(pathMice, full.names=T, recursive=T, pattern = "*.fcs") 

store.allMice1 <- sapply(1:length(allMice), function(x){pathMice})
store.allMice1 <- cbind(store.allMice1, sapply(1:length(allMice), function(x){unlist(strsplit(allMice[x], split = "/"))[length(unlist(strsplit(allMice[x], split = "/")))-1]}))
store.allMice1 <- cbind(store.allMice1, sapply(1:length(allMice), function(x){unlist(strsplit(allMice[x], split = "/"))[length(unlist(strsplit(allMice[x], split = "/")))]}))
store.allMice1 <- cbind(store.allMice1, str_extract(store.allMice1[,3],"L[0-9]+"))
store.allMice1 <- cbind(store.allMice1, NA) # Column for Assay date
store.allMice1 <- cbind(store.allMice1, NA) # Column of Gender

index.Remove <- which(is.na(store.allMice1[,4]))
if(length(index.Remove != 0)){
  store.allMice1 <- store.allMice1[-index.Remove,]
}

# ## Reading the spreadsheet from August 2015 dataset for the MLN Organ
# CSVfile1 <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2015-08-27/attachments/download_phnMLNImmuno.csv")

## Reading the spreadsheet from October 2016 dataset for the MLN Organ
CSVfile1 <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2016-10-11/attachments/download_phnMLNImmuno.csv")
CSVfile1 <- as.matrix(CSVfile1)
Genotype1 <- CSVfile1[,c(11)]
Genotype1 <- sub("/","_", Genotype1)
Mouse_Label1 <- CSVfile1[,c(20)]
Assay_Date1 <- CSVfile1[,c(15)]
Gender1 <- CSVfile1[,c(12)]

# # path to the files sent on March 08, 2016-Spleen Organ
# pathMice <- paste("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2016-03-08/SPLN T 169 +")

# path to the files sent on March 08, 2016-MLN Organ
pathMice <- paste("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2016-03-08/MLN T 169+")
# Reads all folders and files in current path folder and makes a list of all of their paths
allMice <- dir(pathMice, full.names=T, recursive=T, pattern = "*.fcs")

store.allMice2 <- sapply(1:length(allMice), function(x){pathMice})
store.allMice2 <- cbind(store.allMice2, sapply(1:length(allMice), function(x){unlist(strsplit(allMice[x], split = "/"))[length(unlist(strsplit(allMice[x], split = "/")))-1]}))
store.allMice2 <- cbind(store.allMice2, sapply(1:length(allMice), function(x){unlist(strsplit(allMice[x], split = "/"))[length(unlist(strsplit(allMice[x], split = "/")))]}))
store.allMice2 <- cbind(store.allMice2, str_extract(store.allMice2[,3],"L[0-9]+"))
store.allMice2 <- cbind(store.allMice2, NA) # Column for Assay date
store.allMice2 <- cbind(store.allMice2, NA) # Column for Gender

index.Remove <- which(is.na(store.allMice2[,4]))
if(length(index.Remove != 0)){
  store.allMice2<-store.allMice2[-index.Remove,]
}

# ## Reading the spreadsheet from March dataset for Spleen Organ
# CSVfile2 <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2016-03-08/attachments/SPLN data RB.csv")
## Reading the spreadsheet from March dataset for MLN Organ
CSVfile2 <- read.csv("/mnt/f/FCS data/IMPC/IMPC_2016_Pilot_Data/KCL_WTSI_2016-03-08/attachments/MLN data RB.csv")
CSVfile2 <- as.matrix(CSVfile2)
Genotype2 <- CSVfile2[,c(11)]
Genotype2 <- sub("/","_", Genotype2)
Mouse_Label2 <- CSVfile2[,c(20)]
Assay_Date2 <- CSVfile2[,c(15)]
Gender2 <- CSVfile2[,c(12)]


Genotype <- c(Genotype1, Genotype2)
Mouse_Label <- c(Mouse_Label1, Mouse_Label2)
uniqueGT <- unique(Genotype)

## Combining the spreadsheets from August and March.  
CSVfile <- rbind(CSVfile1, CSVfile2) 

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
      store.allMice1[x,5] <- Assay_Date1[temp]
      store.allMice1[x,6] <- Gender1[temp]
      file.copy(paste(store.allMice1[x,1], "/", store.allMice1[x,3], sep=""), paste(FCS.Groups, "/", Genotype1[temp], "/", store.allMice1[x,3], sep=""))
    } else {
      print(paste("Skipping", store.allMice1[x,3]))
      store.allMice1[x,5] <- Assay_Date1[temp]
      store.allMice1[x,6] <- Gender1[temp]
      index.Skip1 <- c(index.Skip1,x)
    }
  } else{print(paste("Did not move", store.allMice1[x,3]))
    index.Remove1 <- c(index.Remove1,x)
    countX <- countX+1}
}
index.Remove1 <-index.Remove1[index.Remove1 !=0]
if(length(index.Remove1 != 0)){
  store.allMice1<-store.allMice1[-index.Remove1,]
}

index.Remove2 <- 0
index.Skip2 <- 0
for(x in 1:nrow(store.allMice2)){
  temp <- grep(store.allMice2[x,4], Mouse_Label2)
  if(length(temp) == 1){
    print(paste("Move:", store.allMice2[x,3]))
    if (!file.exists(paste(FCS.Groups, "/", Genotype2[temp], "/", store.allMice2[x,3], sep="")) ) {
      print(paste("Copying", store.allMice2[x,3]))
      count.FCSmoved <- count.FCSmoved+1
      store.allMice2[x,5] <- Assay_Date2[temp]
      store.allMice2[x,6] <- Gender2[temp]
      file.copy(paste(store.allMice2[x,1], "/", store.allMice2[x,3], sep=""), paste(FCS.Groups, "/", Genotype2[temp], "/", store.allMice2[x,3], sep=""))
    } else {
      print(paste("Skipping", store.allMice2[x,3]))
      store.allMice2[x,5] <- Assay_Date2[temp]
      store.allMice2[x,6] <- Gender2[temp]
      index.Skip2 <- c(index.Skip2,x)
    }
  } else{print(paste("Did not move", store.allMice2[x,3]))
    index.Remove2 <- c(index.Remove2,x)
    countX <- countX+1}
}
index.Remove2 <-index.Remove2[index.Remove2 !=0]
if(length(index.Remove2 != 0)){
    store.allMice2 <- store.allMice2[-index.Remove2,]
}
## Combine information of Datasets from August & March in one matrix.
store.allMice <- rbind(store.allMice1, store.allMice2)
remove(store.allMice1, store.allMice2)
colnames(store.allMice) <- c("Path", "Folder", "FCS file", "Label Barcode", "Assay Date", "Gender")


#count.FCSmoved <- length(intersect(Mouse_Label, store.allMice[,3]))
print(paste("Copied and moved", count.FCSmoved, "files"))

suppressWarnings(dir.create ( "Results/"))
save ( store.allMice , file =  paste("Results/store.allMice.Rdata",sep="") )
save(Genotype, file = paste("Results/Genotype.Rdata", sep = ""))
save(uniqueGT, file = paste("Results/uniqueGT.Rdata", sep = ""))

cat("Total time is: ",TimeOutput(start),sep="")




