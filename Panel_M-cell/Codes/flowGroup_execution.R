# Written by Albina Rahim on January 2019
# Updated by Marjan Najafabadipour on November 2019

#Execution of the flowGroup function

remove(list=ls())
setwd("~/code/Projects/3i/Panel_M-cell/")

#Loding the required packages for executing the flowGroup algorithm
library("flowCore")
library("flowDensity")
library("foreach")
library("Cairo")
library("MASS") 
library("doMC")
library("stringr") 
library("OpenImageR")
library("flowType")
library("Rtsne")
library("oce")
library("plot3D")
library("stringr")
library("flowViz")
library("dplyr")

#Reading flowGroup algorithm and functions
source("Codes/flowGroup.R")
source("Codes/flowGroup_functions.R")
#source("Codes/flowOutlier.R")

##########################################################################

library("plyr")
library("doMC")

##  Extra step for saving CD45 .Rdata files as .fcs files
no_cores <- detectCores() - 4
registerDoMC(no_cores)

input_path <- "/home/rstudio/data/IMPC/3i/Panel_M-cell/SPLEEN/Results"
output_path <- "/home/rstudio/results/IMPC/3i/Panel_BM-cell/flowGroup02"

load(paste0(input_path,"/store.allFCS.Rdata"))
load(paste0(input_path, "/res.clean.Rdata"))
load(paste0(input_path,"/ind.marg.neg.clean.all.Rdata"))
#load(paste0(input_path, "/failedGating.files.Rdata"))



# Remove files which failed TvsF
rem.failed.clean <- which(res.clean["Has the file passed",] == "F")
print(paste0("Number of files flagged by TvsF: ", length(rem.failed.clean)))

# Since I want to manually pass so many I will fail only specific barcodes
rem.failed.clean.barcodes <- c('L000073120', 'L000088655')
rem.failed.clean <- unlist(lapply(rem.failed.clean.barcodes, function(x){ grep(x, store.allFCS[,5])}))
if(length(rem.failed.clean) != 0){
  store.allFCS <- store.allFCS[-rem.failed.clean, ]
  res.clean <- res.clean[, -rem.failed.clean]
  ind.marg.neg.clean.all[rem.failed.clean] <- NULL
}
print(paste0("Number of TvsF-flagged FCS files failed after manual inspection: ", length(rem.failed.clean)))

# According to Adam Laing these files should be disregarded/failed as one of the lineage antibodies wasn’t added
# to the panel & this has a major impact on the downstream data.
manually.failed.barcodes <- c(43170:43182, 43215:43226, 43228:43231, 44310:44323)
manually.failed.idx <- unlist(lapply(manually.failed.barcodes, function(x){ grep(x, store.allFCS[,5])}))
if(length(manually.failed.idx) > 0){
  store.allFCS <- store.allFCS[-manually.failed.idx, ]
  res.clean <- res.clean[, -manually.failed.idx]
  ind.marg.neg.clean.all[manually.failed.idx] <- NULL
}
print(paste0("Number of files manually failed because lineage antibody wasn't added (confirmed by Adam Laing): ", length(manually.failed.idx)))

# Not enough live cells, in addition CD11b vs F4/80 distribution looks niO
manually.failed.barcodes <- c(95648, 90915)
manually.failed.idx <- unlist(lapply(manually.failed.barcodes, function(x){ grep(x, store.allFCS[,5])}))
if(length(manually.failed.idx) > 0){
  store.allFCS <- store.allFCS[-manually.failed.idx, ]
  res.clean <- res.clean[, -manually.failed.idx]
  ind.marg.neg.clean.all[manually.failed.idx] <- NULL
}
print(paste0("Number of files manually failed because not many live cells & CD11b/F4/80 distribution abnormal: ", length(manually.failed.idx)))



rownames(store.allFCS) <- 1:nrow(store.allFCS)

file.names <- data.frame(store.allFCS, stringsAsFactors = F)
file.names[,1] <- "/home/rstudio/data/IMPC/3i/Panel_M-cell/semisupervised/Spleen/ParentPopulation"
## Comment: There were NO files which failed the gating in the first run



props.events <- ldply(1:nrow(file.names), function(i){ 
  
  x <- file.names[i,]
  
  # Load .Rdata
  attach(paste0(x$Path,"/", x$FCS.file, ".Rdata"))
  ls(pos=2)
 
  write.FCS(cd45clean, filename = paste0(input_path,"/CD45_Populations/",x$FCS.files))
  
  data.frame(x$FCS.files)
  
}, .parallel = TRUE) # end ldply


################################################
#Input from the user
input_path <- paste0(input_path,"/CD45_Populations/")
# output_path <- "/result/..."
#Enter the names of the channels (Do not enter the names of the markers)
x_mark <- "PE-Cy7-A"
y_mark <- "BV786-A"

#Calling floWGroup function - You should always keep flwT = 1
res <- flowGroup(dir_path = input_path, xMark = x_mark, yMark = y_mark, 
                 experiment_name = output_path, HOG_cex = 1, 
                 partitions_flowType = 5, vuc = 0, vuc_d = 0, hog = 0, flwT = 1, 
                 DoOverlay=F, DoContour = T, plot_groups = T, verbose=T, 
                 vert_line = NULL, horiz_line = NULL, boundaries = NULL)
