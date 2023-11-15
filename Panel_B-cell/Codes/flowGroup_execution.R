# Written by Albina Rahim on January 2019
# Updated by Marjan Najafabadipour on November 2019

#Execution of the flowGroup function

remove(list=ls())
setwd("~/code/Projects/3i/Panel_B-cell/")

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

input_path <- "/home/rstudio/data/IMPC/3i/Panel_B-cell/SPLEEN/Results"
#output_path <- "/home/rstudio/results/IMPC/3i/Panel_BM-cell/flowGroup02"

load(paste0(input_path,"/store.allFCS.Rdata"))
load(paste0(input_path, "/res.clean.Rdata"))
load(paste0(input_path,"/ind.marg.neg.clean.all.Rdata"))

## There are files which failed the Gating in our previous step. So we exclude those files

# rem.failed.clean <- rem.failed.clean[setdiff(names(rem.failed.clean), c("L000061295","L000097273","L000097279","L000137217","L000133899","L000133905"))] # manually pass these
files.removed <- c("L000058909","L000067085","L000067091","L000137220","L000137221","L000129428","L000135011","L000090915","L000125814",
                   "L000084521", "L000125820")
rem.failed.clean <- sort(match(files.removed,store.allFCS[,"Barcodes"])) # this also contains the TvsF removed files too.

if ( length(rem.failed.clean) != 0){
  store.allFCS <- store.allFCS[-rem.failed.clean,]
  res.clean <- res.clean[,-rem.failed.clean]
  ind.marg.neg.clean.all[rem.failed.clean] <- NULL
}
rownames(store.allFCS) <- 1:nrow(store.allFCS)

file.names <- data.frame(store.allFCS, stringsAsFactors = F)
file.names[,1] <- "/home/rstudio/data/IMPC/3i/Panel_B-cell/semisupervised/SPLEEN/ParentPopulation"




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
