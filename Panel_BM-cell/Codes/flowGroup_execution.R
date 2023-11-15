# Written by Albina Rahim on January 2019
# Updated by Marjan Najafabadipour on November 2019

#Execution of the flowGroup function

remove(list=ls())
setwd("~/code/Projects/3i/Panel_BM-cell/")

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
source("Codes/flowOutlier.R")

##########################################################################

library("plyr")
library("doMC")

##  Extra step for saving CD45 .Rdata files as .fcs files
no_cores <- detectCores() - 4
registerDoMC(no_cores)

input_path <- "/home/rstudio/data/IMPC/3i/Panel_BM-cell/Results"
output_path <- "/home/rstudio/results/IMPC/3i/Panel_BM-cell/flowGroup02"

load(paste0(input_path,"/store.allFCS.Updated.Rdata"))
load(paste0(input_path, "/res.clean.Rdata"))
load(paste0(input_path,"/ind.marg.neg.clean.all.Rdata"))


## Reading the TvsF flagged csv file which was sent to Adam for feedback
flagged.FCS <- as.matrix(read.csv(paste0(input_path, "/flagged.FCS.BoneMarrow.Feedback.csv"), sep = ","))
## List the indices of the flagged FCS files based on TvsF algorithm
flagged.FCS.index <- which(res.clean[c('Has the file passed'),] == "F")

rownames(flagged.FCS) <- flagged.FCS.index

## Manually including some flagged files for further analysis
manualExclude.flagged.FCS <- flagged.FCS[-which(flagged.FCS[,c('Final.Action')] == "Keep"),]
manualExclude.flagged.FCS.index <- as.integer(rownames(manualExclude.flagged.FCS))

## Removing flagged files based on Adam's feedback from further analysis after recieving confirmation from Adam (Jan 04, 2017)
if ( length(manualExclude.flagged.FCS.index) != 0){
  #store.allFCS <- store.allFCS[-manualExclude.flagged.FCS.index,]
  res.clean <- res.clean[,-manualExclude.flagged.FCS.index]
  ind.marg.neg.clean.all[manualExclude.flagged.FCS.index] <- NULL
}

file.names <- data.frame(store.allFCS.Updated, stringsAsFactors = F)
file.names[,1] <- "/home/rstudio/data/IMPC/3i/Panel_BM-cell/semisupervised/ParentPopulation"
## Comment: There were NO files which failed the gating in the first run



props.events <- ldply(1:nrow(file.names), function(i){ 
  
  x <- file.names[i,]
  
  # Load .Rdata
  attach(paste0(x$Path,"/", x$FCS.file, ".Rdata"))
  ls(pos=2)
 
  write.FCS(cd45, filename = paste0(input_path,"/CD45_Populations/",x$FCS.files))
  
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
