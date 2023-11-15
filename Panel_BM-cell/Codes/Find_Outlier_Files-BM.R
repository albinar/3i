# Written by Albina Rahim
# Date: March 21, 2016

remove(list=ls())

setwd("/code/Projects/3i/Panel_BM-cell/")

library("png")

######################################################################

outliersJM <- function(x, SDs = 3){
  mean <- mean(x)
  sd <- sd(x)
  indices <- union(which(x >= (mean + SDs*sd)),which(x <= (mean - SDs*sd)))
  outliers <- x[indices]
  return(list(out=outliers, ind=indices))
}

######################################################################


results.dir <- "/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results"
load(paste0(results.dir,"/store.allFCS.Rdata"))
load(paste0(results.dir,"/ind.marg.neg.clean.all.Rdata"))
load(paste0(results.dir, "/res.clean.Rdata"))
load(paste0(results.dir, "/all.gthres.Rdata"))
allProportions_Table_Initial <- read.csv(paste0(results.dir, "/allProportions_Table_Initial.csv"))

## Reading the TvsF flagged csv file which was sent to Adam for feedback
flagged.FCS <- as.matrix(read.csv(paste0(results.dir, "/flagged.FCS.BoneMarrow.Feedback.csv"), sep = ","))
## List the indices of the flagged FCS files based on TvsF algorithm
flagged.FCS.index <- which(res.clean[c('Has the file passed'),] == "F")

rownames(flagged.FCS) <- flagged.FCS.index

## Manually including some flagged files for further analysis
manualExclude.flagged.FCS <- flagged.FCS[-which(flagged.FCS[,c('Final.Action')] == "Keep"),]
manualExclude.flagged.FCS.index <- as.integer(rownames(manualExclude.flagged.FCS))

## Removing flagged files based on Adam's feedback from further analysis after recieving confirmation from Adam (Jan 04, 2017)
if ( length(manualExclude.flagged.FCS.index) != 0){
  store.allFCS <- store.allFCS[-manualExclude.flagged.FCS.index,]
  res.clean <- res.clean[,-manualExclude.flagged.FCS.index]
  ind.marg.neg.clean.all[manualExclude.flagged.FCS.index] <- NULL
}
rownames(store.allFCS) <- 1:nrow(store.allFCS)

## Comment: There were NO files which failed the gating in the first run

#boxplot(all.gthres[,c('singlets.gate.h','fsc.a.gate.low','cd45.gate.low','b220.gate','cd138.gate','igm.gate', 'igd.gate')], pch=19, cex=0.2)


#############################################################################################################

## Finding outliers and fixing gates for 'singlets.gate.h'
boxplot(all.gthres[,c('singlets.gate.h')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('singlets.gate.h')], SDs = 3.5)$ind, outliersJM(all.gthres[,c('singlets.gate.h')], SDs = 3.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
if(nrow(indices.check) != 0){
  fixed.Gthres <- mean(all.gthres[,c('singlets.gate.h')])
  all.gthres[indices.check[,1],c('singlets.gate.h')] <- fixed.Gthres 
} else{
  print("No gating threshold outliers detected")
}

###########################################################################################################

## Finding outliers and fixing gates for 'singlets.gate.l'
boxplot(all.gthres[,c('singlets.gate.l')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('singlets.gate.l')], SDs = 3.35)$ind, outliersJM(all.gthres[,c('singlets.gate.l')], SDs = 3.35)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
## It was noticed that for singlets.gate.l the gates were fine; however, the outliers detected were mostly of the fsc.a.gate.low
## Therefore, I am saving these indices manually, so that they can be fixed in the section where we find outliers for the fsc.a.gate.low gate
indices.check.fsc.gate.low <- indices.check[-3,]

# if(nrow(indices.check) != 0){
#   fixed.Gthres <- mean(all.gthres[,c('singlets.gate.l')])
#   all.gthres[ind.with.store,c('singlets.gate.l')] <- fixed.Gthres 
# }


###########################################################################################################

## Finding outliers and fixing gates for 'live.gate'
boxplot(all.gthres[,c('live.gate')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('live.gate')], SDs = 3.75)$ind, outliersJM(all.gthres[,c('live.gate')], SDs = 3.75)$out)
#indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}

## Manually adding the files sent by Ania, which had issues with the live.gate. 
manual.index.live <- c('L000046772')
manual.index.live <- sapply(1:length(manual.index.live), function(x){which(store.allFCS[,c('Barcodes')] == manual.index.live[x])})



if(nrow(indices.check) != 0){
  fixed.Gthres.live.gate <- mean(all.gthres[,c('live.gate')])
  all.gthres[indices.check[,1],c('live.gate')] <- fixed.Gthres.live.gate 
  all.gthres[manual.index.live,c('live.gate')] <- 2.9
} else{
  print("No gating threshold outliers detected")
}

###########################################################################################################

## Finding outliers and fixing gates for 'fsc.a.gate.high'
boxplot(all.gthres[,c('fsc.a.gate.high')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('fsc.a.gate.high')], SDs = 2.5)$ind, outliersJM(all.gthres[,c('fsc.a.gate.high')], SDs = 2.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
# if(nrow(indices.check) != 0){
#   fixed.Gthres <- mean(all.gthres[,c('fsc.a.gate.high')])
#   all.gthres[indices.check[,1],c('fsc.a.gate.high')] <- fixed.Gthres 
# }

###########################################################################################################

## Finding outliers and fixing gates for 'fsc.a.gate.low'
boxplot(all.gthres[,c('fsc.a.gate.low')], pch=19, cex=0.2)

indices.check <- cbind(outliersJM(all.gthres[,c('fsc.a.gate.low')], SDs = 2)$ind, outliersJM(all.gthres[,c('fsc.a.gate.low')], SDs = 2)$out)
indices.check <- indices.check[order(indices.check[,1]),]
indices.check.temp.low <- indices.check[which(indices.check[,2] < 29000),]
## Manually excluding the file: BM_L000138004_A12_001.labelled.fcs (index 2652) since the singlets population looks odd
indices.check.temp.low <- indices.check.temp.low[-22,]
#ind.with.store <- unique(sort(indices.check.temp.low[,1]))

indices.check.temp.high <- indices.check[which(indices.check[,2] > 90000),]
indices.check.temp <- indices.check[which(indices.check[,2] >= 29000),]
indices.check.temp <- indices.check.temp[which(indices.check.temp[,2] <= 90000),]
ind.with.store <- unique(sort(indices.check.temp[,1]))
#for (k in 151:length(ind.with.store)){
for (k in 151:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
if(nrow(indices.check) != 0){
  fixed.Gthres.fsc.a.gate.low <- mean(all.gthres[,c('fsc.a.gate.low')])
  all.gthres[indices.check.temp.low[,1],c('fsc.a.gate.low')] <- fixed.Gthres.fsc.a.gate.low 
  all.gthres[indices.check.temp.high[,1],c('fsc.a.gate.low')] <- fixed.Gthres.fsc.a.gate.low 
  ## manually fixing 14 gating thresholds based on how the plots looked:
  all.gthres[indices.check.temp[c(1:5,8,73:74,97:101,104),1],c('fsc.a.gate.low')] <- fixed.Gthres.fsc.a.gate.low
  ## Fixing the fsc.a.gate.low outliers detected while checking for singlets.gate.l outliers
  all.gthres[indices.check.fsc.gate.low[,1],c('fsc.a.gate.low')] <- fixed.Gthres.fsc.a.gate.low
} else{
  print("No gating threshold outliers detected")
}


###########################################################################################################

## Finding outliers and fixing gates for 'ssc.a.gate'
boxplot(all.gthres[,c('ssc.a.gate')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('ssc.a.gate')], SDs = 3)$ind, outliersJM(all.gthres[,c('ssc.a.gate')], SDs = 3)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
# if(nrow(indices.check) != 0){
#   fixed.Gthres <- mean(all.gthres[,c('ssc.a.gate')])
#   all.gthres[indices.check[,1],c('ssc.a.gate')] <- fixed.Gthres 
# }

###########################################################################################################

## Finding outliers and fixing gates for 'cd45.gate.low'
boxplot(all.gthres[,c('cd45.gate.low')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('cd45.gate.low')], SDs = 2.75)$ind, outliersJM(all.gthres[,c('cd45.gate.low')], SDs = 2.75)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])

#for (k in 31:75){
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
#ind.with.store.temp <- ind.with.store[c(1,2,4:7,10:16,18,19,21:26,28:30,33,35,36,39,48)]
# if(nrow(indices.check) != 0){
#   #all.gthres[indices.check[c(32,38:42),1],c('fsc.a.gate.low')] <- fixed.Gthres.fsc.a.gate.low
#   fixed.Gthres <- mean(all.gthres[,c('cd45.gate.low')])
#   all.gthres[ind.with.store,c('cd45.gate.low')] <- fixed.Gthres 
# } else{
#   print("No gating threshold outliers detected")
# }


###########################################################################################################

## Finding outliers and fixing gates for 'cd45.gate.high'
boxplot(all.gthres[,c('cd45.gate.high')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('cd45.gate.high')], SDs = 3.25)$ind, outliersJM(all.gthres[,c('cd45.gate.high')], SDs = 3.25)$out)
#indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])

#for (k in 31:57){
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}

if(nrow(indices.check) != 0){
  fixed.Gthres.cd45.gate.high <- mean(all.gthres[,c('cd45.gate.high')])
  all.gthres[indices.check[6:7,1],c('cd45.gate.high')] <- fixed.Gthres.cd45.gate.high
} else{
  print("No gating threshold outliers detected")
}

###########################################################################################################

## Finding outliers and fixing gates for 'gr1.gate'
boxplot(all.gthres[,c('gr1.gate')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('gr1.gate')], SDs = 4)$ind, outliersJM(all.gthres[,c('gr1.gate')], SDs = 4)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])

for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
# if(nrow(indices.check) != 0){
#   fixed.Gthres <- mean(all.gthres[,c('gr1.gate')])
#   all.gthres[indices.check[,1],c('gr1.gate')] <- fixed.Gthres 
# }


###########################################################################################################

## Finding outliers and fixing gates for 'b220.gate'
boxplot(all.gthres[,c('b220.gate')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('b220.gate')], SDs = 4)$ind, outliersJM(all.gthres[,c('b220.gate')], SDs = 4)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])

for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
if(length(ind.with.store) != 0){
  #fixed.Gthres.b220.gate <- mean(all.gthres[,c('b220.gate')])
  all.gthres[ind.with.store,c('b220.gate')] <- 1.5
} else{
  print("No gating threshold outliers detected")
}

#############################################################################################################

## Finding outliers and fixing gates for 'cd138.gate'
boxplot(all.gthres[,c('cd138.gate')], pch=19, cex=0.2)
#indices.check <- cbind(outliersJM(all.gthres[,c('cd138.gate')], SDs = 2)$ind, outliersJM(all.gthres[,c('cd138.gate')], SDs = 2)$out)
#indices.check <- indices.check[order(indices.check[,1]),]
#ind.with.store <- unique(indices.check[,1])
ind.with.store <- which(allProportions_Table_Initial[,c("Plasma..Parent")] > 1.25) 
prop.temp <- cbind(ind.with.store, allProportions_Table_Initial[ind.with.store,c("Barcodes","Plasma..Parent")])
#ind.with.store.temp <- which(all.gthres[ind.with.store,c('cd138.gate')] < 1)
for (k in 1:40){
#for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
if(length(ind.with.store) != 0){
  fixed.Gthres.cd138.gate <- mean(all.gthres[,c('cd138.gate')])
  all.gthres[ind.with.store,c('cd138.gate')] <- fixed.Gthres.cd138.gate
} else{
  print("No gating threshold outliers detected")
}

#############################################################################################################

## Finding outliers and fixing gates for 'cd43.gate'
boxplot(all.gthres[,c('cd43.gate')], pch=19, cex=0.2)
# indices.check <- cbind(outliersJM(all.gthres[,c('cd43.gate')], SDs = 3.5)$ind, outliersJM(all.gthres[,c('cd43.gate')], SDs = 3.5)$out)
# indices.check <- indices.check[order(indices.check[,1]),]
# ind.with.store <- unique(indices.check[,1])

## Manually adding the files sent by Ania, which had issues with the cd43.gate. 
manual.index.cd43 <- c('L000046772', 'L000046787')
manual.index.cd43 <- sapply(1:length(manual.index.cd43), function(x){which(store.allFCS[,c('Barcodes')] == manual.index.cd43[x])})

all.gthres[manual.index.cd43,c('cd43.gate')] <- 2.5
# for (k in 1:length(ind.with.store)){
#   name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results/Figures/ScatterPlots/",
#                  store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
#   x11(width = 16, height = 9)
#   img <- readPNG(name)
#   grid::grid.raster(img)
# }

#############################################################################################################

## Finding outliers and fixing gates for 'cd24.gate'
boxplot(all.gthres[,c('cd24.gate')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('cd24.gate')], SDs = 3)$ind, outliersJM(all.gthres[,c('cd24.gate')], SDs = 3)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
all.gthres[indices.check[c(8:6,2),1],c('fsc.a.gate.low')] <- fixed.Gthres.fsc.a.gate.low

## Manually adding the files sent by Ania, which had issues with the cd24.gate. 
manual.index.cd24 <- c('L000046772')
manual.index.cd24 <- sapply(1:length(manual.index.cd24), function(x){which(store.allFCS[,c('Barcodes')] == manual.index.cd24[x])})

if(nrow(indices.check) != 0){
  fixed.Gthres.cd24.gate <- mean(all.gthres[,c('cd24.gate')])
  all.gthres[indices.check[c(1:3,5,7:8,11),1],c('cd24.gate')] <- fixed.Gthres.cd24.gate 
  all.gthres[manual.index.cd24,c('cd24.gate')] <- 3.3
} else{
  print("No gating threshold outliers detected")
}

#############################################################################################################

## Finding outliers and fixing gates for 'bp1.gate'
boxplot(all.gthres[,c('bp1.gate')], pch=19, cex=0.2)
# indices.check <- cbind(outliersJM(all.gthres[,c('bp1.gate')], SDs = 3)$ind, outliersJM(all.gthres[,c('bp1.gate')], SDs = 3)$out)
# indices.check <- indices.check[order(indices.check[,1]),]
# ind.with.store <- unique(indices.check[,1])
# for (k in 1:length(ind.with.store)){
#   name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results/Figures/ScatterPlots/",
#                  store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
#   x11(width = 16, height = 9)
#   img <- readPNG(name)
#   grid::grid.raster(img)
# }

## Manually adding the files sent by Ania, which had issues with the bp1.gate. 
manual.index.bp1 <- c('L000046772', 'L000046787')
manual.index.bp1 <- sapply(1:length(manual.index.bp1), function(x){which(store.allFCS[,c('Barcodes')] == manual.index.bp1[x])})

all.gthres[manual.index.bp1,c('bp1.gate')] <- 1


# if(nrow(indices.check) != 0){
#   fixed.Gthres.bp1.gate <- mean(all.gthres[,c('bp1.gate')])
#   all.gthres[indices.check[,1],c('bp1.gate')] <- fixed.Gthres.bp1.gate
# } else{
#   print("No gating threshold outliers detected")
# }
#############################################################################################################

## Finding outliers and fixing gates for 'igm.gate'
boxplot(all.gthres[,c('igm.gate')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('igm.gate')], SDs = 3)$ind, outliersJM(all.gthres[,c('igm.gate')], SDs = 3)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
#all.gthres[indices.check[2,1],c('live.gate')] <- fixed.Gthres.live.gate
all.gthres[indices.check[c(13:20,22:25, 27:31),1],c('fsc.a.gate.low')] <- fixed.Gthres.fsc.a.gate.low

if(nrow(indices.check) != 0){
  fixed.Gthres.igm.gate <- mean(all.gthres[,c('igm.gate')])
  all.gthres[indices.check[c(1,13:nrow(indices.check)),1],c('igm.gate')] <- fixed.Gthres.igm.gate 
} else{
  print("No gating threshold outliers detected")
}

#############################################################################################################

## Finding outliers and fixing gates for 'igd.gate'
boxplot(all.gthres[,c('igd.gate')], pch=19, cex=0.2)
# indices.check <- cbind(outliersJM(all.gthres[,c('igd.gate')], SDs = 3)$ind, outliersJM(all.gthres[,c('igd.gate')], SDs = 3)$out)
# indices.check <- indices.check[order(indices.check[,1]),]
# ind.with.store <- unique(indices.check[,1])

ind.with.store <- which(all.gthres[,c('igd.gate')] > 2)
for (k in 41:length(ind.with.store)){
#for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_BM-cell/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}

all.gthres[ind.with.store[c(52,50,24,23,21,16,77,76,75,74,72,57,55)],c('fsc.a.gate.low')] <- fixed.Gthres.fsc.a.gate.low
all.gthres[ind.with.store[c(1,2,3,16,19, 20,21,23,24,32,83:86,71:77,79,59:66,52,49,46)],c('cd138.gate')] <- fixed.Gthres.cd138.gate
all.gthres[ind.with.store[c(16,19,20)],c('cd24.gate')] <- fixed.Gthres.cd24.gate 
all.gthres[ind.with.store[c(52)],c('igm.gate')] <- fixed.Gthres.igm.gate

## Manually adding the files sent by Ania, which had issues with the igd.gate. 
manual.index.igd <- c('L000046772', 'L000046776', 'L000046777', 'L000046778', 'L000046785')
manual.index.igd <- sapply(1:length(manual.index.igd), function(x){which(store.allFCS[,c('Barcodes')] == manual.index.igd[x])})
ind.with.store <- c(ind.with.store, manual.index.igd)

if(length(ind.with.store) != 0){
  fixed.Gthres <- mean(all.gthres[,c('igd.gate')])
  all.gthres[ind.with.store[c(34:40, 25:32, 11:13, 1:3, 94:82, 80:58, 56,51:54, 48, 41:44)],c('igd.gate')] <- fixed.Gthres 
} else{
  print("No gating threshold outliers detected")
}

#############################################################################################################

all.gthres.Updated <- all.gthres

## Saving the Updated Gating Thresholds as numerical data, since we will need this for running the gating again with updated thresholds
save(all.gthres.Updated, file = paste0(results.dir, "/all.gthres.Updated.Rdata"))
all.gthres.Updated <- cbind(store.allFCS[,c('Panel/Organ', 'Genotype', 'FCS files', 'Barcodes', 'Assay Date', 'Gender', 'Number of Channels', 'Number of Cells')], all.gthres.Updated)

write.csv(all.gthres.Updated, file =  paste0(results.dir, "/allGthres_Table_Updated.csv"), row.names = FALSE)

store.allFCS.Updated <- store.allFCS
save(store.allFCS.Updated, file = paste0(results.dir, "/store.allFCS.Updated.Rdata"))

