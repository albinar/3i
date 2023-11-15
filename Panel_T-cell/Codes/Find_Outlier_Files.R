# Originally written by Justin Meskas 
# Modified & Updated by Albina Rahim
# Date: March 21, 2016

remove(list=ls())

setwd("/code/Projects/3i/Panel_T-cell_MLN/")


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


results.dir <-  "/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results"
load(paste0(results.dir,"/store.allFCS.Rdata"))
load(paste0(results.dir,"/ind.marg.neg.clean.all.Rdata"))
load(paste0(results.dir, "/res.clean.Rdata"))
load(paste0(results.dir, "/all.gthres.Rdata"))
load(paste0(results.dir, "/failedGating.files.Rdata"))

## Reading the TvsF flagged csv file which was sent to Adam for feedback
flagged.FCS <- as.matrix(read.csv(paste0(results.dir, "/flagged.FCS.TcellMLN.Feedback.csv"), sep = ","))
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
## There are files which failed the Gating in our previous step
## So we exclude the thresholds of these files
if (nrow(failedGating.files) != 0){
  failedGating.files.index <- as.integer(rownames(failedGating.files))
  store.allFCS <- store.allFCS[-failedGating.files.index,]
  failedGating.files.gthres <- all.gthres[failedGating.files.index,] # creating a matrix which contains the gating thresholds of the failed files, so that they can be fixed later
  res.clean <- res.clean[,-failedGating.files.index]
  ind.marg.neg.clean.all[failedGating.files.index] <- NULL
  all.gthres <- all.gthres[-failedGating.files.index,]
}

rownames(store.allFCS) <- 1:nrow(store.allFCS)
rownames(all.gthres) <- 1:nrow(all.gthres)


# boxplot(all.gthres, pch=19, cex=0.2)

#############################################################################################################

## Finding outliers and fixing gates for 'singlets.gate.h'
boxplot(all.gthres[,c('singlets.gate.h')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('singlets.gate.h')], SDs = 2.75)$ind, outliersJM(all.gthres[,c('singlets.gate.h')], SDs = 2.75)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
if(nrow(indices.check) != 0){
  fixed.Gthres <- mean(all.gthres[,c('singlets.gate.h')])
  all.gthres[indices.check[,1],c('singlets.gate.h')] <- fixed.Gthres 
}else{
  print("No gating threshold outliers detected")
}


###########################################################################################################

## Finding outliers and fixing gates for 'singlets.gate.l'
boxplot(all.gthres[,c('singlets.gate.l')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('singlets.gate.l')], SDs = 3.75)$ind, outliersJM(all.gthres[,c('singlets.gate.l')], SDs = 3.75)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
if(nrow(indices.check) != 0){
  fixed.Gthres <- mean(all.gthres[,c('singlets.gate.l')])
  all.gthres[indices.check[,1],c('singlets.gate.l')] <- fixed.Gthres 
}else{
  print("No gating threshold outliers detected")
}



###########################################################################################################

## Finding outliers and fixing gates for 'live.gate'
boxplot(all.gthres[,c('live.gate')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('live.gate')], SDs = 3.75)$ind, outliersJM(all.gthres[,c('live.gate')], SDs = 3.75)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 3:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
if(nrow(indices.check) != 0){
  fixed.Gthres <- mean(all.gthres[,c('live.gate')])
  all.gthres[indices.check[,1],c('live.gate')] <- fixed.Gthres 
}else{
  print("No gating threshold outliers detected")
}


###########################################################################################################

## Finding outliers and fixing gates for 'fsc.a.gate.high'
boxplot(all.gthres[,c('fsc.a.gate.high')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('fsc.a.gate.high')], SDs = 3)$ind, outliersJM(all.gthres[,c('fsc.a.gate.high')], SDs = 3)$out)
#ind.with.store <- unique(sort(indices.check[,1]))
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
if(nrow(indices.check) != 0){
  fixed.Gthres <- mean(all.gthres[,c('fsc.a.gate.high')])
  all.gthres[indices.check[,1],c('fsc.a.gate.high')] <- fixed.Gthres 
}else{
  print("No gating threshold outliers detected")
}

###########################################################################################################

## Finding outliers and fixing gates for 'fsc.a.gate.low'
boxplot(all.gthres[,c('fsc.a.gate.low')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('fsc.a.gate.low')], SDs = 1.25)$ind, outliersJM(all.gthres[,c('fsc.a.gate.low')], SDs = 1.25)$out)
indices.check.temp <- which(indices.check[,2] < 50000)
indices.check <- indices.check[indices.check.temp,]
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- sort(indices.check[,1])
#for (k in 1:length(ind.with.store)){
for (k in 1:20){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
if(nrow(indices.check) != 0){
  fixed.Gthres <- mean(all.gthres[,c('fsc.a.gate.low')])
  all.gthres[indices.check[,1],c('fsc.a.gate.low')] <- fixed.Gthres 
}else{
  print("No gating threshold outliers detected")
}


###########################################################################################################

## Finding outliers and fixing gates for 'ssc.a.gate'
boxplot(all.gthres[,c('ssc.a.gate')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('ssc.a.gate')], SDs = 3)$ind, outliersJM(all.gthres[,c('ssc.a.gate')], SDs = 3)$out)
ind.with.store <- unique(sort(indices.check[,1]))
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
if(nrow(indices.check) != 0){
  fixed.Gthres <- mean(all.gthres[,c('ssc.a.gate')])
  all.gthres[indices.check[,1],c('ssc.a.gate')] <- fixed.Gthres 
}else{
  print("No gating threshold outliers detected")
}


###########################################################################################################

## Finding outliers and fixing gates for 'cd45.gate.low'
boxplot(all.gthres[,c('cd45.gate.low')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('cd45.gate.low')], SDs = 3)$ind, outliersJM(all.gthres[,c('cd45.gate.low')], SDs = 3)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:20){
#for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
if(nrow(indices.check) != 0){
  fixed.Gthres <- mean(all.gthres[,c('cd45.gate.low')])
  all.gthres[indices.check[,1],c('cd45.gate.low')] <- fixed.Gthres
}else{
  print("No gating threshold outliers detected")
}


###########################################################################################################

## Finding outliers and fixing gates for 'cd45.gate.high'
boxplot(all.gthres[,c('cd45.gate.high')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('cd45.gate.high')], SDs = 2.5)$ind, outliersJM(all.gthres[,c('cd45.gate.high')], SDs = 2.5)$out)
ind.with.store <- unique(sort(indices.check[,1]))
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
# if(nrow(indices.check) != 0){
#   fixed.Gthres <- mean(all.gthres[,c('cd45.gate.high')])
#   all.gthres[indices.check[,1],c('cd45.gate.high')] <- fixed.Gthres 
# }else{
#   print("No gating threshold outliers detected")
# }


###########################################################################################################

## Finding outliers and fixing gates for 'cd161.gate.high'
boxplot(all.gthres[,c('cd161.gate.high')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('cd161.gate.high')], SDs = 3)$ind, outliersJM(all.gthres[,c('cd161.gate.high')], SDs = 3)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 2:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
# if(nrow(indices.check) != 0){
#   fixed.Gthres <- mean(all.gthres[,c('cd161.gate.high')])
#   all.gthres[indices.check[,1],c('cd161.gate.high')] <- fixed.Gthres
# }else{
#   print("No gating threshold outliers detected")
# }

###########################################################################################################

## Finding outliers and fixing gates for 'cd161.gate'
boxplot(all.gthres[,c('cd161.gate')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('cd161.gate')], SDs = 3)$ind, outliersJM(all.gthres[,c('cd161.gate')], SDs = 3)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
# if(nrow(indices.check) != 0){
#   fixed.Gthres <- mean(all.gthres[,c('cd161.gate')])
#   all.gthres[indices.check[,1],c('cd161.gate')] <- fixed.Gthres
# }else{
#   print("No gating threshold outliers detected")
# }

###########################################################################################################

## Finding outliers and fixing gates for 'cd4.gate.slant'
boxplot(all.gthres[,c('cd4.gate.slant')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('cd4.gate.slant')], SDs = 1.75)$ind, outliersJM(all.gthres[,c('cd4.gate.slant')], SDs = 1.75)$out)
indices.check.temp <- which(indices.check[,2] > 3)
indices.check <- indices.check[indices.check.temp,]
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
fcs.files <- store.allFCS[indices.check[,1],4]
for (k in 50:72){
#for (k in 4:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}

## Manually adding the files sent by Ania, which had issues with the Autofluorescence part
manual.index.autoFluor <- c('C4L000128593', 'L000130249', 'L000131380', 'L000131381', 'L000131386', 'L000131896', 'L000131907', 'L000132406',
                        'L000133190', 'L000133191', 'L000133194', 'L000133200', 'L000133201', 'L000135061', 'L000135062', 'L000135490',
                        'L000136634', 'L000136636', 'L000137233', 'L000137234', 'L000136631', 'L000136634', 'L000136636', 'L000137233',
                        'L000137234', 'L000133848')

manual.index.autoFluor <- unlist(sapply(1:length(manual.index.autoFluor), function(x){which(store.allFCS[,c('Barcodes')] == manual.index.autoFluor[x])}))
ind.with.store <- c(ind.with.store, manual.index.autoFluor)
indices.check <- cbind(ind.with.store,all.gthres[ind.with.store,c('cd4.gate.slant')])
if(nrow(indices.check) != 0){
  fixed.Gthres <- mean(all.gthres[,c('cd4.gate.slant')])
  all.gthres[indices.check[,1],c('cd4.gate.slant')] <- fixed.Gthres 
}else{
  print("No gating threshold outliers detected")
}

###########################################################################################################

## Finding outliers and fixing gates for 'tcrd.gate.slanta'
boxplot(all.gthres[,c('tcrd.gate.slanta')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('tcrd.gate.slanta')], SDs = 1.5)$ind, outliersJM(all.gthres[,c('tcrd.gate.slanta')], SDs = 1.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
fcs.files <- store.allFCS[indices.check[,1],4]
for (k in 1:20){
#for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}

## Manually adding the files sent by Ania, which had issues with the Autofluorescence part
manual.index.autoFluor <- c('C4L000128593', 'L000130249', 'L000131380', 'L000131381', 'L000131386', 'L000131896', 'L000131907', 'L000132406',
                                 'L000133190', 'L000133191', 'L000133194', 'L000133200', 'L000133201', 'L000135061', 'L000135062', 'L000135490',
                                 'L000136634', 'L000136636', 'L000137233', 'L000137234', 'L000136631', 'L000136634', 'L000136636', 'L000137233',
                                 'L000137234', 'L000133848')

manual.index.autoFluor <- unlist(sapply(1:length(manual.index.autoFluor), function(x){which(store.allFCS[,c('Barcodes')] == manual.index.autoFluor[x])}))
ind.with.store <- c(ind.with.store, manual.index.autoFluor)
indices.check <- cbind(ind.with.store,all.gthres[ind.with.store,c('tcrd.gate.slanta')])
if(nrow(indices.check) != 0){
  fixed.Gthres <- mean(all.gthres[,c('tcrd.gate.slanta')])
  all.gthres[indices.check[,1],c('tcrd.gate.slanta')] <- fixed.Gthres 
}else{
  print("No gating threshold outliers detected")
}


###########################################################################################################

## Finding outliers and fixing gates for 'tcrd.gate.slantb'
boxplot(all.gthres[,c('tcrd.gate.slantb')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('tcrd.gate.slantb')], SDs = 1.5)$ind, outliersJM(all.gthres[,c('tcrd.gate.slantb')], SDs = 1.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
fcs.files <- store.allFCS[indices.check[,1],4]
for (k in 1:20){
  #for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}

## Manually adding the files sent by Ania, which had issues with the Autofluorescence part
manual.index.autoFluor <- c('C4L000128593', 'L000130249', 'L000131380', 'L000131381', 'L000131386', 'L000131896', 'L000131907', 'L000132406',
                            'L000133190', 'L000133191', 'L000133194', 'L000133200', 'L000133201', 'L000135061', 'L000135062', 'L000135490',
                            'L000136634', 'L000136636', 'L000137233', 'L000137234', 'L000136631', 'L000136634', 'L000136636', 'L000137233',
                            'L000137234', 'L000133848')

manual.index.autoFluor <- unlist(sapply(1:length(manual.index.autoFluor), function(x){which(store.allFCS[,c('Barcodes')] == manual.index.autoFluor[x])}))
ind.with.store <- c(ind.with.store, manual.index.autoFluor)
indices.check <- cbind(ind.with.store,all.gthres[ind.with.store,c('tcrd.gate.slantb')])
if(nrow(indices.check) != 0){
  fixed.Gthres <- mean(all.gthres[,c('tcrd.gate.slantb')])
  all.gthres[indices.check[,1],c('tcrd.gate.slantb')] <- fixed.Gthres 
}else{
  print("No gating threshold outliers detected")
}

###########################################################################################################


## Finding outliers and fixing gates for 'klrg1.gate'
boxplot(all.gthres[,c('klrg1.gate')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('klrg1.gate')], SDs = 3)$ind, outliersJM(all.gthres[,c('klrg1.gate')], SDs = 3)$out)
indices.check.temp <- which(indices.check[,2] >= 2.8)
indices.check <- indices.check[indices.check.temp,]
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
# ## Manually adding the files sent by Ania, which had issues with the klrg1 part
# manual.index.klrg1 <- c('L000077024', 'L000136632', 'L000137229')
# 
# manual.index.klrg1 <- unlist(sapply(1:length(manual.index.klrg1), function(x){which(store.allFCS[,c('Barcodes')] == manual.index.klrg1[x])}))
# ind.with.store <- c(ind.with.store, manual.index.klrg1)
# indices.check <- cbind(ind.with.store,all.gthres[ind.with.store,c('klrg1.gate')])

if(length(ind.with.store) != 0){
  fixed.Gthres <- mean(all.gthres[,c('klrg1.gate')])
  all.gthres[indices.check[,1],c('klrg1.gate')] <- fixed.Gthres
  # #all.gthres[ind.with.store[47],c('klrg1.gate')] <- 1.5
  # all.gthres[indices.check[,1],c('klrg1.gate')] <- 1.5
} else{
  print("No gating threshold outliers detected")
}

###########################################################################################################


## Finding outliers and fixing gates for 'klrg1.gate.new'
boxplot(all.gthres[,c('klrg1.gate.new')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('klrg1.gate.new')], SDs = 4.5)$ind, outliersJM(all.gthres[,c('klrg1.gate.new')], SDs = 4.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  #for (k in 1:10){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}


if(length(ind.with.store) != 0){
  fixed.Gthres <- mean(all.gthres[,c('klrg1.gate.new')])
  all.gthres[indices.check[,1],c('klrg1.gate.new')] <- fixed.Gthres
} else{
  print("No gating threshold outliers detected")
}

###########################################################################################################


## Finding outliers and fixing gates for 'klrg1.gate.cd4posNKT'
boxplot(all.gthres[,c('klrg1.gate.cd4posNKT')], pch=19, cex=0.2)
#indices.check <- cbind(outliersJM(all.gthres[,c('klrg1.gate.cd4posNKT')], SDs = 2.25)$ind, outliersJM(all.gthres[,c('klrg1.gate.cd4posNKT')], SDs = 2.25)$out)
index <- which(all.gthres[,c('klrg1.gate.cd4posNKT')] > mean(all.gthres[,c('klrg1.gate.cd4posNKT')]))
indices.check <- cbind(index,all.gthres[index,c('klrg1.gate.cd4posNKT')])
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  #for (k in 1:10){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}


if(length(ind.with.store) != 0){
  fixed.Gthres <- mean(all.gthres[,c('klrg1.gate.cd4posNKT')])
  all.gthres[indices.check[,1],c('klrg1.gate.cd4posNKT')] <- fixed.Gthres
} else{
  print("No gating threshold outliers detected")
}



###########################################################################################################

## Finding outliers and fixing gates for 'cd25.gate.low'
boxplot(all.gthres[,c('cd25.gate.low')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('cd25.gate.low')], SDs = 4.75)$ind, outliersJM(all.gthres[,c('cd25.gate.low')], SDs = 4.75)$out)
indices.check.temp <- which(indices.check[,2] > 1)
indices.check <- indices.check[indices.check.temp,]
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
if(nrow(indices.check) != 0){
  fixed.Gthres <- mean(all.gthres[,c('cd25.gate.low')])
  all.gthres[indices.check[,1],c('cd25.gate.low')] <- fixed.Gthres 
} else{
  print("No gating threshold outliers detected")
}

###########################################################################################################

## Finding outliers and fixing gates for 'cd25.gate.high'
boxplot(all.gthres[,c('cd25.gate.high')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('cd25.gate.high')], SDs = 3)$ind, outliersJM(all.gthres[,c('cd25.gate.high')], SDs = 3)$out)
indices.check.temp <- which(indices.check[,2] > 2.5)
indices.check <- indices.check[indices.check.temp,]
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
#for (k in 1:30){
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
# ## Manually adding the files sent by Ania, which had issues with the cd25.high.gate part
# manual.index.cd25.high <- c('L000043315', 'L000043319', 'L000043320', 'L000044419', 'L000045487', 'L000076313', 'L000076314',
#                             'L000076315', 'L000077022', 'L000077028', 'L000077887', 'L000077888', 'L000077890', 'L000078874',
#                             'L000078875', 'L000078880', 'L000078881', 'L000078885', 'L000079581', 'L000079582', 'L000079595',
#                             'L000080323', 'L000076307', 'L000076308', 'L000136635', 'L000136637')
#
# manual.index.cd25.high <- unlist(sapply(1:length(manual.index.cd25.high), function(x){which(store.allFCS[,c('Barcodes')] == manual.index.cd25.high[x])}))
# ind.with.store <- c(ind.with.store, manual.index.cd25.high)
# indices.check <- cbind(ind.with.store,all.gthres[ind.with.store,c('cd25.gate.high')])

if(nrow(indices.check) != 0){
  fixed.Gthres <- mean(all.gthres[,c('cd25.gate.high')])
  all.gthres[indices.check[c(9,11,25,28,39:nrow(indices.check)),1],c('cd25.gate.high')] <- fixed.Gthres
  all.gthres[indices.check[c(5,33,35),1],c('cd25.gate.high')] <- 2.25
} else{
  print("No gating threshold outliers detected")
}


############################################################################################################

## Finding outliers and fixing gates for 'cd161.gate.low'
boxplot(all.gthres[,c('cd161.gate.low')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('cd161.gate.low')], SDs = 4)$ind, outliersJM(all.gthres[,c('cd161.gate.low')], SDs = 4)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}

# if(length(ind.with.store) != 0){
#   # fixed.Gthres <- mean(all.gthres[,c('cd161.gate.low')])
#   # all.gthres[ind.with.store[1:3],c('cd161.gate.low')] <- 2.50
#   all.gthres[ind.with.store,c('cd161.gate.low')] <- 2.15
# } else{
#   print("No gating threshold outliers detected")
# }

# ############################################################################################################

## Finding outliers and fixing gates for 'cd62.gate'
boxplot(all.gthres[,c('cd62.gate')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('cd62.gate')], SDs = 4.5)$ind, outliersJM(all.gthres[,c('cd62.gate')], SDs = 4.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for (k in 1:30){
#for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
# if(nrow(indices.check) != 0){
#   fixed.Gthres <- mean(all.gthres[,c('cd62.gate')])
#   all.gthres[indices.check[,1],c('cd62.gate')] <- fixed.Gthres
# } else{
#   print("No gating threshold outliers detected")
# }


#############################################################################################################

## Finding outliers and fixing gates for 'cd44.gate.high'
boxplot(all.gthres[,c('cd44.gate.high')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('cd44.gate.high')], SDs = 2)$ind, outliersJM(all.gthres[,c('cd44.gate.high')], SDs = 2)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])

for (k in 7:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
if(nrow(indices.check) != 0){
  fixed.Gthres <- mean(all.gthres[,c('cd44.gate.high')])
  all.gthres[indices.check[,1],c('cd44.gate.high')] <- fixed.Gthres
} else{
  print("No gating threshold outliers detected")
}


#############################################################################################################


## Finding outliers and fixing gates for 'tcrd.gate.low'
boxplot(all.gthres[,c('tcrd.gate.low')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('tcrd.gate.low')], SDs = 2.25)$ind, outliersJM(all.gthres[,c('tcrd.gate.low')], SDs = 2.25)$out)
indices.check.temp <- which(indices.check[,2] > 2)
indices.check <- indices.check[indices.check.temp,]
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
for(k in 1:20){
  #for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}

if(nrow(indices.check) != 0){
  fixed.Gthres <- mean(all.gthres[,c('tcrd.gate.low')])
  all.gthres[indices.check[,1],c('tcrd.gate.low')] <- fixed.Gthres 
  #all.gthres[manual.add.tcrd, c('tcrd.gate.low')] <- 1.5
} else{
  print("No gating threshold outliers detected")
}


#############################################################################################################


## Finding outliers and fixing gates for 'cd4.gate'
boxplot(all.gthres[,c('cd4.gate')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('cd4.gate')], SDs = 3)$ind, outliersJM(all.gthres[,c('cd4.gate')], SDs = 3)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])
#for(k in 1:11){
for (k in 13:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}

## Manually adding the files sent by Ania, which had issues with the cd25.high.gate part
manual.index.cd4 <- c('L000133194', 'L000133200', 'L000133201')
manual.index.cd4 <- unlist(sapply(1:length(manual.index.cd4), function(x){which(store.allFCS[,c('Barcodes')] == manual.index.cd4[x])}))
ind.with.store <- c(ind.with.store, manual.index.cd4)
indices.check <- cbind(ind.with.store,all.gthres[ind.with.store,c('cd4.gate')])


if(nrow(indices.check) != 0){
  fixed.Gthres <- mean(all.gthres[,c('cd4.gate')])
  all.gthres[indices.check[c(2:3, 11:nrow(indices.check)),1],c('cd4.gate')] <- fixed.Gthres 
  
} else{
  print("No gating threshold outliers detected")
}


#############################################################################################################

## Finding outliers and fixing gates for 'cd5.gate'
boxplot(all.gthres[,c('cd5.gate')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('cd5.gate')], SDs = 3)$ind, outliersJM(all.gthres[,c('cd5.gate')], SDs = 3)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])

#for (k in 1:10){
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
if(nrow(indices.check) != 0){
  fixed.Gthres <- mean(all.gthres[,c('cd5.gate')])
  all.gthres[indices.check[c(1,20,23),1],c('cd5.gate')] <- fixed.Gthres
  all.gthres[indices.check[c(2:9,26,27,31),1],c('cd5.gate')] <- 2.5
} else{
  print("No gating threshold outliers detected")
}


#############################################################################################################

## Finding outliers and fixing gates for 'cd5.gate.new'
boxplot(all.gthres[,c('cd5.gate.new')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('cd5.gate.new')], SDs = 2.5)$ind, outliersJM(all.gthres[,c('cd5.gate.new')], SDs = 2.5)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])

#for (k in 1:35){
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
if(nrow(indices.check) != 0){
  fixed.Gthres <- mean(all.gthres[,c('cd5.gate.new')])
  all.gthres[indices.check[c(3,4,16,17,23:29,35,40:42,51 ),1],c('cd5.gate.new')] <- fixed.Gthres
} else{
  print("No gating threshold outliers detected")
}


#############################################################################################################

## Finding outliers and fixing gates for 'cd8.gate'
boxplot(all.gthres[,c('cd8.gate')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('cd8.gate')], SDs = 3)$ind, outliersJM(all.gthres[,c('cd8.gate')], SDs = 3)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])

#for (k in 1:35){
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
# if(nrow(indices.check) != 0){
#   fixed.Gthres <- mean(all.gthres[,c('cd8.gate')])
#   all.gthres[indices.check[,1],c('cd8.gate')] <- fixed.Gthres
# } else{
#   print("No gating threshold outliers detected")
# }

#############################################################################################################

## Finding outliers and fixing gates for 'gitr.gate'
boxplot(all.gthres[,c('gitr.gate')], pch=19, cex=0.2)
indices.check <- cbind(outliersJM(all.gthres[,c('gitr.gate')], SDs = 3)$ind, outliersJM(all.gthres[,c('gitr.gate')], SDs = 3)$out)
indices.check <- indices.check[order(indices.check[,1]),]
ind.with.store <- unique(indices.check[,1])

#for (k in 1:35){
for (k in 1:length(ind.with.store)){
  name <- paste0("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_T-cell/MLN/Results/Figures/ScatterPlots/",
                 store.allFCS[ind.with.store[k],"Genotype"], "/Total_", store.allFCS[ind.with.store[k],"FCS files"], ".png")
  x11(width = 16, height = 9)
  img <- readPNG(name)
  grid::grid.raster(img)
}
# if(nrow(indices.check) != 0){
#   fixed.Gthres <- mean(all.gthres[,c('gitr.gate')])
#   all.gthres[indices.check[,1],c('gitr.gate')] <- fixed.Gthres
# } else{
#   print("No gating threshold outliers detected")
# }

###########################################################################################################
###########################################################################################################

all.gthres.Updated <- all.gthres

## Saving the Updated Gating Thresholds as numerical data, since we will need this for running the gating again with updated thresholds
save(all.gthres.Updated, file = paste0(results.dir, "/all.gthres.Updated.Rdata"))
all.gthres.Updated <- cbind(store.allFCS[,c('Panel/Organ', 'Genotype', 'FCS files', 'Barcodes', 'Assay Date', 'Gender', 'Number of Channels', 'Number of Cells')], all.gthres.Updated)

write.csv(all.gthres.Updated, file =  paste0(results.dir, "/allGthres_Table_Updated.csv"), row.names = FALSE)



######################################################################

