# Originally written by Nima Aghaeepour (UFO package) -------------------
# Re-edited by Albina Rahim (August 15, 2016)
## In this version of the UFO script instead of using flowType results on 1000 samples and for method "kmeans", we use the flowType objects, which we created as part of our unsupervised analysis.
##  flowType.res <- flowType(Frame = cd45.NOTP2 , PropMarkers = as.vector(channels.ind.NoLiveNoCD45), MaxMarkersPerPop = 8, PartitionsPerMarker=2, Methods='Thresholds', Thresholds= listGthres.NoLiveNoCD45, verbose=F, MemLimit=400)
## cd45 population is used as input, 10 markers are used (CD45 and Live marker is excluded) and we use out gating thresholds.

remove(list=ls())

# This code uses a flowType .Rdata file for each .fcs file and then runs isomap.

#setwd("/code/Projects/3i/Panel_T-cell/")
setwd("/data/Panel_M-cell")

library(flowCore)
library(flowType)

# ##Load the meta_data and initialize stuff
# NumberOfCPUs = 5 ##for snow
# load("Results/meta_data.Rdata")
# load("Results/File_Names.Rdata")
# load("Results/WT_dates.Rdata")
# load("Results/WT_dates3.Rdata")
# load("Results/channels.ind.NoLiveNoCD45.Rdata")
# 
# # metadata <- read.csv('BatchDayData48.csv')
# # files=metadata[,1]
# # datapath='~/Desktop/UFO/GvHD_Res'
# # datapath <- '/code/Projects/3i/Panel_T-cell/FCS_Groups'
# #datapath='/code/Projects/IMPC/data/H/FR-FCM-ZYCB-WildType/'
# 
# # temp.fcs <- read.FCS("/code/Projects/IMPC/data/CIPHE/1-1-1-B6/Panel1/14-Aug-19_IMPC1_01_labelled.fcs")
# #temp.fcs <- read.FCS("/code/Projects/3i/Panel_T-cell/FCS_Groups/1700001C02Rik_1700001C02Rik/SPLN_L000060557_D01_008.labelled.fcs")
# 
load(file = paste("Results_SPLEEN/After_Clean/", props.events.gates$Genotype[1], "/AftClean_", props.events.gates$FCS.file[1], ".Rdata", sep=""))
temp.fcs <- f
marker.names <- temp.fcs@parameters@data$name
prop.markers <- c(7:18) #as.integer(channels.ind.NoLiveNoCD45)
# 
# N=nrow(meta_data)
# CSN=1
# load(paste0("Results/FlowType/",meta_data[CSN,1],"/FT_",File_Names[CSN],".Rdata"))
# M = length(flowType.res@CellFreqs)
# 
# ##Load the flowType results
# NotAvailable <- vector();
# All.Proportions.Modified <- matrix(0, M, N)
# for (CSN in 1:N){
#   print(CSN)
#   if (file.exists(paste0("Results/FlowType/",meta_data[CSN,1],"/FT_",File_Names[CSN],".Rdata"))){
#     load(paste0("Results/FlowType/",meta_data[CSN,1],"/FT_",File_Names[CSN],".Rdata"))
#   }else{
#     print(paste0("Results/FlowType/",meta_data[CSN,1],"/FT_",File_Names[CSN],".Rdata doesnot exist"))
#     NotAvailable <- c(NotAvailable, CSN)
#     next;
#   }
#   All.Proportions.Modified[,CSN] <- flowType.res@CellFreqs/max(flowType.res@CellFreqs);
# }
# 
# Phenotypes <- unlist(lapply(flowType.res@PhenoCodes, function(x){return(decodePhenotype(
#   x,flowType.res@MarkerNames[prop.markers], flowType.res@PartitionsPerMarker))}))
# 
# rownames(All.Proportions.Modified) <- Phenotypes

# library(RDRToolbox)
# is=Isomap(t(All.Proportions.Modified),2,10)$dim2
# td=is

# SD test 12-12-2016 
# See what Isomap returns when I input only the gating strategy files
load(file = "Results_SPLEEN/PropsEventsGates.Rdata")  
# Don't include singlets/live/lymphocytes/CD45 proportions, also remove NA
# Actually did not make too much difference if I leave out singlets/live/lymphocytes/CD45 or only singlets/live.
All.Proportions.Modified <- props.events.gates[-which(is.na(props.events.gates[,17])), 10:24] 
is=Isomap(as.matrix(All.Proportions.Modified), 2, 10)$dim2
png('IsoMap.png',  height=900, width=1200)
plot(is, pch = 20, xlab='UFO Dimension 1',ylab='UFO Dimension 2', xlim = c(-55, 200))
idx <- union(which(is[,1] > 90), which(is[,2] > 52.5))
idx <- idx[-which(idx == '110')]
fcsfiles <- props.events.gates[-which(is.na(props.events.gates[,17])),1:2]
labels <- paste0(fcsfiles$Genotype[idx], "/", fcsfiles$FCS.file[idx])
text(x = is[idx, 1], y = is[idx, 2], labels = labels, pos = 4, cex = 0.8)
# Barcodes for which the Lin- cell prop is > 3 std deviations from mean 
# The Lin- population (Myeloid) seems probelematic based on Ania's FP test
idx2 <-  grep('47005|49258|51710|64377|116053|116052', fcsfiles$FCS.file)
points(is[idx2,], pch = 20, col = 'red')
dev.off()

# Write some of the "outliers" to file
idx <- which(is[,1] > 75)
idx <- idx[-which(idx == '110')]
Barcodes <- str_extract(fcsfiles$FCS.file[idx],"L[0-9]+")
df <- data.frame(Barcodes, Genotype = fcsfiles$Genotype[idx], FCS.file = fcsfiles$FCS.file[idx],
                  Isomap.d1.value = is[idx,1], Isomap.d2.value = is[idx,2])
df <- df[order(-df$Isomap.d1.value),]
write.table(df, file = "Results_SPLEEN/Isomap_outliers.csv", sep ='\t', col.names = TRUE, row.names = FALSE)

# Compare with tSNE
library(Rtsne)
tsne.results <- Rtsne(All.Proportions.Modified, k = 2, initial_dims = 15, theta = 0, check_duplicates = FALSE, perplexity = 30)
png('tSNE.png',  height=900, width=1200)
plot(tsne.results$Y, pch = 20)
#idx <- which((tsne.results$Y[,2] > 6)  & (tsne.results$Y[,1] < -31))
idx <- which((tsne.results$Y[,1] < -46))
points(tsne.results$Y[idx,], pch = 20, col = 'gray')
idx2 <-  grep('47005|49258|51710|64377|116053|116052', fcsfiles$FCS.file)
points(tsne.results$Y[idx2,], pch = 20, col = 'red')
dev.off()
# So now what does PC analysis look like?
#library(ggfortify)
png('PCA.png',  height=1200, width=1200)
rownames(All.Proportions.Modified) <- fcsfiles$FCS.file
prcomp.results <- prcomp(All.Proportions.Modified, scale. = TRUE)
#autoplot(prcomp.results, label = FALSE, shape = TRUE)
plot(prcomp.results$x, pch = 20, xlim = c(-6, 20))
idx <- union(which(prcomp.results$x[,'PC1'] > 12), which(prcomp.results$x[,'PC2'] > 8))
labels <- paste0(fcsfiles$Genotype[idx], "/", fcsfiles$FCS.file[idx])
text(x = prcomp.results$x[idx, 1], y = prcomp.results$x[idx, 2], labels = labels, pos = 4, cex = 0.8)
idx2 <-  grep('47005|49258|51710|64377|116053|116052', fcsfiles$FCS.file)
points(prcomp.results$x[idx2,], pch = 20, col = 'red')
dev.off()
# end of SD test

td=is

load("Results_SPLEEN/CSVfile.Rdata")
library('stringr')
Barcodes <- str_extract(fcsfiles$FCS.file,"L[0-9]+")
Mouse_Label <-  CSVfile[, 'Label.Barcode'] 
# Thus when I get the indices for the barcodes I just arbitrarily pick the largest because as far as I can tell
# the rows in CSVfile pertaining to one barcode are the same
indices <- sapply(1:length(Barcodes), function(xdat){current.index <- max(grep(Barcodes[xdat], Mouse_Label))})
fcsfiles$Assay.Date <- CSVfile[indices,'Assay.Date']
fcsfiles$Gender <- CSVfile[indices,'Gender']

days <- as.Date(as.Date(fcsfiles$Assay.Date, format = "%d-%b-%Y"))
# Need to correct this for times when the year is not 4 digits, but two instead
for(i in 1:length(fcsfiles$Assay.Date)){
  if(nchar(strsplit(fcsfiles$Assay.Date[i], '-')[[1]][3]) == 2){
    days[i] <- as.Date(fcsfiles$Assay.Date[i], format = "%d-%b-%y")
  }
}
days <- days - min(days)
days <- as.numeric(days)
ts = days

##Make the surface marker annotations. If you want to highlight other populations add them to the "highlights" vector below.
legend.size=1.25
library(gplots)
highlights=sprintf('%s+',marker.names[prop.markers])

suppressWarnings ( dir.create ( "Results_SPLEEN/UFO_Modified_Figures") )
# suppressWarnings ( dir.create ( "Results_SPLEEN/UFO_Modified_Figures/Annotations") )
# 
# xlim=NULL#c(-101.78373,   82.10356)
# ylim=NULL#c(-43.00999,  53.04831)
# for (i in 1:length(highlights)){
#     ylab=sprintf('%% %s',highlights[i])
#     # pdf(sprintf('figs/annotations/%s.pdf',marker.names[prop.markers[i]]))
#     png(sprintf('Results_SPLEEN/UFO_Modified_Figures/Annotations/%s.png',gsub("/", "-", marker.names[prop.markers[i]])), height=900, width=900)
#     palette((rich.colors(1000)))
#     ##par(bg='black')
#     delta.sizex=legend.size/dev.size()[1]
#     split.screen(t(cbind(c(0,1-delta.sizex,0,1),c(1-delta.sizex,1,0,1))))
#     screen(1)
#     ##highlights[i]=Names[which.max(ratsauc[,i])]
#     library(gplots)
#     ts=(All.Proportions.Modified[which(Phenotypes==highlights[i]),])
#     ts=ts-min(ts)
#     temp=quantile(ts,0.99,na.rm=TRUE)
#     ts[which(ts>temp)]=temp
#     temp=quantile(ts,0.01,na.rm=TRUE)
#     ts[which(ts<temp)]=temp
#     cols=round(ts*1000/max(ts,na.rm=TRUE))
#     chulls=list()
#     plot(td,pch=20,col=cols,cex=1.5,axes=FALSE,col.lab=par('fg'),xlim=xlim,ylim=ylim,xlab='UFO Dimension 1',ylab='UFO Dimension 2')
#     axis(1,col=par('fg'));axis(2,col=par('fg'));
#     title(main=highlights[i],col.main=par('fg'))
#     palette('default')
#     screen(2)
#     palette((rich.colors(1000)))
#     par(mar=c(1,4,1,0.2))
#     image(matrix(1:2500, 25), col = 1:1000, xaxt='n', ylab='', yaxt='n',col.lab=8)
#     par(mgp=c(2.5,1,0))
#     title(ylab=ylab,col.lab=par('fg'));
#     myvar=(All.Proportions.Modified[which(Phenotypes==highlights[i]),])
#     axis(2, at=c((1/(length(pretty(myvar))-1))*(0:(length(pretty(myvar))-1))), labels=100*pretty(myvar),col=par('fg'),col.axis=par('fg'))
#     palette('default')
#     close.screen(all.screens=TRUE)
#     dev.off()
# }



## Plot based on days
# days <- as.numeric(meta_data[,4])
# days <- days - min(days)
# pdf(sprintf('figs/days.pdf'))
png(sprintf('Results_SPLEEN/UFO_Modified_Figures/days.png'), height=900, width=900)
ylab='days'
legend.size=1.25
delta.sizex=legend.size/dev.size()[1]
# xlim=NULL
# ylim=NULL
# library(gplots)
split.screen(t(cbind(c(0,1-delta.sizex,0,1),c(1-delta.sizex,1,0,1))))
screen(1)
palette((rich.colors(1000)))
##par(bg='black')
##highlights[i]=Names[which.max(ratsauc[,i])]
#library(gplots)
#days <- as.numeric(days)
#ts = days
#ts=ts-min(ts)
temp=quantile(ts,0.99,na.rm=TRUE)
ts[which(ts>temp)]=temp
temp=quantile(ts,0.01,na.rm=TRUE)
ts[which(ts<temp)]=temp
cols=round(ts*1000/max(ts,na.rm=TRUE))
chulls=list()
plot(td,pch=20,col=cols+1,cex=1.5,axes=FALSE,xlab='UFO Dimension 1',ylab='UFO Dimension 2',col.lab='black') #,xlim=xlim,ylim=ylim)
axis(1,col='black');axis(2,col='black');
title(main='Days since first assay')
screen(2)
par(mar=c(1,4,1,0.2))
image(matrix(1:2500, 25), col = 1:1000, xaxt='n', ylab='', yaxt='n',col.lab=1)
par(mgp=c(2.5,1,0))
title(ylab=ylab,col.lab=1);
axis(2, at=c((1/(length(pretty(days))-1))*(0:(length(pretty(days))-1))), labels=pretty(days),col=1,col.axis=1)
palette('default')
close.screen()
dev.off()

# UFO_days <- cbind(WT_dates, days, td, meta_data[,5], meta_data[,6])
# colnames(UFO_days) <- c("Date", "Days since first assay", "X", "Y", "WT=0, KO=1", "F=0, M=1")
# write.table(UFO_days, file="Results/UFO_Modified_Figures/UFO_days.csv", sep=",", row.names = F)

## Plot based on gender
# groups=as.numeric(meta_data[,6])
# groups <- groups - min(groups)
png(sprintf('Results_SPLEEN/UFO_Modified_Figures/gender.png'), height=900, width=900)
ylab='gender'
legend.size=1.25
delta.sizex=legend.size/dev.size()[1]
xlim=NULL
ylim=NULL
library(gplots)
split.screen(t(cbind(c(0,1-delta.sizex,0,1),c(1-delta.sizex,1,0,1))))
screen(1)
palette((rich.colors(1000)))
##par(bg='black')
##highlights[i]=Names[which.max(ratsauc[,i])]
#library(gplots)
#ts=(groups)
groups <- as.numeric(factor(fcsfiles$Gender)) - 1
ts <- groups
#ts=ts-min(ts)
temp=quantile(ts,0.99,na.rm=TRUE)
ts[which(ts>temp)]=temp
temp=quantile(ts,0.01,na.rm=TRUE)
ts[which(ts<temp)]=temp
cols=round(ts*1000/max(ts,na.rm=TRUE))
chulls=list()
cols[which(cols == 1000)] <- cols[which(cols == 1000)] - 1
plot(td,pch=20,col=cols+1,cex=1.5,axes=FALSE,xlab='UFO Dimension 1',ylab='UFO Dimension 2',col.lab='black',xlim=xlim,ylim=ylim)
axis(1,col='black');axis(2,col='black');
title(main='Gender: F 0 and M 1')
screen(2)
par(mar=c(1,4,1,0.2))
image(matrix(1:2500, 25), col = 1:1000, xaxt='n', ylab='', yaxt='n',col.lab=1)
par(mgp=c(2.5,1,0))
title(ylab=ylab,col.lab=1);
axis(2, at=c((1/(length(pretty(groups))-1))*(0:(length(pretty(groups))-1))), labels=pretty(groups),col=1,col.axis=1)
palette('default')
close.screen()
dev.off()


## Plot based on Wild-type and Knockout group
# groups=as.numeric(meta_data[,5])
# groups <- groups - min(groups)
groups <- as.numeric((fcsfiles$Genotype != '+_+')&(fcsfiles$Genotype != '+_Y')) 
ts <- groups
png(sprintf('Results_SPLEEN/UFO_Modified_Figures/WT_KO.png'), height=900, width=900)
ylab='group'
legend.size=1.25
delta.sizex=legend.size/dev.size()[1]
xlim=NULL
ylim=NULL
library(gplots)
split.screen(t(cbind(c(0,1-delta.sizex,0,1),c(1-delta.sizex,1,0,1))))
screen(1)
palette((rich.colors(1000)))
##par(bg='black')
##highlights[i]=Names[which.max(ratsauc[,i])]
library(gplots)
# ts=(groups)
# ts=ts-min(ts)
temp=quantile(ts,0.99,na.rm=TRUE)
ts[which(ts>temp)]=temp
temp=quantile(ts,0.01,na.rm=TRUE)
ts[which(ts<temp)]=temp
cols=round(ts*1000/max(ts,na.rm=TRUE))
chulls=list()
cols[which(cols == 1000)] <- cols[which(cols == 1000)] - 1
plot(td,pch=20,col=cols+1,cex=1.5,axes=FALSE,xlab='UFO Dimension 1',ylab='UFO Dimension 2',col.lab='black',xlim=xlim,ylim=ylim)
axis(1,col='black');axis(2,col='black');
title(main='Group: WT 0 and KO 1')
screen(2)
par(mar=c(1,4,1,0.2))
image(matrix(1:2500, 25), col = 1:1000, xaxt='n', ylab='', yaxt='n',col.lab=1)
par(mgp=c(2.5,1,0))
title(ylab=ylab,col.lab=1);
axis(2, at=c((1/(length(pretty(groups))-1))*(0:(length(pretty(groups))-1))), labels=pretty(groups),col=1,col.axis=1)
palette('default')
close.screen()
dev.off()

# 
# ###################################################################################################
# ## Finding the outliers in the UFO plots. This part of the script is to be carried out after looking into the UFO plots for days.
# ## Then use the td matrix (N x 2), which is the output from the Isomap() function to determine which of the files are outliers.
# ## Which dimension (column in the td matrix) you are using and which value you are using to check for outliers will depend on the plot.
# meta_data.Unique <- meta_data
# rownames(meta_data.Unique) <- 1:dim(meta_data.Unique)[1]
# outlier.Index <- c(which(td[,2] < -5), which(td[,2] > 5)) 
# #outlier.Index <- which(td[,2] > 5) 
# ## Make sure after finding the outlier.Index look into the td matrix bak to look into the dimension values.
# outlier.List <- meta_data.Unique[outlier.Index,]
# write.table(outlier.List, file="Results/UFO_Modified_Figures/outlierList.csv", sep=",", row.names = F)
# ###################################################################################################
# 
# ##  UFO Analysis Plots for WTs only
# remove(list=ls())
# 
# # This code uses a flowType .Rdata file for each .fcs file and then runs isomap.
# 
# setwd("/code/Projects/3i/Panel_T-cell/")
# library(flowCore)
# library(flowType)
# 
# ## Load the meta_data and initialize stuff
# NumberOfCPUs = 5 ##for snow
# load("Results/meta_dataWT.Rdata")
# load("Results/File_NamesWT.Rdata")
# load("Results/WT_dates.Rdata")
# load("Results/WT_dates3.Rdata")
# load("Results/channels.ind.NoLiveNoCD45.Rdata")
# 
# 
# # temp.fcs <- read.FCS("/code/Projects/IMPC/data/CIPHE/1-1-1-B6/Panel1/14-Aug-19_IMPC1_01_labelled.fcs")
# temp.fcs <- read.FCS("/code/Projects/3i/Panel_T-cell/FCS_Groups/+_+/SPLN 7,2f,20,2f,15_L000092223_024.labelled.fcs")
# marker.names <- temp.fcs@parameters@data$name
# prop.markers <- as.integer(channels.ind.NoLiveNoCD45)
# 
# N=nrow(meta_dataWT)
# # N <- 44
# CSN=1
# load(paste0("Results/FlowType/",meta_dataWT[CSN,1],"/FT_",File_NamesWT[CSN],".Rdata"))
# M = length(flowType.res@CellFreqs)
# 
# ##Load the flowType results
# NotAvailable <- vector();
# All.Proportions.Modified.WT <- matrix(0, M, N)
# for (CSN in 1:N){
#   print(CSN)
#   if (file.exists(paste0("Results/FlowType/",meta_dataWT[CSN,1],"/FT_",File_NamesWT[CSN],".Rdata"))){
#     load(paste0("Results/FlowType/",meta_dataWT[CSN,1],"/FT_",File_NamesWT[CSN],".Rdata"))
#   }else{
#     print(paste0("Results/FlowType/",meta_dataWT[CSN,1],"/FT_",File_NamesWT[CSN],".Rdata does not exist"))
#     NotAvailable <- c(NotAvailable, CSN)
#     next;
#   }
#   All.Proportions.Modified.WT[,CSN] <- flowType.res@CellFreqs/max(flowType.res@CellFreqs);
# }
# 
# PhenotypesWT <- unlist(lapply(flowType.res@PhenoCodes, function(x){return(decodePhenotype(
#   x,flowType.res@MarkerNames[prop.markers], flowType.res@PartitionsPerMarker))}))
# 
# rownames(All.Proportions.Modified.WT) <- PhenotypesWT
# 
# 
# library(RDRToolbox)
# is=Isomap(t(All.Proportions.Modified.WT),2,10)$dim2
# td=is
# 
# 
# ##Make the surface marker annotations. If you want to highlight other populations add them to the "highlights" vector below.
# legend.size=1.25
# library(gplots)
# highlights=sprintf('%s+',marker.names[prop.markers])
# 
# suppressWarnings ( dir.create ( "Results/UFO_Modified_WT_Figures") )
# suppressWarnings ( dir.create ( "Results/UFO_Modified_WT_Figures/Annotations") )
# 
# xlim=NULL#c(-101.78373,   82.10356)
# ylim=NULL#c(-43.00999,  53.04831)
# for (i in 1:length(highlights)){
#   ylab=sprintf('%% %s',highlights[i])
#   # pdf(sprintf('figs/annotations/%s.pdf',marker.names[prop.markers[i]]))
#   png(sprintf('Results/UFO_Modified_WT_Figures/Annotations/%s.png',gsub("/", "-", marker.names[prop.markers[i]])), height=900, width=900)
#   palette((rich.colors(1000)))
#   ##par(bg='black')
#   delta.sizex=legend.size/dev.size()[1]
#   split.screen(t(cbind(c(0,1-delta.sizex,0,1),c(1-delta.sizex,1,0,1))))
#   screen(1)
#   ##highlights[i]=Names[which.max(ratsauc[,i])]
#   library(gplots)
#   ts=(All.Proportions.Modified.WT[which(PhenotypesWT==highlights[i]),])
#   ts=ts-min(ts)
#   temp=quantile(ts,0.99,na.rm=TRUE)
#   ts[which(ts>temp)]=temp
#   temp=quantile(ts,0.01,na.rm=TRUE)
#   ts[which(ts<temp)]=temp
#   cols=round(ts*1000/max(ts,na.rm=TRUE))
#   chulls=list()
#   plot(td,pch=20,col=cols,cex=1.5,axes=FALSE,col.lab=par('fg'),xlim=xlim,ylim=ylim,xlab='UFO Dimension 1',ylab='UFO Dimension 2')
#   axis(1,col=par('fg'));axis(2,col=par('fg'));
#   title(main=highlights[i],col.main=par('fg'))
#   palette('default')
#   screen(2)
#   palette((rich.colors(1000)))
#   par(mar=c(1,4,1,0.2))
#   image(matrix(1:2500, 25), col = 1:1000, xaxt='n', ylab='', yaxt='n',col.lab=8)
#   par(mgp=c(2.5,1,0))
#   title(ylab=ylab,col.lab=par('fg'));
#   myvar=(All.Proportions.Modified.WT[which(PhenotypesWT==highlights[i]),])
#   axis(2, at=c((1/(length(pretty(myvar))-1))*(0:(length(pretty(myvar))-1))), labels=100*pretty(myvar),col=par('fg'),col.axis=par('fg'))
#   palette('default')
#   close.screen(all.screens=TRUE)
#   dev.off()
# }
# 
# 
# 
# ## Plot based on days
# days <- as.numeric(meta_dataWT[,4])
# days <- days - min(days)
# # pdf(sprintf('figs/days.pdf'))
# png(sprintf('Results/UFO_Modified_WT_Figures/daysWT.png'), height=900, width=900)
# ylab='days'
# legend.size=1.25
# delta.sizex=legend.size/dev.size()[1]
# xlim=NULL
# ylim=NULL
# library(gplots)
# split.screen(t(cbind(c(0,1-delta.sizex,0,1),c(1-delta.sizex,1,0,1))))
# screen(1)
# palette((rich.colors(1000)))
# library(gplots)
# ts = days
# #ts=ts-min(ts)
# temp=quantile(ts,0.99,na.rm=TRUE)
# ts[which(ts>temp)]=temp
# temp=quantile(ts,0.01,na.rm=TRUE)
# ts[which(ts<temp)]=temp
# cols=round(ts*1000/max(ts,na.rm=TRUE))
# chulls=list()
# plot(td,pch=20,col=cols+1,cex=1.5,axes=FALSE,xlab='UFO Dimension 1',ylab='UFO Dimension 2',col.lab='black',xlim=xlim,ylim=ylim)
# axis(1,col='black');axis(2,col='black');
# title(main='Days since first assay')
# screen(2)
# par(mar=c(1,4,1,0.2))
# image(matrix(1:2500, 25), col = 1:1000, xaxt='n', ylab='', yaxt='n',col.lab=1)
# par(mgp=c(2.5,1,0))
# title(ylab=ylab,col.lab=1);
# axis(2, at=c((1/(length(pretty(days))-1))*(0:(length(pretty(days))-1))), labels=pretty(days),col=1,col.axis=1)
# palette('default')
# close.screen()
# dev.off()
# 
# UFO_WT_days <- cbind(meta_dataWT[,3], days, td, meta_dataWT[,5], meta_dataWT[,6])
# colnames(UFO_WT_days) <- c("Date", "Days since first assay", "X", "Y", "WT=0, KO=1", "F=0, M=1")
# write.table(UFO_WT_days, file="Results/UFO_Modified_WT_Figures/UFO_WT_days.csv", sep=",", row.names = F)
# 
# ## Plot based on gender
# groups=as.numeric(meta_dataWT[,6])
# groups <- groups - min(groups)
# png(sprintf('Results/UFO_Modified_WT_Figures/genderWT.png'), height=900, width=900)
# ylab='gender'
# legend.size=1.25
# delta.sizex=legend.size/dev.size()[1]
# xlim=NULL
# ylim=NULL
# library(gplots)
# split.screen(t(cbind(c(0,1-delta.sizex,0,1),c(1-delta.sizex,1,0,1))))
# screen(1)
# palette((rich.colors(1000)))
# ##par(bg='black')
# ##highlights[i]=Names[which.max(ratsauc[,i])]
# library(gplots)
# ts=(groups)
# #ts=ts-min(ts)
# temp=quantile(ts,0.99,na.rm=TRUE)
# ts[which(ts>temp)]=temp
# temp=quantile(ts,0.01,na.rm=TRUE)
# ts[which(ts<temp)]=temp
# cols=round(ts*1000/max(ts,na.rm=TRUE))
# chulls=list()
# cols[which(cols == 1000)] <- cols[which(cols == 1000)] - 1
# plot(td,pch=20,col=cols+1,cex=1.5,axes=FALSE,xlab='UFO Dimension 1',ylab='UFO Dimension 2',col.lab='black',xlim=xlim,ylim=ylim)
# axis(1,col='black');axis(2,col='black');
# title(main='Gender: F 0 and M 1')
# screen(2)
# par(mar=c(1,4,1,0.2))
# image(matrix(1:2500, 25), col = 1:1000, xaxt='n', ylab='', yaxt='n',col.lab=1)
# par(mgp=c(2.5,1,0))
# title(ylab=ylab,col.lab=1);
# axis(2, at=c((1/(length(pretty(groups))-1))*(0:(length(pretty(groups))-1))), labels=pretty(groups),col=1,col.axis=1)
# palette('default')
# close.screen()
# dev.off()
# 
# ###################################################################################################
# ## Finding the outliers in the UFO plots for Wild types. This part of the script is to be carried out after looking into the WT UFO plots for days.
# ## Then use the td matrix (N x 2), which is the output from the Isomap() function to determine which of the files are outliers.
# ## Which dimension (column in the td matrix) you are using and which value you are using to check for outliers will depend on the plot.
# meta_dataWT.Unique <- meta_dataWT
# rownames(meta_dataWT.Unique) <- 1:dim(meta_dataWT.Unique)[1]
# outlierWT.Index <- c(which(td[,2] < -3), which(td[,2] > 3))
# ## Make sure after finding the outlier.Index look into the td matrix bak to look into the dimension values.
# outlierWT.List <- meta_dataWT.Unique[outlierWT.Index,]
# write.table(outlierWT.List, file="Results/UFO_Modified_WT_Figures/outlierWTList.csv", sep=",", row.names = F)
# 
