# Originally written by Nima Aghaeepour (UFO package) -------------------
# Re-edited by Albina Rahim (July 11, 2016)

remove(list=ls())

# This code uses a flowType .Rdata file for each .fcs file and then runs isomap.

setwd("/code/Projects/3i/Panel_T-cell_MLN/")

##Load the meta_data and initialize stuff
NumberOfCPUs = 5 ##for snow
load("Results/meta_data.Rdata")
load("Results/File_Names.Rdata")
load("Results/WT_dates.Rdata")
load("Results/WT_dates3.Rdata")

datapath <- '/code/Projects/3i/Panel_T-cell_MLN/FCS_Groups'


library(flowCore)
# temp.fcs <- read.FCS("/code/Projects/IMPC/data/CIPHE/1-1-1-B6/Panel1/14-Aug-19_IMPC1_01_labelled.fcs")
temp.fcs <- read.FCS("/code/Projects/3i/Panel_T-cell_MLN/FCS_Groups/1700001C02Rik_1700001C02Rik/MLN_L000060574_D04_025.labelled.fcs")
marker.names <- temp.fcs@parameters@data$name
prop.markers=c(7:18)


N=nrow(meta_data)
CSN=1
load(paste0("Results/UFO/",File_Names[CSN],".Rdata"))
M = length(FT@CellFreqs)

##Load the flowType results
NotAvailable <- vector();
All.Proportions <- matrix(0, M, N)
for (CSN in 1:N){
  print(CSN)
  if (file.exists(paste0("Results/UFO/",File_Names[CSN],".Rdata"))){
    load(paste0("Results/UFO/",File_Names[CSN],".Rdata"))
  }else{
    print(paste0("Results/UFO/",File_Names[CSN],".Rdata Does Not Exist"))
    NotAvailable <- c(NotAvailable, CSN)
    next;
  }
  All.Proportions[,CSN] <- FT@CellFreqs/max(FT@CellFreqs);
}

library(flowType)
library(Rmpi)
##Calculate phenotype names just in case we need them later.
tempfunc <- function(x){return(decodePhenotype(x,FT@MarkerNames[prop.markers],FT@PartitionsPerMarker))}
library(snow)
cl <- makeCluster(NumberOfCPUs, type = "SOCK",outfile='/tmp/snowout.log')
clusterExport(cl, list=c('decodePhenotype','FT','prop.markers'))
Names=unlist(clusterApplyLB(cl, FT@PhenoCodes, tempfunc));
stopCluster(cl)
rownames(All.Proportions) <- Names

library(RDRToolbox)
is=Isomap(t(All.Proportions),2,10)$dim2
td=is


##Make the surface marker annotations. If you want to highlight other populations add them to the "highlights" vector below.
legend.size=1.25
library(gplots)
highlights=sprintf('%s+',marker.names[prop.markers])

suppressWarnings ( dir.create ( "Results/UFO_Figures") )
suppressWarnings ( dir.create ( "Results/UFO_Figures/Annotations") )

xlim=NULL#c(-101.78373,   82.10356)
ylim=NULL#c(-43.00999,  53.04831)
for (i in 1:length(highlights)){
    ylab=sprintf('%% %s',highlights[i])
    # pdf(sprintf('figs/annotations/%s.pdf',marker.names[prop.markers[i]]))
    png(sprintf('Results/UFO_Figures/Annotations/%s.png',gsub("/", "-", marker.names[prop.markers[i]])), height=900, width=900)
    palette((rich.colors(1000)))
    ##par(bg='black')
    delta.sizex=legend.size/dev.size()[1]
    split.screen(t(cbind(c(0,1-delta.sizex,0,1),c(1-delta.sizex,1,0,1))))
    screen(1)
    ##highlights[i]=Names[which.max(ratsauc[,i])]
    library(gplots)
    ts=(All.Proportions[which(Names==highlights[i]),])
    ts=ts-min(ts)
    temp=quantile(ts,0.99,na.rm=TRUE)
    ts[which(ts>temp)]=temp
    temp=quantile(ts,0.01,na.rm=TRUE)
    ts[which(ts<temp)]=temp
    cols=round(ts*1000/max(ts,na.rm=TRUE))
    chulls=list()
    plot(td,pch=20,col=cols,cex=1.5,axes=FALSE,col.lab=par('fg'),xlim=xlim,ylim=ylim,xlab='UFO Dimension 1',ylab='UFO Dimension 2')
    axis(1,col=par('fg'));axis(2,col=par('fg'));
    title(main=highlights[i],col.main=par('fg'))
    palette('default')
    screen(2)
    palette((rich.colors(1000)))
    par(mar=c(1,4,1,0.2))
    image(matrix(1:2500, 25), col = 1:1000, xaxt='n', ylab='', yaxt='n',col.lab=8)
    par(mgp=c(2.5,1,0))
    title(ylab=ylab,col.lab=par('fg'));
    myvar=(All.Proportions[which(Names==highlights[i]),])
    axis(2, at=c((1/(length(pretty(myvar))-1))*(0:(length(pretty(myvar))-1))), labels=100*pretty(myvar),col=par('fg'),col.axis=par('fg'))
    palette('default')
    close.screen(all.screens=TRUE)
    dev.off()
}


# # pdf('figs/Groups.pdf')
# png('figs/Groups.png', height=900, width=900)
# categorical=metadata[,4]
# ##palette(c('white',palette()[c(2:8)]))
# ##par(bg='black')
# ##par(fg='gray')
# plot(data.frame(td),col=categorical,pch=20, cex=3,axes=FALSE,col.lab=par('fg'),xlab='UFO Dimension 1',ylab='USO Dimension 2')
# axis(1);axis(2);
# legend('topright',col=seq(max(categorical)),pt.cex=3,cex=2,legend=seq(max(categorical)), pch=20, title='Groups')
# ##palette('default')
# dev.off()


## Plot based on days
days <- as.numeric(meta_data[,4])
days <- days - min(days)
# pdf(sprintf('figs/days.pdf'))
png(sprintf('Results/UFO_Figures/days.png'), height=900, width=900)
ylab='days'
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
ts = days
#ts=ts-min(ts)
temp=quantile(ts,0.99,na.rm=TRUE)
ts[which(ts>temp)]=temp
temp=quantile(ts,0.01,na.rm=TRUE)
ts[which(ts<temp)]=temp
cols=round(ts*1000/max(ts,na.rm=TRUE))
chulls=list()
plot(td,pch=20,col=cols+1,cex=1.5,axes=FALSE,xlab='UFO Dimension 1',ylab='UFO Dimension 2',col.lab='black',xlim=xlim,ylim=ylim)
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

UFO_days <- cbind(WT_dates, days, td, meta_data[,5], meta_data[,6])
colnames(UFO_days) <- c("Date", "Days since first assay", "X", "Y", "WT=0, KO=1", "F=0, M=1")
write.table(UFO_days, file="Results/UFO/UFO_days.csv", sep=",", row.names = F)

## Plot based on gender
groups=as.numeric(meta_data[,6])
groups <- groups - min(groups)
png(sprintf('Results/UFO_Figures/gender.png'), height=900, width=900)
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
library(gplots)
ts=(groups)
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
groups=as.numeric(meta_data[,5])
groups <- groups - min(groups)
png(sprintf('Results/UFO_Figures/WT_KO.png'), height=900, width=900)
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
ts=(groups)
ts=ts-min(ts)
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

###################################################################################################
## Finding the outliers in the UFO plots. This part of the script is to be carried out after looking into the UFO plots for days.
## Then use the td matrix (N x 2), which is the output from the Isomap() function to determine which of the files are outliers.
## Which dimension (column in the td matrix) you are using and which value you are using to check for outliers will depend on the plot.
meta_data.Unique <- meta_data
rownames(meta_data.Unique) <- 1:dim(meta_data.Unique)[1]
outlier.Index <- which(td[,1] > 40)
## Make sure after finding the outlier.Index look into the td matrix bak to look into the dimension values.
outlier.List <- meta_data.Unique[outlier.Index,]
write.table(outlier.List, file="Results/UFO_Figures/outlierList.csv", sep=",", row.names = F)

###################################################################################################

## UFO Analysis Plots for WTs only
## Load the meta_data and initialize stuff
NumberOfCPUs = 5 ##for snow
load("Results/meta_dataWT.Rdata")
load("Results/File_NamesWT.Rdata")


# metadata <- read.csv('BatchDayData48.csv')
# files=metadata[,1]

datapath <- '/code/Projects/3i/Panel_T-cell_MLN/FCS_Groups'


library(flowCore)
# temp.fcs <- read.FCS("/code/Projects/IMPC/data/CIPHE/1-1-1-B6/Panel1/14-Aug-19_IMPC1_01_labelled.fcs")
temp.fcs <- read.FCS("/code/Projects/3i/Panel_T-cell_MLN/FCS_Groups/+_+/MLN 7,2f,20,2f,15_L000092105_008.labelled.fcs")
marker.names <- temp.fcs@parameters@data$name
prop.markers=c(7:18)

N=nrow(meta_dataWT)
# N <- 44
CSN=1
load(paste0("Results/UFO_WT/",File_NamesWT[CSN],".Rdata"))
M = length(FT@CellFreqs)

##Load the flowType results
NotAvailable <- vector();
All.Proportions.WT <- matrix(0, M, N)
for (CSN in 1:N){
  print(CSN)
  if (file.exists(paste0("Results/UFO_WT/",File_NamesWT[CSN],".Rdata"))){
    load(paste0("Results/UFO_WT/",File_NamesWT[CSN],".Rdata"))
  }else{
    print(paste0("Results/UFO_WT/",File_NamesWT[CSN],".Rdata Does Not Exist"))
    NotAvailable <- c(NotAvailable, CSN)
    next;
  }
  All.Proportions.WT[,CSN] <- FT@CellFreqs/max(FT@CellFreqs);
}

library(flowType)
library(Rmpi)
##Calculate phenotype names just in case we need them later.
tempfuncWT <- function(x){return(decodePhenotype(x,FT@MarkerNames[prop.markers],FT@PartitionsPerMarker))}
library(snow)
cl <- makeCluster(NumberOfCPUs, type = "SOCK",outfile='/tmp/snowout.log')
clusterExport(cl, list=c('decodePhenotype','FT','prop.markers'))
Names=unlist(clusterApplyLB(cl, FT@PhenoCodes, tempfuncWT));
stopCluster(cl)
rownames(All.Proportions.WT) <- Names

library(RDRToolbox)
is=Isomap(t(All.Proportions.WT),2,10)$dim2
td=is




##Make the surface marker annotations. If you want to highlight other populations add them to the "highlights" vector below.
legend.size=1.25
library(gplots)
highlights=sprintf('%s+',marker.names[prop.markers])

suppressWarnings ( dir.create ( "Results/UFO_WT_Figures") )
suppressWarnings ( dir.create ( "Results/UFO_WT_Figures/Annotations") )

xlim=NULL#c(-101.78373,   82.10356)
ylim=NULL#c(-43.00999,  53.04831)
for (i in 1:length(highlights)){
  ylab=sprintf('%% %s',highlights[i])
  # pdf(sprintf('figs/annotations/%s.pdf',marker.names[prop.markers[i]]))
  png(sprintf('Results/UFO_WT_Figures/Annotations/%s.png',gsub("/", "-", marker.names[prop.markers[i]])), height=900, width=900)
  palette((rich.colors(1000)))
  ##par(bg='black')
  delta.sizex=legend.size/dev.size()[1]
  split.screen(t(cbind(c(0,1-delta.sizex,0,1),c(1-delta.sizex,1,0,1))))
  screen(1)
  ##highlights[i]=Names[which.max(ratsauc[,i])]
  library(gplots)
  ts=(All.Proportions.WT[which(Names==highlights[i]),])
  ts=ts-min(ts)
  temp=quantile(ts,0.99,na.rm=TRUE)
  ts[which(ts>temp)]=temp
  temp=quantile(ts,0.01,na.rm=TRUE)
  ts[which(ts<temp)]=temp
  cols=round(ts*1000/max(ts,na.rm=TRUE))
  chulls=list()
  plot(td,pch=20,col=cols,cex=1.5,axes=FALSE,col.lab=par('fg'),xlim=xlim,ylim=ylim,xlab='UFO Dimension 1',ylab='UFO Dimension 2')
  axis(1,col=par('fg'));axis(2,col=par('fg'));
  title(main=highlights[i],col.main=par('fg'))
  palette('default')
  screen(2)
  palette((rich.colors(1000)))
  par(mar=c(1,4,1,0.2))
  image(matrix(1:2500, 25), col = 1:1000, xaxt='n', ylab='', yaxt='n',col.lab=8)
  par(mgp=c(2.5,1,0))
  title(ylab=ylab,col.lab=par('fg'));
  myvar=(All.Proportions.WT[which(Names==highlights[i]),])
  axis(2, at=c((1/(length(pretty(myvar))-1))*(0:(length(pretty(myvar))-1))), labels=100*pretty(myvar),col=par('fg'),col.axis=par('fg'))
  palette('default')
  close.screen(all.screens=TRUE)
  dev.off()
}



## Plot based on days
days <- as.numeric(meta_dataWT[,4])
days <- days - min(days)
# pdf(sprintf('figs/days.pdf'))
png(sprintf('Results/UFO_WT_Figures/daysWT.png'), height=900, width=900)
ylab='days'
legend.size=1.25
delta.sizex=legend.size/dev.size()[1]
xlim=NULL
ylim=NULL
library(gplots)
split.screen(t(cbind(c(0,1-delta.sizex,0,1),c(1-delta.sizex,1,0,1))))
screen(1)
palette((rich.colors(1000)))
library(gplots)
ts = days
#ts=ts-min(ts)
temp=quantile(ts,0.99,na.rm=TRUE)
ts[which(ts>temp)]=temp
temp=quantile(ts,0.01,na.rm=TRUE)
ts[which(ts<temp)]=temp
cols=round(ts*1000/max(ts,na.rm=TRUE))
chulls=list()
plot(td,pch=20,col=cols+1,cex=1.5,axes=FALSE,xlab='UFO Dimension 1',ylab='UFO Dimension 2',col.lab='black',xlim=xlim,ylim=ylim)
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

UFO_WT_days <- cbind(meta_dataWT[,3], days, td, meta_dataWT[,5], meta_dataWT[,6])
colnames(UFO_WT_days) <- c("Date", "Days since first assay", "X", "Y", "WT=0, KO=1", "F=0, M=1")
write.table(UFO_WT_days, file="Results/UFO_WT/UFO_WT_days.csv", sep=",", row.names = F)

## Plot based on gender
groups=as.numeric(meta_dataWT[,6])
groups <- groups - min(groups)
png(sprintf('Results/UFO_WT_Figures/genderWT.png'), height=900, width=900)
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
library(gplots)
ts=(groups)
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

###################################################################################################
## Finding the outliers in the UFO plots for Wild types. This part of the script is to be carried out after looking into the WT UFO plots for days.
## Then use the td matrix (N x 2), which is the output from the Isomap() function to determine which of the files are outliers.
## Which dimension (column in the td matrix) you are using and which value you are using to check for outliers will depend on the plot.
meta_dataWT.Unique <- meta_dataWT
rownames(meta_dataWT.Unique) <- 1:dim(meta_dataWT.Unique)[1]
outlierWT.Index <- which(td[,1] > 40)
## Make sure after finding the outlier.Index look into the td matrix bak to look into the dimension values.
outlierWT.List <- meta_dataWT.Unique[outlierWT.Index,]
write.table(outlierWT.List, file="Results/UFO_WT_Figures/outlierWTList.csv", sep=",", row.names = F)
