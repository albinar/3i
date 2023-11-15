# ############################################################################################
# # removes all cells that exist on the edges of FSC-A, FSC-H, SSC-A, and SSC-H. 
# removeMargins <-function(f.temp){
#     
# 
# for ( q1 in 1:length(colnames(f.temp))){ # 6 because there may be 3 FSC and 3 SSC columns 
#         if (length(which(f.temp@exprs[,q1]==262143))>0){f.temp <- f.temp[-which(f.temp@exprs[,q1]==262143)]}
#         if (length(which(f.temp@exprs[,q1]==0))>0)     {f.temp <- f.temp[-which(f.temp@exprs[,q1]==0)]}
# }
# return(f.temp)
# 
# }

removeMargins<- function(f,chans,sens=1, debris=FALSE,return.ind=F,neg=500, verbose = T)
{
  neg <-cbind(chans,neg)[,2]
  #Size is a vector of size 2, to be passed to mfrow in case of plotting
  data <- exprs(f)
  margins <- c()
  marg.list <-list()
  if(!debris)
  {
    for(chan in chans)
    {
      stain.max <-max(data[,chan])
      margins <- which ( data[, chan] >= stain.max*sens)
      marg.list <- c(marg.list, list(margins))
      data <- data[ -margins, ]
      if(verbose == T){print(paste(length(margins), "margin events in",colnames(f)[chan], "will be removed.",sep =" "))}
    }
    
  }else
  {
    for(i in 1:length(chans))
    {
      stain.min <-min(data[,chans[i]])
      margins <- which ( data[, chans[i]] <= stain.min*sens)
      if (neg[i]<500)
      {
        negs <- which ( data[, chans[i]] < neg[i])
        margins <- negs
      }
      marg.list <- c(marg.list, list(margins))
      if (length(margins)!=0){
        data <- data[ -margins, ]
      }
      if(verbose == T){print(paste(length(margins), "debris events in",colnames(f)[chans[i]], "will be removed.",sep =" "))}
    }
  }
  exprs(f) <- data
  if (!return.ind)
    return(f)
  else
    return(list(frame=f,ind=marg.list))
  
}

############################################################################################
# removes all cells that exist on the edges of FSC-A, FSC-H, SSC-A, and SSC-H. 
removeMarginsAndNegatives <-function(f.temp){
    

for ( q1 in 1:length(colnames(f.temp))){ # 6 because there may be 3 FSC and 3 SSC columns 
        if (length(which(f.temp@exprs[,q1]==262143))>0){f.temp <- f.temp[-which(f.temp@exprs[,q1]==262143)]}
        if (length(which(f.temp@exprs[,q1]<=0))>0)     {f.temp <- f.temp[-which(f.temp@exprs[,q1]<=0)]}
}
return(f.temp)

}


############################################################################################
#Remove all cell types on the outside of an elipse from the FSC-A and FSC-H plot.
removeDoublets <-function(f.temp, temp.nameU, q){
if (length(which(f.temp@parameters@data$name=="FSC-H"))>=1) {   
  png ( file = paste("Plots/Figures/",temp.nameU, "_", q, "_B4DoubMarg.png", sep=""))
    plotDens(f.temp, c(1,2))  
  dev.off()
  
  singlet <- flowDensity(f.temp,c(which(f.temp@parameters@data$name=="FSC-A"),which(f.temp@parameters@data$name=="FSC-H")),
                                position=c(T,F),percentile=c(.01,.99),use.percentile=c(T,T),ellip.gate=T)
  f.temp <- getflowFrame(singlet)
  
  png ( file = paste("Plots/Figures/",temp.nameU, "_", q, "_AftDoubMarg.png", sep=""))
    plotDens(f.temp, c(1,2))  
  dev.off()
}
return(f.temp)
}



############################################################################################
# Remove all cell types on the outside of an elipse from the FSC-A and FSC-H plot.
removeDoubletsSSC <-function(f.temp, name1, name2, directory, Plt){
   
if (length(which(f.temp@parameters@data$name=="SSC-H"))>=1) {   
  singlet <- flowDensity(f.temp,c(which(f.temp@parameters@data$name=="SSC-A"),which(f.temp@parameters@data$name=="SSC-H")),
                                position=c(T,F),percentile=c(0.01,.99),use.percentile=c(T,T),ellip.gate=T,scale = 0.999999)

  f.temp <- getflowFrame(singlet)
}
return(f.temp)
}

#################################################################################
# rearranges the order of the flowFrame
reorderfSet <- function(f,commonMarker,columnPlacement){

NumOfChannels <- length(f@parameters@data$desc)    
temp.number <- which(f@parameters@data$desc==commonMarker)
if ( length(temp.number)>1){ temp.number <- temp.number[length(temp.number)]}
row.order <- NULL
for(q in 1:NumOfChannels) {row.order <- c(row.order,q)}
row.order[temp.number]    <- columnPlacement
row.order[columnPlacement]<- temp.number
f <- f[,row.order]
return(f)
}


#################################################################################
# Prints out the time since start_time. Used for optimizing code.
TimeOutput <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}
TimeOutput(Sys.Date())


#################################################################################
#SVD reduction for flowType to be used for grouping patients based on reduced matrix
#Kmeans used for grouping patients, any other method can be used

#Author: Mehrnoush Malek
#Date: April 2012

svd.reduction <- function(cell.prop,kmean.centres=3,kmeans.start=1000)
  #cell.prop is a matrix of size 3^k by m, where k is number of markers used in flowType and m is number of samples
{
  
svdr.Result<-svd(t(cell.prop))
#Find a threshold where samples get far from others
x11();plot(sort(svdr.Result$d))

inds<-which(svdr.Result$d > locator()$y)
acf<-t(cell.prop) %*% (svdr.Result$v[,inds])
kmeans.cluster<-kmeans(acf, centers=kmean.centres, nstart=kmeans.start)
x11();plot(data.frame(acf), pch=16,col=kmeans.cluster$cluster)
return(kmeans.cluster)
}

# length(flowType.res@CellFreqs)
# which(names(flowType.res@CellFreqs)=="")
# flowType.res@CellFreqs[1093]
# flowType.res@CellFreqs[1:2]
# ?heatmap
# d<-plot(flowType.res,flowFrame_fB)
# d
#######################################################################################
## This part was taken from the helperfunc.R script
##Finds markers in the FCS file
Find.markers <- function(frame,marker.list)
{
  #Parameters:
  #*frame: a flowFrame in the flowSet
  #**marker.list: A vector of characters
  #Output:
  #*channels.ind: a vector of channels numbers in the frame  corresponding to marker.list
  channels.ind <- unlist(lapply(marker.list, function(x) {
    ind <- grep(x, frame@parameters@data[,2], ignore.case=T)
    ind_store <- ind
    if(length(ind)==0){
      warning(paste (x, "not found, check markers!"))
      return(NA)
    } else {
      if(length(ind)>1) {
        cnt <- 0
        repeat{
          cnt <- cnt + 1
          fs.markers<-unlist(lapply(frame@parameters@data[,2], function(x) unlist(strsplit(x," "))[cnt]))
          ind<-match(x,fs.markers)
          if (is.na(ind))
          {
            fs.markers<-unlist(lapply(frame@parameters@data[,2], function(x) unlist(strsplit(x,"-"))[cnt]))
            ind<-match(x,fs.markers)
            if(!is.na(ind))
              break;
          } else {
            break;
          }
          if(cnt >= 10) {
            
            if (length(ind_store) >= 2){
              ind <- ind_store[1]
              warning(paste (x, "found more than one, choosing first. Check markers!"))
            } else {
              warning(paste (x, "not found, check markers!"))
            }
            break;
          }
        }
      }
    }
    return(ind)
  }))
  names(channels.ind)<-marker.list
  #Removing NAs in the channel vector, as not all centres have all of these markers
  #Note that most centres should have Live/CD4/CD8/CD44/CD62/CD5/CD161/CD25 in their channels
  ind <- which (is.na(channels.ind))
  if (length(ind)!=0)
    channels.ind <- channels.ind[-ind]
  return(channels.ind)
}