### R code from vignette source 'RchyOptimyx.Rnw'

###################################################
### code chunk number 1: c10
###################################################
library(flowType)
data(HIVData)
data(HIVMetaData)
HIVMetaData <- HIVMetaData[which(HIVMetaData[,'Tube']==2),];


###################################################
### code chunk number 2: c11
###################################################
Labels=(HIVMetaData[,2]=='+')+1;


###################################################
### code chunk number 3: c12
###################################################
library(flowCore)
library(RchyOptimyx)

##Markers for which cell proportions will be measured.
PropMarkers <- 5:10

##Markers for which MFIs will be measured.
MFIMarkers <- PropMarkers

##Marker Names
MarkerNames <- c('Time', 'FSC-A','FSC-H','SSC-A',
                 'IgG','CD38','CD19','CD3',
                 'CD27','CD20', 'NA', 'NA')

##Apply flowType
ResList <- fsApply(HIVData, 'flowType', PropMarkers, 
                   MFIMarkers, 'kmeans', MarkerNames);

##Extract phenotype names
phenotype.names=unlist(lapply(ResList[[1]]@PhenoCodes,function(x){return(decodePhenotype(x,MarkerNames[PropMarkers],ResList[[1]]@PartitionsPerMarker))}))
names(ResList[[1]]@PhenoCodes)=phenotype.names


###################################################
### code chunk number 4: c13
###################################################
all.proportions <- matrix(0,length(ResList[[1]]@CellFreqs),length(HIVMetaData[,1]))
for (i in 1:length(ResList))
  all.proportions[,i] = ResList[[i]]@CellFreqs / ResList[[i]]@CellFreqs[1]


###################################################
### code chunk number 5: c14
###################################################
Pvals <- vector();
EffectSize <- vector();
for (i in 1:dim(all.proportions)[1]){

  ##If all of the cell proportions are 1 (i.e., the phenotype 
  ##with no gates) the p-value is 1.
  if (length(which(all.proportions[i,]!=1))==0){
    Pvals[i]=1;
    EffectSize[i]=0;
    next;
  }
  temp=t.test(all.proportions[i, Labels==1], 
    all.proportions[i, Labels==2])
  Pvals[i] <- temp$p.value
  EffectSize[i] <- abs(temp$statistic)  
}

Pvals[is.nan(Pvals)]=1
names(Pvals)=phenotype.names

##Bonferroni's correction
## TODO: no good results shows up and there is no adjusted p-value <0.6
selected <- which(p.adjust(Pvals)<0.65);

print(names(selected))


###################################################
### code chunk number 6: c16
###################################################
res<-RchyOptimyx(ResList[[1]]@PhenoCodes, -log10(Pvals), 
		ResList[[1]]@PhenoCodes[selected[38]], factorial(6),FALSE)
plot(res, phenotypeScores=-log10(Pvals), phenotypeCodes=ResList[[1]]@PhenoCodes, marker.names=MarkerNames[PropMarkers], ylab='-log10(Pvalue)')


###################################################
### code chunk number 7: c17
###################################################
res<-RchyOptimyx(pheno.codes=ResList[[1]]@PhenoCodes, phenotypeScores=-log10(Pvals), 
		startPhenotype=ResList[[1]]@PhenoCodes[selected[38]], 2, FALSE)
plot(res, phenotypeScores=-log10(Pvals), phenotypeCodes=ResList[[1]]@PhenoCodes, marker.names=MarkerNames[PropMarkers], ylab='-log10(Pvalue)')


###################################################
### code chunk number 8: c18
###################################################
res<-RchyOptimyx(pheno.codes=ResList[[1]]@PhenoCodes, phenotypeScores=-log10(Pvals), 
  	startPhenotype=ResList[[1]]@PhenoCodes[selected[1]], 1, FALSE,trim.level=4)
for (i in 2:length(selected)){
temp<-RchyOptimyx(pheno.codes=ResList[[1]]@PhenoCodes, phenotypeScores=-log10(Pvals), 
    startPhenotype=ResList[[1]]@PhenoCodes[selected[i]], 1, FALSE,trim.level=4)  
res=merge(res,temp)
}
plot(res, phenotypeScores=-log10(Pvals), phenotypeCodes=ResList[[1]]@PhenoCodes, marker.names=MarkerNames[PropMarkers], ylab='-log10(Pvalue)')


