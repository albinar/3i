###3i semi-supervised pipeline written for BM panel
##Written by : Albina Rahim, Dec,2017
#Updated: July 3, 2018
#**************************************
#              Notes
#**************************************
#Files are extracted by Albina
#flowTypeFilter package is not public.
rm(list=ls())

## Setting the path based on docker
setwd("/home/rstudio/code/Projects/3i/semisupervised")


library(fields)
library(flowCore)
library(flowTypeFilter)
library(flowDensity)
library(ggplot2)
#library(cowplot)
library(data.table)
library(Cairo)
library(RSVGTipsDevice)


##########################
##########################

##########################
##########################
##########################
source("3i-Helper.R")
source("Ania_helper.R")
#Input.path needs to be changed based on user's directory
input.path <- "/home/rstudio/data/3i/"
##This is for BM:
panel <- c("Panel_B-cell", "Panel_BM-cell","Panel_M-cell", "Panel_P2-cell","Panel_T-cell" )[3]
organ <- c("","MLN","Spleen")[3]

## Loading Gating thresholds and filters for each of the FCS files analysed in supervised analysis step
temp <- load(paste0(input.path,panel,"/", "semisupervised/",organ,"/Gates.list.Rdata"))
gates.list <- get(temp)

## Loading the FCS files - starting from CD45 population - saved as .Rdata
load(list.files(paste0(input.path,panel,"/semisupervised/",organ,"/flowType/"),pattern = ".Rdata",full.names = T)[1])
#load(list.files(paste0(input.path,panel,"/semisupervised/",organ,"/flowType/"),pattern = ".Rdata",full.names = T)[1034])

## Loading proportions of all populations generated using flowType for each FCS file 
load(paste0(input.path,panel,"/semisupervised/",organ,"/Allproportions-FT.RData"))
colnames(ft.props)[1]<-"CD45+"
pop.names <- colnames(ft.props)



man.pops<- as.matrix(read.csv(paste(input.path,panel,"/semisupervised/",organ,"/Matchingpopulation.csv",sep="")))
variable<- rep(NA,length(pop.names))
variable[match(man.pops[,1],pop.names)]<- as.vector(man.pops[,2])

allinfo <- as.matrix(read.csv(file = paste(input.path,panel,"semisupervised",organ,"Mnumber_Lnumber.csv",sep="/"),check.names = F))
whole_dt <- as.data.frame(read.csv(paste(input.path,panel,"/semisupervised/",organ,"/allProportions_Table.csv",sep=""),check.names = F))
ind<-match(rownames(ft.props), whole_dt$`FCS files`)
whole_dt<-whole_dt[ind,]
colnames(whole_dt)[which(colnames(whole_dt)=="Assay Date")]<-"Assay_Date"
whole_dt[,"Assay_Date"]<- as.Date(whole_dt[,"Assay_Date"],format = "%d-%b-%y")
temp<-c(); for ( i in 1:nrow(whole_dt)){
  ind <-which(allinfo[,"Label_Barcode"]==whole_dt$Barcodes[i])
  if(length(ind)>1)
    print(length(ind))
  temp <- c(temp,paste(allinfo[ind, "Colony_Prefix"],allinfo[ind, "zygosity"],sep="_"))
}
whole_dt <- as.data.table(whole_dt)

whole_dt<-whole_dt[,strain:=temp]

if (panel=="Panel_BM-cell")
{
  outliers <- read.csv(paste0(input.path,"Panel_BM-cell/semisupervised/suspected_dates_2med_bothsex_2animals_BM.csv"))
  for (i1 in as.vector(outliers$Assay_Date)){
    ind <- grep(i1,whole_dt[,Assay_Date])
    if (length(ind)>0)
    {
      
      whole_dt<- whole_dt[-ind,]  
    }else{
      print(i1)
    }
  }
}else{
  source("detect_outliers_proportions.R")
  for (i1 in as.vector(suspected_dates_BM$`Assay Date`)){
    ind <- grep(i1,whole_dt[,Assay_Date])
    if (length(ind)>0)
    {
      
      whole_dt<- whole_dt[-ind,]  
    }else{
      print(i1)
    }
  }
}
allko <- whole_dt[Genotype!="WT",]
strains<- unique(allko[,strain])
updated.genes <-c()
#Remove genotypes that have less than 4 mice, or 0 samples for either sexes.
for (i in strains)
{
  d<- which(allko[,strain]==i) 
  
  if (allko[d,.N]<4 | (length(which(allko[d,"Gender"]=="Female"))<1| length(which(allko[d,"Gender"]=="Male"))<1 )){
    print("skipping, less than 2 mice.")
  }else{
    updated.genes <- c(updated.genes, i)
  }
}




#Based on Ania's code of choosing WT samples for each KO sample
#Fixed order issue in Ania's code, and fixed misising--> missing

source("find_RR_dates.R")
take_only_WTs <- function(all_data_dt){
  #for given data dt only WT individuals, requires the column "Genotype"
  return(all_data_dt[Genotype%in%c("+/+","+/Y","WT")])
  
}


param<-ifelse(panel=="Panel_M-cell",yes="Lineage neg-%Parent",no=ifelse(panel=="Panel_B-cell",no="CD45+-%Parent",yes="CD45-%Parent")) # I don't understand the need of this, though
all_res <- find_RR_dates(whole_dt, param, by_sex=TRUE,min_ind = 70)
#Finding WT samples for each KO samples within each genotype
size<-70
counter2<-1
all.props<-c()
hits<-c()
all.pvals<-c()
all.effects<-c()
all.files<-c()
rm(genotype)
for ( g1 in updated.genes)
{
  print(counter2)
  ko.dates<-allko[strain==g1,]
  all.ranks <- lapply( 1:nrow(ko.dates), function(i1){
    print(i1)
    dd <- ko.dates[i1,Assay_Date]
    which.sex<- as.vector(ko.dates[i1,Gender])
    st.date <- all_res[Assay_Date==dd & sex==tolower(which.sex) ,st_date]
    end.date <- all_res[Assay_Date==dd & sex==tolower(which.sex) ,end_date]
    WTs<-as.vector(subset(whole_dt,Genotype=="WT"& Gender==which.sex&Assay_Date>=st.date& Assay_Date<=end.date)[,"FCS files"][[1]])
    m<-match(WTs, rownames(ft.props))
    if (length(is.na(m))>0)
      WTs<-WTs[which(!is.na(m))]
    if (length(WTs)>size)
    {
      set.seed(250)
      WTs<- sample(WTs,size)
    }
    ranks<- apply(ft.props[c(as.vector(ko.dates[i1,"FCS files"][[1]]),WTs),],2,rank)
    return(list(sex=rep(which.sex,length(WTs)),ranks=ranks,wt=WTs))
  })
  
  
  ko.ranks <- sapply(all.ranks,function(x) return(x$ranks[1,]))
  genotype<-as.vector(unique(allko$Genotype[which(allko$strain==g1)]))
  #  This is to keep the file selection information, later on for flowRepository file upload
  
  all.files <- rbind(all.files,cbind(rbind(cbind(unique(unlist(lapply(all.ranks , function(x) x$wt))), "WT",genotype[1]), cbind(as.vector(ko.dates$`FCS files`),"KO",genotype[1])),
                                     unlist(strsplit(g1,"_"))[1],unlist(strsplit(g1,"_"))[2]))
  #  counter2 <- (counter2+1)
  #}
  #####################################################################
  ######P-values: Using Ania's get_subset_wt_ranks function
  #####################################################################
  pvals<- apply(ko.ranks,1, function(rank1){
    rank2<-get_subset_wt_ranks(1:71, rank1)
    return(wilcox.test(x = rank1,y=rank2, var.equal = varaince)$p.value)
  })
  
  pvals.adjust <- p.adjust(pvals,method = "fdr")
  
  #####################################################################
  # Adding effect size
  #####################################################################
  columnMeans <- function(mat)
  {
    if (is.vector(mat)){
      return(mat)
    }else{
      return(colMeans(mat))
    }
  }
  # calculate.effect <- function(all.ranks, ko)
  # {
  #   all.wt <- unlist(lapply(all.ranks, function(x) x$wt))
  #   unique.samples <- which(duplicated(all.wt)==FALSE)
  #   wt.sex <- unlist(lapply(all.ranks, function(x) x$sex))
  #   wt.sex<- wt.sex[unique.samples]
  #   all.wt <- all.wt[unique.samples]
  #   sd_WT_male <-apply(ft.props[all.wt[which(wt.sex=="Male")],],2,sd)
  #   sd_WT_female <-apply(ft.props[all.wt[which(wt.sex=="Female")],],2,sd)
  #   
  #   
  #   mean_KO_male<-columnMeans(ft.props[as.vector(ko$`FCS files`)[which(ko$Gender=="Male")],])
  #   mean_WT_male<-colMeans(ft.props[all.wt[which(wt.sex=="Male")],])
  #   
  #   Effsize_male <- ( mean_KO_male- mean_WT_male)/sd_WT_male
  #   
  #   mean_KO_female <-columnMeans(ft.props[as.vector(ko$`FCS files`)[which(ko$Gender=="Female")],])
  #   mean_WT_female<-colMeans(ft.props[all.wt[which(wt.sex=="Female")],])
  #   
  #   Effsize_female <- ( mean_KO_female- mean_WT_female)/sd_WT_female
  #   
  #   return(colMeans(rbind(Effsize_male,  Effsize_female)))
  #   
  # }
  calculate.effect <- function(all.ranks, ko)
  {
    all.wt <- unlist(lapply(all.ranks, function(x) x$wt))
    unique.samples <- which(duplicated(all.wt)==FALSE)
    wt.sex <- unlist(lapply(all.ranks, function(x) x$sex))
    wt.sex<- wt.sex[unique.samples]
    all.wt <- all.wt[unique.samples]
    mean_WT_male <-colMeans(ft.props[all.wt[which(wt.sex=="Male")],])
    sd_WT_male <-apply(ft.props[all.wt[which(wt.sex=="Male")],],2,sd)
    mean_WT_female <-colMeans(ft.props[all.wt[which(wt.sex=="Female")],])
    sd_WT_female <-apply(ft.props[all.wt[which(wt.sex=="Female")],],2,sd)
    
    
    mean_KO_male<-columnMeans(ft.props[as.vector(ko$`FCS files`)[which(ko$Gender=="Male")],])
    mean_KO_male_scaled <- mean_KO_male/sd_WT_male
    mean_KO_female <-columnMeans(ft.props[as.vector(ko$`FCS files`)[which(ko$Gender=="Female")],])
    mean_KO_female_scaled <- mean_KO_female/sd_WT_female
    return(colMeans(rbind(mean_KO_male_scaled , mean_KO_female_scaled)))
    
  }
  genotype<-as.vector(unique(allko$Genotype[which(allko$strain==g1)]))
  if (length(genotype)>1){
    warning("Something is not right.")
    genotype<-unlist(strsplit(as.vector(unique(allko$Genotype[which(allko$strain==g1)]))[1],":"))[1]
  }
  effect.size <- calculate.effect(all.ranks = all.ranks,ko = ko.dates)
  #All effect sizes for all KO line and samples
  all.effects<- rbind(all.effects,cbind(variable,unlist(strsplit(g1,"_"))[1],unlist(strsplit(g1,"_"))[2],genotype,effect.size))
  all.pvals<- rbind(all.pvals,cbind(variable,unlist(strsplit(g1,"_"))[1],unlist(strsplit(g1,"_"))[2],genotype,pvals.adjust))
  
  
  #*************************************************************************************************
  #*************************************************************************************************
  #*************************************************************************************************
  #*************************************************************************************************
  
  ##Commenting the cheunk below until we agree on the threshold for Effect size and p-values.
  
  ##########################################################
  # Cutoff threshold for effect size and p-values
  ##########################################################
  #Comment it out for now till Ania gives threshold for both effect,size and p.values
  
  # significant.pvals <-which(pvals.adjust<0.02)
  # if (length(significant.pvals)>0)
  # {
  #   print(g1)
  #   print("found phenodeviants.")
  #   wt.files<-unlist(lapply(all.ranks,function(x) x$wt))
  #   dup.wt<-duplicated(wt.files)
  #   ko.files<- as.vector(ko.dates$`FCS files`)
  #   lbls<-c(rep("KO",nrow(ko.dates)),rep("WT",length(unique(wt.files))))
  #   
  #   wt.files<- wt.files[-which(dup.wt)]
  #   genders<-c(unlist(lapply(all.ranks,function(x) unique(x$sex))),unlist(lapply(all.ranks,function(x) x$sex))[-which(dup.wt)])
  #   all.props[[counter2]]<-list(props=ft.props[c(ko.files,wt.files),], labels=lbls,gender=genders, genotype=g1,pvals=pvals.adjust)
  #   if (length(significant.pvals)==1){
  #     hits <-rbind(hits, c(paste(length(which(lbls=="KO")), "vs.",length(which(lbls=="WT")),sep=" "),g1,pvals.adjust[significant.pvals,]))
  #   }else{
  #     hits <-rbind(hits, cbind(rep(paste(length(which(lbls=="KO")), "vs.",length(which(lbls=="WT")),sep=" "),length(significant.pvals)),rep(g1,length(significant.pvals)),pvals.adjust[significant.pvals]))
  #   }
  #   counter2 <- (counter2+1)
  # }else{
  #   print("Nothing significant")
  # }
  
  #*************************************************************************************************
  #*************************************************************************************************
  #*************************************************************************************************
  #*************************************************************************************************
  counter2 <- (counter2+1)
  
}
#colnames(hits)<- c("WT-KO counts","Gene","Adj pvalue_WRank")
#save(hits,file = paste0(input.path,panel,"/",organ,"/semisupervised/Hits-WRank.RData"))
#save(all.props,file = paste0(input.path,panel,"/",organ,"/semisupervised/PropsHits-WRank.RData"))


all.pvals<-cbind(rownames(all.pvals),all.pvals)
all.effects <- cbind(rownames(all.effects),all.effects)
rownames(all.pvals)<-rownames(all.effects) <- NULL
colnames(all.pvals)<-c("Cell type","Variable","Colony_Prefix","zygosity","Genotype","Adjusted P-value")
colnames(all.effects)<-c("Cell type","Variable","Colony_Prefix","zygosity","Genotype","Effect size")
colnames(all.files)<-c("FCS File","Type","Genotype","Colony_Prefix","zygosity")
write.csv(x = all.pvals,file=paste(input.path,panel,"/semisupervised/",organ,"/All_genes_Pvalues-old.csv",sep=""))
save(all.pvals,file=paste(input.path,panel,"/semisupervised/",organ,"/All_genes_Pvalues-old.RData",sep=""))
save(all.effects,file=paste(input.path,panel,"/semisupervised/",organ,"/All_genes_EffectSize-old.RData",sep=""))
write.csv(x = all.files,file=paste(input.path,panel,"/semisupervised/",organ,"/All_FilesInfo.csv",sep=""))

##################################
## This line added by Albina:
load(paste0(input.path,panel,"/semisupervised/", organ, "/Hits-WRank.RData"))
reduced.genes <- unique(hits[,2])

tab <- c()
for ( u1  in reduced.genes)
{
  tab <- rbind(tab, c(u1, names(table(hits[which(hits[,2]==u1)[1],1])), length(which(hits[,2]==u1)),
                      round(median(as.numeric(hits[which(hits[,2]==u1),3])),3)))
}
colnames (tab) <- c("genes","Number of mice","Number of hits","Average of adjusted p-value: Wilcoxon Rank")
write.csv(tab, file=paste0(input.path,panel,"/",organ,"/semisupervised/ListOfGeneswithPhenodeviants-WRank.csv"))


####Boxplot
if ( panel== "Panel_BM-cell"){
  markers<-c("GR1","Tcells","Plasma","B220","CD43","CD24","BP1","IgM","IgD")
}else if ( panel== "Panel_M-cell"){
  markers<-c("Lin","F4/80","?","Ly6G","Ly6C","CD11b","CD137","CD11c","MHCII","CD86","CD103") 
}else if ( panel== "Panel_B-cell") {
  markers <- c("b220","Plasma","B1Bcells", "CD95","GL7","IgG1","MZP","Transitional cells",
               "CD23","IgM","IgD" ,"IgG MemoryBcells")
}
path <- paste("/mnt/data/IMPC-flowType",panel,organ,"V-Bplots/",sep="/")
#reduced.genes<- sapply(reduced.genes,function(x) paste(unlist(strsplit(x, "_")),"__",sep="",collapse=""))
draw.violinplot(reduced.genes = reduced.genes,ft = ft,pop.names,all.props,plot.path=path,markers=markers)


##########
#Rchy generation
#paste(input.path,panel,"/semisupervised/",organ,"RchyOptimyx/",sep="")
path <- paste(input.path,panel,"/semisupervised/",organ,"/RchyOptimyx/",sep="")
#path <- paste("/mnt/data/IMPC-flowType",panel,organ,"Rchy/",sep="/")
if (panel=="Panel_BM-cell")
{
  filt <- c(F,T,T,F,F,F,F,F,F)
}else if (panel=="Panel_M-cell"){
  filt <-  c(T,T,T,T,F,F,T,F,F,T, F)
}else if (panel=="Panel_B-cell")
{
  filt <- c(F,T,T,F,F,F,T,T,F,F,F,T)
}

tmp<- generate.rchy.plots(output.path = path, marker.names = markers,all.pvals = selected.pval,reduced.genes = selected.reduced.genes,ft = ft,filters = filt )



#########
#setwd("/home/rstudio/code/Projects/IMPC-Universal/Panel2")
load("/home/rstudio/results/3i/Panel_M-cell/semisupervised/Spleen/All_genes_EffectSize.RData")
ind <- which(all.effects[,5]== "Vps13a(b):Hom") 
#ind <- which(all.effects[,5]== "Pclaf:Hom")
ind2 <- which(as.numeric(all.effects[ind,6]) > 1.6)
load("/home/rstudio/results/3i/Panel_M-cell/semisupervised/Spleen/All_genes_Pvalues.RData")
ind3 <- which(all.pvals[,5]=="Vps13a(b):Hom")
#ind3 <- which(all.pvals[,5]=="Pclaf:Hom")
ind4 <- which(as.numeric(all.pvals[ind3,6]) < 0.2)
ind.Pval.Effect <- intersect(ind[ind2],ind3[ind4])
#ind.Pval.Effect <- ind4


all.pvals.original <- all.pvals
selected.reduced.genes <- "MUDS_Hom"
#selected.reduced.genes <- "MGRK_Hom"
selected.pval <- all.pvals[ind.Pval.Effect,]
selected.pval <- cbind(selected.pval,NA)
selected.pval[,ncol(selected.pval)] <- selected.reduced.genes
colnames(selected.pval)[7] <- "Gene"

min.pvalue <- min(all.pvals.original[,6])

path <- paste(input.path,panel,"/semisupervised/",organ,"/RchyOptimyx/",sep="")
suppressWarnings(dir.create(path))

tmp<- generate.rchy.plots(output.path = path, marker.names = markers,all.pvals = selected.pval,reduced.genes = selected.reduced.genes,ft = ft,filters = filt )

####################################################################

## Written by Mehrnoush so I can test
set.seed(100)
library(flowTypeFilter)
library(RchyOptimyx)
# load("~/data/Specimen_001_L000073886_C02_006.labelled.fcs_FT.Rdata")   
# load("~/data/All_genes_EffectSize.RData")

path <- paste(input.path,panel,"/semisupervised/",organ,"/RchyOptimyx/",sep="")
load(paste(input.path,panel,"/semisupervised/", organ,"/flowType/Specimen_001_L000073886_C02_006.labelled.fcs_FT.Rdata", sep = ""))   

## Effect size is the Group mean divided by the SD of the controls
load(paste(input.path,panel,"/semisupervised/", organ,"/All_genes_EffectSize.RData",sep=""))

## p-values
load(paste(input.path,panel,"/semisupervised/", organ,"/All_genes_Pvalues.RData", sep = ""))

ind <- which(all.effects[,"Genotype"]=="Vps13a(b):Hom")
ind2 <- which(as.numeric(all.effects[ind,6])>1.6)

ind3 <- which(all.pvals[,"Genotype"]=="Vps13a(b):Hom")
ind4 <- which(as.numeric(all.pvals[ind3,6])<0.2)
filt <-  c(T,T,T,T,F,F,T,F,F,T, F)
markers<-c("Lin","F4/80","?","Ly6G","Ly6C","CD11b","CD137","CD11c","MHCII","CD86","CD103")
pvals <- all.pvals[ind3,6]
pvals <-as.numeric(pvals)
marker.names = markers
output.path = path
filters = filt
pheno.names <- all.pvals[ind3,1]
selected.pheno <- ft@PhenoCodes[match(all.pvals[intersect(ind2,ind4),1],pheno.names)]
library(fields)
if (length(selected.pheno)>10)
  selected.pheno <- find.farthest.pheno(pheno = selected.pheno,min.limit = 8)

rchy.res <- list()

for (i in 1:length(selected.pheno))
{
  rchy.res[[i]]<-RchyOptimyx(pheno.codes=ft@PhenoCodes, phenotypeScores=-log10(pvals),
                             startPhenotype=selected.pheno[i], pathCount=3,trimPaths=FALSE,trim.level = 6) 
}
merged.res <- merge(rchy.res[[1]],rchy.res[[2]])
if (length(rchy.res)>2)
{
  for ( i in 3:length(rchy.res))
    merged.res<-merge(merged.res,rchy.res[[i]])
}


# Phenotypes <- unlist(lapply(ft@PhenoCodes, function(x){return(decodePhenotype(
#   x, ft@MarkerNames[ft@PropMarkers], ft@PartitionsPerMarker, filters=filters))}))
# 
# selectedPhenocodes.index <- NULL
# selectedPhenocodes <- NULL
# for(i in 1:ncol(merged.res@nodes)){
#   selectedPhenocodes.index[i] <- which(merged.res@nodes[1,i] == ft@PhenoCodes)
#   selectedPhenocodes[i] <- Phenotypes[selectedPhenocodes.index[i]]
# }



panel<-"Panel_M-cell"
gene <- reduced.genes <- "Vps13a(b):Hom"
#png(paste(path,gene,"_RchyOptimyx.png",sep=""),width = 1900,height = 1600,res=150,pointsize = 12 )

svg(filename = paste(path,gene,"_RchyOptimyx.svg",sep=""),width = 12,height = 8,pointsize = 12)
temp<- plot(merged.res,phenotypeCodes=ft@PhenoCodes,phenotypeScores=-log10(pvals),root.name=paste(gene,ifelse(panel=="Panel_M-cell",yes = "Lin-",no = ": CD45+")),legend.size=1.2,
     partitions.per.marker=2,filters=filters,
     node.lwd = 1,marker.names=marker.names)
dev.off()


# merged.res.original <- merged.res
# index.edges <- which(merged.res@edges[2,] == "0")
# 
# for(i in 1:length(index.edges)){
#   merged.res@edges[2,index.edges[i]] <- "0.05"
# }
# 
# index.edges <- which(merged.res@edges[2,] < 0.006)
# 
# for(i in 1:length(index.edges)){
#   merged.res@edges[2,index.edges[i]] <- 0.05
# }
# 




# store.edges <- merged.res@edges[2,]
# 
# store.edges <- as.numeric(store.edges)+1
# 
# index.edges <- which(store.edges < 1)
# for(i in 1:length(index.edges)){
#   store.edges[index.edges[i]] <- 0.99
# }
# 
# for(i in 1:length(merged.res@edges[2,])){
#   merged.res@edges[2,i] <- store.edges[i]
# }

for(i in 1:length(merged.res@edges[2,])){
  merged.res@edges[2,i] <- 1
}

set.seed(101)
#--------
devSVGTips(paste(path,gene,"_RchyOptimyx_Tooltip.svg",sep=""), toolTipMode=1, width=18,height=12, toolTipFontSize = 12)
plot(merged.res,phenotypeCodes=ft@PhenoCodes,phenotypeScores=-log10(pvals),root.name=paste(gene,ifelse(panel=="Panel_M-cell",yes = "Lin-",no = ": CD45+")),legend.size=1.2,
     partitions.per.marker=2,filters=filters,
     node.lwd = 1,marker.names=marker.names)
setSVGShapeToolTip(title=NULL, desc="Vps13a(b):Hom")
text(x=-0.48, y=1.01, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="Ly6G_")
text(x=-0.8, y=0.81, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="CD103_")
text(x=-0.65, y=0.81, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="CD86_")
text(x=-0.54, y=0.81, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="CD11b+_")
text(x=-0.41, y=0.81, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="NotCD137_")
text(x=-0.28, y=0.81, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="NotLin_")
text(x=-0.14, y=0.81, lab=".", col="black")


setSVGShapeToolTip(title=NULL, desc="Ly6G_CD86_")
text(x=-0.89, y=0.6, lab=".", col="black")


setSVGShapeToolTip(title=NULL, desc="CD86_CD103-_")
text(x=-0.67, y=0.6, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="CD11b+_CD86_")
text(x=-0.47, y=0.6, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="NotCD137_CD86_")
text(x=-0.24, y=0.6, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="NotLin_CD86_")
text(x=-0.04, y=0.6, lab=".", col="black")


setSVGShapeToolTip(title=NULL, desc="Ly6G_CD11b+_CD86_")
text(x=-1.02, y=0.4, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="CD11b+_CD86_CD103-_")
text(x=-0.8, y=0.4, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="Ly6C-_CD11b+_CD86_")
text(x=-0.49, y=0.4, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="CD11b+_MHCII+_CD86_")
text(x=-0.25, y=0.4, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="CD11b+_NotCD137_CD86_")
text(x=0.31, y=0.4, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="NotLin_CD11b+_CD86_")
text(x=0.67, y=0.4, lab=".", col="black")


setSVGShapeToolTip(title=NULL, desc="Ly6G_CD11b+_CD11c-_CD86_")
text(x=-2.12, y=0.19, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="Ly6G_CD11b+_CD86_CD103-_")
text(x=-1.60, y=0.19, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="Ly6G_Ly6C-_CD11b+_CD86_")
text(x=-1.20, y=0.19, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="CD11b+_CD11c-_CD86_CD103-")
text(x=-0.85, y=0.19, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="Ly6C-_CD11b+_CD11c-_CD86_")
text(x=-0.5, y=0.19, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="CD11b+_CD11c-_MHCII+_CD86_")
text(x=-0.2, y=0.19, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="CD11b+_NotCD137_CD11c-_CD86_")
text(x=0.32, y=0.19, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="Ly6C-_CD11b+_NotCD137_CD86_")
text(x=0.72, y=0.19, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="CD11b+_NotCD137_CD86_CD103-_")
text(x=1.15, y=0.19, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="NotLin_CD11b+_CD11c-_CD86_")
text(x=1.52, y=0.19, lab=".", col="black")



setSVGShapeToolTip(title=NULL, desc="Ly6G_CD11b+_NotCD137_CD11c-_CD86_")
text(x=-2.12, y=-0.01, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="Ly6G_CD11b+_CD11c-_CD86_CD103-_")
text(x=-1.7, y=-0.01, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="Ly6G_Ly6C-_CD11b+_CD11c-_CD86_")
text(x=-1.3, y=-0.01, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="CD11b+_CD11c-_MHCII+_CD86_CD103-")
text(x=-0.9, y=-0.01, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="Ly6C-_CD11b+_CD11c-_MHCII+_CD86_")
text(x=-0.5, y=-0.01, lab=".", col="black")



setSVGShapeToolTip(title=NULL, desc="Ly6G_CD11b+_CD11c-_MHCII+_CD86_")
text(x=-0.1, y=-0.01, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="CD11b+_NotCD137_CD11c-_MHCII+_CD86_")
text(x=0.32, y=-0.01, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="Ly6C-_CD11b+_NotCD137_CD11c-_CD86_")
text(x=0.76, y=-0.01, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="CD11b+_NotCD137_CD11c-_CD86_CD103-")
text(x=1.19, y=-0.01, lab=".", col="black")

setSVGShapeToolTip(title=NULL, desc="NotLin_Ly6G_CD11b+_CD11c-_CD86_")
text(x=1.61, y=-0.01, lab=".", col="black")

dev.off()
