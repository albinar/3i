# Developed by Albina Rahim

# This code creates RchyOptimyx plots from Matrix3iTcell.csv.
remove(list=ls())

# setwd("/home/arahim/Desktop/3iTcell_AFC_Oct_2014/")
setwd("/data/Panel_M-cell")

library("flowCore")
library("flowBin")
library("flowDensity")
library("flowType")
library("flowClean")
library("RchyOptimyx")
source("Codes/3iTcellfunctions-M.R")

start <- Sys.time()

load(file =  paste("Results/Phenotypes.Rdata",sep=""))
load(file =  paste("Results/Genotypes.Rdata",sep=""))
load(file =  paste("Results/GenotypesLong.Rdata",sep=""))
load(file =  paste("Results/allFT.Rdata",sep=""))
load(file =  paste("Results/flowType.res_Sample.Rdata",sep=""))
load("Results/GlobalFrame657.Rdata")
load(file =  paste("Results/Matrix3iTcell.Rdata",sep=""))


dir.create ( "Results/Figures/AdjPvalues")
dir.create ( "Results/Figures/SigPvalues")
dir.create ( "Results/Figures/RchyOptimyx")
dir.create ( "Results/Figures/RchyOptimyxTrim4")
dir.create ( "Results/Figures/RchyOptimyxTrim3")
dir.create ( "Results/Figures/RchyOptimyx_High_Proportion")
dir.create ( "Results/Figures/PropVsProp")
dir.create ( "Results/Figures/BoxPlotProportions")



# dir.create ( "Results/Figures/RchyOptimyx_Non_Adj")
dir.create ( "Results/Figures/RchyOptimyx_Gating_Strategy_Pop")
# dir.create ( "Results/Figures/RchyOptimyx_Gating_Strategy_Pop_Non_Adj")

# if(!exists("Matrix3iTcell")) {
#     print("loading")
#     Matrix3iTcell <- read.csv("Results/Matrix3iTcell.csv")
# }

TotalAdjPvals <- NULL
TotalPvals <- NULL

AverageCellProportions <- rowSums(Matrix3iTcell)/dim(Matrix3iTcell)[2]

uniqueWT <- unique(Genotypes)

nameGenotypes <- (Genotypes[duplicated(Genotypes)])
nameGenotypes <- nameGenotypes[which(nameGenotypes == uniqueWT[4])[1]:length(nameGenotypes)]
nameGenotypes <- unique(nameGenotypes)

WTindices <- c(which(Genotypes == uniqueWT[1]), which(Genotypes == uniqueWT[2]) )

flowType.resWT <- NULL
flowType.resWT_Avg <- matrix(0, dim(Matrix3iTcell)[1], 1 )

for ( r1 in 1:length(WTindices)){
    load (file =  allFT[WTindices[r1]] )
    flowType.resWT <- rbind(flowType.resWT, flowType.res@CellFreqs/(flowType.res@CellFreqs[1]))
    flowType.resWT_Avg <- flowType.resWT_Avg + flowType.res@CellFreqs
}
flowType.resWT_Avg <- flowType.resWT_Avg / length(WTindices)


SignifNames <- list()
count <- 0
countSkip <- 0
countZero <- 0

###############################################################################################################
# Phenotypes names for populations in the gating strategy

# Get names of populations which arise in flowType data
Ly6Gneg <- paste0(Phenotypes[3], Phenotypes[4])
eosin3 <- paste0(Phenotypes[3],Phenotypes[5])    
neu2 <- Phenotypes[2]

mono.Ly6c.hi <- paste0(Ly6Gneg, Phenotypes[7], Phenotypes[9])
pDC <- paste0(Ly6Gneg, Phenotypes[6], Phenotypes[9], Phenotypes[11])
# by plotting populations I found all cDC comes from the Ly6C- portion of the Ly6Gneg population
cDC <- paste0(Ly6Gneg, Phenotypes[8], Phenotypes[13], Phenotypes[15])   
CD8A.TypeDC <- paste0(cDC, Phenotypes[17])
CD11bp.CD86lo <- paste0(cDC, Phenotypes[16])
CD103p.dc <- paste0(CD8A.TypeDC, Phenotypes[19])

# Note: the NOT(mono.Ly6c.hi) and NOT(cDC) populations come from multiple flowType populations
# But if either of these populations change as a function of the starting Lin-Mac-, so will their inverse
# (mono.Ly6chi or cDC) which only comes from one flowType population

# Contains all indices for the gating strategy
GatingStrategyPop0 <- NULL
GatingStrategyPop0[1] <- as.numeric( which(Phenotypes == "") )  # Lin-Mac-
GatingStrategyPop0[2]  <- as.numeric(  which(Phenotypes == Ly6Gneg) ) 
GatingStrategyPop0[3]  <- as.numeric(  which(Phenotypes == eosin3) ) 
GatingStrategyPop0[4]  <- as.numeric(  which(Phenotypes == neu2) ) 
GatingStrategyPop0[5]  <- as.numeric(  which(Phenotypes == mono.Ly6c.hi) ) 
GatingStrategyPop0[6] <- as.numeric(  which(Phenotypes == pDC) ) 
GatingStrategyPop0[7] <- as.numeric(  which(Phenotypes == cDC) ) 
GatingStrategyPop0[8] <- as.numeric(  which(Phenotypes == CD8A.TypeDC) ) 
GatingStrategyPop0[9] <- as.numeric(  which(Phenotypes == CD11bp.CD86lo) ) 
GatingStrategyPop0[10] <- as.numeric(  which(Phenotypes == CD103p.dc) ) 

################################################################################################################


ForLoopIndices <- 3:length(uniqueWT)
#ForLoopIndices <- 9:50

for ( k1 in ForLoopIndices){

    print(paste(k1,": Starting", uniqueWT[k1]))
    KOindices <- c(which(Genotypes == uniqueWT[k1]))
#     if(length(KOindices) <= 1){
    if(length(KOindices) <= 2){
#         print("Skipped, only one file")
        print("Skipped, 2 or less files")
        countSkip <- countSkip + 1
        next
    }
    
    Labels <- c(WTindices, KOindices)
    
    WCpValues   <- vector()


    load(file=paste0("Results/WilcoxTest/", uniqueWT[k1], ".Rdata", sep=""))
#     load(file=paste0("Results/T_Test/", uniqueWT[k1], ".Rdata", sep=""))
#     
#     WCpValues <- Pvals
    
#     WCpValues <- WCpValues[1:19683]  # this was for when I added p-values which I no longer want
#     WCpValues[is.nan(WCpValues)] <- 1

    
   # KOindices <- c(which(Genotypes == uniqueWT[k1]))
    
    flowType.resKO <- NULL
    flowType.resKO_Avg <- matrix(0, dim(Matrix3iTcell)[1], 1 )
#     flowType.resKO_Avg <- matrix(0, calcNumPops(c(3,2,2,2,2,2,2,2,3,2), 10), 1 )
    
    for ( r1 in 1:length(KOindices)){
        load (file =  allFT[KOindices[r1]] )
#         flowType.resKO[[r1]] <- flowType.res@CellFreqs
        flowType.resKO <- rbind(flowType.resKO, flowType.res@CellFreqs/(flowType.res@CellFreqs[1]))
        flowType.resKO_Avg <- flowType.resKO_Avg + flowType.res@CellFreqs
    }
    flowType.resKO_Avg <- flowType.resKO_Avg / length(KOindices)
  
        
        GatingStrategyPop <- NULL
        GatingStrategyPop <- GatingStrategyPop0


#         GatingStrategyPop[39] <- as.numeric(  which(PhenotypesB == "TCRb+CD161-TCRd-") )

        Parental <- NULL
        Parental[122]  <- as.numeric(  which(Phenotypes == "" ) )
        
#         TotalAdjPvals <- rbind(TotalAdjPvals, )
    
        selected1 <- intersect ( which(flowType.resWT_Avg <= 500 ), which(flowType.resKO_Avg <= 500 ))
#         selected2 <- union ( which(flowType.resWT_Avg <= 500 ), which(flowType.resKO_Avg <= 500 ))
#         selected3 <- which(AverageCellProportions <= 500/300000 )

#         selected1 <- selected1[setdiff(selected1, GatingStrategyPop[c(3:38)])]
        
        # SD: I need this to get rid of some populations that I added manually to the flowType populations
#         upper.index <-calcNumPops(c(2,2,2,2,2,2,2,2,2), 9)
#         selected1 <- union(selected1, ((upper.index + 1):length(flowType.resWT_Avg)))
#         
        selected1 <- setdiff(selected1, c(GatingStrategyPop[c(1:10)],Parental[k1]) )
        WCpValuesB <- WCpValues[-selected1]
        flowType.res_PhenoCodes_Reduced <- flowType.res@PhenoCodes[-selected1]
        PhenotypesB <- Phenotypes[-selected1]
    
 #      AdjWCpValuesBH <- p.adjust(WCpValuesB, method= "BH")
        #AdjWCpValuesBH <- p.adjust(WCpValuesB, method= "bonferroni")
        AdjWCpValuesBH <- p.adjust(WCpValuesB, method= "fdr")
#         print(min(p.adjust(WCpValuesB, method= "bonferroni")))
#         print(min(p.adjust(WCpValuesB, method= "fdr")))
#         print(min(p.adjust(WCpValuesB, method= "BH")))
        #selectedSignif <- which(AdjWCpValuesBN < 0.05)
        selectedSignif <- which(AdjWCpValuesBH < 0.05)
        
    
        if (length(selectedSignif) == 0 ) {
            print("Zero selected: Skipping RchyOptimyx Plots")
            countZero <- countZero + 1
            next;
        }
        cat("\n")
    
        count <- count + 1
        SignifNames[[count]] <- names(AdjWCpValuesBH[selectedSignif])
    
        write.table(cbind(names(AdjWCpValuesBH[selectedSignif]), AdjWCpValuesBH[selectedSignif]), paste0("Results/Figures/SigPvalues/SigPvalues_", uniqueWT[k1], ".csv"), sep=",", col.names=F, row.names=F)
    
        #trim 4
        NumberOfTags <- NULL
        for ( x1 in 1:length(selectedSignif)){
            NumberOfTags[x1] <- length( gregexpr("[+-]", names(AdjWCpValuesBH[selectedSignif[x1]]))[[1]] )
        }
        selectedTrim4 <- selectedSignif[1:length(which(NumberOfTags <= 3))]
        #trim 3
        selectedTrim3 <- selectedSignif[1:length(which(NumberOfTags <= 2))]
    
        AverageCellProportionsB <- AverageCellProportions[-selected1]
        # find the indices for the 15 significant populations that have the highest proportion
        if ( length(selectedSignif) >= 15 ) {
            temp <- sort(AdjWCpValuesBH, index.return=T)
            selectedTrimProportion <- temp$ix[1:length(selectedSignif)]    
            temp <- sort(AverageCellProportionsB[selectedTrimProportion], decreasing=T, index.return=T)
            selectedTrimProportion <- selectedTrimProportion[temp$ix][1:15] 
            
            # find the indices for the 15 populations with the smallest P-values. It they are tied then the ones with less markers will be automatically used.
            temp <- sort(AdjWCpValuesBH, index.return=T)
            selected15 <- temp$ix[1:15]    
#             selected <- selected[which(AdjWCpValuesBH[selected] <= sort(AdjWCpValuesBH[selected], decreasing=F)[10])]
            
        } else {
            selected15 <- selectedSignif
            selectedTrimProportion <- selectedSignif
        }
    
        flowType.resKO_AvgTemp <- flowType.resKO_Avg[-selected1]/max(flowType.resKO_Avg[-selected1])
        flowType.resWT_AvgTemp <- flowType.resWT_Avg[-selected1]/max(flowType.resWT_Avg[-selected1])
        pdf ( file = paste("Results/Figures/PropVsProp/PropVsProp_", uniqueWT[k1], ".pdf", sep = "" ) )
        plot(flowType.resKO_AvgTemp, flowType.resWT_AvgTemp, cex=0.2, pch=19, 
             xlab="Proportion of Knock-Out", ylab="Proportion of Wild-Type", main=paste0("Proportion vs Proportion: ", uniqueWT[k1]))
        points(flowType.resKO_AvgTemp[selectedSignif],flowType.resWT_AvgTemp[selectedSignif], col="blue" , cex=0.6, pch=19)
        points(flowType.resKO_AvgTemp[selected15]    ,flowType.resWT_AvgTemp[selected15]    , col="green", cex=0.8, pch=1)
        points(flowType.resKO_AvgTemp[selectedTrimProportion]    ,flowType.resWT_AvgTemp[selectedTrimProportion]    , col="orange" , cex=0.8, pch=1)
        selectedTemp <- intersect(selectedTrimProportion, selected15)
        if ( length(selectedTemp) != 0 ){
            points(flowType.resKO_AvgTemp[selectedTemp]    ,flowType.resWT_AvgTemp[selectedTemp], col="red" , cex=0.8, pch=1)
        }
        abline(0,1)
        dev.off()
    
        dude1 <- as.numeric(AdjWCpValuesBH)
#         dude1B <- dude1[-selected1]
#         AverageCellProportionsB <- AverageCellProportions[-selected1]
        pdf ( file = paste("Results/Figures/AdjPvalues/AdjPvalues_", uniqueWT[k1], ".pdf", sep = "" ) )
        
        plot(AverageCellProportionsB, log10(dude1),  col=1, main=paste(uniqueWT[k1], ": ", length(which(dude1<=0.05))," Signif Pops", sep=""), 
             xlab="Proportion", ylab="log10(AdjPvalues)", cex=0.2, pch=19)
    
        abline ( h = log10(0.05), col="black")
#         abline ( v = 500/300000, col="green", cex=40)
        points(AverageCellProportionsB[selectedSignif],log10(dude1)[selectedSignif], col="blue", cex=0.2, pch=19)
        points(AverageCellProportionsB[selected15],log10(dude1)[selected15], col="green", cex=0.4, pch=1)
        points(AverageCellProportionsB[selectedTrimProportion],log10(dude1)[selectedTrimProportion], col="orange", cex=0.4, pch=1)
        selectedTemp <- intersect(selectedTrimProportion, selected15)
        if ( length(selectedTemp) != 0 ){
            points(AverageCellProportionsB[selectedTemp], log10(dude1)[selectedTemp], col="red" , cex=0.4, pch=1)
        }

#         text(1000,log10(min(dude1))/4, length(which(dude1<=0.05)), col="red")
    
        dev.off()
    
    
    
        AdjPvals <- AdjWCpValuesBH
    
        #trim 10 (all 9 markers used)
        res <-RchyOptimyx(pheno.codes=flowType.res_PhenoCodes_Reduced, phenotypeScores=-log10(AdjPvals), startPhenotype=flowType.res_PhenoCodes_Reduced[selected15[1]], 1, FALSE,trim.level=8)
        
        if(length(selected15) > 1){
          for (i in 2:length(selected15)){
            temp<-RchyOptimyx(pheno.codes=flowType.res_PhenoCodes_Reduced, phenotypeScores=-log10(AdjPvals), startPhenotype=flowType.res_PhenoCodes_Reduced[selected15[i]], 1, FALSE,trim.level=8)
            res=merge(res,temp)
          }
        }
        pdf ( file = paste("Results/Figures/RchyOptimyx/Default_", uniqueWT[k1], ".pdf", sep=""), height= 18, width = 28)
        PropMarkers <- c(4,9,13,10,15,7,16,18,17)
        plot(res, phenotypeScores=-log10(AdjPvals), phenotypeCodes=flowType.res_PhenoCodes_Reduced, marker.names=as.vector(g@parameters@data$desc[PropMarkers]), ylab='-log10(Pvalue)', max.score = max(-log10(AdjPvals)))
        dev.off();
    
        #trim high proportion
        res <-RchyOptimyx(pheno.codes=flowType.res_PhenoCodes_Reduced, phenotypeScores=-log10(AdjPvals), startPhenotype=flowType.res_PhenoCodes_Reduced[selectedTrimProportion[1]], 1, FALSE,trim.level=9)
        if(length(selectedTrimProportion) > 1){
          for (i in 2:length(selectedTrimProportion)){
            temp<-RchyOptimyx(pheno.codes=flowType.res_PhenoCodes_Reduced, phenotypeScores=-log10(AdjPvals), startPhenotype=flowType.res_PhenoCodes_Reduced[selectedTrimProportion[i]], 1, FALSE,trim.level=9)
            res=merge(res,temp)
          }
        }
        pdf ( file = paste("Results/Figures/RchyOptimyx_High_Proportion/Hi_Prop_", uniqueWT[k1], ".pdf", sep=""), height= 18, width = 28)
        plot(res, phenotypeScores=-log10(AdjPvals), phenotypeCodes=flowType.res_PhenoCodes_Reduced, marker.names=as.vector(g@parameters@data$desc[PropMarkers]), ylab='-log10(Pvalue)', max.score = max(-log10(AdjPvals)))
        dev.off();

        #trim high proportion - part 2 - using proportions as value
        pdf ( file = paste("Results/Figures/RchyOptimyx_High_Proportion/Hi_Prop_", uniqueWT[k1], "Proportion.pdf", sep=""), height= 18, width = 28)
        plot(res, phenotypeScores=AverageCellProportionsB, phenotypeCodes=flowType.res_PhenoCodes_Reduced, marker.names=as.vector(g@parameters@data$desc[PropMarkers]), ylab='Proportion', max.score = 1)
        dev.off();
    
        #trim 4
        res <-RchyOptimyx(pheno.codes=flowType.res_PhenoCodes_Reduced, phenotypeScores=-log10(AdjPvals), startPhenotype=flowType.res_PhenoCodes_Reduced[selectedTrim4[1]], 1, FALSE,trim.level=4)
    
        if(length(selectedTrim3) >1){
          for (i in 2:length(selectedTrim4)){
            temp<-RchyOptimyx(pheno.codes=flowType.res_PhenoCodes_Reduced, phenotypeScores=-log10(AdjPvals), startPhenotype=flowType.res_PhenoCodes_Reduced[selectedTrim4[i]], 1, FALSE,trim.level=4)
            res=merge(res,temp)
          }
        }
        pdf ( file = paste("Results/Figures/RchyOptimyxTrim4/Trim4_", uniqueWT[k1], ".pdf", sep=""), height= 18, width = 28)
        plot(res, phenotypeScores=-log10(AdjPvals), phenotypeCodes=flowType.res_PhenoCodes_Reduced, marker.names=as.vector(g@parameters@data$desc[PropMarkers]), ylab='-log10(Pvalue)')#, max.score = max(-log10(AdjPvals[selectedOrg])))
        dev.off();
    
        #trim 3
        res <-RchyOptimyx(pheno.codes=flowType.res_PhenoCodes_Reduced, phenotypeScores=-log10(AdjPvals), startPhenotype=flowType.res_PhenoCodes_Reduced[selectedTrim3[1]], 1, FALSE,trim.level=4)
    
        if(length(selectedTrim3) >1){
          for (i in 2:length(selectedTrim3)){
            temp<-RchyOptimyx(pheno.codes=flowType.res_PhenoCodes_Reduced, phenotypeScores=-log10(AdjPvals), startPhenotype=flowType.res_PhenoCodes_Reduced[selectedTrim3[i]], 1, FALSE,trim.level=4)
            res=merge(res,temp)
          }
        }
        pdf ( file = paste("Results/Figures/RchyOptimyxTrim3/Trim3_", uniqueWT[k1], ".pdf", sep=""), height= 18, width = 28)
        plot(res, phenotypeScores=-log10(AdjPvals), phenotypeCodes=flowType.res_PhenoCodes_Reduced, marker.names=as.vector(g@parameters@data$desc[PropMarkers]), ylab='-log10(Pvalue)')#, max.score = max(-log10(AdjPvals[selectedOrg])))
        dev.off();
    
        GatingStrategyPop <- NULL
        GatingStrategyPop[1] <- as.numeric( which(Phenotypes == "") )  # Lin-Mac-
        GatingStrategyPop[2]  <- as.numeric(  which(PhenotypesB == Ly6Gneg) ) 
        GatingStrategyPop[3]  <- as.numeric(  which(PhenotypesB == eosin3) ) 
        GatingStrategyPop[4]  <- as.numeric(  which(PhenotypesB == neu2) ) 
        GatingStrategyPop[5]  <- as.numeric(  which(PhenotypesB == mono.Ly6c.hi) ) 
        GatingStrategyPop[6] <- as.numeric(  which(PhenotypesB == pDC) ) 
        GatingStrategyPop[7] <- as.numeric(  which(PhenotypesB == cDC) ) 
        GatingStrategyPop[8] <- as.numeric(  which(PhenotypesB == CD8A.TypeDC) ) 
        GatingStrategyPop[9] <- as.numeric(  which(PhenotypesB == CD11bp.CD86lo) ) 
        GatingStrategyPop[10] <- as.numeric(  which(PhenotypesB == CD103p.dc) ) 

        AdjPvals <- AdjWCpValuesBH
#     AdjPvals <- AdjHighest
    
        # Populations for the Gating Strategy part 1
        GSPselected <- c(1:10)

        res <-RchyOptimyx(pheno.codes=flowType.res_PhenoCodes_Reduced, phenotypeScores=-log10(AdjPvals), startPhenotype=flowType.res_PhenoCodes_Reduced[ GatingStrategyPop[GSPselected[1]] ], 1, FALSE,trim.level=9)
        
        for (i in 2:length(GatingStrategyPop[GSPselected])){
            temp<-RchyOptimyx(pheno.codes=flowType.res_PhenoCodes_Reduced, phenotypeScores=-log10(AdjPvals), startPhenotype=flowType.res_PhenoCodes_Reduced[ GatingStrategyPop[GSPselected[i]] ], 1, FALSE,trim.level=9)
            res=merge(res,temp)
        }
    
        maxScoreIndices <- NULL
        for ( p1 in 1:length(res@nodes[1,])){
            if ( length(which(res@nodes[1,p1]==flowType.res_PhenoCodes_Reduced)) == 0){
                 maxScoreIndices[p1] <-  maxScoreIndices[p1-1]
            } else {
                maxScoreIndices[p1] <- which(res@nodes[1,p1]==flowType.res_PhenoCodes_Reduced)
            }
        }
    
        pdf ( file = paste("Results/Figures/RchyOptimyx_Gating_Strategy_Pop/GatingSP_", uniqueWT[k1], "_1st_Half.pdf", sep=""), height= 18, width = 28)
        plot(res, phenotypeScores=-log10(AdjPvals), phenotypeCodes=flowType.res_PhenoCodes_Reduced, marker.names=as.vector(g@parameters@data$desc[PropMarkers]), 
             ylab='-log10(Pvalue)', max.score = c(0,max(-log10(AdjPvals[maxScoreIndices] ))) )
#              ylab='-log10(Pvalue)', max.score = c(1.4) )
        dev.off();
    

#         # Populations for the Gating Strategy part 2
#         GSPselected <- c(4:21)
# 
#         res <-RchyOptimyx(pheno.codes=flowType.res_PhenoCodes_Reduced, phenotypeScores=-log10(AdjPvals), startPhenotype=flowType.res_PhenoCodes_Reduced[ GatingStrategyPop[GSPselected[1]] ], 1, FALSE,trim.level=9)
#         
#         for (i in 2:length(GatingStrategyPop[GSPselected])){
#             temp<-RchyOptimyx(pheno.codes=flowType.res_PhenoCodes_Reduced, phenotypeScores=-log10(AdjPvals), startPhenotype=flowType.res_PhenoCodes_Reduced[ GatingStrategyPop[GSPselected[i]] ], 1, FALSE,trim.level=9)
#             res=merge(res,temp)
#         }
#         pdf ( file = paste("Results/Figures/RchyOptimyx_Gating_Strategy_Pop/GatingSP_", uniqueWT[k1], "_2nd_Half.pdf", sep=""), height= 18, width = 28)
#         plot(res, phenotypeScores=-log10(AdjPvals), phenotypeCodes=flowType.res_PhenoCodes_Reduced, marker.names=as.vector(g@parameters@data$desc[PropMarkers]), 
#              ylab='-log10(Pvalue)', max.score = c(0,max(-log10(AdjPvals[maxScoreIndices] ))) )
# #              ylab='-log10(Pvalue)', max.score = max(-log10(AdjPvals)))
#         dev.off();
    
    
    allFT2 <- allFT
    for ( p4 in 1:length(allFT2)) {
      temp.name <- strsplit(allFT2[p4],split="/")
      temp.name <- temp.name[[1]][length(unlist(temp.name[[1]]))]
      temp.name <- strsplit(temp.name,split="[.]")[[1]][1]
      allFT2[p4] <- strsplit(temp.name,split="_")[[1]][3]
    }
    GenderKOindices <- allFT2[KOindices]
    
#     CSVfile <- read.csv("/home/arahim/Desktop/3iTcell/data_output-99.csv", check.names = FALSE)
#    CSVfile <- read.csv("/home/arahim/Desktop/3iTcell_AFC_Oct_2014/FCS_Oct_2014/attachments/Spleen_Data_III.csv", check.names = FALSE)

    load("Results/CSVfile.Rdata")
    #CSVfile <- as.matrix(CSVfile)
    
    #GenderData <- CSVfile[,c(12,20)]
    GenderData <- CSVfile[,c('Gender','Label.Barcode')]
    GenderCrop <- NULL
    for ( p5 in 1:length(GenderKOindices)){
        #GenderCrop[p5] <- which(GenderData[,2] == GenderKOindices[p5])
      GenderCrop[p5] <- min(grep(GenderKOindices[p5], GenderData[,2]))      # I added the min because CSVfile contains seemingly duplicate entries (same barcodes, gender, assay date ect)
    }
    GenderData <- GenderData[GenderCrop,]
    
    GenderFemaleIndices <- which(tolower(GenderData[,1]) == "female")
    GenderMaleIndices   <- which(tolower(GenderData[,1]) == "male")
    
    
    selectedBoxPlot <- selectedTrimProportion
    AdjPvalsBoxPlot <- AdjWCpValuesBH
    
    flowType.resWTB <- flowType.resWT[,-selected1]
    flowType.resKOB <- flowType.resKO[,-selected1]
    
#     TempIndices <- c(1,2,3)
#     selectedTrimProportion

    Parental <- NULL
#     Parental[23]   <- as.numeric(  which(PhenotypesB == "KLRG1-" ) )
#     Parental[26]   <- as.numeric(  which(PhenotypesB == "TCRd+" ) )
#     Parental[29]   <- as.numeric(  which(PhenotypesB == "CD62L-" ) )
#     Parental[33]   <- as.numeric(  which(PhenotypesB == "CD8a-" ) )
#     Parental[41]   <- as.numeric(  which(PhenotypesB == "TCRb+CD8a-CD4-" ) )
#     Parental[67]   <- as.numeric(  which(PhenotypesB == "CD62L+CD4-" ) )
#     Parental[69]   <- as.numeric(  which(PhenotypesB == "TCRd+" ) )
#     Parental[76]   <- as.numeric(  which(PhenotypesB == "CD62L+TCRb+" ) )
#     Parental[88]   <- as.numeric(  which(PhenotypesB == "CD25-CD4-" ) )
#     Parental[188]  <- as.numeric(  which(PhenotypesB == "CD4+" ) )         # 9 alternative
#     Parental[122]  <- as.numeric(  which(PhenotypesB == "CD62L-TCRb+" ) )

    Parental[140]  <- as.numeric(  which(PhenotypesB == "" ) )
    

    
 #   c(13,27,33,38,54,69,72,78,87,98,99,100,107,112,113,123,140)
    
#     if( q2 == 123 ) {
#         selectedBoxPlot <- c(Parental[k1], Parental[k1+1], Parental[k1],Parental[k1], selectedBoxPlot)
#     } else {
#         selectedBoxPlot <- c(Parental[k1], Parental[k1],Parental[k1], selectedBoxPlot)
#     }
#     
    boxplotDataWT <- flowType.resWTB[,selectedBoxPlot]
    boxplotDataKO <- flowType.resKOB[,selectedBoxPlot]
     
    ylimMin <- min(c(boxplotDataWT,boxplotDataKO))
    ylimMax <- max(c(boxplotDataWT,boxplotDataKO))

    
    png ( file = paste("Results/Figures/BoxPlotProportions/Box_Plot_", uniqueWT[k1],".png", sep = "" ), height=1000, width=1600 )
    par(mar=c(20,5,8,2)) # margins
#         boxplot(dude[,], cex=0.2, pch=19,  col=rep(c("Blue","red"),33), main=paste("Box Plots: ", PlotNames[w2], sep=""), ylab=PlotYaxis[w2], axes=FALSE, frame.plot=FALSE,las=2,na.action=na.omit)
        boxplot(boxplotDataWT, cex=0.2, pch=19,  col=rep(c("gray83"),length(selectedBoxPlot)), axes=FALSE, frame.plot=FALSE,las=2,na.action=na.omit, ylim=c(ylimMin,ylimMax))
        if(length(selectedBoxPlot) > 1){
          for ( p6 in 1:length(selectedBoxPlot)) {
            if (length(GenderFemaleIndices) >=1) {
              stripchart(as.list(boxplotDataKO[GenderFemaleIndices,p6]), vertical = TRUE, method = "jitter", pch = 21, 
                         col = "black", bg = "green", add = TRUE, at=c(rep(p6,length(GenderFemaleIndices))) )    
            }
            if (length(GenderMaleIndices) >=1) {
              stripchart(as.list(boxplotDataKO[GenderMaleIndices,p6]), vertical = TRUE, method = "jitter", pch = 21, 
                         col = "black", bg = "deepskyblue", add = TRUE, at=c(rep(p6,length(GenderMaleIndices))) )     
            }
          }
        }else{
            if (length(GenderFemaleIndices) >=1) {
              stripchart(as.list(boxplotDataKO[GenderFemaleIndices]), vertical = TRUE, method = "jitter", pch = 21, 
                         col = "black", bg = "green", add = TRUE, at=c(rep(1,length(GenderFemaleIndices))) )    
            }
            if (length(GenderMaleIndices) >=1) {
              stripchart(as.list(boxplotDataKO[GenderMaleIndices]), vertical = TRUE, method = "jitter", pch = 21, 
                         col = "black", bg = "deepskyblue", add = TRUE, at=c(rep(1,length(GenderMaleIndices))) )     
            }
        }
        tempPhenotypes <- PhenotypesB[selectedBoxPlot]

        axis(1,at=1:length(selectedBoxPlot), labels=tempPhenotypes,las=2)
        axis(2,at=seq(0, 1, by=0.01) )
        tempPvals <- round(AdjPvalsBoxPlot[selectedBoxPlot], digits=4)

        axis(3,at=1:length(selectedBoxPlot), labels=tempPvals, las=2)
        mtext(side=1, line=18, "Phenotypes")
        mtext(side=2, line=2.5, "Proportion")
        mtext(side=3, line=5, "P-values")
#         if( q2 == 123 ) {
#             legend(x=2.0, y=mean(c(ylimMin,ylimMax)), legend=c("Wild-type","Female Knock-out","Male Knock-out"), fill=c("gray83","green","deepskyblue"))
#         } else {
#             legend(x=1.55, y=mean(c(ylimMin,ylimMax)), legend=c("Wild-type","Female Knock-out","Male Knock-out"), fill=c("gray83","green","deepskyblue"))
#         }
#         legend(x=1.0,y=c(ylimMin+0.5), legend=c("Wild-type","Female Knock-out","Male Knock-out"), fill=c("gray83","green","deepskyblue"))
    dev.off()

#     Chart_JM <- flowType.res@CellFreqs
#     load(file =  paste("Results/Phenotypes.Rdata",sep=""))
#     Chart_JM <- rbind(Phenotypes, Chart_JM)

    
    
    cat("Time is: ",TimeOutput(start),"\n",sep="")

}

cat("Time is: ",TimeOutput(start),"\n",sep="")

# write.table(t(nameGenotypes), paste("Results/TotalAdjPvals_nameGenotypes.csv", sep=""), sep=",", row.names=F, col.names=F)


# pValue vs CellFreq
# plot(flowType.res@CellFreqs, AdjPvals, pch=19, cex=0.02, xlim=c(0,15000))
# dude <- AdjPvals[which(flowType.res@CellFreqs>=10000)]
# dude[which(dude<=0.5)]

SignifNamesTemp <- unlist(SignifNames)

SignifNamesTemp <- SignifNamesTemp[duplicated(SignifNamesTemp)]; print(length(SignifNamesTemp))
SignifNamesTemp <- SignifNamesTemp[duplicated(SignifNamesTemp)]; print(length(SignifNamesTemp))
SignifNamesTemp <- SignifNamesTemp[duplicated(SignifNamesTemp)]; print(length(SignifNamesTemp))
SignifNamesTemp <- SignifNamesTemp[duplicated(SignifNamesTemp)]; print(length(SignifNamesTemp))

NumberOfCommon <- matrix(0, length(SignifNamesTemp), length(SignifNames))

for ( t1 in 1: length(SignifNamesTemp) ) {
    for ( t2 in 1:length(SignifNames) ) {
        NumberOfCommon[t1,t2] <- length(which(SignifNames[[t2]] == SignifNamesTemp[t1]))
    }
}

# rbind(uniqueWT[c(23,26,29,33,41,67,69,76,88,122)], colSums(NumberOfCommon))

# Common_Signif_Phenotypes <- cbind(c("Significant Knock-out","All Phenotypes",SignifNamesTemp),rbind(uniqueWT[c(23,26,29,33,41,67,69,76,88,122)], colSums(NumberOfCommon),NumberOfCommon) )
# Common_Signif_Phenotypes <- cbind(c("Significant Knock-out","All Phenotypes",SignifNamesTemp),rbind(uniqueWT[c(13,27,33,38,54,69,72,78,87,98,99,100,107,112,113,123,140)], colSums(NumberOfCommon),NumberOfCommon) )
# write.table(Common_Signif_Phenotypes, paste0("Results/Figures/Common_Signif_Phenotypes.csv"), sep=",", col.names=F, row.names=F)
