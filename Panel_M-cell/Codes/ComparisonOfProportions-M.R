# Originally written by Albina Rahim 
# Modified by Sibyl Drissler
# Last modified: June 14, 2016


# Returns: /Results/FinalTableOfProportions.csv. 
# Eight rows per file with the following data for each cell population (genotype)
# %cd45 (fD), %cd45 (fT), %total (fD), %total (fT), %parent (fD), %parent (fT), cell count (fD), cell count (fT)
# fD = supervised flowDensity gating, fT = flowType gating

setwd("/data/Panel_M-cell")

load( file =  paste("Results/Phenotypes.Rdata",sep=""))
load( file =  paste("Results/Genotypes.Rdata",sep=""))
load( file =  paste("Results/GenotypesLong.Rdata",sep=""))
load( file =  paste("Results/allFT.Rdata",sep=""))

FinalTable <- NULL

# load the supervised flowDensity gating data
CSVfile <- read.csv("Results/Events_Proportions_Table_Updated.csv", check.names = FALSE, stringsAsFactors = F)
#colnames(CSVfile)[c(7, 15, 17, 19, seq(23, ncol(CSVfile), 2))]
cell.count.fD <- CSVfile[, c(7, 15, 17, 19, seq(23, ncol(CSVfile), 2))]
percent.parent.fD <- CSVfile[, c(8, 16, 18, 20, seq(24, ncol(CSVfile), 2))]

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
GatingStrategyPop <- NULL
GatingStrategyPop[1] <- as.numeric( which(Phenotypes == "") )  # Lin-Mac-
GatingStrategyPop[2]  <- as.numeric(  which(Phenotypes == Ly6Gneg) ) 
GatingStrategyPop[3]  <- as.numeric(  which(Phenotypes == neu2) ) 
GatingStrategyPop[4]  <- as.numeric(  which(Phenotypes == eosin3) ) 
GatingStrategyPop[5]  <- as.numeric(  which(Phenotypes == mono.Ly6c.hi) ) 
GatingStrategyPop[6] <- as.numeric(  which(Phenotypes == pDC) ) 
GatingStrategyPop[7] <- as.numeric(  which(Phenotypes == cDC) ) 
GatingStrategyPop[8] <- as.numeric(  which(Phenotypes == CD8A.TypeDC) ) 
GatingStrategyPop[9] <- as.numeric(  which(Phenotypes == CD11bp.CD86lo) ) 
GatingStrategyPop[10] <- as.numeric(  which(Phenotypes == CD103p.dc) ) 
################################################################################################################

#for ( q2 in 1:2){
for ( q2 in 1: length(allFT)){

  load(file = allFT[q2])
  
  # row of CSVfile corresponding to allFT[q2]
  csv.ind <- grep(sub(".*?SPLN (.*?)(.Rdata).*", "\\1", allFT[q2]), CSVfile[, 3]) 
  
  # calculate cell counts for NOT(x) populations, where x is a cell population calculated in flowType
  not.mono.Ly6c.hi.cell.count <- flowType.res@CellFreqs[GatingStrategyPop[2]] - flowType.res@CellFreqs[GatingStrategyPop[5]]
  not.pDC.cell.count <- not.mono.Ly6c.hi.cell.count - flowType.res@CellFreqs[GatingStrategyPop[6]] 
  misc.cell.count <- not.pDC.cell.count -  flowType.res@CellFreqs[GatingStrategyPop[7]] 
  
  cell.count.fT <- c(flowType.res@CellFreqs[GatingStrategyPop[1:5]], not.mono.Ly6c.hi.cell.count,
                     flowType.res@CellFreqs[GatingStrategyPop[6]], not.pDC.cell.count, 
                     flowType.res@CellFreqs[GatingStrategyPop[7]], misc.cell.count,
                     flowType.res@CellFreqs[GatingStrategyPop[8:10]])
  percent.total.fT <- cell.count.fT/cell.count.fD[csv.ind, 1]*100
  percent.cd45.fT <- cell.count.fT/cell.count.fD[csv.ind, 2]*100
  
  # Parents are (Lin-, Lin-Mac-, Lin-Mac-, Lin-Mac-, Ly6G-, Ly6G-, not.mono.Ly6c.hi, not.mono.Ly6c.hi, not.pDC, not.pDC, cDC, cDC, CD8A.TypeDC
  parent.cell.count <- c(cell.count.fD[csv.ind, 3], cell.count.fT[c(1, 1, 1, 2, 2, 6, 6, 8, 8, 9, 9, 11)])
  percent.parent.fT <- cell.count.fT/parent.cell.count*100
  
  cell.count.fD.temp <- unlist(cell.count.fD[csv.ind, 4:16])
  names(cell.count.fD.temp) <- NULL
  percent.total.fD <- cell.count.fD.temp/cell.count.fD[csv.ind, 1]*100
  percent.cd45.fD <- cell.count.fD.temp/cell.count.fD[csv.ind ,2]*100
  percent.parent.fD.temp <- unlist(percent.parent.fD[csv.ind, 4:16])
  names(percent.parent.fD.temp) <- NULL
  
  current.barcode <- sub(".*(L000.*?)(_).*", "\\1", allFT[q2])
  d <- rbind(percent.total.fD, percent.total.fT, percent.cd45.fD, percent.cd45.fT, percent.parent.fD.temp, percent.parent.fT, cell.count.fD.temp, cell.count.fT)
  data.type <- rownames(d)
  rownames(d) <- NULL
  d <- cbind(rep(current.barcode, 8), data.type, d)
  FinalTable <- rbind(FinalTable, d)

  print(paste("Done", allFT[q2]))
    
}

colnames(FinalTable) <- c('barcode', 'data.type','Lin-Mac-', 'Ly6G-', 'Neutrophils 2', 'Eosinophils 3', 'Monocytes Ly6c hi',
                          'NOT(Monocytes Ly6c hi)', 'pDC', 'NOT(pDC)', 'cDC', 'Misc', 'CD8A Type DC', 'CD11b+ CD86lo',
                          'CD103+ DC')

write.table(FinalTable, "Results/FinalTableOfProportions.csv", sep=",", col.names=F, row.names=F)

