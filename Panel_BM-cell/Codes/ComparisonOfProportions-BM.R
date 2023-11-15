# CSVfile <- read.csv("/home/arahim/Desktop/3iTcell/3i T-cell Data Spleen 04-2014.csv", check.names = FALSE)
# CSVfile <- as.matrix(CSVfile)
# 
# Chart <- CSVfile[,c(21:(26*5+1))]
# CellCount <- CSVfile[,c((26*5+2):(7*26+14))]
# 
# SolCellCount <- CSVfile[,c(20,(26*5+2), (26*5+6):(6*26+12))]
# # SolCellCount <- CSVfile[,c(20,(26*5+2):(6*26+12))]
# 
# PercentageOfCD45 <- CSVfile[,c(20,(26*3+1):(4*26+4))]
# 
# PercentageOfParent <- CSVfile[,c(20,(21):(2*26+1))]

#############################################################################3
CSVfile <- read.csv("/home/arahim/Desktop/3iTcell/data_output-99.csv", check.names = FALSE)
CSVfile <- as.matrix(CSVfile)

# Chart <- CSVfile[,c(21:(26*5+1))]
# CellCount <- CSVfile[,c((26*5+2):(7*26+14))]

SolCellCount <- CSVfile[,c(20,(26*5+2+10),(26*5+6+10):(26*5+22+10),(26*5+22+12):(6*26+12+10+1))]
# SolCellCount <- CSVfile[,c(20,(26*5+2):(6*26+12))]

PercentageOfCD45 <- CSVfile[,c(20,(26*3+1+5):(26*3+21),(26*3+23):(4*26+4+5+1), 42,51,21)]

PercentageOfParent <- CSVfile[,c(20,21:37, 39:(2*26+2))]
#############################################################################3
CSVfile <- read.csv("/home/arahim/Desktop/3iTcell_AFC_Oct_2014/FCS_Oct_2014/attachments/Spleen_Data_III.csv", check.names = FALSE)
CSVfile <- as.matrix(CSVfile)

# Chart <- CSVfile[,c(21:(26*5+1))]
# CellCount <- CSVfile[,c((26*5+2):(7*26+14))]

SolCellCount <- CSVfile[,c(20,(26*5+2+10+10),(26*5+6+10+10):(26*5+22+10+10),(26*5+22+12+10):(6*26+12+10+1+10))]
# SolCellCount <- CSVfile[,c(20,(26*5+2):(6*26+12))]

PercentageOfCD45 <- CSVfile[,c(20,(26*3+1+5):(26*3+21),(26*3+23):(4*26+4+5+1), 42,51,21)]

PercentageOfParent <- CSVfile[,c(20,21:37, 39:(2*26+2))]
#############################################################################3

pathFiles <- paste("/home/arahim/Desktop/3iTcell_Results/FlowType")
# Reads all folders and files in current path folder and makes a list of all of their paths
allFiles  <- dir(pathFiles,  full.names=T, recursive=T) 

allFiles2 <- allFiles
allFiles3 <- allFiles
Genotype  <- allFiles
for ( k1 in 1:length(allFiles)) {
  temp.name <- strsplit(allFiles2[k1],split="/")
  temp.name2 <- temp.name[[1]][length(unlist(temp.name[[1]]))]
  Genotype[k1] <- temp.name[[1]][length(unlist(temp.name[[1]]))-1]
  temp.name2 <- strsplit(temp.name2,split="[.]")
  allFiles2[k1] <- temp.name2[[1]][1]
  allFiles3[k1] <- strsplit(allFiles2[k1],split="_")[[1]][3]
  allFiles2[k1] <- paste(Genotype[k1], "/", allFiles2[k1], sep="")
  allFiles2[k1] <- sub ( "/SPLN", "/FT_SPLN", allFiles2[k1])
  allFiles2[k1] <- sub ( "/SLPN", "/FT_SLPN", allFiles2[k1])
  allFiles2[k1] <- sub ( "/Spleen", "/FT_Spleen", allFiles2[k1])
}

FinalTable <- NULL
TotalSaveIt <- NULL

for ( q2 in 1: length(allFiles)){
# for ( q2 in 43: 44){
    
    load(paste("/home/arahim/Desktop/3iTcell_Results/FlowType/", allFiles2[q2], ".Rdata", sep=""))
    # load("/home/arahim/Desktop/3iTcell_Results/FlowType/+_+/FT_SLPN_L000017420_Size_007.Rdata")
    # f <- read.FCS( "FCS_Groups/+_+/SLPN_L000017420_Size_007.fcs")
#     f <- read.FCS( paste( allFiles[q2], sep=""))
    
    
    Chart_JM <- flowType.res@CellFreqs
    
    # PropNames <- unlist(lapply(flowType.res@PhenoCodes, function(x){return(decodePhenotype(
    #                           x,as.vector(f@parameters@data$desc[7:16]), flowType.res@PartitionsPerMarker))}))
    
    # save(PropNames, file =  paste("../3iTcell_Results/PropNames.Rdata",sep=""))
    
#     load(file =  paste("../3iTcell_Results/PropNames.Rdata",sep=""))
    load(file =  paste("../3iTcell_Results/Phenotypes.Rdata",sep=""))
    
    # colnames(Chart_JM) <- PropNames 
    Chart_JM <- rbind(Phenotypes, Chart_JM)
    
    # Chart"CD25-TCRb+KLRG1+CD161-CD4+"
    
    Who <- NULL
    # Who[1]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD25+TCRb+CD161-CD4+GITR+") ]) / as.numeric(Chart_JM[2,1])
    # Who[2]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD25+TCRb+CD161-CD4+GITR+") ]) / as.numeric(Chart_JM[2,1])
    # Who[3]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD25+TCRb+CD161-CD4+GITR+") ]) / as.numeric(Chart_JM[2,1])
    
    Who[4]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-CD8a-CD161+TCRd-") ]) / as.numeric(Chart_JM[2,1]) 
    Who[5]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-TCRb-CD8a-CD161+TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[6]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+TCRb-CD8a-CD161+TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[7]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+TCRb-KLRG1+CD8a-CD161+TCRd-") ]) / as.numeric(Chart_JM[2,1])
    
    # Who[8]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-CD161-TCRd+") ]) / as.numeric(Chart_JM[2,1])
    # Who[9]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-TCRb-CD161-TCRd+") ]) / as.numeric(Chart_JM[2,1])
    # Who[10] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+TCRb-CD161-TCRd+") ]) / as.numeric(Chart_JM[2,1])
    # Who[11] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+TCRb-KLRG1+CD161-TCRd+") ]) / as.numeric(Chart_JM[2,1])
    
    Who[8]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-TCRd+") ]) / as.numeric(Chart_JM[2,1])
    Who[9]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-KLRG1+GITR-TCRd+") ]) / as.numeric(Chart_JM[2,1])
    Who[10] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-TCRb-TCRd+") ]) / as.numeric(Chart_JM[2,1])
    Who[11] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+TCRb-TCRd+") ]) / as.numeric(Chart_JM[2,1])
        
    Who[12] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+CD161+CD4-TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[14] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-TCRb+CD161+CD4-TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[15] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+TCRb+CD161+CD4-TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[16] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+KLRG1+CD161+CD4-TCRd-") ]) / as.numeric(Chart_JM[2,1])
    
    Who[13] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+CD161+CD4+TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[17] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-TCRb+CD161+CD4+TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[18] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+TCRb+CD161+CD4+TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[19] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+KLRG1+CD161+CD4+TCRd-") ]) / as.numeric(Chart_JM[2,1])
    
    Who[20] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD25+TCRb+CD8a-CD161-CD4+GITR+TCRd-") ]) / as.numeric(Chart_JM[2,1]) 
    Who[21] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-CD25+TCRb+CD8a-CD161-CD4+GITR+TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[22] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+CD25+TCRb+CD8a-CD161-CD4+GITR+TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[23] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD25+TCRb+KLRG1+CD8a-CD161-CD4+GITR+TCRd-") ]) / as.numeric(Chart_JM[2,1])
    
    Who[24] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+CD8a-CD161-CD4+TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[26] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD25-TCRb+CD8a-CD161-CD4+TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[27] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-CD25-TCRb+CD8a-CD161-CD4+TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[28] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD62L+CD25-TCRb+CD8a-CD161-CD4+TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[29] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD25-TCRb+KLRG1+CD8a-CD161-CD4+TCRd-") ]) / as.numeric(Chart_JM[2,1])
    
    Who[25] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+CD8a+CD161-CD4-TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[30] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-TCRb+CD8a+CD161-CD4-TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[31] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+TCRb+CD8a+CD161-CD4-TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[32] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+TCRb+KLRG1+CD8a+CD161-CD4-TCRd-") ]) / as.numeric(Chart_JM[2,1])
    
    Who[33] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[34] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+CD161+TCRd-") ]) / as.numeric(Chart_JM[2,1])
#     Who[34] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+CD161-TCRd-") ]) / as.numeric(Chart_JM[2,1]) # negative to be consistent with manual data

    Who[35] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-KLRG1-GITR+TCRd+") ]) / as.numeric(Chart_JM[2,1])
    Who[36] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44-CD62L+TCRb-TCRd+") ]) / as.numeric(Chart_JM[2,1])
    Who[37] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-TCRd-") ]) / as.numeric(Chart_JM[2,1])
    Who[38] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44-CD62L+TCRb+CD8a+CD161-CD4-TCRd-") ]) / as.numeric(Chart_JM[2,1])

    
    
    WhoCount <- NULL

    WhoCount[3] <- as.numeric(Chart_JM[2,1])
    
    WhoCount[4]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-CD8a-CD161+TCRd-") ]) 
    WhoCount[5]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-TCRb-CD8a-CD161+TCRd-") ]) 
    WhoCount[6]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+TCRb-CD8a-CD161+TCRd-") ]) 
    WhoCount[7]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+TCRb-KLRG1+CD8a-CD161+TCRd-") ]) 

    WhoCount[8]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-TCRd+") ]) 
    WhoCount[9]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-KLRG1+GITR-TCRd+") ]) 
    WhoCount[10] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-TCRb-TCRd+") ]) 
    WhoCount[11] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+TCRb-TCRd+") ]) 
    

    
    WhoCount[12] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+CD161+CD4-TCRd-") ]) 
    WhoCount[14] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-TCRb+CD161+CD4-TCRd-") ]) 
    WhoCount[15] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+TCRb+CD161+CD4-TCRd-") ]) 
    WhoCount[16] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+KLRG1+CD161+CD4-TCRd-") ]) 
    
    WhoCount[13] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+CD161+CD4+TCRd-") ]) 
    WhoCount[17] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-TCRb+CD161+CD4+TCRd-") ]) 
    WhoCount[18] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+TCRb+CD161+CD4+TCRd-") ]) 
    WhoCount[19] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+KLRG1+CD161+CD4+TCRd-") ]) 
    
    WhoCount[24] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+CD8a-CD161-CD4+TCRd-") ]) 
    WhoCount[20] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD25+TCRb+CD8a-CD161-CD4+GITR+TCRd-") ]) 
    WhoCount[21] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-CD25+TCRb+CD8a-CD161-CD4+GITR+TCRd-") ]) 
    WhoCount[22] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+CD25+TCRb+CD8a-CD161-CD4+GITR+TCRd-") ]) 
    WhoCount[23] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD25+TCRb+KLRG1+CD8a-CD161-CD4+GITR+TCRd-") ]) 
    
    WhoCount[26] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD25-TCRb+CD8a-CD161-CD4+TCRd-") ]) 
    WhoCount[27] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-CD25-TCRb+CD8a-CD161-CD4+TCRd-") ])
    WhoCount[28] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD62L+CD25-TCRb+CD8a-CD161-CD4+TCRd-") ])
    WhoCount[29] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD25-TCRb+KLRG1+CD8a-CD161-CD4+TCRd-") ])
    
    WhoCount[25] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+CD8a+CD161-CD4-TCRd-") ])
    WhoCount[30] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-TCRb+CD8a+CD161-CD4-TCRd-") ])
    WhoCount[31] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+TCRb+CD8a+CD161-CD4-TCRd-") ])
    WhoCount[32] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+TCRb+KLRG1+CD8a+CD161-CD4-TCRd-") ]) 
    
    WhoCount[33] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+TCRd-") ]) 
    WhoCount[34] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+CD161+TCRd-") ])
    
    WhoCount[35] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-KLRG1-GITR+TCRd+") ])
    WhoCount[36] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44-CD62L+TCRb-TCRd+") ])
    WhoCount[37] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-TCRd-") ])
    WhoCount[38] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44-CD62L+TCRb+CD8a+CD161-CD4-TCRd-") ])

    
    
    WhoPercentParent <- NULL

#     WhoPercentParent[4]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-CD8a-CD161+TCRd-") ]) / as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-TCRd-") ])
    WhoPercentParent[4]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-CD8a-CD161+TCRd-") ]) / as.numeric(Chart_JM[2,1])
    
    WhoPercentParent[5]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-TCRb-CD8a-CD161+TCRd-") ]) / WhoCount[4]
    WhoPercentParent[6]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+TCRb-CD8a-CD161+TCRd-") ]) / WhoCount[4]
    WhoPercentParent[7]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+TCRb-KLRG1+CD8a-CD161+TCRd-") ]) / WhoCount[4]

#     WhoPercentParent[8]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-TCRd+") ]) / as.numeric(Chart_JM[2,1])
    WhoPercentParent[8]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-TCRd+") ]) / (as.numeric(Chart_JM[2,1]) - as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-TCRd-") ]) )
    WhoPercentParent[9]  <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-KLRG1+GITR-TCRd+") ]) / WhoCount[8]
    WhoPercentParent[10] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-TCRb-TCRd+") ]) / WhoCount[8]
    WhoPercentParent[11] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+TCRb-TCRd+") ]) / WhoCount[8]
    
#     WhoPercentParent[33] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+TCRd-") ]) / as.numeric(Chart_JM[2,1])
    WhoPercentParent[33] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+TCRd-") ]) / (as.numeric(Chart_JM[2,1]) - as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-TCRd-") ]) )

    WhoCountB    <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+CD161-TCRd-") ])
#     WhoPercentParent[34] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+CD161+TCRd-") ]) / WhoCount[33]
    WhoPercentParent[34] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+CD161+TCRd-") ]) / as.numeric(Chart_JM[2,1])

    
    WhoPercentParent[12] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+CD161+CD4-TCRd-") ]) / WhoCount[34]
    WhoPercentParent[14] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-TCRb+CD161+CD4-TCRd-") ]) / WhoCount[12]
    WhoPercentParent[15] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+TCRb+CD161+CD4-TCRd-") ]) / WhoCount[12]
    WhoPercentParent[16] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+KLRG1+CD161+CD4-TCRd-") ]) / WhoCount[12]
    
    WhoPercentParent[13] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+CD161+CD4+TCRd-") ]) / WhoCount[34]
    WhoPercentParent[17] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-TCRb+CD161+CD4+TCRd-") ]) / WhoCount[13]
    WhoPercentParent[18] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+TCRb+CD161+CD4+TCRd-") ]) / WhoCount[13]
    WhoPercentParent[19] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+KLRG1+CD161+CD4+TCRd-") ]) / WhoCount[13]
    
    WhoPercentParent[24] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+CD8a-CD161-CD4+TCRd-") ]) / WhoCountB
    
    WhoPercentParent[20] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD25+TCRb+CD8a-CD161-CD4+GITR+TCRd-") ]) / WhoCount[24]
    WhoPercentParent[21] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-CD25+TCRb+CD8a-CD161-CD4+GITR+TCRd-") ]) / WhoCount[20]
    WhoPercentParent[22] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+CD25+TCRb+CD8a-CD161-CD4+GITR+TCRd-") ]) / WhoCount[20]
    WhoPercentParent[23] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD25+TCRb+KLRG1+CD8a-CD161-CD4+GITR+TCRd-") ]) / WhoCount[20]
    
    WhoPercentParent[26] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD25-TCRb+CD8a-CD161-CD4+TCRd-") ]) / WhoCount[24]
    WhoPercentParent[27] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-CD25-TCRb+CD8a-CD161-CD4+TCRd-") ]) / WhoCount[26]
    WhoPercentParent[28] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD62L+CD25-TCRb+CD8a-CD161-CD4+TCRd-") ]) / WhoCount[26]
    WhoPercentParent[29] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD25-TCRb+KLRG1+CD8a-CD161-CD4+TCRd-") ]) / WhoCount[26]
    
    WhoPercentParent[25] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb+CD8a+CD161-CD4-TCRd-") ]) / WhoCountB
    WhoPercentParent[30] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-TCRb+CD8a+CD161-CD4-TCRd-") ]) / WhoCount[25]
    WhoPercentParent[31] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L+TCRb+CD8a+CD161-CD4-TCRd-") ]) / WhoCount[25]
    WhoPercentParent[32] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+TCRb+KLRG1+CD8a+CD161-CD4-TCRd-") ]) / WhoCount[25]
    
    WhoPercentParent[35] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-KLRG1-GITR+TCRd+") ]) / WhoCount[8]
    WhoPercentParent[36] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44-CD62L+TCRb-TCRd+") ]) / WhoCount[8]
    WhoPercentParent[37] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "TCRb-TCRd-") ]) / WhoCount[8]
    WhoPercentParent[38] <- as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44-CD62L+TCRb+CD8a+CD161-CD4-TCRd-") ]) / WhoCount[25]

    
    # 
    # 
    # Chart_JM[2,6072]
    # 
    # Chart_JM[2,1]
    # 
    # 461/347704
    # 
    # What <- NULL
    # What[20] <- prod(as.numeric(Chart[5,c(1,2,3,8)])) / 100^4
    # What[21] <- prod(as.numeric(Chart[5,c(1,2,3,8,9)])) / 100^5
    # What[22] <- prod(as.numeric(Chart[5,c(1,2,3,8,10)])) / 100^5
    # What[23] <- prod(as.numeric(Chart[5,c(1,2,3,8,11)])) / 100^5
    # What[24] <- prod(as.numeric(Chart[5,c(1,2,3)])) / 100^3
    # What[26] <- prod(as.numeric(Chart[5,c(1,2,3,4)])) / 100^4
    # What[27] <- prod(as.numeric(Chart[5,c(1,2,3,4,5)])) / 100^5
    # What[28] <- prod(as.numeric(Chart[5,c(1,2,3,4,6)])) / 100^5
    # What[29] <- prod(as.numeric(Chart[5,c(1,2,3,4,7)])) / 100^5
    # 
    # What[25] <- prod(as.numeric(Chart[5,c(1,2,12)])) / 100^3
    # What[30] <- prod(as.numeric(Chart[5,c(1,2,12,13)])) / 100^4
    # What[31] <- prod(as.numeric(Chart[5,c(1,2,12,14)])) / 100^4
    # What[32] <- prod(as.numeric(Chart[5,c(1,2,12,16)])) / 100^4
    # 
    #  as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD44+CD62L-CD25-TCRb+CD161-CD4+") ]) 
    #  as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD62L+CD25-TCRb+CD161-CD4+") ]) 
    #  as.numeric(Chart_JM[2, which(Chart_JM[1,] == "CD25-TCRb+KLRG1+CD161-CD4+") ]) 
    # 
    # 11026
    # 35522
    # 499
    
    # dude <- rbind(What, Who)
    # saveIt <- dude[1,c(20:24,26:29)]/dude[2,c(20:24,26:29)]
    # saveIt <- dude[1,c(20:32)]/dude[2,c(20:32)]
    # print(saveIt)
    # 1/mean(saveIt)
    
    
    
    # which(PercentageOfCD45[,1]=="L000017420")
#     temp <- PercentageOfCD45[which(PercentageOfCD45[,1]=="L000017420"), ]
    temp <- PercentageOfCD45[which(PercentageOfCD45[,1]== allFiles3[q2] ), ]
#     temp <- temp[c(29,29,30,31, 17,20,18,19,21,25,22,23,24,26,27,28, 8,9,10,11,3,12,4,5,6,7,13,14,16,2)]
#     temp <- temp[c(17,17,20,18,19, 1,29,30,31, 2,  21,22,23,24,25,26,27,28,29, 3,8,10,9,11, 4,5,6,7, 12,14,13,16)]
    temp <- temp[c(34,17,20,18,19, 33,29,30,31, 2,32,21,22,23,24,25,26,27,28,    3,8,10,9,11, 4,5,6,7, 12,14,13,15,16)] # 15 Naive CD8
    namesTemp <- names(temp)
    temp <- as.numeric(temp)
  
    
    Who2             <- Who             [c(1, 8,9,10,11, 4,5,6,7, 33,34,12,14,15,16,13,17,18,19,24,20,21,22,23,26,27,28,29,25,31,30,38,32)]*100
    WhoCount         <- WhoCount        [c(3, 8,9,10,11, 4,5,6,7, 33,34,12,14,15,16,13,17,18,19,24,20,21,22,23,26,27,28,29,25,31,30,38,32)]
    WhoPercentParent <- WhoPercentParent[c(1, 8,9,10,11, 4,5,6,7, 33,34,12,14,15,16,13,17,18,19,24,20,21,22,23,26,27,28,29,25,31,30,38,32)]*100

    temp2 <- PercentageOfParent[which(PercentageOfParent[,1]== allFiles3[q2] ), ]
    temp2 <- temp2[c(2, 18,21,19,20, 31,32,33,34, 3,22,23,24,25,26,27,28,29,30, 4,9,11,10,12, 5,6,7,8, 13,15,14,16,17)] # 16 Naive CD8
    namesTemp2 <- names(temp2)
    temp2 <- as.numeric(temp2)
    
    temp3 <- SolCellCount[which(SolCellCount[,1]== allFiles3[q2] ), ]
    temp3 <- temp3[c(2, 19,22,20,21, 32,33,34,35, 4,23,24,25,26,27,28,29,30,31, 5,10,12,11,13, 6,7,8,9, 14,16,15,17,18)] # 17 Naive CD8
    namesTemp3 <- names(temp3)
    temp3 <- as.numeric(temp3)
    
#     dude <- rbind( as.numeric(temp) , Who2)
    
    # saveIt <- dude[1,c(20:24,26:29)]/dude[2,c(20:24,26:29)]
    # saveIt <- dude[1,c(20:32)]/dude[2,c(20:32)]
    # print(saveIt)
    # 1/mean(saveIt)
#     dude2 <-dude[,2:32]
#     colnames(dude2) <- names(temp)[2:32]
#     saveIt <- dude2[1,]/dude2[2,]
#     print(saveIt)
#     saveIt <- c(0, saveIt)
#     dude2A <- c(0, dude2[1,])
#     dude2B <- c(0, dude2[2,])
    
#     temp[1] <- 0
    Who2[1] <- 0


#     1/mean(saveIt)
#     FinalTable <- rbind( FinalTable, saveIt )
#     FinalTable <- rbind( FinalTable, dude2A, dude2B, temp2, WhoPercentParent*100, temp3, WhoCount )
    FinalTable <- rbind( FinalTable, temp, Who2, temp2, WhoPercentParent, temp3, WhoCount )

#     TotalSaveIt <- rbind(TotalSaveIt, saveIt)
#     rownames(FinalTable)[q2] <- allFiles3[q2]
    print(paste("Done", allFiles3[q2]))
}

#     FinalTable <- cbind(allFiles3[1:4], FinalTable)
# FinalTable <- cbind(allFiles3[1:length(allFiles)], FinalTable)

allMiceF <- allFiles3[1:length(allFiles)]
# allMiceF[seq(1,(2*length(allFiles)-1),by=2)] <- allFiles3[1:length(allFiles)]
# allMiceF[seq(2, (2*length(allFiles)), by=2)] <- allFiles3[1:length(allFiles)]

# allMiceF[seq(1, (7*length(allFiles)-6), by=7)] <- allFiles3[1:length(allFiles)]
# allMiceF[seq(2, (7*length(allFiles)-5), by=7)] <- allFiles3[1:length(allFiles)]
# allMiceF[seq(3, (7*length(allFiles)-4), by=7)] <- allFiles3[1:length(allFiles)]
# allMiceF[seq(4, (7*length(allFiles)-3), by=7)] <- allFiles3[1:length(allFiles)]
# allMiceF[seq(5, (7*length(allFiles)-2), by=7)] <- allFiles3[1:length(allFiles)]
# allMiceF[seq(6, (7*length(allFiles)-1), by=7)] <- allFiles3[1:length(allFiles)]
# allMiceF[seq(7, (7*length(allFiles))  , by=7)] <- allFiles3[1:length(allFiles)]

allMiceF[seq(1, (6*length(allFiles)-5), by=6)] <- allFiles3[1:length(allFiles)]
allMiceF[seq(2, (6*length(allFiles)-4), by=6)] <- allFiles3[1:length(allFiles)]
allMiceF[seq(3, (6*length(allFiles)-3), by=6)] <- allFiles3[1:length(allFiles)]
allMiceF[seq(4, (6*length(allFiles)-2), by=6)] <- allFiles3[1:length(allFiles)]
allMiceF[seq(5, (6*length(allFiles)-1), by=6)] <- allFiles3[1:length(allFiles)]
allMiceF[seq(6, (6*length(allFiles))  , by=6)] <- allFiles3[1:length(allFiles)]


FinalTableSave <- FinalTable
FinalTable <- FinalTableSave

FinalTable <- cbind(allMiceF, FinalTable)


#     colnames(FinalTable)[1] <- "Mice"
    FinalTable <- rbind(c("", namesTemp),
                        c("", namesTemp),
                        c("", namesTemp2),
                        c("", namesTemp2),
                        c("", namesTemp3), 
                        c("", namesTemp3), FinalTable)
FinalTable[1:6,1] <- c("UK Res", "JM Res", "UK Res", "JM Res", "UK Res", "JM Res")
#     colnames(FinalTable)[6:35] <- colnames(FinalTable)[5:34]
write.table(FinalTable, "../3iTcell_Results/FinalTableOfProportions.csv", sep=",", col.names=F, row.names=F)