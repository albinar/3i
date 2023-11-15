
for (MLNorSPLN in c("SPLN", "MLN")){

    # load locally
    source("/data/IMPC_Sanger/IMPC_Functions.R")
    source("/data/IMPC_Sanger/timeVsFluorescence.R")
    source("/data/IMPC_Sanger/SangerCleanFunc.R")
    #source("/data/Panel_M-cell/Codes/SangerCleanFunc.R")

    if(MLNorSPLN == "MLN"  ){
        setwd("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_M-cell/MLN/")
        results.dir <- "/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_M-cell/MLN/Results/"
        #load("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_M-cell/MLN/Results/MLNdesc.Rdata") # storeDesc
    }
    if(MLNorSPLN == "SPLN" ){
        setwd("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_M-cell/SPLEEN/")
        results.dir <- "/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_M-cell/SPLEEN/Results/"
        #load("/mnt/f/FCS data/IMPC/IMPC-Results/3i/Panel_M-cell/SPLEEN/Results/SLPEENdesc.Rdata") # storeDesc
    }

    marker.names <- c('MHCII', 'F4/80', 'Ly6G', 'Ly6C', 'Live/Dead', '(CD3,CD19,CD161)', 'CD11b', 'CD45', 'CD317', 'CD11c', 'CD103', 'CD86')

    res.clean <- SangerCleanFunc(results.dir, marker.names) #, overwrite.channel.names = storeDesc)

}