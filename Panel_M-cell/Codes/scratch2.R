# This is for plotting the gthres outliers

temp <-   data.frame(x$FCS.files, all.props, all.events, gthres)
load(file = "Results/PropsEventsGates.Rdata")
props.events.gates[i, 10:ncol(props.events.gates)] <- temp[1, 2:ncol(temp)]
save(props.events.gates, file = "Results/PropsEventsGates.Rdata")

library(reshape2)
library(ggplot2)
library(dplyr)

gthres.temp <- data.frame(gthres[, 4:ncol(gthres)])
gthres.temp <- na.omit(gthres.temp)
gthres.temp <- melt(gthres.temp, id.vars = c('FCS.files','Barcodes'))
normf <- function(m){ # z-score
  (m - mean(m))/sd(m)
}
gthres.temp <- gthres.temp %>%
  group_by(variable) %>%
  mutate(value.norm = normf((value)))

boxplot = ggplot(gthres.temp, mapping=aes_string(y = "value.norm", x = "variable")) + theme_bw()
boxplot = boxplot + geom_boxplot() 
boxplot


library('dplyr')
test <- data.frame(Matrix3iTcell.index = 1:length(allFT), allFT, Genotypes)
test <- mutate(test, Label.Barcode = sub(".*(L000.*?)(_).*", "\\1", allFT))

load("Results_MLN/CSVfile.Rdata")

# unfortunately the Label.Barcode data in the file below can have multiple barcodes on a line
GenderData <- CSVfile[,c('Gender','Label.Barcode','Assay.Date')]

GenderCrop <- NULL
for(idx in 1:nrow(test)){
  GenderCrop[idx] <- min(grep(test$Label.Barcode[idx], GenderData[,2])) # I added the min because CSVfile contains seemingly duplicate entries (same barcodes, gender, assay date ect)
}
GenderData <- GenderData[GenderCrop,]

test <- cbind(test, GenderData[ ,c('Gender', 'Assay.Date')]) # I cannot merge because the GenderData still contains the multiple barcodes per row

WTfemale.indices <- which(tolower(test[which(test$Genotypes %in% c("+_+","+_Y")), "Gender"]) == "female")
WTmale.indices <- which(tolower(test[which(test$Genotypes %in% c("+_+","+_Y")), "Gender"]) == "male")                       

test$Assay.Date <- as.Date( as.character(test$Assay.Date), format = "%d-%b-%Y")

uniqueWT <- unique(test$Genotypes)
k1 <- 9
KOindices <- c(which(test$Genotypes == uniqueWT[k1]))
test$Assay.Date[KOindices]

min(test$Assay.Date[KOindices])
max(test$Assay.Date[KOindices])


sw.WTfemale.indices <- which(test$Assay.Date[WTfemale.indices] >= min(test$Assay.Date[KOindices]) & 
                               test$Assay.Date[WTfemale.indices] <= max(test$Assay.Date[KOindices]) )
sw.WTmale.indices <- which(test$Assay.Date[WTmale.indices] >= min(test$Assay.Date[KOindices]) & 
                               test$Assay.Date[WTmale.indices] <= max(test$Assay.Date[KOindices]) )
# now I need to check the length. If these have lengths of < 70, I would have to add more dates or files
# add closest dates first until go over 70 files


uniqueDates <- unique(test$Assay.Date)
uniqueDates <- sort(uniqueDates)

DateRange <- c(which(uniqueDates == min(test$Assay.Date[KOindices])), which(uniqueDates == max(test$Assay.Date[KOindices])))
while(length(sw.WTfemale.indices) < 70){
  start.idx <- max(1, DateRange[1] - 1)
  stop.idx <- min(length(uniqueDates), DateRange[2] + 1)
  sw.WTfemale.indices <- which(test$Assay.Date[WTfemale.indices] >= min(uniqueDates[start.idx]) & 
                                 test$Assay.Date[WTfemale.indices] <= max(uniqueDates[stop.idx]))
  DateRange <- c(start.idx, stop.idx)
}
if(length(sw.WTfemale.indices) > 70){ # cut down length to 70
  q <- (length(sw.WTfemale.indices) -70)%/%2
  r <- (length(sw.WTfemale.indices) -70)%%2
  sw.WTfemale.indices <- sw.WTfemale.indices[(1+q+r):(70+q+r)]
}



DateRange <- c(which(uniqueDates == min(test$Assay.Date[KOindices])), which(uniqueDates == max(test$Assay.Date[KOindices])))
while(length(sw.WTmale.indices) < 70){
  start.idx <- max(1, DateRange[1] - 1)
  stop.idx <- min(length(uniqueDates), DateRange[2] + 1)
  sw.WTfemale.indices <- which(test$Assay.Date[WTmale.indices] >= min(uniqueDates[start.idx]) & 
                                 test$Assay.Date[WTmale.indices] <= max(uniqueDates[stop.idx]))
  DateRange <- c(start.idx, stop.idx)
}
if(length(sw.WTmale.indices) > 70){ # cut down length to 70
  q <- (length(sw.WTmale.indices) -70)%/%2
  r <- (length(sw.WTmale.indices) -70)%%2
  sw.WTmale.indices <- sw.WTmale.indices[(1+q+r):(70+q+r)]
}

# WTfemale.indices[sw.WTfemale.indices]
# WTmale.indices[sw.WTmale.indices]


# Scratch for SangerCleanFunc.R
for(i in 1:nrow(store.allFCS)){
  if(row.names(results[[i]]$data)[1] == "error message"){
    print(results[[i]]$data)
  }
}


