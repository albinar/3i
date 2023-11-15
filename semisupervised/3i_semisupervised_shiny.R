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


##########################
##########################

##########################
##########################
##########################
source("3i-Helper.R")
source("Ania_helper.R")
#Input.path needs to be changed based on user's directory
input.path <- "/home/rstudio/results/3i/"
##This is for BM:
panel <- c("Panel_B-cell", "Panel_BM-cell","Panel_M-cell", "Panel_P2-cell","Panel_T-cell" )[3]
organ <- c("","MLN","Spleen")[3]
temp <- load(paste0(input.path,panel,"/", "semisupervised/",organ,"/Gates.list.Rdata"))
gates.list <- get(temp)
load(list.files(paste0(input.path,panel,"/semisupervised/",organ,"/flowType/"),pattern = ".Rdata",full.names = T)[1])
#load(list.files(paste0(input.path,panel,"/semisupervised/",organ,"/flowType/"),pattern = ".Rdata",full.names = T)[1034])


####################################################################

## Written by Mehrnoush so I can test
set.seed(100)
library(flowTypeFilter)
library(RchyOptimyx)
# load("~/data/Specimen_001_L000073886_C02_006.labelled.fcs_FT.Rdata")   
# load("~/data/All_genes_EffectSize.RData")

path <- paste(input.path,panel,"/semisupervised/",organ,"/RchyOptimyx/",sep="")
load(paste(input.path,panel,"/semisupervised/", organ,"/flowType/Specimen_001_L000073886_C02_006.labelled.fcs_FT.Rdata", sep = ""))   
load(paste(input.path,panel,"/semisupervised/", organ,"/All_genes_EffectSize.RData",sep=""))

ind <- which(all.effects[,"Genotype"]=="Vps13a(b):Hom")
ind2 <- which(as.numeric(all.effects[ind,6])>1.6)
load(paste(input.path,panel,"/semisupervised/", organ,"/All_genes_Pvalues.RData", sep = ""))
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


library(RchyOptimyx)
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
panel<-"Panel_M-cell"
gene <- reduced.genes <- "Vps13a(b):Hom"


# library("RSVGTipsDevice")
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
# merged.res@edges[2,1] <- "0.05"
# merged.res@edges[2,3] <- "0.05"
# merged.res@edges[2,4] <- "0.05"
# merged.res@edges[2,5] <- "0.05"
# merged.res@edges[2,6] <- "0.05"
# merged.res@edges[2,7] <- "0.05"
# merged.res@edges[2,8] <- "0.05"
# merged.res@edges[2,9] <- "0.05"
# 
# merged.res@edges[2,10] <- "0.05"
# merged.res@edges[2,12] <- "0.05"
# merged.res@edges[2,16] <- "0.05"
# merged.res@edges[2,17] <- "0.05"
# merged.res@edges[2,22] <- "0.05"
# merged.res@edges[2,23] <- "0.05"
# merged.res@edges[2,24] <- "0.05"
# merged.res@edges[2,25] <- "0.05"
# 
# merged.res@edges[2,27] <- "0.05"
# merged.res@edges[2,28] <- "0.05"
# merged.res@edges[2,31] <- "0.05"
# merged.res@edges[2,32] <- "0.05"
# merged.res@edges[2,34] <- "0.05"
# merged.res@edges[2,35] <- "0.05"
# merged.res@edges[2,39] <- "0.05"
# merged.res@edges[2,43] <- "0.05"
# merged.res@edges[2,44] <- "0.05"
# 
# merged.res@edges[2,14] <- "0.05"
# merged.res@edges[2,18] <- "0.05"
# merged.res@edges[2,20] <- "0.05"
# merged.res@edges[2,29] <- "0.05"
# merged.res@edges[2,36] <- "0.05"
# merged.res@edges[2,40] <- "0.05"
# merged.res@edges[2,41] <- "0.05"
# merged.res@edges[2,42] <- "0.05"
# 
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
# 
# 
# set.seed(101)
# devSVGTips(paste(path,gene,"_RchyOptimyx_Tooltip.svg",sep=""), toolTipMode=1, width=12,height=10)
# plot(merged.res,phenotypeCodes=ft@PhenoCodes,phenotypeScores=-log10(pvals),root.name=paste(gene,ifelse(panel=="Panel_M-cell",yes = "Lin-",no = ": CD45+")),legend.size=1.2,
#      partitions.per.marker=2,filters=filters,
#      node.lwd = 1,marker.names=marker.names)
# setSVGShapeToolTip(title=NULL, desc="Vps13a(b):Hom")
# #points(x=-0.5, y=1.01, cex=3, pch=19, col=NULL)
# text(x=0, y=1.01, lab="Test", col="black")
# dev.off()
# 
# devSVGTips(filename = paste(path,gene,"_RchyOptimyxTest.svg",sep=""), toolTipMode=1)
# plot(0, type="n")
# setSVGShapeToolTip("A rectangle","it is red")
# rect(.1,.1,.4,.6, col='red')
# dev.off()
# 
# #cat('<embed src="plot.svg" type="image/svg+xml" /'>)
# 
# library("RSVGTipsDevice")
# devSVGTips("svgplot1.svg", toolTipMode=1,
#            title="SVG example plot 1: shapes and points, tooltips are title + 1 line")
# plot(c(0,10),c(0,10), type="n", xlab="x", ylab="y",
#      main="Example SVG plot with title+ 1 line tips (mode=1)")
# setSVGShapeToolTip(title="A rectangle", desc="that is yellow")
# rect(1,1,4,6, col='yellow')
# setSVGShapeToolTip(title="1st circle with title only")
# points(5.5,7.5,cex=20,pch=19,col='red')
# setSVGShapeToolTip(title="A triangle", desc="big and green")
# polygon(c(3,6,8), c(3,6,3), col='green')
# # no tooltips on these points
# points(2:8, 8:2, cex=3, pch=19, col='black')
# # tooltips on each these points
# invisible(sapply(1:7, function(x) {
#   setSVGShapeToolTip(title=paste("point", x))
#   points(x+1, 8-x, cex=3, pch=1, col='black')
# }))
# setSVGShapeToolTip(title="Text", desc="can have a tool tip too!")
# text(x=4, y=9, lab="Poke me!", col="blue")
# dev.off()
# 
# 


##########################################################################


ui <- basicPage(
  plotOutput("plot1",
             #click = "plot_click",
             #dblclick = "plot_dblclick",
             hover = "plot_hover"
             #brush = "plot_brush"
  ),
  verbatimTextOutput("info")
)

server <- function(input, output) {
  output$plot1 <- renderPlot({
    #plot(mtcars$wt, mtcars$mpg)
    plot(merged.res,phenotypeCodes=ft@PhenoCodes,phenotypeScores=-log10(pvals),root.name=paste(gene,ifelse(panel=="Panel_M-cell",yes = "Lin-",no = ": CD45+")),legend.size=1.2,
         partitions.per.marker=2,filters=filters,
         node.lwd = 1,marker.names=marker.names)
  })
  
  output$info <- renderText({
    paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)
  })
  
  # output$info <- renderText({
  #   
  #   xy_str <- function(e) {
  #     if(is.null(e)) return("NULL\n")
  #     paste0("x=", round(e$x, 1), " y=", round(e$y, 1), "\n")
  #   }
  #   # xy_range_str <- function(e) {
  #   #   if(is.null(e)) return("NULL\n")
  #   #   paste0("xmin=", round(e$xmin, 1), " xmax=", round(e$xmax, 1), 
  #   #          " ymin=", round(e$ymin, 1), " ymax=", round(e$ymax, 1))
  #   # }
  #   
  #   paste0(
  #     #"click: ", xy_str(input$plot_click),
  #     #"dblclick: ", xy_str(input$plot_dblclick),
  #     "hover: ", xy_str(input$plot_hover)
  #     #"brush: ", xy_range_str(input$plot_brush)
  #   )
  # })
}

shinyApp(ui, server)
