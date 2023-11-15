setClass("OptimizedHierarchy",
         representation(nodes="matrix", edges="matrix", pathScores="vector"))


setMethod("summary", signature(object="OptimizedHierarchy"), function(object){
  cat(sprintf("An optimized hierarchy with %d nodes, %d edges, and %d paths\n", dim(object@nodes)[2], dim(object@edges)[2], length(object@pathScores)))
})

setMethod("plot", signature(x="OptimizedHierarchy", y="ANY"),
          function(x, 
                   phenotypeScores,
                   phenotypeCodes,
                   marker.names,
                   filters=rep(F,length(marker.names)),
                   partitions.per.marker = NULL,
                   uniformColors=FALSE, 
                   ylab=NULL, 
                   xlab=NULL, 
                   colors=c('blue','cyan','yellow','red'), 
                   edgeWeights=TRUE, 
                   edgeLabels=TRUE, 
                   nodeLabels=TRUE, 
                   min.score=NA, 
                   max.score=NA, 
                   cell.proportions=NULL, 
                   cell.counts=NULL,
                   min.proportion=NA, 
                   max.proportion=NA, 
                   proportion.colors=c("black", "white"),
                   node.lwd=5,
                   hide.small.prop=NA,
                   hide.small.count=NA,
                   root.name='All Cells',
                   legend.size=1.25,
                   plot.legend=TRUE,
                   textcol=par('fg')){            
            if (!is.vector(phenotypeScores))
              stop("phenptypeScores must be a numeric vector.")
            if (!is.vector(marker.names))
              stop("a character string must be presented as marker.names.")
            if (!is.logical(uniformColors))
              stop("uniformColors must be a logical value.")
            if (!(is.null(ylab) | is.character(ylab)))
              stop("ylab must be a character string or NULL.")
            if (!(is.null(xlab) | is.character(xlab)))
              stop("xlab must be a character string or NULL.")
            if (!is.vector(colors))
              stop("colors must be a color vector.")
            if (!is.vector(proportion.colors))
              stop("proportion.colors must be a color vector.")
            if (!is.logical(edgeWeights))
              stop("edgeWeights must be a logical value.")
            if (!is.logical(edgeLabels))
              stop("edgeLabels must be a logical value.")
            if (!is.logical(nodeLabels))
              stop("nodeLabels must be a logical value.")
            if (!(is.na(min.score) | is.numeric(min.score)))
              stop("min.score must be a numeric value or NA.")
            if (!(is.na(max.score) | is.numeric(max.score)))
              stop("max.score must be a numeric value or NA.")
            if (!(is.vector(cell.proportions) | is.null(cell.proportions)))
              stop("cell.proportions must be a numeric vector.")
            if (!(is.na(min.proportion) | is.numeric(min.proportion)))
              stop("min.proportion must be a numeric value or NA.")
            if (!(is.na(max.proportion) | is.numeric(max.proportion)))
              stop("max.proportion must be a numeric value or NA.")
            if (!is.numeric(node.lwd))
              stop("node.lwd must be a numeric value.")
            if (!is.character(root.name))
              stop("root.name must be a character string.")
            
            par(fg=textcol)
            SetTextContrastColor <- function(color){
              ifelse( color=="#FFFFFF", yes="white", no=ifelse(mean(col2rgb(color)) > 127,  yes="black", no="white"))
            }
            
            if (is.null(partitions.per.marker))
              partitions.per.marker <- rep(2, length(marker.names))
            if (is.null(filters))
              {
               warning("If filters are used in flowType make sure to pass it to plot(.), otherwise the population names will be incorrect!")
               rep(FALSE,length(marker.names))
              }
            OptH=x
            ab=OptH
            Scores=phenotypeScores
            names(Scores) <- phenotypeCodes
            ed <- vector("list", length=dim(ab@nodes)[2])
            V=vector();
            for (i in 1:dim(ab@nodes)[2]){  
            if (any(filters)){
              V[i]= flowTypeFilter::decodePhenotype(pheno.code=ab@nodes[1,i], marker.names=marker.names, partitions.per.marker=partitions.per.marker,filters=filters)
            }else{ 
                  V[i]= flowType::decodePhenotype(pheno.code=ab@nodes[1,i], marker.names=marker.names, partitions.per.marker=partitions.per.marker)
                 }
            }

            
            V[1]=root.name
            names(ed) <- V
            for (i in 1:dim(ab@nodes)[2]){  
              ed[[i]] <- list(edges=c(), weights=c(), labels=c())
            }
            for (i in 1:dim(ab@edges)[2]){
              temp=getNodeNumber(ab@edges[1,i], ab@nodes[1,])
              ed[[temp[1]]]$edges=c(ed[[temp[1]]]$edges,temp[2])
              if(edgeWeights){
                ab@edges[2,i] <- max(0, as.numeric(ab@edges[2,i]))
                ed[[temp[1]]]$weights=c(ed[[temp[1]]]$weights,as.numeric(ab@edges[2,i]))
                ##     Ves <- unlist(strsplit(ab@edges[1,i],'~'))
                ##      Ves[1]=ab@nodes[2,which(ab@nodes[1,]==Ves[1])]
                ##      Ves[2]=ab@nodes[2,which(ab@nodes[1,]==Ves[2])]
                ed[[temp[1]]]$labels=c(ed[[temp[1]]]$labels,ab@edges[3,i])      
              }
            }

            for (i in 1:length(ed)){
              if (length(ed[[i]]$edges)<2)
                next;
              yechiz<-order(ed[[i]]$edges,decreasing=FALSE)
              ed[[i]]$edges=ed[[i]]$edges[yechiz]
              ed[[i]]$weights=ed[[i]]$weights[yechiz]
              ed[[i]]$labels=ed[[i]]$labels[yechiz]
            }
            
            
            g <- new("graphNEL", nodes=V, edgeL=ed, edgemode='directed')
            attrs <- list(node = list(shape = "ellipse", fixedsize = FALSE))
            
            colorFunc <- colorRampPalette(colors)
            colorFunc2 <- colorRampPalette(proportion.colors)
            scores=unlist(lapply(1:length(ab@nodes[1,]), 
              function(i) {return(Scores[ab@nodes[1,i]])}))
            
            
            if (is.na(min.score))
              min.score<-min(scores)
            if (is.na(max.score))
              max.score<-max(scores)
            
            proportions=vector();
            if(!is.null(cell.proportions)){
              proportions=unlist(lapply(1:length(ab@nodes[1,]), 
                function(i) {
                          if (any(filters)) {
                           code <-flowTypeFilter::decodePhenotype(pheno.code =ab@nodes[1,i],partitions.per.marker = 2,marker.names = markers,filters=filters)
                           return(ifelse(code=="", yes=cell.proportions[1],no=cell.proportions[code]))

                         }else{
                                code <-flowType::decodePhenotype(pheno.code =ab@nodes[1,i],partitions.per.marker = 2,marker.names = markers)
                           return(ifelse(code=="", yes=cell.proportions[1],no=cell.proportions[code])) 
                         }
                        }))
            }
            
      counts=vector();
            if(!is.null(cell.counts)){
              counts=unlist(lapply(1:length(ab@nodes[1,]), 
                function(i) {
                          if (any(filters)) {
                           code <-flowTypeFilter::decodePhenotype(pheno.code =ab@nodes[1,i],partitions.per.marker = 2,marker.names = markers,filters=filters)
                             return(ifelse(code=="", yes=cell.counts[1],no=cell.counts[code]))
                         }else{
                               code <-flowType::decodePhenotype(pheno.code =ab@nodes[1,i],partitions.per.marker = 2,marker.names = markers)
                              return(ifelse(code=="", yes=cell.counts[1],no=cell.counts[code]))
                         }
                        }))
            }




            if(!is.null(cell.proportions)){
              if (is.na(min.proportion))
                min.proportion<-min(proportions)
              if (is.na(max.proportion))
                max.proportion<-max(proportions)
            }
            
            if(!uniformColors){
              z=pretty(c(min.score, max.score))
              colinds=unlist(lapply(1:length(scores), function(i){max(1,ceiling((scores-min(z))*1000/(max(z)-min(z)))[i])}))
              cols=colorFunc(1000)[colinds]
          if (!is.na(hide.small.prop) & !is.null(cell.proportions))
          { 
#Above vector is of same length as score, you can add arguments to make indices of elements with cell.prop less than some number to be white
               white.inds<- which(proportions < hide.small.prop)
               if (length(white.inds)>0)
                  cols[white.inds]<-"#FFFFFF"
          }
          if (!is.na(hide.small.count) & !is.null(cell.counts))
          { 
            print("Thresholds for cell count will be used. if you prefer cell proportion, set hide.small.count to NA and hide.small.prop to desired number.")
               #Above vector is of same length as score, you can add arguments to make indices of elements with cell.prop less than some number to be white
               white.inds<- which(counts < hide.small.count)
               if (length(white.inds)>0)
                  cols[white.inds]<-"#FFFFFF"
          }
          if ((!is.na(hide.small.count) & is.null(cell.counts)) |(!is.na(hide.small.prop) & is.null(cell.proportions)) )
              warnings("You need to provide either cell counts or cell proportion when you'd like to hide nodes below some thresholds.")
              
              if (!is.null(cell.proportions)){
                z2=pretty(c(min.proportion, max.proportion))
                colinds2=unlist(lapply(1:length(proportions), function(i){max(1,ceiling((proportions-min(z2))*1000/(max(z2)-min(z2)))[i])}))
                cols2=colorFunc2(1000)[colinds2]
                if (!is.na(hide.small.prop)){
                   white.inds<- which(proportions < hide.small.prop)
                   if (length(white.inds)>0)
                    cols2[white.inds]<-"#FFFFFF"
                 }
                if (!is.na(hide.small.count)){
                   white.inds<- which(counts < hide.small.count)
                   if (length(white.inds)>0)
                    cols2[white.inds]<-"#FFFFFF"
                 }
              }else if (!is.null(cell.counts)){
                 z2=pretty(c(min(counts,na.rm=T), max(counts, na.rm=T)))
                colinds2=unlist(lapply(1:length(counts), function(i){max(1,ceiling((counts-min(z2))*1000/(max(z2)-min(z2)))[i])}))
                cols2=colorFunc2(1000)[colinds2]
               if (!is.na(hide.small.count)){
                   white.inds<- which(counts < hide.small.count)
                   if (length(white.inds)>0)
                    cols2[white.inds]<-"#FFFFFF"
                 }

              }
          }
            
            nAttrs <- list()
            nAttrs$fillcolor=rep('white', length(ab@nodes[3,]))
            nAttrs$border.color=rep('darkgrey', length(ab@nodes[3,]))
            nAttrs$fontcolor=rep('black', length(ab@nodes[3,]))
            nAttrs$width=rep(0.75, length(ab@nodes[3,]))
            nAttrs$style=rep('bold', length(ab@nodes[3,]))
            ##nAttrs$shape=rep('ellipse', length(ab@nodes[3,]))
            if (uniformColors==FALSE){
              if (!is.null(cell.proportions)){
                nAttrs$border.color=cols2
                nAttrs$border.lwd=rep(node.lwd, length(ab@nodes[3,]))
              }
              nAttrs$fillcolor=cols
              nAttrs$fontcolor=unlist(sapply(cols, SetTextContrastColor))
            }
            nAttrs$label=rep('',length(V))
            if(nodeLabels)
              nAttrs$label=V
            ##nAttrs$fillcolor=ab$nodes[3,]
            
            names(nAttrs$fillcolor)=V
            names(nAttrs$border.color)=V
            names(nAttrs$fontcolor)=V
            names(nAttrs$width)=V
            names(nAttrs$label)=V
            names(nAttrs$style)=V
            ##names(nAttrs$shape)=V
            if (!is.null(cell.proportions)){
              names(nAttrs$border.color)=V
              names(nAttrs$border.lwd)=V
            }  
            eAttrs <- list()
            eAttrs$lwd=rep(1,length(ab@edges[2,]))
            edgevalues <- as.numeric(ab@edges[2,])
            if(uniformColors==FALSE)
              if (edgeWeights)
                {
                  eAttrs$lwd=as.numeric(unlist(lapply(1:length(ed), function(x) {return(0.75+ 15*(ed[[x]]$weights - min(edgevalues))/(max(edgevalues) - min(edgevalues)))})))##[unlist(lapply(1:length(ed), function(x) {return(ed[[x]]$weights)}))]
                }
            eAttrs$label=rep('',length(ab@edges))
            if(edgeLabels)
              {
            	raw.labels <- unlist(lapply(1:length(ed), function(x) { return(ed[[x]]$labels)}))
              	eAttrs$label <- sapply(raw.labels, decodePhenotype, marker.names, partitions.per.marker,filters)
              }
            eAttrs$color=rep('gray', length(ab@edges[3,]))
            eAttrs$arrowhead=rep('empty', length(ab@edges[3,]))
            names(eAttrs$lwd)=edgeNames(g)
            names(eAttrs$arrowhead)=edgeNames(g)
            names(eAttrs$color)=edgeNames(g)
            names(eAttrs$label)=edgeNames(g)
            
            if ((!uniformColors)){
              if (plot.legend){
                delta.sizex=legend.size/dev.size()[1]
                delta.sizey=legend.size/dev.size()[2]
                ##delta.sizex=0.2
                ##delta.sizey=0.2   
                if (is.null(cell.proportions) | is.null(cell.counts))
                  split.screen(t(cbind(c(0,1-delta.sizex,0,1),c(1-delta.sizex,1,0,1))))
                if (!is.null(cell.proportions)| is.null(cell.counts))
                  split.screen(t(cbind(c(0,1-delta.sizex,delta.sizey,1),c(1-delta.sizex,1,delta.sizey,1),c(0,1-delta.sizex,0,delta.sizey))))
                screen(1)
                par(mar=c(0,0,0,0))
              }
            }
            plot(g, 'dot', attrs=attrs,nodeAttrs=nAttrs,edgeAttrs=eAttrs)
            if (!uniformColors){
              if (plot.legend){
                screen(2)
                par(mar=c(1,4,1,0.2))
                image(matrix(1:2500, 50), col = colorFunc(50), xaxt='n', ylab='', yaxt='n',col.lab=textcol)
                par(mgp=c(2.5,1,0))    
                title(ylab=ylab,col.lab=textcol);
                ##z=z[which(z<max(scores))]
                axis(2, at=c((1/(length(z)-1))*(0:(length(z)-1))), labels=z,col=textcol,col.axis=textcol)
                ##axis(2, at=z*max(scores), labels=z)
                #Commented out legend for proportion legend for now
                #if (!is.null(cell.proportions)){
                #  screen(3);
                 # par(mar=c(4,1,0.2,1))
                  #image(t(matrix(1:2500,50)), col = colorFunc2(50), yaxt='n', xlab='', xaxt='n')
                  #par(mgp=c(2.5,1,0))
                  #title(xlab=xlab,col.lab=textcol);
                  #axis(1, at=c((1/(length(z2)-1))*(0:(length(z2)-1))), labels=z2,col=textcol,col.axis=textcol)
                }
                close.screen(all.screens = TRUE)
              }              
            })
