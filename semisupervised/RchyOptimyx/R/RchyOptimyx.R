setGeneric("RchyOptimyx",function(pheno.codes, 
                                  phenotypeScores, 
                                  startPhenotype, 
                                  pathCount=1, 
                                  trimPaths=FALSE, 
                                  trim.tolerance=0, 
                                  trim.level=0) {standardGeneric("RchyOptimyx")})
setMethod("RchyOptimyx", signature=signature(pheno.codes = "character", 
                           phenotypeScores = "numeric", 
                           startPhenotype = "character", 
                           pathCount = "numeric", 
                           trimPaths = "logical", 
                           trim.tolerance="ANY", 
                           trim.level="ANY"),
          function(pheno.codes, 
                   phenotypeScores, 
                   startPhenotype, 
                   pathCount, 
                   trimPaths, 
                   trim.tolerance=0, 
                   trim.level=0){
            subChar <- function(the.str, old.val, new.val)
              {
                split.string <- strsplit(the.str, '')[[1]]
                split.string[which(split.string==old.val)] <- new.val
                paste(split.string, sep='', collapse="")
              }
#startPhenotype=subChar(startPhenotype,'0','3')
#startPhenotype=subChar(startPhenotype,'1','0')
#startPhenotype=subChar(startPhenotype,'3','1')
            ##startPhenotype[which(startPhenotype==0)] <- 3
            ##startPhenotype[which(startPhenotype==1)] <- 0
            ##startPhenotype[which(startPhenotype==3)] <- 1
            ##startPhenotype=paste(startPhenotype, collapse='')
            
            phenotype.count <- nchar(startPhenotype) + 1
            node.attr.count <- 3
            edge.attr.count <- 5
            ab <- .C("c_analyze", 
                     pheno.codes, 
                     as.integer(length(pheno.codes)), 
                     as.double(phenotypeScores), 
                     as.double(max(phenotypeScores)),
                     as.double(min(phenotypeScores)),
                     as.character(startPhenotype), 
                     as.integer(pathCount),
                     as.integer(trimPaths),
                     as.integer(trim.tolerance),
                     as.integer(trim.level),
                     ## outputs
                     best.paths = rep(paste(rep(' ', 255), collapse=''), pathCount * phenotype.count),
                     node.pvals = double(pathCount * phenotype.count),
                     gnodes = rep(paste(rep(' ', 255), collapse=''), pathCount * phenotype.count * node.attr.count), 
                     gedges = rep(paste(rep(' ', 255), collapse=''), pathCount * phenotype.count * edge.attr.count),
                     node.count = integer(1),
                     edge.count = integer(1))
            nodes <- matrix(ab[[13]][1:(ab[[15]] * node.attr.count)], node.attr.count, ab[[15]])
            edges <- matrix(ab[[14]][1:(ab[[16]] * edge.attr.count)], edge.attr.count, ab[[16]])
            res = list()
            res$nodes <- nodes
            res$edges <- edges
            res$scores <- rowSums(matrix(ab[[12]], pathCount, phenotype.count))
            ##names(res$scores)=nodes[2,]
            return (OptH=new("OptimizedHierarchy", nodes=nodes, edges=edges, pathScores=res$scores));
          })
