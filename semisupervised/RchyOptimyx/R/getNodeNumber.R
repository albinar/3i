setGeneric("getNodeNumber",function(edges, nodes){standardGeneric("getNodeNumber")})
setMethod("getNodeNumber", signature=signature(edges="vector", nodes="vector"),
          function(edges, nodes){
            return(c(which(nodes==strsplit(edges,'~')[[1]][1]),which(nodes==strsplit(edges,'~')[[1]][2])))
          })

