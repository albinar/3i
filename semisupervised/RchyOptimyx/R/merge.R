##Merges optimized hierarchies x and y
setGeneric("merge",function(x, y) {standardGeneric("merge")})
setMethod("merge", signature=signature(x="OptimizedHierarchy", y="OptimizedHierarchy"),
          function(x,y)
          {
            res1<-x
            res2<-y
            res<-res1
            node.count<-dim(res@nodes)[2]
            for (i in 1:dim(res2@nodes)[2])
              {
                if (!res2@nodes[1,i] %in% res1@nodes[1,])
                  {
                    res@nodes<-cbind(res@nodes, res2@nodes[,i])
                    node.count<-node.count+1
                  }
              }
            
            edge.count<-dim(res@edges)[2]
            for (i in 1:dim(res2@edges)[2])
              {
                if (!res2@edges[1,i] %in% res1@edges[1,])
                  {
                    res@edges<-cbind(res@edges, res2@edges[,i])
                    edge.count<-edge.count+1
                  }
              }
            res@pathScores<-c(res@pathScores,res2@pathScores)
            return(res)
          })
