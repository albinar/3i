##########################
##########################
########################## 

#Naive way of generating WT rank based on what I understood from the call
rank.generator <- function(ranks)
{
  dup<-duplicated(ranks)
  wt.ranks<-setdiff(1:140, unique(ranks))
  for ( i in 1:length(ranks))
  {
    if(dup[i])
      wt.ranks <-setdiff(wt.ranks,wt.ranks[wt.ranks-ranks[i]>0][1])
    
  }
  return(wt.ranks)
}



find.farthest.pheno <- function(pheno,max.limit=300, min.limit=5)
{
  pheno.mat <-sapply(pheno,function(x) as.numeric(unlist(strsplit(x, ""))))
  rownames(pheno.mat)<-colnames(pheno.mat) <- NULL
  distance <- rdist(t(pheno.mat))
  dmax <-apply(distance,2,which.max)
  dmax <- unique(dmax)
  mat <- pheno.mat[,dmax]
  prev.ind <- length(dmax)
  while(ncol(mat)>min.limit)
  {
    distance <- rdist(t(mat))
    dmax <-unique(apply(distance,2,which.max))
    ind <-sort(dmax,decreasing = T,index.return=T)$ix[1:min(max.limit,length(dmax))]
    ind2<-unique(apply(distance[ind,],1, which.max))
    mat <-mat[,ind2]
    if (length(ind2)==1)
      return(paste(mat, sep="",collapse=""))
    if ( ncol(mat)==prev.ind){
      break
    }else{
      prev.ind <- ncol(mat)
    }
  }
  return(apply(mat,2,paste,sep="",collapse=""))
}
##########################
##########################
########################## 
Get.radii<- function(pop,ft,all.pops,which.filter.neg=NA)
  #filter.to.threshold)
  #match.marker is a vector that correponds to marker used for flowType. For filters, just grab on of the marker for now
{

  #color<- colors()[c(26,47,555,114,88,503,431,103,370)]
  code <- ft@PhenoCodes[which(all.pops==pop)]
  if (length(code)==0)
    stop("Population can not be found")
  code.string <-unlist(lapply(code, function(x) unlist(strsplit(x, ""))))
  # for ( i in filter.to.threshold)
  # {
  #   
  # }
  if (!is.na(which.filter.neg)){
    
    radii<- sapply(1:length(code.string), function(x){
              if (which.filter.neg[x]==T)
                return( ifelse(code.string[x]=="0",yes = .1,no = ifelse(code.string[x]=="2",yes =1.3,no = 5)))
              else
                return( ifelse(code.string[x]=="0",yes = .1,no = ifelse(code.string[x]=="2",yes =5,no = 1.3)))
      
      })
  }else{
  radii<- sapply(code.string, function(x) return( ifelse(x=="0",yes = .1,no = ifelse(x=="2",yes =5,no = 1.3))))
  }
  return(radii)
  
}
##########################
##########################
########################## 
pval.calc <- function(props, label,geneders=NULL,varaince, WT=F,which.sex=NULL,wilcox=F, rank.wilcox=F,method="BH")
{
  if (WT)
  {
    print("Comparing wildtype to KO")
    label.1  <- which(label==which.sex)
    label.2 <- which(label=="WT") 
  }else{
    label.1  <- which(label=="Female")
    label.2 <- which(label=="Male")
  }
  pvals <- apply(props,2, function(prop.x) {
    if (all(prop.x==1))
      return(1)
    if (wilcox)
    { 
      p<-tryCatch(wilcox.test(prop.x[label.1],prop.x[label.2],var.equal = varaince)$p.value, error=function(x) {return("wrong")})
      
    }else{
      p<-tryCatch(t.test(prop.x[label.1],prop.x[label.2],var.equal = varaince)$p.value, error=function(x) {return("wrong")})
    }
    if (mode(p)=="character" | is.nan(p))
    {
      #print("t.test failed, as data was almost constant.")
      p<-1
    }
    if (rank.wilcox)
    { 
      m_rank<- rank(c(prop.x[which(genders$ko=="Male")],prop.x[which(genders$wt=="Male")]))
      f_rank <- rank(c(prop.x[which(genders$ko=="Female")],prop.x[which(genders$wt=="Female")]))
      
      #And this is significant: 
      p2<-tryCatch(wilcox.test(x = c( m_rank[1:length(which(genders$ko=="Male"))], f_rank[1:length(which(genders$ko=="Female"))] ),
                               y= c( m_rank[(length(which(genders$ko=="Male"))+1): length(m_rank)], f_rank[(length(which(genders$ko=="Female"))+1): length(f_rank)] ),
                               var.equal = varaince)$p.value, error=function(x) {return("wrong")})
    }
    
    if (mode(p2)=="character"| is.nan(p2))
    {
      #print("t.test failed, as data was almost constant.")
      p2<-1
    }
    return(c(p,p2))
  })
  p.w<- unlist(lapply(pvals, function(x) x[1]))
  p.r<- unlist(lapply(pvals, function(x) x[2]))
  mean.1 <- apply(props[label.1,],2,median)
  mean.2 <- apply(props[label.2,],2,median)
  names(pvals)<- colnames(props)
  pvalsw.adjust <- p.adjust(p.w, method=method)
  pvalsr.adjust <- p.adjust(p.r, method=method)
  return(cbind(p.w,pvalsw.adjust, mean.1,mean.2,p.r,pvalsr.adjust))
}

##########################
##########################
##########################
Find.WTsamples<-function(g1, allProportions,ko.ind,wt.ind, sex)
{
  ind1 <- which(allProportions[ko.ind,"Gender"]==sex)
  group <- c("Female","Male")
  if (length(ind1)==0)
  {
    warning(paste("There were not any", sex, "KO for this gene. Taking the range from the opposite sex instead.", sep=" "))
    ind1 <- which(allProportions[ko.ind,"Gender"]==group[which(group!=sex)])
  }  
  ind2 <- which(allProportions[wt.ind,"Gender"]==sex)
  assay.date<-as.Date(allProportions[ko.ind[ind1],"Assay Date"],format = "%d-%b-%Y")
  wt.date <- as.Date(allProportions[wt.ind[ind2],"Assay Date"],format = "%d-%b-%Y")
  t1<-assay.date[order(assay.date)][1]
  t2 <- tail(assay.date[order(assay.date)],1)
  # on.range1 <-which( wt.date ==t1 )
  within.range <-  which( wt.date >=t1 & wt.date  <=t2)
  #on.range2 <-  which(wt.date  ==t2)
  out.range.left <- which(wt.date  <t1)
  out.range.right <- which(wt.date  >t2)
  
  if (length(within.range)>=70)
  {
    set.seed(240)
    wt.range <- sample(size =70,within.range)
  }else{
    needed <- 70-length(within.range)
    # if (needed%%2==0)
    # {
    if (length(out.range.left)>=needed/2)
    {
      wt.l.range <-tail(out.range.left,floor(needed/2))
      if (length(out.range.right)<needed/2)
      {
        wt.r.range <-out.range.right
        wt.l.range <- tail(out.range.left,needed-length(out.range.right))
      }else{
        wt.r.range <-tail(out.range.right,ceiling(needed/2))
      }
    }else
    {
      wt.l.range <-out.range.left
      if (length(out.range.right)<(needed/2-length(out.range.left)))
      {
        print("Not enough WT sample")
        next;
      }
      wt.r.range <- tail(out.range.right,needed-length(out.range.left))
    }
    wt.range <- c( wt.l.range,within.range,wt.r.range)
  }
  
  return(wt.range)
}
###################################################################################
###################################################################################
###################################################################################

draw.violinplot <- function(reduced.genes,ft,pop.names,all.props,plot.path,markers,which.filter.neg=NA)
{
  counter <- 1
  
  while(counter<length(reduced.genes)+1)
  {
    gene <- reduced.genes[counter]
    selected.pheno <- ft@PhenoCodes[match(names(which(hits[,2]==gene)),pop.names)]
    if (length(selected.pheno)>10)
      selected.pheno <- find.farthest.pheno(pheno = selected.pheno,min.limit = 8)
    temp <-all.props[[which(unlist(lapply(all.props, function(x) return(ifelse(length(which(x$genotype==gene))>0, yes =1,no=0 ))))==1)]]
    for ( i in 1:length(selected.pheno))
    {
      
      
      prop.hits <-as.vector(temp$props[,match(selected.pheno[i],ft@PhenoCodes)])
      pheno <- as.vector(sapply( colnames(temp$props[,c(2,match(selected.pheno[i],ft@PhenoCodes))])[2],rep, length(temp$labels)))
      g<- rep(gene,length(temp$labels))
      labels <- paste(temp$gender, temp$labels)
      groups<-temp$labels
      sex <- temp$gender
      size <- as.vector(sapply(groups, function(x) return(ifelse(x=="KO",yes = 3,no = 2.8))))
      dat <- data.frame(Phenotypes=pheno,Labels=labels ,Proportions=prop.hits, Gene=g,Sex=sex, Groups=groups,size=size)  
      
      dat$Labels <- factor(dat$Labels,
                           levels =c("Female WT", "Male WT","Female KO", "Male KO"),ordered = TRUE)
      dat$Groups <- factor(dat$Groups, levels=c("WT","KO"))
      plot.name <-paste(unlist(strsplit(gene, "_")),"__",sep="",collapse="")
      png(paste(plot.path,plot.name,unique(pheno),"_KOWT.png",sep=""),pointsize = 20,width = 1000,height = 700)
      p1 <- ggplot()+#fill=Groups))+
        geom_violin(data=dat[which(dat$Groups=="WT"),],width=.5, mapping=aes(x=Labels,y=Proportions,fill=Labels),position = position_dodge(width = 0.2),alpha=.7)+
        scale_fill_manual(values = c("Female KO"="#0978A1", "Male KO"="#0978A1","Female WT"="#F26300", "Male WT"="#F26300"),breaks=NULL,labels=NULL,name=NULL)+
        guides(fill=guide_legend(title=NULL))+
        geom_jitter(data=dat,mapping=aes(x=Labels,y=Proportions,shape=Sex,color=Groups),cex=3,position = position_jitter(width = .07))+ 
        scale_color_manual(values = c("KO"="#0978A1","WT"="#F26300"))+ labs(y="%-Proportion of CD45")+
        theme(#title=element_text( hjust = 1, size=16,color=1),
          axis.text.x = element_text( hjust = 1, size=12,color=1),
          axis.text.y = element_text( hjust = 1, size=12,color=1),
          #legend.text = element_text( hjust = 1, size=12,color=1),
          strip.text.x = element_text(size = 6.5))
      
      DF<- data.frame(variable=as.factor(1:length(markers)), 
                      markers
                      ,values=as.vector(Get.radii(pop = unique(pheno),ft = ft,all.pops = pop.names,which.filter.neg = which.filter.neg)))
      
      plot <- ggplot(DF, aes(markers, values, fill = markers,label=markers)) +
        geom_bar(width = 1, stat = "identity", color = "white") +geom_text(aes(y = 5.5,label = markers),size=2.8)+
        scale_y_continuous(breaks = 0:nlevels(DF$variable)) +guides(fill=FALSE)+
        coord_polar()+theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          panel.grid  = element_blank())
      plot <- plot+scale_fill_brewer(palette="Paired")
      pfinal<-plot_grid(plot, p1, align = "v", nrow = 2, rel_heights = c(.7, 1),scale = 0.99)
      title <- ggdraw() + draw_label(paste(gene,unique(pheno),sep=": "), fontface='bold')
      print(plot_grid(title, pfinal, ncol=1, rel_heights=c(0.07, 1)) )
      dev.off()
    }
    counter <- counter+1
    
  }
}

##################################################################################
#Rchy
##################################################################################
generate.rchy.plots<- function(output.path,marker.names, all.pvals,reduced.genes,ft,filters,svg=F)
{
library(RchyOptimyx)
counter <- 1
prop.hits<-c()
pheno <-c()
labels <- c()
g<-c()

pheno.names<-all.pvals[which(all.pvals[,3]==unlist(strsplit(reduced.genes[1],"_"))[1] & all.pvals[,4]==unlist(strsplit(reduced.genes[1],"_"))[2]),"Cell type"]

while(counter<length(reduced.genes)+1)
{
  gene <- reduced.genes[counter]
  selected.pheno <- ft@PhenoCodes[match(names(which(hits[,2]==gene)),pheno.names)]
  if (length(selected.pheno)>10)
    selected.pheno <- find.farthest.pheno(pheno = selected.pheno,min.limit = 8)
  rchy.res <- list()
  if (is.na(unlist(strsplit(gene,"_"))[2]))
  {
    pvals <- as.numeric(all.pvals[which(all.pvals[,3]==unlist(strsplit(gene,"_"))[1] ),"Adjusted P-value"])   
  }else{
 pvals <- as.numeric(all.pvals[which(all.pvals[,3]==unlist(strsplit(gene,"_"))[1] & all.pvals[,4]==unlist(strsplit(gene,"_"))[2]),"Adjusted P-value"])
}
  library(RchyOptimyx)
  for (i in 1:length(selected.pheno))
  {
    rchy.res[[i]]<-RchyOptimyx(pheno.codes=ft@PhenoCodes, phenotypeScores=-log10(pvals),
                               startPhenotype=selected.pheno[i], pathCount=3,trimPaths=FALSE, trim.level = 6) 
  }
  merged.res <- merge(rchy.res[[1]],rchy.res[[2]])
  if (length(rchy.res)>2)
  {
  for ( i in 3:length(rchy.res))
    merged.res<-merge(merged.res,rchy.res[[i]])
  }
  path <- paste(input.path,panel,"/semisupervised/",organ,"/RchyOptimyx/",sep="")
  suppressWarnings(dir.create(path))
  #if(!svg)
    png(paste(path,gene,"_RchyOptimyx.png",sep=""),width = 1900,height = 1600,res=150,pointsize = 12 )
  #else
    #svg(filename = paste(path,gene,"_RchyOptimyx.svg",sep=""),width = 12,height = 8,pointsize = 12)
  plot(merged.res,phenotypeCodes=ft@PhenoCodes[ind.Pval.Effect],phenotypeScores=-log10(pvals),root.name=paste(gene,ifelse(panel=="Panel_M-cell",yes = "Lin-",no = ": CD45+")),legend.size=1.2,
       partitions.per.marker=2,filters=filters,
       node.lwd = 1,marker.names=marker.names, uniformColors = FALSE)
  dev.off()
  counter<-counter+1
}
}

