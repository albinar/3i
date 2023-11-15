# Written by Albina Rahim on December 2018
# Updated by Marjan Najafabadipour on November 2019

#FlowGroup algorithm code:
#The flowGroup algorithm uses the sum of 4 distance measures to perform 
#heirarchical clustring in order to find homogeneous samples in a dataset. 
#These mesures are: 1. Value Under the Curve (VUC); 2. VUC derivative; 
#3. Histogram of Oriented Gradients (HOG) glotting; and 4. flowType partitioning.'

#Satrting flowGroup function 
flowGroup <- function(dir_path, xMark, yMark, experiment_name = NULL, 
                      plot_groups = F, HOG_cex = 0.2, grid_width = 50, 
                      partitions_flowType = 4, vuc = 1, vuc_d = 1, hog = 1, 
                      flwT = 1, vert_line = NULL, horiz_line = NULL, 
                      DoOverlay = F, DoContour = F, boundaries = NULL, 
                      verbose = F){

  #Fetching the system date and time
  start <- Sys.time()
  #Converting the date object to the character to be used for recording the excution time at each step
  date_time <- strftime(Sys.time(), "%y%m%d_%H%M")
  
  #Detecting the number of CPU cores for performing parallelization
  no_cores <- detectCores() - 2
  if (no_cores <= 0){
    no_cores <- 1
  }
  no_cores<-8
  registerDoMC(no_cores)
  
  #Creating a folder to store all the figures in it
  if (is.null(experiment_name)){
    experiment_name <- date_time
    if(verbose == T){ 
      cat("Experiment name set to: ", date_time, "\n", sep = "")}
  }
  suppressWarnings(dir.create(experiment_name, recursive = T) )
  #Creating 2 other folders to store the primilinary results of distant measures in them
  suppressWarnings(dir.create(paste0(experiment_name,"/HogImages/"), recursive = T))
  suppressWarnings(dir.create(paste0(experiment_name,"/VucDerivatives/"), recursive = T))
  
  #Counting the numebr of input files
  num_of_files <- 0
  for (i in 1:length(dir_path)){
    if(dir.exists(dir_path[i])){
      if(file.exists(dir_path[i])){
        f_temp <- data.frame(list.files(dir_path[i], pattern = "\\.fcs$"))
        n_files <- length(f_temp$list.files.dir_path.i..)
        num_of_files <- num_of_files + n_files
        } else {
          print("No files found in directory", dir_path[i])
          }
      } else {
        print("The directory", dir_path[i], "you have enterd is not found. Please, enter a valid directory")
      }
  }
  
  #if no files found in the path, exit from the flowGroup function
  if (num_of_files == 0){
    return("No files is found in the directory you have entered. Please, enter a valid directory")
  }
  
  #list the file names in the path
  full_file_names <- list.files(dir_path, pattern = "\\.fcs$", full.names = TRUE)
  file_names <- list.files(dir_path, pattern = "\\.fcs$", full.names = FALSE)
  #Finding begin and end of X and Y markers limits
  start_kde2d <- Sys.time()
  xyMark_values <- foreach(q1 = 1:num_of_files) %dopar% {
    f_temp <- load_flowFrame(Fcs_file_name = full_file_names[q1])
    x <- f_temp@exprs[, xMark]
    y <- f_temp@exprs[, yMark]
    xyList <- list(x,y)
    return(xyList) 
  }

  markers_limit_df <- data.frame(marker = integer(), begin = integer()
                            , end = integer())
  for (i in 1:2){
    temp_var <- lapply(xyMark_values, function(x){return(x[[i]])})
    temp_var <- unlist(temp_var, use.names = FALSE)
    begin_value <- min(temp_var)
    end_value <- max(temp_var)
    result_df <- data.frame("marker" = i, "begin" =  begin_value, "end" = end_value)
    markers_limit_df <- rbind( markers_limit_df, result_df)
  }
  
  xlim = c(markers_limit_df$begin[1], markers_limit_df$end[1])
  ylim = c(markers_limit_df$begin[2], markers_limit_df$end[2])
  if(verbose == T){ cat("Time for finding the begin and end limits of markers: ", TimeOutput(start_kde2d),"\n",sep="") }
  
'#Performing Quantile function for getting the begin and the end limit for markers
  quantile_df <- data.frame(marker = integer(), begin = integer()
                            , end = integer())
  for (i in 1:2){
    temp_var <- lapply(xyMark_values, function(x){return(x[[i]])})
    temp_var <- unlist(temp_var, use.names = FALSE)
    result  <- quantile(temp_var, probs = c(0.02, 0.98))
    result_df <- data.frame("marker" = i, "begin" = result[[1]], "end" = result[[2]])
    quantile_df <- rbind(quantile_df, result_df)
    }
  
  xlim = c(quantile_df$begin[1], quantile_df$end[1])
  ylim = c(quantile_df$begin[2], quantile_df$end[2])
  if(verbose == T){ cat("Time for the quantile function: ", TimeOutput(start_kde2d),"\n",sep="") }'
  
  ########################################################################
  if (vuc_d != 0 || vuc != 0) {
        # density representation
        start_kde2d <- Sys.time()
        kde2d_res <- foreach(q1 = 1:num_of_files) %dopar% {
            f_temp <- load_flowFrame(Fcs_file_name = full_file_names[q1])
            f1 <- kde2d(f_temp@exprs[,xMark], f_temp@exprs[,yMark], n = grid_width, lims = c(xlim,ylim))
            return(f1) # this returns roughly the same value. On one test it can a +- 1.6% difference.
        }
    }
  
  #VUC - calculate distant matrix
  if (vuc != 0) {
    volume_perc_under_both_curves <- foreach(q1 = 1:num_of_files, .combine = rbind, .maxcombine = num_of_files, .multicombine = T) %dopar% {
      f1 <- kde2d_res[[q1]]
      ratio_value <- vector(length = num_of_files)
      for (q2 in 1:num_of_files){
        if (q2 >= q1) {
          ratio_value[q2] <- 0
          } else {
            f2 <- kde2d_res[[q2]]
            trap_value <- sapply(1:grid_width^2, function(x) { abs(f1$z[x] - f2$z[x]) })
            max_value <- sapply(1:grid_width^2, function(x) { max(c(f1$z[x], f2$z[x])) })
            ratio_value[q2] <- sum(trap_value)/sum(max_value)
          }
        }
      return(ratio_value)
      }
      
    row.names(volume_perc_under_both_curves) <- 1:num_of_files
    # distance is not euclidean or manhattan, it is my own custom distance.
    dist_vuc <- stats::as.dist(volume_perc_under_both_curves) # I use as.dist because I want my values to be in the distance matrix. I dont want my points wo be thought of as coordinates and for dist() to calculate something I dont want / is incorrect.

    # hclust and dendrogram
    hc_vuc <- hclust(dist_vuc)
    png(paste0(experiment_name, "/", date_time,"_VUC_Dendrogram.png"), width = 1500, height = 1000)
    plot(hc_vuc, main = "Volume under the curve")
    dev.off()

    dist_vuc_norm <- dist_vuc/max(dist_vuc)

    if(verbose == T){ cat("Time for the volume under both curves method: ", TimeOutput(start_kde2d), "\n", sep = "") }
    } else {
      dist_vuc <- stats::dist(matrix(0, num_of_files, num_of_files))
      dist_vuc_norm <- stats::dist(matrix(0, num_of_files, num_of_files))
      if(verbose == T){ cat("Skipped the volume under both curves method.\n", sep = "") }
    }
  
    ########################################################################
    # VUC derivative
    if (vuc_d != 0) {
      start_kde2d <- Sys.time()
      derivative <- foreach(q1 = 1:num_of_files) %dopar% {
        derivative = grad(h = kde2d_res[[q1]]$z, x = kde2d_res[[q1]]$x, y = kde2d_res[[q1]]$y)
        return(derivative)
        }
      doubleDerivative <- foreach(q1 = 1:num_of_files) %dopar% {
        doubleDerivative <- sqrt(derivative[[q1]]$gx^2 + derivative[[q1]]$gy^2)
        return(doubleDerivative)
      }
      
      Do_plots <- T
      if (Do_plots == T){
        figs_temp <- foreach(q2 = 1:num_of_files, .combine = list, .maxcombine = num_of_files, .multicombine = T) %dopar% {
          CairoPNG (file = paste0(experiment_name, "/VucDerivatives/",str_pad(q2, 3, pad = "0"),".png"), width = 2000, height = 2000)
          par(mfrow = c(2,2), mar = c(3, 3, 1, 1), mgp = c(2, 0.7, 0))
          image2D(kde2d_res[[q2]]$z, xlab = "", ylab = "", frame = FALSE, axes = F, main = "")
          image2D(derivative[[q2]]$gx, xlab = "", ylab = "", frame = FALSE, axes = F, main = "")
          image2D(derivative[[q2]]$gy, xlab = "", ylab = "", frame = FALSE, axes = F, main = "")
          image2D(doubleDerivative[[q2]], xlab = "", ylab = "", frame = FALSE, axes = F, main = "")
          dev.off()
          return(1)
        }
        }

        # calculate distant matrix
        volume_perc_under_both_curves_dervative <- foreach(q1 = 1:num_of_files, .combine = rbind, .maxcombine = num_of_files, .multicombine = T) %dopar% {
            f1 <- doubleDerivative[[q1]]
            ratio_value <- vector(length = num_of_files)
            for ( q2 in 1:num_of_files){
                if (q2 >= q1) {
                    ratio_value[q2] <- 0
                } else {
                    f2 <- doubleDerivative[[q2]]
                    trap_value <- sapply(1:grid_width^2, function(x) { abs(f1[[x]] - f2[[x]]) })
                    max_value <- sapply(1:grid_width^2, function(x) { max(c(f1[[x]], f2[[x]])) })
                    ratio_value[q2] <- sum(trap_value)/sum(max_value)
                }
            }
            return(ratio_value)
        }
        
        row.names(volume_perc_under_both_curves_dervative) <- 1:num_of_files
        # distance is not euclidean or manhattan, it is my own custom distance.
        dist_vuc_derv <- stats::as.dist(volume_perc_under_both_curves_dervative) # I use as.dist because I want my values to be in the distance matrix. I dont want my points wo be thought of as coordinates and for dist() to calculate something I dont want / is incorrect.

        # hclust and dendrogram
        hc_vuc <- hclust(dist_vuc_derv)
        png(paste0(experiment_name,"/", date_time,"_VUC_Derv_Dendrogram.png"), width = 1500, height = 1000)
        plot(hc_vuc, main = "Volume under the curve of the derivatives")
        dev.off()

        dist_vuc_derv_norm <- dist_vuc_derv/max(dist_vuc_derv)

        if(verbose == T){ cat("Time for the volume under both curves of the derivatives method: ", TimeOutput(start_kde2d), "\n",sep = "") }
        } else {
          dist_vuc_derv <- stats::dist(matrix(0, num_of_files, num_of_files))
          dist_vuc_derv_norm <- stats::dist(matrix(0, num_of_files, num_of_files))
          if(verbose == T){ cat("Skipped the volume under both curves of the derivative method.\n",sep="") }
          }

    ########################################################################
    # HOG plotting
    if (hog != 0) {
      start_HOG <- Sys.time()
      Do_plots <- T
      if (Do_plots == T){
        figs_temp <- foreach(q1 = 1:num_of_files, .combine = list, .maxcombine = num_of_files, .multicombine = T) %dopar% {
          CairoPNG (file = paste0(experiment_name, "/HogImages/", str_pad(q1, 3, pad = "0"),".png"), width = 4000, height = 4000)
          par(mar=c(1,1,1,1)) # margins
          f_temp <- load_flowFrame(Fcs_file_name = full_file_names[q1])
          plot(f_temp@exprs[,xMark],f_temp@exprs[,yMark], xlab= "", ylab = "", frame = FALSE,axes = F,
               cex.main = 2, cex.lab = 2, cex.axis = 2, cex = HOG_cex, pch = 19, main = "", xlim = xlim, ylim = ylim)
          dev.off()
        }
      }
      
      res_hog <- HOG_apply(paste0(experiment_name, "/HogImages/"), cells =  100, orientations = 12, threads = no_cores)
      
      dist_hog <- stats::dist(res_hog$hog, method = "manhattan") # with cells being 100 and orientations being 12, there are 120000 data values for each file/picutre. Manhattan is better than euclidean for large dimensions.
      hc_hog <- hclust(dist_hog, method = "complete") # default
      png(paste0(experiment_name, "/", date_time,"_HOG_Dendrogram_100_4k.png"), width = 1500, height = 1000)
      plot(hc_hog, main = "HOG", xlim = c(0,4), ylim = c(0,4))
      dev.off()
      
      dist_hog_norm <- dist_hog/max(dist_hog)
      
      if(verbose==T){ cat("Time for the HOG method: ", TimeOutput(start_HOG),"\n",sep="") }
    } else {
      dist_hog <- stats::dist(matrix(0, num_of_files, num_of_files))
      dist_hog_norm <- stats::dist(matrix(0, num_of_files, num_of_files))
      if(verbose == T){ cat("Skipped the HOG method.\n", sep = "") }
    }
  

    ########################################################################
    # flowType partitioning
    if (flwT != 0) {
      start_flowType <- Sys.time()
      flowType_res <- foreach(q1 = 1:num_of_files ) %dopar% {
        f_temp <- load_flowFrame(Fcs_file_name = full_file_names[q1])
        if(!is.na(xlim) && !is.na(ylim)){
          xThres <- seq(from = xlim[1], to = xlim[2], length.out = partitions_flowType+1)
          yThres <- seq(from = ylim[1], to = ylim[2], length.out = partitions_flowType+1)
          } else {
            xThres <- seq(from = min(f_temp@exprs[,xMark]), to = max(f_temp@exprs[,xMark]), length.out = partitions_flowType+1)
            yThres <- seq(from = min(f_temp@exprs[,yMark]), to = max(f_temp@exprs[,yMark]), length.out = partitions_flowType+1)
            }
        xThres <- xThres[-c(1,length(xThres))] # do not want partitions of nothing on the borders
        yThres <- yThres[-c(1,length(yThres))] # do not want partitions of nothing on the borders
        res <- flowType(Frame = f_temp, PropMarkers = c(xMark, yMark), MFIMarkers = c(xMark, yMark), Methods = "Thresholds",
                        MarkerNames = c(xMark, yMark), Thresholds = list(xThres, yThres), MemLimit = 10, verbose = T,
                        PartitionsPerMarker = c(length(xThres)+1,length(yThres)+1) );
        rem.ind <- grep(x = res@PhenoCodes, pattern = "0")
        CellFreqs <- res@CellFreqs[-rem.ind]/res@CellFreqs[1]
        return(CellFreqs)
        }

        flowType_matrix <- t(cbind(sapply(1:length(flowType_res), function(x){as.vector(flowType_res[[x]])})))
        dist_flwT <- stats::dist(flowType_matrix, method = "manhattan") # choose manhattan over euclidean because there are 16 dimensions, and manhattan performs better with large amount of dimensions.
        hc_flwT <- hclust(dist_flwT)
        
        # png(paste0(experiment_name, "/", date_time,"_flowtype_Dendrogram.png"),width=1500,height=1000)
        # plot(hc_flwT, main = "flowType partitions")
        # dev.off()

        dist_flwT_norm <- dist_flwT/max(dist_flwT)

        if(verbose == T){ cat("Time for the flowType method: ", TimeOutput(start_flowType), "\n", sep= "") }
        } else {
          dist_flwT <- stats::dist(matrix(0, num_of_files, num_of_files))
          dist_flwT_norm <- stats::dist(matrix(0, num_of_files, num_of_files))
          if(verbose == T){ cat("Skipped the flowType method.\n", sep = "") }
        }
    
    ########################################################################
    # combine 4 distance measures
    dist_all <- vuc*dist_vuc_norm + vuc_d*dist_vuc_derv_norm + hog*dist_hog_norm + flwT*dist_flwT_norm
    hclust_group <- "complete"
    suppressWarnings( dir.create(paste0(experiment_name,"/", hclust_group), recursive = T) )
    dist_all_norm <- (dist_all-min(dist_all))/(max(dist_all)-min(dist_all))
    hc_all <- hclust(dist_all_norm, method = hclust_group)
    png(paste0(experiment_name, "/", hclust_group, "/", date_time,"_combined_Dendrogram_", hclust_group, ".png"),width=1500,height=1000)
    plot(hc_all, main = "combined")
    dev.off()

    ########################################################################
    # best k values
    k_count <- NULL
    seq_t <- seq(min(hc_all$height), max(hc_all$height), length.out = 100)
    for(t1 in 1:100){
        k_count[t1] <- length(unique(cutree(hc_all, h=seq_t[t1])))
    }

    find_best_k <- sort(table(k_count),decreasing=TRUE)
    ind_rem <- sort(unique(c(which(as.numeric(names(find_best_k)) > 50), which(find_best_k <= 1)))) # dont want clusters of size 1 or k too large where the largeness is overweighting the score
    find_best_k <- find_best_k[-ind_rem]
    temp <- as.numeric(names(find_best_k))
    temp[temp>=15] <- 15
    scores <- temp*find_best_k # random scoring system. I want a medium sized k with many points. (i.e. large vertical gab on dendrogram)
    best_k <- as.numeric(names(which(scores >= 0.85*max(scores))))

    # make the table nice enough to display
    scores_vis <- t(as.matrix(scores))
    scores_vis <- rbind(scores_vis, scores_vis)
    scores_vis[1,] <- as.numeric(colnames(scores_vis))
    row.names(scores_vis) <- c("Best k","Score")
    colnames(scores_vis) <- rep("", ncol(scores_vis))
    scores_vis <- scores_vis[,sort(scores_vis[2,], index.return = T, decreasing = T)$ix]
    if(verbose==T){print(scores_vis[,1])}

    ###############################################################
    #plot groups
    #Preperation for generating plots
    set.seed(4342)
    res.tsne <- Rtsne(dist_all_norm, theta = 0.0, perplexity = 1)
    prin_comp <- prcomp(dist_all_norm) # PCA
    return_cutree_order <- list()
    order_of_plotting_total <- list()
    k_to_run <- sort(c(unique(c(scores_vis[1,1]))))

    cutree_order <- cutree(hc_all,k=k_to_run)
    names(cutree_order) <- NULL
    printing_order <- unique(cutree_order[hc_all$order])
    group_ind <- cutree_order[hc_all$order]
    width <- sort(table(group_ind),decreasing=TRUE)[1]
    order_of_plotting <- lapply(1:k_to_run, function(x){temp <- hc_all$order[which(group_ind==printing_order[x])]
    c(temp,rep(num_of_files+1,width-length(temp)))})
    order_of_plotting_total[[length(order_of_plotting_total)+1]] <- order_of_plotting
    colours <- as.numeric(sapply(1:num_of_files, function(x) {names(find_number(x, order_of_plotting))}))
    colours <- rainbow(k_to_run+2)[colours]
    return_cutree_order[[length(return_cutree_order)+1]] <- cutree_order
    
    # tSNE plots
    png(paste0(experiment_name, "/", hclust_group, "/tsne_", date_time,"_", k_to_run, "_", hclust_group, ".png"), width = 1000, height = 1000, pointsize = 18)
    plot(res.tsne$Y, col=1:num_of_files, pch = '.', xlab = '', ylab = '')
    text(res.tsne$Y, labels = c(1:num_of_files), cex = 0.8, col = colours)
    dev.off()
    
    # PCA plots
    png(paste0(experiment_name, "/", hclust_group, "/pca_", date_time,"_", k_to_run, "_", hclust_group, ".png"), width = 1600, height = 1600, pointsize = 18)
    par(mfrow = c(2,2), mar = c(3, 3, 1, 1), mgp = c(2, 0.7, 0))
    plot(prin_comp$x[,c(1,2)], pch = '.', col = return_cutree_order[k_to_run], xlab = '', ylab = '')
    text(prin_comp$x[,c(1,2)], labels = c(1:num_of_files), cex=0.8, col=colours)
    plot(prin_comp$x[,c(1,3)], pch = '.', col = return_cutree_order[k_to_run], xlab = '', ylab = '')
    text(prin_comp$x[,c(1,3)], labels = c(1:num_of_files), cex = 0.8, col = colours)
    plot(prin_comp$x[,c(2,3)], pch='.', col = return_cutree_order[k_to_run], xlab = '', ylab ='')
    text(prin_comp$x[,c(2,3)], labels = c(1:num_of_files), cex = 0.8, col = colours)
    dev.off()
    
    #Creating a data farme with results about groups
    df <- data.frame(file_names, full_file_names, groups = cutree_order)
    df <- df[order(df[,3]),]
    
    #Counting number of samples in each group
    n_samples_in_groups <- data.frame(count(df, groups))
    n_samples_in_groups <- n_samples_in_groups[order(n_samples_in_groups[,2], decreasing = TRUE),]
    n_samples_in_groups <- transform(n_samples_in_groups, id = match(groups, unique(groups)))
    n_samples_in_groups <- select(n_samples_in_groups, -c(groups))
    names(n_samples_in_groups) <- c("count", "groups")
    
    #Changing the order of groups by setting the group with largest number of sample as group 1
    df <- data.frame(df %>% group_by(groups) %>% mutate(count = n()))
    df <- df[order(df[,4], decreasing = TRUE),]
    df <- transform(df, id = match(groups, unique(groups)))
    df <- select(df, -c(groups))
    names(df) <- c("file_names", "full_file_names", "count", "groups")
    
    #Setting directory for stroing clustering results in .csv format
    result_path <- paste(as.character(experiment_name),"/", as.character(hclust_group)
    ,"/csv_results")
    result_path <- gsub(x = result_path, pattern = " ", replacement = "")
    if(!dir.exists(result_path)){
      dir.create(paste0(experiment_name,"/", hclust_group, "/csv_results/"))
    }
    if(!file.exists(result_path)){
      file.remove(file.path(result_path, "result.csv"))
      write.csv(df, file.path(result_path, "result.csv"))
    } else {
      write.csv(df, file.path(result_path, "result.csv"))
    }
    
    #Plotting samples of each group
    for (n_of_groups in 1:k_to_run){
      start_groups <- Sys.time()
      if(plot_groups){
        df_temp <- data.frame(df$full_file_names[which(df$groups == n_samples_in_groups$groups[n_of_groups])])
        df_temp1 <- data.frame(df$file_names[which(df$groups == n_samples_in_groups$groups[n_of_groups])])
        colnames(df_temp) <- "fileNames"
        colnames(df_temp1) <- "fileNames"
        suppressWarnings(dir.create(paste0(experiment_name,"/", hclust_group, "/" , "group", n_of_groups, "/"), recursive = T))
        foreach(q1 = 1:n_samples_in_groups$count[n_of_groups])  %dopar% {
          f_temp <- load_flowFrame(Fcs_file_name = as.character(df_temp$fileNames[q1]))
          png(paste0(experiment_name, "/", hclust_group, "/" , "group", n_of_groups, "/", as.character(df_temp1$fileNames[q1]), ".png"), height = 800, width = 800, pointsize = 18)
          par(mar = c(rep(0.25,4)))
          plotDens(f_temp, c(xMark,yMark), cex.main=2, cex.lab=4, cex.axis=4, main="", axes=T, xlab="", ylab="", xaxt='n', yaxt='n', xlim = xlim, ylim = ylim)
          
          # Add overlay 1D plot
          if(DoOverlay){
            overlay.dens <- density(f_temp@exprs[,xMark])
            overlay.dens$y <- overlay.dens$y / max(overlay.dens$y)
            if (identical(ylim, NA)) {
              overlay.scale <- max(f_temp@exprs[,yMark])
              overlay.translation <- min(f_temp@exprs[,yMark])
            } else {
              overlay.scale <- (ylim[2]-ylim[1])
              overlay.translation <- ylim[1]
            }
            par(new = T)
            lines(overlay.dens$x, overlay.dens$y*overlay.scale + overlay.translation, col="grey48")
            
            #horiz
            overlay.dens <- density(f_temp@exprs[,yMark])
            overlay.dens$y <- overlay.dens$y / max(overlay.dens$y)
            if (identical(xlim, NA)) {
              overlay.scale <- max(f_temp@exprs[,xMark])
              overlay.translation <- min(f_temp@exprs[,xMark])
            } else {
              overlay.scale <- (xlim[2]-xlim[1])
              overlay.translation <- xlim[1]
            }
            par(new = T)
            lines(overlay.dens$y*overlay.scale + overlay.translation, overlay.dens$x, col="grey48")#)
          }
          
          # Add contour lines
          if(DoContour){
            z <- kde2d(f_temp@exprs[,xMark], f_temp@exprs[,yMark], n = 50)
            z$z <- (z$z)^(1/4)
            contour(z, drawlabels = FALSE, add = TRUE, nlevels = 10, lty = 2, col="grey49")
          }
          
          if (!is.null(vert_line)){   
            abline(v=vert_line[q1])                    
          }
          if (!is.null(horiz_line)){
            abline(h=horiz_line[q1])
          }
          
          # should replace with labels parameter and with upper left corner instead of relying on the xlim and ylim
          if (!is.na(xlim) && !is.na(ylim)){
            text(x = (xlim[1]), y = (ylim[2]), labels = as.character(df_temp1$fileNames[q1]), adj = c(0, 1), cex=1.75)
          }
          
          if (!is.null(boundaries)){
            lines(boundaries[[q1]]@boundaries, lwd = 2)
          }
          
          box(lwd = 2)
          dev.off()
        }
      }
      if(verbose==T){ cat("Time to plot for ", "group ", n_of_groups, ": ", TimeOutput(start_groups),"\n",sep="")}
    }

    cat("Total time: ", TimeOutput(start),"\n", sep = "")
    
    return(list(Best_k = best_k, k_scores = scores_vis, 
                Groups = return_cutree_order, dist_all = dist_all, 
                dist_all_norm = dist_all_norm, 
                order_of_plotting = order_of_plotting_total))
}
