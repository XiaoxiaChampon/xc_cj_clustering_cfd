

######################################################################
#
# Purpose: Functions to cluster Twitter Categorical Functional Data 
# Author:  
# Date: 4/6/2023
# Revised : Date: 12/18/2024
##############################################################

#############################################################################################################
#############################################################################################################
#Part I: Simulation 
########################################################################################################

# For: profiling and visualization of profiling
library(profvis)

# For: gam 
library(mgcv)

# For: cubicspline
library(pracma)

# For: fpca.face
library(refund)

# For: FADPclust
library(FADPclust)

# For: kNNdist
library(dbscan)

# For: elbow
devtools::install_github("ahasverus/elbow")
library(elbow)

# For: rand.index
library(fossil)

# For: cfda method
library(cfda)

# For: gather method
library(tidyverse)

# ---- For: parallelization ----
# For: foreach loop
library(foreach)

run_parallel <- TRUE
time_elapsed <- list()
if(run_parallel)
{
  print("RUNNING PARALLEL")
  
  # For: makeCluster
  library(doParallel)
  
  # For: %dorng% or registerDoRNG for reproducable parallel random number generation
  library(doRNG)
  
  if(exists("initialized_parallel") && initialized_parallel == TRUE)
  {
    parallel::stopCluster(cl = my.cluster)
  }
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  doParallel::registerDoParallel(cl = my.cluster)
  cat("Parellel Registered: ", foreach::getDoParRegistered(), "\n")
  initialized_parallel <- TRUE
  
  # registerDoRNG(123) # ///<<<< THIS CREATES THE ERROR FOR FADPClust !!!
}

#' Create directories 
if (!dir.exists("outputs")){
  dir.create("outputs")
}
if (!dir.exists("outputs/clustersims")){
  dir.create("outputs/clustersims")
}

# profvis({


#**** evaluate cluster accuracy using RI and ARI
# INPUT
# VECTORS with true memberships (true_cluster)
#        with estimated memberships (coming from cluster_method$label)
#noise: scalar, 0 or 3, if dbscan clustering, use 0, else, use 3.
# OUTPUT: RI and ARI, CPN
#
evaluate_cluster <- function(true_cluster, new_cluster,noise){
  ri=rand.index(true_cluster,new_cluster)
  ari=adj.rand.index(true_cluster, new_cluster)
  
  c3=which(true_cluster==noise)
  cc3=which(new_cluster==noise)
  cpn=sum(cc3 %in% c3)/length(c3)
  return(list(ri=ri, ari=ari,cpn=cpn))
}

#**** clustering the scores using KMEANS
#INPUT
# n by K scores - or features extracted from the multivariate FD
#OUTPUT
# list with 2 elements: nclust (number of clusters), label (vector with cluster membership)
#
kmeans_cluster <- function(data=scores_K){
  out_kmeans = NbClust::NbClust(data = data, diss = NULL,
                                distance = "euclidean", min.nc = 2, max.nc = 5,
                                method = "kmeans",index="silhouette")
  
  return(list(nclust=as.numeric(out_kmeans$Best.nc[1]), label = out_kmeans$Best.partition))
}


#**** clustering the scores using KMEANS
#INPUT
# n by K scores - or features extracted from the multivariate FD
# scale_eps, scalar, the scale of the epsilon
#OUTPUT
# list with 2 elements: nclust (number of clusters), label (vector with cluster membership)
#
dbscan_cluster <- function(data=scores_K, scale_eps){
  
  dimz=dim(data)[2]
  if (dimz<=2){minPts=4}
  if (dimz>2){minPts=2*dimz+1}
  dist=kNNdist(data, k = minPts-1)
  #########change to max increase
  distdataelbow=data.frame(sort(dist))
  distdataelbow$index=1:(dim(data)[1])
  ipoint <- elbow(data = distdataelbow)
  epsoptimal=(ipoint$sort.dist._selected)*scale_eps
  
  out_dbscan <- dbscan(data, eps =epsoptimal , minPts = minPts)
  return(list(nclust=dim(table(out_dbscan$cluster)), label = out_dbscan$cluster))
}


#*** clustering the curves directly using FADP
#INPUT
# Z1,Z2 - m by n latent processes for the categ FD with 3 categories
# tt - grid of points
# PVE - percentage of explained variance (default value 0.95)
# OUTPUT
# list with 2 elements: nclust (number of clusters), label (vector with cluster membership)
#
fadp_cluster <- function(mZ1, mZ2, tt=tt, PVE=0.95){
  
  fdabasis = fda::create.bspline.basis(c(0,1), 25, 4)
  fdatime = tt
  fdafd1 = fda::smooth.basis(tt, mZ1, fdabasis)$fd
  fdafd2 = fda::smooth.basis(tt, mZ2, fdabasis)$fd
  
  FADlist=list(fdafd1,fdafd2)
  
  FADP <- FADPclust(fdata = FADlist, cluster = 2:5, method = "FADP1",
                    proportion = seq(0.02, 0.2, 0.02), f.cut = 0.15,
                    pve = PVE, stats = "Avg.silhouette")
  
  return(list(nclust=FADP$nclust, label = FADP$clust))
}



#*** extract scores using UNIVARIAFE fpca, FAST COV estimation.  SUPER FAST
#INPUT
# Z1,Z2 - m by n latent processes for the categ FD with 3 categories
# tt - grid of points
# PVE - percentage of explained variance (default value 0.95)
# OUTPUT
# list with estimated PHI (m by K) with function norm 1
#           projection of data onto PHI (called scores - nby K matrix)
#

extract_scores_UNIVFPCA <- function (mZ1,mZ2, tt , PVE=0.95)
{
  m<- nrow(mZ1)
  n<-ncol(mZ1)
  
  out1 <- fpca.face(Y=t(mZ1), argvals =tt, pve = 0.99)
  out2 <- fpca.face(Y=t(mZ2), argvals =tt, pve = 0.99)
  O_m1 <- matrix(0, nrow=m, ncol=out1$npc)
  O_m2  <- matrix(0, nrow=m, ncol=out2$npc)
  
  #construct PHI
  Phi_est0 <-  rbind(cbind(out1$efunctions*sqrt(m),O_m2 ),
                     cbind(O_m1,out2$efunctions*sqrt(m)))
  Scores0 <- cbind(out1$scores, out2$scores)
  ScoresCov_0 <- cov(Scores0 )
  oute <- eigen(ScoresCov_0)
  
  K<- which(cumsum( oute$values)/sum(oute$values)>=PVE)[1]
  count_iter = 0
  delta=0.01
  while (K<2 && count_iter<100) {
    count_iter = count_iter + 1
    cat("count_iter: ", count_iter, "\n")
    K<- which(cumsum( oute$values)/sum(oute$values)>=(PVE+delta))[1]
    delta=delta+0.01
  }
  
  Phi_est <-  Phi_est0%*% oute$vectors[,1:K] # correct eigenfns
  
  mZ <- rbind(mZ1, mZ2)
  Scores_est <- t(mZ) %*%Phi_est/sqrt(m)  # they are not demeaned
  
  return (list(scores=Scores_est, Phi= Phi_est))
}

#Function to use trapzfnum function and find L2 distance for 2D array, n of them
#' @param  truecurve- 2D array (true curve)
#' @param   estcurve-2D array (estimated curve)
#' @return : scalar-L2 distance

mse_bw_matrix <- function(truecurve,estcurve,timestamps01)
{
  
  n=dim(truecurve)[2]
  # datapoints=dim(truecurve)[1]
  # mseall=c(0)
  ######could probably use apply function here it's also subject level
  mseall <- foreach(i = 1:n, .combine = c, .packages = c("pracma")) %dorng% {
    source("R/acj/trapzfnum_function.R")
    return(rbind(trapzfnum(truecurve[,i], estcurve[,i],timestamps01)))
  }
  
  return(mseall)
}


#Function to use trapzfnump function and find Hellinger distance for 2D array, n of them
mse_bw_matrixp <- function(truecurve,estcurve,timestamps01)
{
  n=dim(truecurve)[2]
  # datapoints= dim(truecurve)[1]
  # mseall=c(0)
  # ######could probably use apply function here it's also subject level
  # for (i in 1:n)
  # {
  #   mseall[i]=trapzfnump(truecurve[,i],estcurve[,i])/sqrt(2)
  # }
  
  sqrt_2 <- sqrt(2)
  mseall <- foreach(i = 1:n, .combine = c, .packages = c("pracma")) %dorng% {
    source("R/acj/trapzfnum_function.R")
    return(rbind(trapzfnump(truecurve[,i], estcurve[,i],timestamps01)/sqrt_2))
  }
  
  return(mseall)
}


#cfda
####
#write a function to produce the scores
#input one  X_nt matrix n*t: 1 of  Nth simulation, n*t* N
#output scores matrix for that specific Nth simulation   n*M     M is the column number of scores
cfda_score_function <- function(cfda_data, nCores,timestamps01,basis_num ){
  cfda_data <- t(cfda_data)
  timeseries_length <- dim(cfda_data)[2]
  times <- seq(1,timeseries_length,by=1)
  colnames(cfda_data) <- paste0("time_",times)
  time <- rep(timestamps01 ,times=dim(cfda_data)[1])
  subjects <- seq(1,dim(cfda_data)[1],by=1)
  cfda_new <- data.frame(cbind(cfda_data,subjects))
  data_long <- gather(cfda_new, time_t, state, paste0("time_",times[1]):paste0("time_",times[timeseries_length]),factor_key=TRUE)
  cfda_d <- data_long%>%arrange(subjects)%>%dplyr::rename(id=subjects)%>%dplyr::select(id,state)
  cfda_final <- data.frame(cbind(cfda_d,time))
  Tmax <- max(cfda_final$time)
  cfda_cut <- cut_data(cfda_final, Tmax = Tmax)
  b <- create.bspline.basis(c(0, Tmax), nbasis = basis_num, norder = 4)
  fmca <- compute_optimal_encoding(cfda_cut, b, nCores =nCores, verbose = FALSE)
  delta <- 0.01
  pve <- 0.95
  nPc90 <-which(cumsum(prop.table(fmca$eigenvalues)) > pve)[1]
  while (nPc90<2 && pve < 1){
    pve <- pve+delta
    nPc90 <-which(cumsum(prop.table(fmca$eigenvalues)) > pve)[1]
  }
  cfda_score <- fmca$pc[, 1:nPc90]
  return(cfda_score)
}




ClusterSimulation <- function(num_indvs, timeseries_length,
                              scenario, num_replicas, est_choice, run_hellinger, temp_folder)
{
  cat("Cluster Simulation\nNum Indvs:\t", num_indvs,
      "\nTimeseries Len:\t", timeseries_length,
      "\nScenario:\t", scenario,
      "\nNum Replicas:\t", num_replicas)
  
  occur_fraction <- GetOccurrenceFractions(scenario)
  cat("\nOccurrence Fractions: ", occur_fraction)
  
  cluster_allocation <- occur_fraction * num_indvs
  cat("\nIndividuals in each cluster (cluster alloc.): ", cluster_allocation)
  
  true_cluster <- rep(c(1,2,3), cluster_allocation)
  cat("\nTrue Cluster:\n", true_cluster)
  
  true_cluster_db <- rep(c(1,2,0), cluster_allocation)
  cat("\nTrue Cluster DB:\n", true_cluster_db)
  
  rmse<- array(0, c(num_replicas, 2, 3))       #2 components and 3 settings
  hellinger <- array(0, c(num_replicas, 3, 3)) #3 events probab and 3 settings
  true_kmeans <- est_kmeans  <- NULL # records clustering membership on the TRUE SCORES/ESTIMATED SCORES
  true_fadp <- est_fadp <-  NULL # records clustering membership on the TRUE Z/ESTIMATED Z
  true_dbscan <- est_dbscan <-  NULL
  
  time_elapsed <<- list()
  # "Xiaoxia"=NULL, "univfpca"=NULL, "kmeans"=NULL, "fadp"=NULL, "dbscan"=NULL, "cfd"=NULL)
  last_time <- 0
  row_name <- NULL
  timeKeeperStart <- function(rn)
  {
    row_name <<- rn
    if(FALSE == row_name %in% names(time_elapsed))
    {
      time_elapsed[[row_name]] <<- NULL
    }
    last_time <<- Sys.time()
  }
  timeKeeperNext <- function()
  {
    this_time <- Sys.time()
    this_section_time <- this_time - last_time
    cat(row_name, "calc time taken:", capture.output(this_section_time), "\n")
    time_elapsed[[row_name]] <<- append(time_elapsed[[row_name]], this_section_time)
    last_time <<- this_time
  }
  
  ######
  Z_true_curve <- Z_est_curve <- list()
  p_true_curve <- p_est_curve <- list()
  W_cfd <- list()
  total_regens <- 0
  for(replica_idx in 1:num_replicas)
  {
    cat("\nNum Indvs:", num_indvs,
        "\tTimeseries Len:", timeseries_length,
        "\tScenario:", scenario,
        "\tNum Replicas:", num_replicas)
    cat("\nReplica: ", replica_idx, "\n")
    
    # generate clusters
    # set.seed(seed_cluster + 100 * replica_idx)
    # cat("\nCluster", replica_idx, " --> seed: ", seed_cluster + 100 * replica_idx, "\n")
    cat("Cluster", replica_idx, "\n")
    
    cluster_f1 <- GenerateClusterData(1, scenario, 3, cluster_allocation[1], timeseries_length)
    cluster_f2 <- GenerateClusterData(2, scenario, 3, cluster_allocation[2], timeseries_length)
    cluster_f3 <- GenerateClusterData(3, scenario, 3, cluster_allocation[3], timeseries_length)
    
    # Recover the latent Gaussian process --> is this always 2 ???
    Z1 <- cbind(cluster_f1$Z1, cluster_f2$Z1, cluster_f3$Z1)
    Z2 <- cbind(cluster_f1$Z2, cluster_f2$Z2, cluster_f3$Z2)
    
    
    # Recover the true probability curves --> could there be more than 3 ???
    p1 <- cbind(cluster_f1$p1, cluster_f2$p1, cluster_f3$p1)
    p2 <- cbind(cluster_f1$p2, cluster_f2$p2, cluster_f3$p2)
    p3 <- cbind(cluster_f1$p3, cluster_f2$p3, cluster_f3$p3)
    prob_curves <- list(p1=p1, p2=p2, p3=p3)
    #########9/11/2023  
    Z_true_curve[[replica_idx]]=array(c(Z1,Z2),dim=c(timeseries_length,num_indvs,2))
    p_true_curve[[replica_idx]]=array(c(p1,p2,p3),dim=c(timeseries_length,num_indvs,3))
    ############
    
    
    # generate categFuncData
    # set.seed(seed_cfd + 100 * replica_idx)
    # cat("\nCategFD", replica_idx, " --> seed: ", seed_cfd + 100 * replica_idx, "\n")
    cat("CategFD", replica_idx, "\n")
    
    categ_func_data_list <- GenerateCategFuncData(prob_curves)
    
    # what is Q vals ? better name???
    Q_vals <- unique(c(categ_func_data_list$W))
    if(is.numeric(Q_vals))
    {
      Q_vals <- sort(Q_vals)
    }
    
    # I need to know what this loop is meant to do !??? maybe there is a better way
    for(indv in 1:num_indvs)
    {
      if(indv %in% 1:cluster_allocation[1])
      {
        setting_choice <- 1
      }
      if(indv %in% (cluster_allocation[1] + 1):(cluster_allocation[1] + cluster_allocation[2]))
      {
        setting_choice <- 2
      }
      if(indv %in% (cluster_allocation[1] + cluster_allocation[2] + 1):num_indvs)
      {
        setting_choice <- 3
      }
      
      # better names for the following variables ???
      # 1. check weather one category only appears 1 time and is it in the end of the timeseries
      # 2. OR is it appearing only one time in the begining
      # 3. OR if the category is less thant Q total category
      # In general, one category only occurs 2 times
      # If timepoints=300, one category only occurs less than 4 times 3/300=0.01
      #If timepoints=750, one category only occurs less than 6 times 5/750=0.0067
      tolcat <- table(categ_func_data_list$W[,indv])
      catorder <- order(tolcat, decreasing = TRUE)
      numcat <- length(catorder)
      refcat <- catorder[numcat]
      count_iter <- 0
      while (count_iter < 100 && 
             ( (numcat < length(Q_vals))
               ||(timeseries_length==300  && min(as.numeric(tolcat)) < 4)
               ||(timeseries_length==750  && min(as.numeric(tolcat)) < 10)
             )
      )
      {
        count_iter <- count_iter + 1
        
        new_cluster_data <- GenerateClusterData(setting_choice, scenario, 3, 5, timeseries_length)
        
        new_prob_curves <- list(p1 = new_cluster_data$p1, p2 = new_cluster_data$p2, p3 = new_cluster_data$p3)
        new_categ_func_data_list <- GenerateCategFuncData( new_prob_curves )
        
        # what is this 3 ?? arbitrarily chosen?
        categ_func_data_list$W[, indv] <- new_categ_func_data_list$W[, 3]
        Z1[, indv] <- new_cluster_data$Z1[, 3] # latent curves Z1 and Z2
        categ_func_data_list$X[indv, , ] <- 0
        Z2[, indv] <- new_cluster_data$Z2[, 3]
        
        for (this_time in 1:timeseries_length)
        {
          categ_func_data_list$X[indv, this_time, which(Q_vals == categ_func_data_list$W[, indv][this_time])] <- 1
        }
        
        tolcat <- table(categ_func_data_list$W[, indv])
        catorder <- order(tolcat, decreasing = TRUE)
        numcat <- length(catorder)
        refcat <- catorder[numcat]
      } # end while
      total_regens <- total_regens + count_iter
    } # end for(indv in 1:num_indvs)
    
    #Estimation
    cat("EstimateCategFuncData", replica_idx, "\n")
    timestamps01 <- seq(from = 0.0001, to = 1, length=timeseries_length)
    timeKeeperStart("Xiaoxia")
    categFD_est <- EstimateCategFuncData(est_choice, timestamps01, categ_func_data_list$W)
    #####################9/11/2023
    Z_est_curve[[replica_idx]]=array(c(categFD_est$Z1_est,categFD_est$Z2_est),dim=c(timeseries_length,num_indvs,2))
    p_est_curve[[replica_idx]]=array(c(categFD_est$p1_est,categFD_est$p2_est,categFD_est$p3_est),dim=c(timeseries_length,num_indvs,3))
    W_cfd[[replica_idx]]=categ_func_data_list$W
    #####################9/11/2023
    timeKeeperNext()
    
    if (run_hellinger)
    {
      # evaluate performance Z and P
      rmse1_temp <- c(by(mse_bw_matrix(Z1, categFD_est$Z1_est, timestamps01) , true_cluster, mean))
      rmse2_temp <- c(by(mse_bw_matrix(Z2, categFD_est$Z2_est, timestamps01), true_cluster, mean))
      rmse[replica_idx, ,] <- rbind(rmse1_temp, rmse2_temp )
      
      error.p1 <- mse_bw_matrixp(p1, categFD_est$p1_est, timestamps01)
      error.p2 <- mse_bw_matrixp(p2, categFD_est$p2_est, timestamps01)
      error.p3 <- mse_bw_matrixp(p3, categFD_est$p3_est, timestamps01)
      
      
      hellinger[replica_idx, ,] <-  rbind( c(by(error.p1, true_cluster, mean)),
                                           c(by(error.p2, true_cluster, mean)),
                                           c(by(error.p3, true_cluster, mean)))
      
    }
    
    # Record clustering performance. USE only Z
    #**** extract_scores_UNIVFPCA is the FASTEST
    
    mfpca_true <-extract_scores_UNIVFPCA(mZ1=Z1, mZ2=Z2,  tt=timestamps01 , PVE=0.95)
    #plot(  scores_true$scores[, 1:2])
    timeKeeperStart("univfpca")
    mfpca_est <- extract_scores_UNIVFPCA(mZ1=categFD_est$Z1_est, mZ2=categFD_est$Z2_est, tt=timestamps01, PVE=0.95)
    timeKeeperNext()
    
    #KMEANS
    true_kmeans_temp <- kmeans_cluster(data=mfpca_true$scores)$label
    timeKeeperStart("kmeans")
    est_kmeans_temp <- kmeans_cluster(data=mfpca_est$scores)$label
    timeKeeperNext()
    
    #FADP
    true_fadp_temp <- fadp_cluster(mZ1=Z1, mZ2=Z2, tt=timestamps01)$label
    timeKeeperStart("fadp")
    est_fadp_temp <- fadp_cluster(mZ1=categFD_est$Z1_est, mZ2=categFD_est$Z2_est, tt=timestamps01)$label
    timeKeeperNext()
    #dbscan
    true_dbscan_temp <- dbscan_cluster(data=mfpca_true$scores,1)$label
    timeKeeperStart("dbscan")
    est_dbscan_temp <- dbscan_cluster(data=mfpca_est$scores,1)$label
    timeKeeperNext()
    ##dbscan cfda
    ##########no cfda
    basis_num <- 10
    nCores= n.cores
    true_dbscan_temp_cfda <- true_cluster_db
    parallel::stopCluster(cl = my.cluster)
    timeKeeperStart("cfd")
    cfd_scores <- cfda_score_function(categ_func_data_list$W, nCores,timestamps01,basis_num )
    my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
    doParallel::registerDoParallel(cl = my.cluster)
    cat("Parellel Registered: ", foreach::getDoParRegistered(), "\n")
    est_dbscan_temp_cfda <- dbscan_cluster(data=cfd_scores,1)$label
    ##########no cfda
    timeKeeperNext()
    #record results
    true_dbscan<- cbind(true_dbscan, true_dbscan_temp)
    est_dbscan <- cbind(est_dbscan, est_dbscan_temp)
    
    true_dbscan_cfda<- cbind(true_dbscan, true_dbscan_temp_cfda)
    est_dbscan_cfda <- cbind(est_dbscan, est_dbscan_temp_cfda)
    
    true_kmeans<- cbind(true_kmeans, true_kmeans_temp)
    est_kmeans <- cbind(est_kmeans,  est_kmeans_temp)
    
    true_fadp<- cbind(true_fadp, true_fadp_temp)
    est_fadp <- cbind(est_fadp, est_fadp_temp)
    
    cat("Done replica:", replica_idx, "\n")
  }# END of "for(replica_idx in 1:num_replicas)'
  
  cat("\n replicas done \n")
  
  if (run_hellinger)
  {
    # Assess accuracy in the simulation study
    mse_sim= apply(rmse, c(2,3), mean)
    hellinger_sim= apply(hellinger, c(2,3), mean)
  }
  
  ### Code could be simplified into a loop???
  results<- list(
    #dbscan
    true_dbscan_ri=apply(true_dbscan, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster_db, new_cluster=cluster,0)$ri),
    true_dbscan_ari=apply(true_dbscan, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster_db, new_cluster=cluster,0)$ari),
    true_dbscan_cpn=apply(true_dbscan, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster_db, new_cluster=cluster,0)$cpn),
    
    
    est_dbscan_ri=apply(est_dbscan, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster_db, new_cluster=cluster,0)$ri),
    est_dbscan_ari=apply(est_dbscan, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster_db, new_cluster=cluster,0)$ari),
    est_dbscan_cpn=apply(est_dbscan, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster_db, new_cluster=cluster,0)$cpn),
    #dbscan cfda
    true_dbscan_ri_cfda=apply(true_dbscan_cfda, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster_db, new_cluster=cluster,0)$ri),
    true_dbscan_ari_cfda=apply(true_dbscan_cfda, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster_db, new_cluster=cluster,0)$ari),
    true_dbscan_cpn_cfda=apply(true_dbscan_cfda, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster_db, new_cluster=cluster,0)$cpn),
    
    
    est_dbscan_ri_cfda=apply(est_dbscan_cfda, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster_db, new_cluster=cluster,0)$ri),
    est_dbscan_ari_cfda=apply(est_dbscan_cfda, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster_db, new_cluster=cluster,0)$ari),
    est_dbscan_cpn_cfda=apply(est_dbscan_cfda, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster_db, new_cluster=cluster,0)$cpn),
    ###kmeans
    true_kmeans_ri=apply(true_kmeans, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster, new_cluster=cluster,3)$ri),
    true_kmeans_ari=apply(true_kmeans, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster, new_cluster=cluster,3)$ari),
    true_kmeans_cpn=apply(true_kmeans, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster, new_cluster=cluster,3)$cpn),
    
    
    
    est_kmeans_ri=apply(est_kmeans, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster, new_cluster=cluster,3)$ri),
    est_kmeans_ari=apply(est_kmeans, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster, new_cluster=cluster,3)$ari),
    est_kmeans_cpn=apply(est_kmeans, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster, new_cluster=cluster,3)$cpn),
    
    #fadp
    true_fadp_ri=apply(true_fadp, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster, new_cluster=cluster,3)$ri),
    true_fadp_ari=apply(true_fadp, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster, new_cluster=cluster,3)$ari),
    true_fadp_cpn=apply(true_fadp, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster, new_cluster=cluster,3)$cpn),
    
    
    est_fadp_ri=apply(est_fadp, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster, new_cluster=cluster,3)$ri),
    est_fadp_ari=apply(est_fadp, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster, new_cluster=cluster,3)$ari),
    est_fadp_cpn=apply(est_fadp, 2, function(cluster)
      evaluate_cluster(true_cluster=true_cluster, new_cluster=cluster,3)$cpn)
    
  )
  
  print("cluster_table")
  cluster_table_true=c(mean(results$true_dbscan_ri_cfda),
                       mean(results$true_dbscan_ari_cfda),
                       mean(results$true_dbscan_cpn_cfda),
                       
                       mean(results$true_fadp_ri),
                       mean(results$true_fadp_ari),
                       mean(results$true_fadp_cpn),
                       
                       mean(results$true_kmeans_ri),
                       mean(results$true_kmeans_ari),
                       mean(results$true_kmeans_cpn),
                       
                       mean(results$true_dbscan_ri),
                       mean(results$true_dbscan_ari),
                       mean(results$true_dbscan_cpn)
  )
  names(cluster_table_true)=c("cfda-db RI","cfda-db ARI","cfda-db cpn",
                              "fadp RI","fadp ARI","fadp cpn",
                              "kmeans RI","kmeans ARI","kmeans cpn",
                              "dbscan RI","dbscan ARI","dbscan cpn"
  )
  
  cluster_table_est=c( mean(results$est_dbscan_ri_cfda),
                       mean(results$est_dbscan_ari_cfda),
                       mean(results$est_dbscan_cpn_cfda),
                       
                       mean(results$est_fadp_ri),
                       mean(results$est_fadp_ari),
                       mean(results$est_fadp_cpn),
                       
                       mean(results$est_kmeans_ri),
                       mean(results$est_kmeans_ari),
                       mean(results$est_kmeans_cpn),
                       
                       mean(results$est_dbscan_ri),
                       mean(results$est_dbscan_ari),
                       mean(results$est_dbscan_cpn)
  )
  names(cluster_table_est)=c("cfda-db RI","cfda-db ARI","cfda-db cpn",
                             "fadp RI","fadp ARI","fadp cpn",
                             "kmeans RI","kmeans ARI","kmeans cpn",
                             "dbscan RI","dbscan ARI","dbscan cpn"
  )
  
  
  cluster_table_est_se=c(sd(results$est_dbscan_ri_cfda)/sqrt(num_replicas),
                         sd(results$est_dbscan_ari_cfda)/sqrt(num_replicas),
                         sd(results$est_dbscan_cpn_cfda)/sqrt(num_replicas),
                         
                         sd(results$est_fadp_ri)/sqrt(num_replicas),
                         sd(results$est_fadp_ari)/sqrt(num_replicas),
                         
                         sd(results$est_fadp_cpn)/sqrt(num_replicas),
                         
                         
                         sd(results$est_kmeans_ri)/sqrt(num_replicas),
                         sd(results$est_kmeans_ari)/sqrt(num_replicas),
                         sd(results$est_kmeans_cpn)/sqrt(num_replicas),
                         
                         sd(results$est_dbscan_ri)/sqrt(num_replicas),
                         sd(results$est_dbscan_ari)/sqrt(num_replicas),
                         sd(results$est_dbscan_cpn)/sqrt(num_replicas))
  names(cluster_table_est_se)=c("cfda-db RI","cfda-db ARI","cfda-db cpn",
                                "fadp RI","fadp ARI","fadp cpn",
                                "kmeans RI","kmeans ARI","kmeans cpn",
                                "dbscan RI","dbscan ARI","dbscan cpn")
  print("returning")
  
  save(time_elapsed, total_regens, file=file.path(temp_folder, paste("time_elapsed_", num_indvs, "_", timeseries_length, "_",
                                                                     scenario, "_", num_replicas, "_", est_choice, "_", run_hellinger, ".RData", sep="")))
  
  return_vals <- list("cluster_table_true"=cluster_table_true,
                      "cluster_table_est"=cluster_table_est,
                      "cluster_table_est_se"=cluster_table_est_se,
                      "Z_true_curves"=Z_true_curve,
                      "Z_est_curves"=Z_est_curve,
                      "p_true_curves"=p_true_curve,
                      "p_est_curves"=p_est_curve,
                      "W_cfd"=W_cfd)
  
  save(return_vals, file=file.path(temp_folder, paste("ClusterSim_", num_indvs, "_", timeseries_length, "_",
                                                      scenario, "_", num_replicas, "_", est_choice, "_", run_hellinger, ".RData", sep="")))
  
  if (run_hellinger)
  {
    return_vals$mse <- mse_sim
    return_vals$hellinger <- hellinger_sim
  }
  return(return_vals)
}

#' Function to produce functional dummy variables X from categorical functional data W
#' @param W 2D array, t*n: t is the timestamp and n is the number of the observation
#' @return X 3D array, n*t*Q, Q: the total number of the category
GetXFromW <- function(W)
{
  num_indv <- ncol(W)
  timeseries_length <-nrow(W)
  category_count<- length(unique(c(W)))
  Q_vals <- unique(c(W))
  if(is.numeric(Q_vals)) Q_vals<- sort(Q_vals)
  
  X<- array(0, c(num_indv,timeseries_length,category_count))
  for(indv in 1:num_indv)
  {
    for(timestamps01 in 1:timeseries_length)
    {
      X[indv, timestamps01, which(Q_vals==W[, indv][timestamps01])] <- 1
    }
  }
  return(X)
}

#' Function to select 
#' @param choice "probit", "binomial",  or "multinormial"
#' @param timestamps01, 1D array, time interval that cfd is observed
#' @param  W: 2D array, t*n, t: the number of time points, n: the number of individuals
#' @param  basis_size=25, the number of basis function used 
#' @param  method="ML"
#' @return list of 2D array: True Z curves , Est Z curves, True p curves, Est p curves
#'                           all have dimension t*n
EstimateCategFuncData <- function(choice, timestamps01, W, basis_size=25, method="ML")
{
  if(choice == "probit"){
    X <- GetXFromW(W)
    return(EstimateCategFuncData_probit(timestamps01, X, basis_size, method, 1/150))
  }else if(choice == "binomial"){
    X <- GetXFromW(W)
    return(EstimateCategFuncData_binorm(timestamps01, X, basis_size, method))
  }else if(choice == "multinomial"){
    return(EstimateCategFuncData_multinormial(timestamps01, W, basis_size, method))
  }
}

#'Function to estimate z and p using wood_multinormial
#' @param timestamps01, 1D array, time interval that cfd is observed
#' @param  W: 2D array, t*n, t: the number of time points, n: the number of individuals
#' @param  basis_size=25, the number of basis function used 
#' @param  method="ML"
#' @return list of 2D array: True Z curves , Est Z curves, True p curves, Est p curves
#'                           all have dimension t*n

EstimateCategFuncData_multinormial <- function(timestamps01, W, basis_size=25, method="ML")
{
  
  num_indv<- ncol(W)
  timeseries_length <-nrow(W)
  
  Z<-NULL
  prob<-array(0, c(num_indv, timeseries_length , 3))
  for (i in 1:num_indv){
    fit_binom<-gam(list(W[,i]-1~s(timestamps01,bs = "cr", m=2, k = basis_size),
                        ~s(timestamps01,bs = "cr", m=2, k = basis_size)),
                   family=multinom(K=2), method = method,
                   control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                   optimizer=c("outer","bfgs")) 
    
    z1<- fit_binom$linear.predictors[,1]
    z2<- fit_binom$linear.predictors[,2]
    Z<- cbind(Z, c(z1,z2))
    ##find probability
    Z_cbind=cbind(z1,z2)
    exp_z=exp(Z_cbind)
    denominator_p=1+exp_z[,1]+exp_z[,2]
    p1 <- exp_z[,1]/denominator_p
    p2 <- exp_z[,2]/denominator_p
    p3=1/denominator_p
    prob[i,,] <- cbind(p1, p2, p3)
    
  }
  
  return(list(Z1_est=Z[1:timeseries_length ,], Z2_est=Z[1:timeseries_length +timeseries_length ,], 
              p1_est=t(prob[,,1]), p2_est=t(prob[,,2]), p3_est=t(prob[,,3]) ))
}


#'Function to estimate z and p using probit
#' @param timestamps01, 1D array, time interval that cfd is observed
#' @param  X: 3D array, t*n*Q, t: the number of time points, n: the number of individuals, Q: the number of categories
#' @param  basis_size=25, the number of basis function used 
#' @param  method="ML"
#' @return list of 2D array: True Z curves , Est Z curves, True p curves, Est p curves
#'                           all have dimension t*n
EstimateCategFuncData_probit <- function(timestamps01, X, basis_size=25, method="ML", threshold_probability=0.004)
{
  num_indv<- dim(X)[1]
  timeseries_length<- dim(X)[2]
  category_count <- dim(X)[3]
  
  Z<-NULL
  p<-array(0, c(num_indv, timeseries_length, category_count))
  
  
  zp <- foreach (indv = 1:num_indv, .combine = cbind, .init = NULL, .packages = c("mgcv")) %dorng%
    {
      source("R/acj/run_gam_function.R")
      
      x1<- X[indv,,1]
      x2<- X[indv,,2]
      x3<- X[indv,,3]
      
      probit_binom <- function(x_binary){
        if (sum(x_binary)/timeseries_length < threshold_probability){
          gam_result_binary <- RunGam(timestamps01, x_binary, "probit", basis_size, method)
          p_binary <- gam_result_binary$prob
          p_binary_linpred <- gam_result_binary$linpred
        }else{
          gam_result_binary <- RunGam(timestamps01, x_binary, "binomial", basis_size, method)
          p_binary <- gam_result_binary$prob
          p_binary_linpred <- gam_result_binary$linpred
        }
        return(list("p_binary"=p_binary,"p_binary_linpred"=p_binary_linpred))
      }
      
      r_1 <- probit_binom(x1)
      p1 <- r_1$p_binary
      p1_linpred <- r_1$p_binary_linpred
      
      r_2 <- probit_binom(x2)
      p2 <- r_2$p_binary
      p2_linpred <- r_2$p_binary_linpred
      
      r_3 <- probit_binom(x3)
      p3 <- r_3$p_binary
      p3_linpred <- r_3$p_binary_linpred
      
      # estimate the latent curves Z
      denominator_p <- 1 + exp(p3_linpred)
      z1 <- (p1_linpred-p3_linpred)-log( (1+exp(p1_linpred))/(denominator_p))
      z2 <- (p2_linpred-p3_linpred)-log( (1+exp(p2_linpred))/(denominator_p))
      
      psum <- p1 + p2 + p3
      return(c(c(z1,z2), cbind(p1/psum, p2/psum, p3/psum)))
    }
  # Unravel the two variables from zp
  z_rows_count <- timeseries_length * 2
  Z <- array(zp[1:z_rows_count, ], c(z_rows_count, num_indv))
  p <- array(t(matrix(zp[(z_rows_count + 1):dim(zp)[1], ], ncol=num_indv)), c(num_indv, timeseries_length, category_count))
  
  
  return(list(Z1_est=Z[1:timeseries_length,],
              Z2_est=Z[1:timeseries_length+timeseries_length,],
              p1_est=t(p[,,1]),
              p2_est=t(p[,,2]),
              p3_est=t(p[,,3]) ))
}

#'Function to estimate z and p using binom
#' @param timestamps01, 1D array, time interval that cfd is observed
#' @param  X: 3D array, t*n*Q, t: the number of time points, n: the number of individuals, Q: the number of categories
#' @param  basis_size=25, the number of basis function used 
#' @param  method="ML"
#' @return list of 2D array: True Z curves , Est Z curves, True p curves, Est p curves
#'                           all have dimension t*n

EstimateCategFuncData_binorm <- function(timestamps01, X, basis_size=25, method="ML")
{
  num_indv<- dim(X)[1]
  timeseries_length<- dim(X)[2]
  category_count <- dim(X)[3]
  
  Z<-NULL
  p<-array(0, c(num_indv, timeseries_length, category_count))
  
  
  zp <- foreach (indv = 1:num_indv, .combine = cbind, .init = NULL, .packages = c("mgcv")) %dorng%
    {
      source("R/acj/run_gam_function.R")
      
      x1<- X[indv,,1]
      x2<- X[indv,,2]
      x3<- X[indv,,3]
      
      gam_result_1 <- RunGam(timestamps01, x1, "binomial", basis_size, method)
      p1 <- gam_result_1$prob
      p1_linpred <- gam_result_1$linpred
      
      gam_result_2 <- RunGam(timestamps01, x2, "binomial", basis_size, method)
      p2 <- gam_result_2$prob
      p2_linpred <- gam_result_2$linpred
      
      gam_result_3 <- RunGam(timestamps01, x3, "binomial", basis_size, method)
      p3 <- gam_result_3$prob
      p3_linpred <- gam_result_3$linpred
      
      # estimate the latent tranjecotries Z
      denominator_p <- 1 + exp(p3_linpred)
      z1 <- (p1_linpred-p3_linpred)-log( (1+exp(p1_linpred))/(denominator_p))
      z2 <- (p2_linpred-p3_linpred)-log( (1+exp(p2_linpred))/(denominator_p))
      
      psum <- p1 + p2 + p3
      return(c(c(z1,z2), cbind(p1/psum, p2/psum, p3/psum)))
    }
  # Unravel the two variables from zp
  z_rows_count <- timeseries_length * 2
  Z <- array(zp[1:z_rows_count, ], c(z_rows_count, num_indv))
  p <- array(t(matrix(zp[(z_rows_count + 1):dim(zp)[1], ], ncol=num_indv)), c(num_indv, timeseries_length, category_count))
  
  return(list(Z1_est=Z[1:timeseries_length,],
              Z2_est=Z[1:timeseries_length+timeseries_length,],
              p1_est=t(p[,,1]),
              p2_est=t(p[,,2]),
              p3_est=t(p[,,3]) ))
}



GenerateCategFuncData <- function(prob_curves)
{
  curve_count <- length(prob_curves);
  
  # we could have just passed these arguments ???
  num_indvs <- ncol(prob_curves$p1)
  timeseries_length <- nrow(prob_curves$p1)
  
  # better names for W and X ???
  W <- matrix(0, ncol=num_indvs, nrow=timeseries_length)
  X_array <- array(0, c(num_indvs, timeseries_length, curve_count))
  
  for(indv in c(1:num_indvs))
  {
    X <- sapply(c(1:timeseries_length),
                function(this_time) rmultinom(n=1,
                                              size=1,
                                              prob = c(prob_curves$p1[this_time,indv],
                                                       prob_curves$p2[this_time,indv],
                                                       prob_curves$p3[this_time,indv]) ))
    W[,indv] <- apply(X, 2, which.max)
    X_array[indv,,] <- t(X)
  }
  
  return(list(X=X_array, W=W)) # X_binary W_catfd
}

#' Get clustered data
#'
#'
GenerateClusterData <- function(setting, scenario, k, num_indvs, timeseries_length)
{
  setting_object <- GetMuAndScore(setting, scenario, k)
  cluster_f <- GenerateClusterDataScenario(num_indvs,
                                           timeseries_length,
                                           k,
                                           mu_1 = setting_object$mu_1,
                                           mu_2 = setting_object$mu_2,
                                           score_vals = setting_object$score_vals)
  return (cluster_f)
}

#' Get fraction of occurrence of each class for a given scenario
#' @param scenario scenario name as a string "A", "B", "C"
#' @return a vector containing the fractions
#'
GetOccurrenceFractions <- function(scenario)
{
  occur_fraction <- switch (scenario,
                            "A" = c(0.75, 0.22, 0.03),
                            "B" = c(0.5, 0.3, 0.2),
                            "C" = c(0.1, 0.6, 0.3)
  )
  
  return (occur_fraction)
}

#' Get mu_1, mu_2 functions, and score_vals objects for a given context.
#' @param setting setting identified as an integer 1,2,3
#' @param scenario scenario name as a string "A", "B", "C"
#' @param k number of points along the score decay axis
#' @return A list that contains mu_1, mu_2, score_vals
#'
# setting <- 2
# scenario <- "B"
# k <- 2
# t <- timestamps01

#abc <- GetMuAndScore(setting, scenario, k)
GetMuAndScore <- function(setting, scenario, k)
{
  all_score_values = rep(0, k)
  
  if(1 == setting)
  {
    mu_1 <- function(t) -1 + 2 * t + 2 * t^2
    
    mu_2 <- switch(scenario,
                   "A" = function(t) -2.5 + exp(t * 2),
                   "B" = function(t) -0.5 + exp(t * 2),
                   "C" = function(t) -2.5 + exp(t * 2)
    )
    score_front <- switch(scenario,
                          "A" = c(1, 1/2, 1/4),
                          "B" = c(1, 1/2, 1/4),
                          "C" = c(50, 25, 5)
    )
  } else if(2 == setting)
  {
    mu_1 <- function(t) 4 * t^2 - 1.2
    
    mu_2 <- function(t) 4 * t^2 - 3.5
    
    score_front <- c(1, 1/2, 1/4)
  } else if(3 == setting)
  {
    mu_1 <- function(t) -2.2 + 4 * t^2
    
    mu_2 <- function(t) -7 + 6 * t^2
    
    score_front <- c(1, 1/4, 1/16)
  }
  
  for(idx in 1:length(score_front))
  {
    all_score_values[idx] <- score_front[idx]
  }
  
  return(list("mu_1" = mu_1, "mu_2" = mu_2, "score_vals" = all_score_values))
}


#' Generate cluster data for a given scenario
#' @param num_indvs number of individuals
#' @param timeseries_length length of time-series as an integer
#' @param k  number of eigen(psi) functions
#' @param mu_1 mean function for the first latent curve
#' @param mu_2 mean function for the second latent curve
#' @param score_vals the variance of the principal component scores
#'
GenerateClusterDataScenario <- function(num_indvs,
                                        timeseries_length,
                                        k = 3,
                                        mu_1,
                                        mu_2,
                                        score_vals)
{
  timestamps01 <- seq(from = 0.0001, to = 1, length=timeseries_length)
  
  # noise octaves
  # cat("octave", num_indvs, k, num_indvs * k, "\n")
  scores_standard <- matrix(rnorm(num_indvs * k), ncol = k)
  scores <- scores_standard %*% diag(sqrt(score_vals))
  
  #
  BIG_mu <- c(mu_1(timestamps01), mu_2(timestamps01))
  BIG_phi <- PsiFunc(k, timestamps01)
  
  Z <- BIG_phi %*% t(scores) + BIG_mu
  Z1 <- Z[1:timeseries_length, ]
  Z2 <- Z[1:timeseries_length + timeseries_length, ]
  expZ1 <- exp(Z1)
  expZ2 <- exp(Z2)
  denom <- 1 + expZ1 + expZ2
  p1 <- expZ1 / denom
  p2 <- expZ2 / denom
  p3 <- 1 / denom
  
  # vectorize for future work!!!
  return(list(Z1 = Z1, Z2 = Z2,
              p1 = p1, p2 = p2, p3 = p3,
              MEAN = BIG_mu, PHI = BIG_phi, MFPC = scores))
}

#' Psi function
#'
PsiFunc <- function(klen, timestamps01)
{
  psi_k1 <- sapply(c(1:klen), function(i) sin((2 * i + 1) * pi * timestamps01))
  psi_k2 <- sapply(c(1:klen), function(i) cos(2 * i * pi * timestamps01))
  return(rbind(psi_k1, psi_k2))
}

RunExperiment <- function(scenario, num_replicas, est_choice, some_identifier="noid")
{
  temp_folder <- file.path("outputs", "clustersims", paste(scenario, "_", num_replicas, "_", est_choice, "_", some_identifier, sep=""))
  # Empty the directory if it exists
  if(dir.exists(temp_folder)){
    unlink(temp_folder, recursive = TRUE)
  }
  dir.create(temp_folder)
  print(temp_folder)
  
  
  n100t300C <- ClusterSimulation(100,300,scenario,num_replicas,est_choice,TRUE,temp_folder)
  n100t750C <- ClusterSimulation(100,750,scenario,num_replicas,est_choice,TRUE,temp_folder)
  n100t2000C <- ClusterSimulation(100,2000,scenario,num_replicas,est_choice,TRUE,temp_folder)
  
  
  n500t300C <- ClusterSimulation(500,300,scenario,num_replicas,est_choice,TRUE,temp_folder)
  n500t750C <- ClusterSimulation(500,750,scenario,num_replicas,est_choice,TRUE,temp_folder)
  n500t2000C <- ClusterSimulation(500,2000,scenario,num_replicas,est_choice,TRUE,temp_folder)
  
  
  n1000t300C <- ClusterSimulation(1000,300,scenario,num_replicas,est_choice,TRUE,temp_folder)
  n1000t750C <- ClusterSimulation(1000,750,scenario,num_replicas,est_choice,TRUE,temp_folder)
  n1000t2000C <- ClusterSimulation(1000,2000,scenario,num_replicas,est_choice,TRUE,temp_folder)
  
  
  true_tableC <- rbind(n100t300C$cluster_table_true,n100t750C$cluster_table_true,n100t2000C$cluster_table_true,
                       n500t300C$cluster_table_true,n500t750C$cluster_table_true,n500t2000C$cluster_table_true,
                       n1000t300C$cluster_table_true,n1000t750C$cluster_table_true,n1000t2000C$cluster_table_true)
  rownames(true_tableC) <- c("n100t300","n100t750","n100t2000",
                             "n500t300","n500t750","n500t2000",
                             "n1000t300","n1000t750","n1000t2000")
  
  est_tableC <- rbind(n100t300C$cluster_table_est,n100t750C$cluster_table_est,n100t2000C$cluster_table_est,
                      n500t300C$cluster_table_est,n500t750C$cluster_table_est,n500t2000C$cluster_table_est,
                      n1000t300C$cluster_table_est,n1000t750C$cluster_table_est,n1000t2000C$cluster_table_est)
  rownames(est_tableC) <- c("n100t300","n100t750","n100t2000",
                            "n500t300","n500t750","n500t2000",
                            "n1000t300","n1000t750","n1000t2000")
  
  
  est_tableC_se <- rbind(n100t300C$cluster_table_est_se,n100t750C$cluster_table_est_se,n100t2000C$cluster_table_est_se,
                         n500t300C$cluster_table_est_se,n500t750C$cluster_table_est_se,n500t2000C$cluster_table_est_se,
                         n1000t300C$cluster_table_est_se,n1000t750C$cluster_table_est_se,n1000t2000C$cluster_table_est_se)
  rownames(est_tableC_se) <- c("n100t300","n100t750","n100t2000",
                               "n500t300","n500t750","n500t2000",
                               "n1000t300","n1000t750","n1000t2000")
  
  mse_tableC100 <- rbind(
    c(n100t300C$mse[1,1],n100t300C$mse[2,1],n100t300C$hellinger[1,1],n100t300C$hellinger[2,1],n100t300C$hellinger[3,1]),
    c(n100t750C$mse[1,1],n100t750C$mse[2,1],n100t750C$hellinger[1,1],n100t750C$hellinger[2,1],n100t750C$hellinger[3,1]),
    c(n100t2000C$mse[1,1],n100t2000C$mse[2,1],n100t2000C$hellinger[1,1],n100t2000C$hellinger[2,1],n100t2000C$hellinger[3,1]),
    
    c(n100t300C$mse[1,2],n100t300C$mse[2,2],n100t300C$hellinger[1,2],n100t300C$hellinger[2,2],n100t300C$hellinger[3,2]),
    c(n100t750C$mse[1,2],n100t750C$mse[2,2],n100t750C$hellinger[1,2],n100t750C$hellinger[2,2],n100t750C$hellinger[3,2]),
    c(n100t2000C$mse[1,2],n100t2000C$mse[2,2],n100t2000C$hellinger[1,2],n100t2000C$hellinger[2,2],n100t2000C$hellinger[3,2]),
    
    c(n100t300C$mse[1,3],n100t300C$mse[2,3],n100t300C$hellinger[1,3],n100t300C$hellinger[2,3],n100t300C$hellinger[3,3]),
    c(n100t750C$mse[1,3],n100t750C$mse[2,3],n100t750C$hellinger[1,3],n100t750C$hellinger[2,3],n100t750C$hellinger[3,3]),
    c(n100t2000C$mse[1,3],n100t2000C$mse[2,3],n100t2000C$hellinger[1,3],n100t2000C$hellinger[2,3],n100t2000C$hellinger[3,3])
  )
  
  rownames(mse_tableC100) <- c("s1n100t300","s1n100t750","s1n100t2000",
                               "s2n100t300","s2n100t750","s2n100t2000",
                               "s3n100t300","s3n100t750","s3n100t2000")
  colnames(mse_tableC100) <- c("z1","z2","p1","p2","p3")
  
  mse_tableC500 <- rbind(
    c(n500t300C$mse[1,1],n500t300C$mse[2,1],n500t300C$hellinger[1,1],n500t300C$hellinger[2,1],n500t300C$hellinger[3,1]),
    c(n500t750C$mse[1,1],n500t750C$mse[2,1],n500t750C$hellinger[1,1],n500t750C$hellinger[2,1],n500t750C$hellinger[3,1]),
    c(n500t2000C$mse[1,1],n500t2000C$mse[2,1],n500t2000C$hellinger[1,1],n500t2000C$hellinger[2,1],n500t2000C$hellinger[3,1]),
    
    c(n500t300C$mse[1,2],n500t300C$mse[2,2],n500t300C$hellinger[1,2],n500t300C$hellinger[2,2],n500t300C$hellinger[3,2]),
    c(n500t750C$mse[1,2],n500t750C$mse[2,2],n500t750C$hellinger[1,2],n500t750C$hellinger[2,2],n500t750C$hellinger[3,2]),
    c(n500t2000C$mse[1,2],n500t2000C$mse[2,2],n500t2000C$hellinger[1,2],n500t2000C$hellinger[2,2],n500t2000C$hellinger[3,2]),
    
    c(n500t300C$mse[1,3],n500t300C$mse[2,3],n500t300C$hellinger[1,3],n500t300C$hellinger[2,3],n500t300C$hellinger[3,3]),
    c(n500t750C$mse[1,3],n500t750C$mse[2,3],n500t750C$hellinger[1,3],n500t750C$hellinger[2,3],n500t750C$hellinger[3,3]),
    c(n500t2000C$mse[1,3],n500t2000C$mse[2,3],n500t2000C$hellinger[1,3],n500t2000C$hellinger[2,3],n500t2000C$hellinger[3,3])
  )
  
  
  rownames(mse_tableC500) <- c("s1n500t300","s1n500t750","s1n500t2000",
                               "s2n500t300","s2n500t750","s2n500t2000",
                               "s3n500t300","s3n500t750","s3n500t2000")
  colnames(mse_tableC500) <- c("z1","z2","p1","p2","p3")
  
  
  mse_tableC1000 <- rbind(
    c(n1000t300C$mse[1,1],n1000t300C$mse[2,1],n1000t300C$hellinger[1,1],n1000t300C$hellinger[2,1],n1000t300C$hellinger[3,1]),
    c(n1000t750C$mse[1,1],n1000t750C$mse[2,1],n1000t750C$hellinger[1,1],n1000t750C$hellinger[2,1],n1000t750C$hellinger[3,1]),
    c(n1000t2000C$mse[1,1],n1000t2000C$mse[2,1],n1000t2000C$hellinger[1,1],n1000t2000C$hellinger[2,1],n1000t2000C$hellinger[3,1]),
    
    c(n1000t300C$mse[1,2],n1000t300C$mse[2,2],n1000t300C$hellinger[1,2],n1000t300C$hellinger[2,2],n1000t300C$hellinger[3,2]),
    c(n1000t750C$mse[1,2],n1000t750C$mse[2,2],n1000t750C$hellinger[1,2],n1000t750C$hellinger[2,2],n1000t750C$hellinger[3,2]),
    c(n1000t2000C$mse[1,2],n1000t2000C$mse[2,2],n1000t2000C$hellinger[1,2],n1000t2000C$hellinger[2,2],n1000t2000C$hellinger[3,2]),
    
    c(n1000t300C$mse[1,3],n1000t300C$mse[2,3],n1000t300C$hellinger[1,3],n1000t300C$hellinger[2,3],n1000t300C$hellinger[3,3]),
    c(n1000t750C$mse[1,3],n1000t750C$mse[2,3],n1000t750C$hellinger[1,3],n1000t750C$hellinger[2,3],n1000t750C$hellinger[3,3]),
    c(n1000t2000C$mse[1,3],n1000t2000C$mse[2,3],n1000t2000C$hellinger[1,3],n1000t2000C$hellinger[2,3],n1000t2000C$hellinger[3,3])
  )
  
  
  rownames(mse_tableC1000) <- c("s1n1000t300","s1n1000t750","s1n1000t2000",
                                "s2n1000t300","s2n1000t750","s2n1000t2000",
                                "s3n1000t300","s3n1000t750","s3n1000t2000")
  colnames(mse_tableC1000) <- c("z1","z2","p1","p2","p3")
  
  
  
  #save(true_tableC,est_tableC,est_tableC_se,mse_tableC100,mse_tableC500,mse_tableC1000,
  #file=paste(scenario,num_replicas,est_choice,"clustering.RData",sep="_"))
  
  true_est_w_data <- (list("true_tableC"=true_tableC,
                           "est_tableC"=est_tableC,
                           "est_tableC_se"=est_tableC_se,
                           "mse_tableC100"=mse_tableC100,
                           "mse_tableC500"=mse_tableC500,
                           "mse_tableC1000"=mse_tableC1000,
                           
                           "Z_true_curves"=list(n100t300C$Z_true_curves,
                                                n100t750C$Z_true_curves,
                                                n100t2000C$Z_true_curves,
                                                n500t300C$Z_true_curves,
                                                n500t750C$Z_true_curves,
                                                n500t2000C$Z_true_curves,
                                                n1000t300C$Z_true_curves,
                                                n1000t750C$Z_true_curves,
                                                n1000t2000C$Z_true_curves),
                           "Z_est_curves"=list(n100t300C$Z_est_curves,
                                               n100t750C$Z_est_curves,
                                               n100t2000C$Z_est_curves,
                                               n500t300C$Z_est_curves,
                                               n500t750C$Z_est_curves,
                                               n500t2000C$Z_est_curves,
                                               n1000t300C$Z_est_curves,
                                               n1000t750C$Z_est_curves,
                                               n1000t2000C$Z_est_curves),
                           "p_true_curves"=list(n100t300C$p_true_curves,
                                                n100t750C$p_true_curves,
                                                n100t2000C$p_true_curves,
                                                n500t300C$p_true_curves,
                                                n500t750C$p_true_curves,
                                                n500t2000C$p_true_curves,
                                                n1000t300C$p_true_curves,
                                                n1000t750C$p_true_curves,
                                                n1000t2000C$p_true_curves),
                           "p_est_curves"=list(n100t300C$p_est_curves,
                                               n100t750C$p_est_curves,
                                               n100t2000C$p_est_curves,
                                               n500t300C$p_est_curves,
                                               n500t750C$p_est_curves,
                                               n500t2000C$p_est_curves,
                                               n1000t300C$p_est_curves,
                                               n1000t750C$p_est_curves,
                                               n1000t2000C$p_est_curves),
                           "W_cfd"=list(n100t300C$W_cfd,
                                        n100t750C$W_cfd,
                                        n100t2000C$W_cfd,
                                        n500t300C$W_cfd,
                                        n500t750C$W_cfd,
                                        n500t2000C$W_cfd,
                                        n1000t300C$W_cfd,
                                        n1000t750C$W_cfd,
                                        n1000t2000C$W_cfd)
  ))
  
  save(true_est_w_data,file=file.path("outputs", paste(scenario,num_replicas,est_choice,some_identifier,"true_est_w_data_clustering.RData",sep="_")))
  
  return(true_est_w_data)
}



set.seed(123)
A_100_binomial <- RunExperiment("A",100,"binomial","paper1")



if(run_parallel)
{
  parallel::stopCluster(cl = my.cluster)
  initialized_parallel <- FALSE
}


#############################################################################################################
#############################################################################################################
#Part II: Functions for estimating the weekend effect
EstimateCategFuncData_multinormial_weekend_parallel <- function(timestamps01, W, basis_size=25, method="ML")
{
  
  
  num_indv<- nrow(W)
  timeseries_length <-ncol(W)
  category_count <- length(unique(c(W)))
  
  weekend_vector <- as.factor(c(rep(c(rep(0,5*24*60/20),rep(1,2*24*60/20)),4)))[1:timeseries_length]
  Z<-NULL
 
  
  
  Z_P_WeekendCoef <- foreach(i = 1:num_indv , .packages = c("mgcv")) %dorng%
    
    #T_rep <- foreach(this_row = 1:5) %dorng%
    { #source("./source_code/R/data_generator.R")
      
      Z_P_WeekendCoef <- numeric(timeseries_length*(category_count+category_count-1)+category_count-1)
      
      print(i)
      if ( length(table(W[i,])) == category_count) {
        
        fit_binom<-gam(list(W[i,]-1~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector,
                            ~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector
        ),
        family=multinom(K=category_count-1), method = method,
        control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
        optimizer=c("outer","bfgs")) 
        #####
        #print("Dino")
        ####to find design matrix
        g_design <- predict(fit_binom,type = "lpmatrix")
        g_mul <- g_design[,c(1,category_count:basis_size)]
        coef_fit <- fit_binom$coefficients[c(1,category_count:basis_size)]
        #extract z
        z1 <- g_mul %*% as.matrix(coef_fit,ncol=1)
        #####
        g_mul_2 <- g_design[,c(1,category_count:basis_size)+basis_size]
        coef_fit_2 <- fit_binom$coefficients[c(1,category_count:basis_size)+basis_size]
        z2 <- g_mul_2 %*% as.matrix(coef_fit_2,ncol=1)
        
        weekend_vector_coef <- fit_binom$coefficients[c(category_count-1,basis_size+category_count-1)]
        pp <- predict(fit_binom,type="response")
        p1 <- pp[,1]
        p2 <- pp[,2]
        p3 <- pp[,3]
      } else {
        if (names(table(W[i,]))[2]=="3"){
          W[i,][W[i,]==3] <- 2
          basis_size_rev <- max(min(round(min(unname(table(W[i,])[2]), sum(1-unname(table(W[i,])[2])))/2), basis_size ), 5)
          fit_binom <- gam(W[i,]-1~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                           family = "binomial", method = method,
                           control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                           optimizer=c("outer","bfgs"))
          # print("GODAK ADERI")
          ######################
          ####to find design matrix
          g_design <- predict(fit_binom,type = "lpmatrix")
          g_mul <- g_design[,c(1,category_count:basis_size_rev)]
          coef_fit <- fit_binom$coefficients[c(1,category_count:basis_size_rev)]
          #extract z
          z2 <- g_mul %*% as.matrix(coef_fit,ncol=1)
          #####
          z1 <- rep(0,timeseries_length)
          
          weekend_vector_coef <- c(0, fit_binom$coefficients[category_count-1])
          p3 <- predict(fit_binom,type="response")
          p1 <- 1-p3
          p2 <- rep(0,timeseries_length)
          ##########################
        }else {
          basis_size_rev <- max(min(round(min(unname(table(W[i,])[2]), sum(1-unname(table(W[i,])[2])))/2), basis_size ), 5)
          fit_binom <- gam(W[i,]-1~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                           family = "binomial", method = method,
                           control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                           optimizer=c("outer","bfgs"))
          #print("MATA OYA NATHWA PALUI")
          ######################
          ####to find design matrix
          g_design <- predict(fit_binom,type = "lpmatrix")
          g_mul <- g_design[,c(1,category_count:basis_size_rev)]
          coef_fit <- fit_binom$coefficients[c(1,category_count:basis_size_rev)]
          #extract z
          z1 <- g_mul %*% as.matrix(coef_fit,ncol=1)
          z2 <- rep(0,timeseries_length)
          weekend_vector_coef <- c(fit_binom$coefficients[category_count-1],0)
          ##########################
          ##########################
          p2 <- predict(fit_binom,type="response")
          p1 <- 1-p2
          p3 <- rep(0,timeseries_length)
        }
      } 
    
      Z_P_WeekendCoef <- c(z1,z2,p1,p2,p3, weekend_vector_coef)
      
      return(Z_P_WeekendCoef)
      
    }
  
  Z_P_WeekendCoef <- do.call(rbind, Z_P_WeekendCoef)
  #t*n
  return(list(Z1_est=t(Z_P_WeekendCoef[,1:timeseries_length]), 
              Z2_est=t(Z_P_WeekendCoef[,(1+timeseries_length):(2*timeseries_length)]),
              p1_est=t(Z_P_WeekendCoef[,(1+timeseries_length*2):(3*timeseries_length)]), 
              p2_est=t(Z_P_WeekendCoef[,(1+timeseries_length*3):(4*timeseries_length)]), 
              p3_est=t(Z_P_WeekendCoef[,(1+timeseries_length*4):(5*timeseries_length)]) ,
              weekend_vector_coef = Z_P_WeekendCoef[,(1+timeseries_length*(category_count+category_count-1)):dim(Z_P_WeekendCoef)[2]]))
}

# library("pracma")

#Function to find the L2 distance between two latent curves
#' @param yy 1D vector (true curve)
#' @param yy2 1D vector (estimated curve)
#' @return scalar-L2 distance
trapzfnum <- function(yy,yy2,timestamps01)
{
  st=timestamps01[1]
  et=tail(timestamps01,n=1)
  x=seq(st,et,length=5000)
  xx=seq(st,et,length=length(yy))
  y1 <- cubicspline(xx, yy,x)
  y2 <- cubicspline(xx, yy2,x)
  out=sqrt(trapz(x, (y1-y2)^2) )
  return(out)
}

#'Function to find the Hellinger distance between two probability curves
#' @param yy 1D vector (true curve)
#' @param yy2 1D vector (estimated curve)
#' @param timestamps01 1D vector 
#' @return scalar-Hellinger distance
trapzfnump <- function(yy,yy2,timestamps01)
{
  st=timestamps01[1]
  et=tail(timestamps01,n=1)
  x=seq(st,et,length=5000)
  xx=seq(st,et,length=length(yy))
  y1 <- cubicspline(xx, yy,x)
  y1[y1<0]=0
  y2 <- cubicspline(xx, yy2,x)
  y2[y2<0]=0
  out=sqrt(trapz(x, (sqrt(y1)-sqrt(y2))^2) )
  return(out)
}


time_elapsed <<- list()
last_time <<- 0
row_name <<- NULL


#' Function to start recording the time for one task
#' @param task_name : task name that needs to track the time
timeKeeperStart <- function(task_name)
{
  row_name <<- task_name
  if(FALSE == row_name %in% names(time_elapsed))
  {
    time_elapsed[[row_name]] <<- NULL
  }
  last_time <<- Sys.time()
}


#' Function to print the time taken by the task 
timeKeeperNext <- function()
{
  this_time <- Sys.time()
  this_section_time <- this_time - last_time
  cat("\n--------------------\n",
      row_name, "\n\ttook:", capture.output(this_section_time), 
      "\n====================\n")
  time_elapsed[[row_name]] <<- append(time_elapsed[[row_name]], this_section_time)
  last_time <<- this_time
}


#############################################################################################################
#############################################################################################################



#############################################################################################################
#############################################################################################################
#Part IV: Twitter Data Simulation
########################################################################################################
ClusterSimulation <- function(num_indvs, timeseries_length,
                              scenario, num_replicas, est_choice, run_hellinger, temp_folder)
{
  cat("Cluster Simulation\nNum Indvs:\t", num_indvs,
      "\nTimeseries Len:\t", timeseries_length,
      "\nScenario:\t", scenario,
      "\nNum Replicas:\t", num_replicas)
  
  occur_fraction <- GetOccurrenceFractions(scenario)
  cat("\nOccurrence Fractions: ", occur_fraction)
  
  cluster_allocation <- occur_fraction * num_indvs
  cat("\nIndividuals in each cluster (cluster alloc.): ", cluster_allocation)
  
  true_cluster <- rep(c(1,2,3), cluster_allocation)
  cat("\nTrue Cluster:\n", true_cluster)
  
  true_cluster_db <- rep(c(1,2,0), cluster_allocation)
  cat("\nTrue Cluster DB:\n", true_cluster_db)
  
  rmse<- array(0, c(num_replicas, 2, 3))       #2 components and 3 settings
  hellinger <- array(0, c(num_replicas, 3, 3)) #3 events probab and 3 settings
  true_kmeans <- est_kmeans  <- NULL # records clustering membership on the TRUE SCORES/ESTIMATED SCORES
  true_fadp <- est_fadp <-  NULL # records clustering membership on the TRUE Z/ESTIMATED Z
  true_dbscan <- est_dbscan <-  NULL
  
  time_elapsed <<- list()
  # "X "=NULL, "univfpca"=NULL, "kmeans"=NULL, "fadp"=NULL, "dbscan"=NULL, "cfd"=NULL)
  last_time <- 0
  row_name <- NULL
  timeKeeperStart <- function(rn)
  {
    row_name <<- rn
    if(FALSE == row_name %in% names(time_elapsed))
    {
      time_elapsed[[row_name]] <<- NULL
    }
    last_time <<- Sys.time()
  }
  timeKeeperNext <- function()
  {
    this_time <- Sys.time()
    this_section_time <- this_time - last_time
    cat(row_name, "calc time taken:", capture.output(this_section_time), "\n")
    time_elapsed[[row_name]] <<- append(time_elapsed[[row_name]], this_section_time)
    last_time <<- this_time
  }
  
  ######
  Z_true_curve <- Z_est_curve <- list()
  p_true_curve <- p_est_curve <- list()
  W_cfd <- list()
  click_evals <- NULL
  total_regens <- 0
  for(replica_idx in 1:num_replicas)
  {
    cat("\nNum Indvs:", num_indvs,
        "\tTimeseries Len:", timeseries_length,
        "\tScenario:", scenario,
        "\tNum Replicas:", num_replicas)
    cat("\nReplica: ", replica_idx, "\n")
    
    # generate clusters
    # set.seed(seed_cluster + 100 * replica_idx)
    # cat("\nCluster", replica_idx, " --> seed: ", seed_cluster + 100 * replica_idx, "\n")
    cat("Cluster", replica_idx, "\n")
    
    cluster_f1 <- GenerateClusterData(1, scenario, 3, cluster_allocation[1], timeseries_length)
    cluster_f2 <- GenerateClusterData(2, scenario, 3, cluster_allocation[2], timeseries_length)
    cluster_f3 <- GenerateClusterData(3, scenario, 3, cluster_allocation[3], timeseries_length)
    
    # Recover the latent Gaussian process --> is this always 2 ???
    Z1 <- cbind(cluster_f1$Z1, cluster_f2$Z1, cluster_f3$Z1)
    Z2 <- cbind(cluster_f1$Z2, cluster_f2$Z2, cluster_f3$Z2)
    
    
    # Recover the true probability curves --> could there be more than 3 ???
    p1 <- cbind(cluster_f1$p1, cluster_f2$p1, cluster_f3$p1)
    p2 <- cbind(cluster_f1$p2, cluster_f2$p2, cluster_f3$p2)
    p3 <- cbind(cluster_f1$p3, cluster_f2$p3, cluster_f3$p3)
    prob_curves <- list(p1=p1, p2=p2, p3=p3)
    #########9/11/2023  
    Z_true_curve[[replica_idx]]=array(c(Z1,Z2),dim=c(timeseries_length,num_indvs,2))
    p_true_curve[[replica_idx]]=array(c(p1,p2,p3),dim=c(timeseries_length,num_indvs,3))
    ############
    
    
    # generate categFuncData
    # set.seed(seed_cfd + 100 * replica_idx)
    # cat("\nCategFD", replica_idx, " --> seed: ", seed_cfd + 100 * replica_idx, "\n")
    cat("CategFD", replica_idx, "\n")
    
    categ_func_data_list <- GenerateCategFuncData(prob_curves)
    
    # what is Q vals ? better name???
    Q_vals <- unique(c(categ_func_data_list$W))
    if(is.numeric(Q_vals))
    {
      Q_vals <- sort(Q_vals)
    }
    
    # I need to know what this loop is meant to do !??? maybe there is a better way
    for(indv in 1:num_indvs)
    {
      if(indv %in% 1:cluster_allocation[1])
      {
        setting_choice <- 1
      }
      if(indv %in% (cluster_allocation[1] + 1):(cluster_allocation[1] + cluster_allocation[2]))
      {
        setting_choice <- 2
      }
      if(indv %in% (cluster_allocation[1] + cluster_allocation[2] + 1):num_indvs)
      {
        setting_choice <- 3
      }
      
      # better names for the following variables ???
      # 1. check weather one category only appears 1 time and is it in the end of the timeseries
      # 2. OR is it appearing only one time in the begining
      # 3. OR if the category is less thant Q total category
      # In general, one category only occurs 2 times
      # If timepoints=300, one category only occurs less than 4 times 3/300=0.01
      #If timepoints=750, one category only occurs less than 6 times 5/750=0.0067
      tolcat <- table(categ_func_data_list$W[,indv])
      catorder <- order(tolcat, decreasing = TRUE)
      numcat <- length(catorder)
      refcat <- catorder[numcat]
      count_iter <- 0
      while (count_iter < 100 && 
             ( (numcat < length(Q_vals))
               ||(timeseries_length==300  && min(as.numeric(tolcat)) < 4)
               ||(timeseries_length==750  && min(as.numeric(tolcat)) < 10)
             )
      )
      {
        count_iter <- count_iter + 1
        
        new_cluster_data <- GenerateClusterData(setting_choice, scenario, 3, 5, timeseries_length)
        
        new_prob_curves <- list(p1 = new_cluster_data$p1, p2 = new_cluster_data$p2, p3 = new_cluster_data$p3)
        new_categ_func_data_list <- GenerateCategFuncData( new_prob_curves )
        
        # what is this 3 ?? arbitrarily chosen?
        categ_func_data_list$W[, indv] <- new_categ_func_data_list$W[, 3]
        Z1[, indv] <- new_cluster_data$Z1[, 3] # latent curves Z1 and Z2
        categ_func_data_list$X[indv, , ] <- 0
        Z2[, indv] <- new_cluster_data$Z2[, 3]
        
        for (this_time in 1:timeseries_length)
        {
          categ_func_data_list$X[indv, this_time, which(Q_vals == categ_func_data_list$W[, indv][this_time])] <- 1
        }
        
        tolcat <- table(categ_func_data_list$W[, indv])
        catorder <- order(tolcat, decreasing = TRUE)
        numcat <- length(catorder)
        refcat <- catorder[numcat]
      } # end while
      total_regens <- total_regens + count_iter
    } # end for(indv in 1:num_indvs)
    
    # browser()
    click_based_clusts <- ComputeClickClustering(t(categ_func_data_list$W))
    
    # browser()
    cs_ri=rand.index(true_cluster,click_based_clusts$cs)
    cs_ari=adj.rand.index(true_cluster, click_based_clusts$cs)
    
    cc_ri=rand.index(true_cluster,click_based_clusts$cc)
    cc_ari=adj.rand.index(true_cluster, click_based_clusts$cc)
    
    cluster_check <- list(cs_ri=cs_ri, cs_ari=cs_ari, cc_ri=cc_ri, cc_ari=cc_ari)
    click_evals <- rbind(click_evals, cluster_check)
    
    cat("Done replica:", replica_idx, "\n")
  }# END of "for(replica_idx in 1:num_replicas)'
  
  cat("\n replicas done \n")
  
  return(click_evals)
}


ComputeClickClustering <- function(W_mat) 
{
  
  
  # browser()
  
  xw <- apply(W_mat, 1, paste, collapse = ",")
  
  clickstream_data <- as.clickstreams(xw, header = FALSE)
  
  # cluster
  min_totss <- Inf
  best_clustering <- NULL
  for (cs_kval in 2:5) {
    cc_clusters <- clusterClickstreams(clickstream_data, order = 0, centers = cs_kval)
    print(cc_clusters$totss)
    if(cc_clusters$totss <= min_totss){
      min_totss <- cc_clusters$totss
      best_clustering <- cc_clusters
    }
  }
  # clusters <- clusterClickstreams(clickstream_data, order = 0, centers = 3)
  
  print(paste("Num clusters:", length(best_clustering$clusters)))
  
  W_lables <- rep(0,dim(W_mat)[1])
  #for (cltr in best_clustering$clusters) {
  for (cltr in 1:length(best_clustering$clusters) ){
    cat("Cluster:" , cltr,"\n")
    #print(paste(names(cltr)))
    cluster_1 <- as.numeric(noquote(names(best_clustering$clusters[[cltr]])))
    W_lables[cluster_1] <- cltr
  }
  print(table(W_lables))
  
  # convert clickstream data to click cluster compatible data
  clickclust_data <- as.ClickClust(clickstream_data)
  
  # run EM on the data
  best_emclust <- NULL
  min_bic <- Inf
  for (kval in 2:5) {
    emclust <- click.EM(clickclust_data$X, K = kval)
    print(emclust$BIC)
    if(!is.na(emclust$BIC) && emclust$BIC <= min_bic){
      best_emclust <- emclust
      min_bic <- emclust$BIC
    }
  }
  
  if(is.null(best_emclust)){
    # browser()
    best_emclust <- click.EM(clickclust_data$X, K = 2)
    min_bic <- best_emclust$BIC
  }
  
  W_label_clickclust <- best_emclust$id
  
  
  return(list(cs=W_lables, cs_totss=min_totss, cc=W_label_clickclust, cc_BIC=min_bic))
}

#############################################################################################################
#############################################################################################################
#Part IV: Twitter data pre processing
########################################################################################################
#step one: sample code to collect Twitter user data directly from Twitter API
#######
#Intro to TwitteR
######


#Need to put in your credentials from the developers account
api_key <- ""
api_secret <- ""
access_token <- ""
access_token_secret <- ""

setup_twitter_oauth(api_key, api_secret, access_token, access_token_secret)


#Get current rate limits
getCurRateLimitInfo()

#Change the key word or number of tweets
# can filter by language or other features
cur_tweets = searchTwitter("anti-vax", n=25)
#convert data to dataframe
cur_tweets_df = twListToDF(cur_tweets)

#Use the 1st account
cur_username = cur_tweets_df$screenName[[1]]

#get user information
cur_user  = getUser(cur_username)

#get the timelines
cur_user_tweets <- userTimeline(cur_username, includeRts=TRUE, n=3200)
cur_user_tweets=twListToDF(cur_user_tweets)

#get location of user
cur_user$location

#######################################################################################################
#######
#Timeline collection
#Input: List of usernames and reference tweets (date time of ref tweet)
#output:  1) List of Tweets from users whose last 3200 tweets contain the reference tweet
######
#Change to your APIs
api_key <- ""
api_secret <- ""
access_token <- ""
access_token_secret <- ""

setup_twitter_oauth(api_key, api_secret, access_token, access_token_secret)

#Current Rate Limit is 900 every 15 minutes

getCurRateLimitInfo()
twt_lim=900 #number of times able to call userTimeline
#time_lim=15*60 #wait time in sec

#All_users is the vector containing the users of interest
#CHANGE to users that you want to collect
#e.g. c("aweishampel", "billrand")
all_users = c()

#variables that you want to keep look
keeps <- c("screenName","text", "created", "id")

for(i in 1:length(all_users)){
  
  cur_user=all_users[i]
  print(paste(i," ", cur_user, " ", Sys.time()))
  skip=FALSE
  
  #Only collect users whose data exist
  #(prevents errors)
  cur_user_tweets <- try(userTimeline(cur_user, includeRts=TRUE, n=3200))
  
  if(!grepl("Error", cur_user_tweets[1]) & length(cur_user_tweets)!=0){
    cur_user_tweets=twListToDF(cur_user_tweets)
  }else{
    skip=TRUE
  }
  #the furtherst I can go is 3200 tweets
  if(!skip){
    
    cur_user_tweets = cur_user_tweets[ , keeps, drop=FALSE]
    
    #Is reference tweet included in current users pulled tweets
    #cur_user_ref_tweet=subset(ref_tweet, ref_tweet$users==cur_user)
    
    #ref_tweet_included=min(cur_user_tweets$created)<min(cur_user_ref_tweet$ref_time)
    
    if(dim(cur_user_tweets)[1]==0){
      temp_hist_data[i, 2]=NA
    }else{
      temp_hist_data[i, 2]=max(cur_user_tweets$created)
    }
    
    if(i==1){
      users_tweets=cur_user_tweets
    }else{
      users_tweets=rbind.data.frame(users_tweets, cur_user_tweets)
    }
  }
  
  if((i %% (twt_lim-15)) == 0){
    sleeptime=as.numeric(difftime(getCurRateLimitInfo()[41,4], Sys.time(), unit="secs"))
    print(paste( "We reached our limit. Resting for ", sleeptime, " seconds", sep=""))
    Sys.sleep(sleeptime)
  }
  
  
}
########################################################################################################
#step two: processing raw Twitter user data to the categorical functional data format
#Due to the amount of the data and the size of the excel file, the following code are written in Python
#Please download all the original csv file in order to proceed to step two
#Please either uncomment and copy code to python or download the python file (Fortmat Data.py) to run this step

#######################################################################################################
#Python code to process the raw Twitter data collected through Twitter API
#Please download all the data from the foler: TwitterRawData and make sure the path is correct
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import os
import re
import glob
import os

#%%

absolute_path = os.path.dirname(__file__)
#%%
profile_df = pd.read_csv(absolute_path + '\\twitter_data\\Users_profile_data.csv')
profile_df = profile_df.dropna()

#%%
file = absolute_path + '\\twitter_data\\Tweets_of_4776users_CSVs\\'
csv_files = glob.glob(os.path.join(file, "*.csv"))
df_tweets = pd.DataFrame()

for f in csv_files:
  temp = pd.read_csv(f, encoding='latin-1')
df_tweets = df_tweets.append(temp, ignore_index=True)

df_tweets = df_tweets.add_suffix('_tweeted_data')
print(df_tweets)
#%%
df_tweets['UserID_tweeted_data'] = df_tweets['UserID_tweeted_data'].str.strip('@')
#%%
def parser(fn):
  txt = open(fn).read()
preparse = re.findall('"\d+",\d+,".+?","[^\t\n\r\f\v]+?","[^\t\n\r\f\v]+?"', txt, re.DOTALL)
parsedRows = []
for line in preparse:
  columns = line.split(',')
output = {}

output['Index'] = columns.pop(0)
output['Tweet_ID'] = columns.pop(0)
output['UserID'] = columns.pop()
output['Company'] = columns.pop()
output['Text'] = ','.join(columns)

parsedRows.append(output)

return parsedRows

#parsed = [t.split(',') for t in preparse]
#print(parsed)

#%%

data = parser(absolute_path + r'\twitter_data\EDITED - Reference_tweet_data.txt')
reference_df = pd.DataFrame(data)

reference_df['Index'] = reference_df['Index'].str.strip('"')
reference_df['Index'] = reference_df['Index'].astype(int)
reference_df['UserID'] = reference_df['UserID'].str.strip('"')
reference_df =reference_df.add_suffix('_ref_tweets')
print ("\nUnique values :  \n",reference_df.nunique())

duplicated = reference_df[reference_df.duplicated(subset='UserID_ref_tweets', keep=False)]
reference_without_duplicates = reference_df.drop_duplicates(
  subset = ['UserID_ref_tweets', 'Company_ref_tweets'],
  keep = 'last').reset_index(drop=True)
#%%
df = pd.merge(reference_without_duplicates, profile_df, left_on='UserID_ref_tweets', right_on='user_id')
print ("\nFeatures : \n" ,df.columns.tolist())
print ("\nMissing values :  ", df.isnull().sum().values.sum())
print ("\nUnique values :  \n",df.nunique())

#%%
finalDF = pd.merge(df_tweets,df, left_on='UserID_tweeted_data', right_on='UserID_ref_tweets')
print ("\nFeatures : \n" ,finalDF.columns.tolist())
print ("\nMissing values :  ", finalDF.isnull().sum().values.sum())
print ("\nUnique values :  \n",finalDF.nunique())

#%%
finalDF['text_tweeted_data'] = finalDF['text_tweeted_data'].astype(str)
finalDF['Company_ref_tweets'] = finalDF['Company_ref_tweets'].astype(str)
finalDF['Company_ref_tweets'] = finalDF['Company_ref_tweets'].str.strip('"')
finalDF['mention'] = finalDF.apply(lambda x: x.Company_ref_tweets in x.text_tweeted_data, axis=1)
print(finalDF['mention'].value_counts())

#%%
finalDF['mention_type'] = 2
finalDF.loc[finalDF['mention']==False,'mention_type'] = 1
print(finalDF['mention_type'].value_counts())

#%%
finalDF['DateTime2_tweeted_data'] = pd.to_datetime(finalDF['DateTime2_tweeted_data'])
finalDF['DateTime_tweeted_data'] = finalDF.apply(lambda x: pd.Timestamp(x.DateTime2_tweeted_data), axis=1)
finalDF['DateTime_tweeted_data'] = finalDF.apply(lambda x: x.DateTime_tweeted_data.timestamp(), axis=1)
#%%
"""WRITE FORMATTED DATA TO CSV"""
finalDF.to_csv(absolute_path + r'\final_output.csv')

#%%
"""IMPORT FORMATTED DATA"""
finalDF = pd.read_csv(r'C:\Users\Rob\OneDrive\NCSU PhD\twitter\final_output.csv')

print(finalDF.columns)

#%%
"""REMOVE EXTRA COLUMNS"""
final_column_list = ['DateTime2_tweeted_data', 'Index_ref_tweets', 'mention_type']
finalDF = finalDF[finalDF.columns.intersection(final_column_list)]
finalDF['mention_type'] = finalDF['mention_type'].astype(int)

#%%
"""lAST MONTH OF ACTIVITY"""

finalDF['DateTime2_tweeted_data'] = pd.to_datetime(finalDF['DateTime2_tweeted_data'])
month_ago = finalDF['DateTime2_tweeted_data'].max() - datetime.timedelta(days=30)

last_30_days_DF = finalDF.query('DateTime2_tweeted_data > @month_ago')
print(last_30_days_DF.head())

#%%

def round_time(dt: datetime, unit=20):
  #seconds = dt - dt.date()
  #unit_seconds = unit.total_seconds()
  #rounded_seconds = seconds - (seconds % unit_seconds)
  #return dt.date() + rounded_seconds
  return dt.replace(second=0, microsecond=0, minute=dt.minute-(dt.minute%unit), hour=dt.hour)

"""Bin Each User"""
each_user = last_30_days_DF.groupby('Index_ref_tweets')
output = pd.DataFrame()
counter = 0
start_time=round_time(dt=last_30_days_DF['DateTime2_tweeted_data'].min())
end_time=round_time(dt=last_30_days_DF['DateTime2_tweeted_data'].max())
print(start_time)
print(end_time)

t_index = pd.DatetimeIndex(pd.date_range(start=start_time, end=end_time, freq="20min"))

for user_name, user_data in each_user:
  print(f"Read {counter} out of {len(each_user)}")
#df_bin=df_subset.set_index('DateTime2_tweeted_data',drop=False)
binuser=user_data.groupby(pd.Grouper(key='DateTime2_tweeted_data', freq='20min')).max()
binuser = binuser.reindex(t_index)
binuser['mention_type'] = binuser['mention_type'].fillna(5)
binuser['Index_ref_tweets'] = binuser['Index_ref_tweets'].fillna(user_name)
binuser = binuser.transpose()
print(binuser)
output = output.append(binuser, ignore_index=True)
counter += 1

#%%
"""SPLIT INTO THREE DATAFRAMES"""
"""Mentions"""
mentions = output.replace(to_replace=[1, 5, 2], value=[0, 0, 1])
mentions.to_csv(absolute_path + r'\mentions_dataset.csv')

no_mentions = output.replace(to_replace=[1, 5, 2], value=[1, 0, 0])
no_mentions.to_csv(absolute_path + r'\no_mentions_dataset.csv')

no_tweets = output.replace(to_replace=[1, 5, 2], value=[0, 1, 0])
no_tweets.to_csv(absolute_path + r'\no_tweet_dataset.csv')
Output the three csv files: One file for each category
#######################################################################################################
#R code to pre process three csv files in the form of the categorical functional data
#Please download all the files from processedTwitterData and make sure the path is correct
X_i1r=read.csv(file='no_tweet_dataset.csv')[,-1]   #3875, 2162  to 1000,2161  1st col originally the id X
X_i2r=read.csv(file='no_mentions_dataset.csv')[,-1]  
X_i3r=read.csv(file='mentions_dataset.csv')[,-1]
timet=colnames(X_i1r)
timetnew=gsub("X","",timet)
timetnewn=gsub('\\.', '-', timetnew)
time1=stringi::stri_replace_last_regex(timetnewn,'-',':')
time2=stringi::stri_replace_last_regex(time1,'-',':')
time3=stringi::stri_replace_last_regex(time2,'-',' ')
startendt=time3[c(1,length(time3))] # "2017-02-14 03:00:00" "2017-03-16 03:00:00"
time4=as.numeric(as.POSIXct(strptime(time3,"%Y-%m-%d %H:%M:%S")))
which(is.na(time4)) #1870 1871 1872
time4[c(1870,1871,1872)]=as.numeric(as.POSIXct(c("2017-03-12 02:00:00", "2017-03-12 02:20:00" ,"2017-03-12 02:40:00"),"%Y-%m-%d %H:%M:%S",tz="EST"))
time5=time4-min(time4)
timefinal=time5/max(time5)
t=timefinal[101:2100]
st=t[1]
et=tail(t,n=1)
datapoints=length(t)
X_i1=X_i1r[101:2100]
X_i2=X_i2r[101:2100]
X_i3=X_i3r[101:2100]
colnames(X_i1)=t
colnames(X_i2)=t
colnames(X_i3)=t

save(X_i1,file="X_i1.RData")
save(X_i2,file="X_i2.RData")
save(X_i3,file="X_i3.RData")

#########################################################################################################
#########################################################################################################
#Part V: Using proposed method to cluster Twitter user data
#Please make sure the previous Rdata file is created and they are under the same file path
#apply to Twitter

load("X_i1.RData")
load("X_i2.RData")
load("X_i3.RData")

datapoints=dim(X_i1)[2]
n=dim(X_i1)[1]
st=as.numeric(colnames(X_i1))[1]
et=tail(as.numeric(colnames(X_i1)),n=1)
t=as.numeric(colnames(X_i1))

startt=Sys.time()
#step1
#get_estimate_latent=function(X_i1,X_i2,X_i3,ps,k,t)
trueest=get_estimate_latent(as.matrix(X_i1),as.matrix(X_i2),as.matrix(X_i3),1,3,t)
#MFPCA

#ramsay MFPCA method
startt=Sys.time()
mfpcaram=mfpca_ram(trueest,2,2000,q=3,1000,0.1,0.9,seq(0.1,0.9,length=2000))
save(trueest,mfpcaram,file="Twitterram.RData")
##cumsum(mfpcaram$eigenvalue)/sum(mfpcaram$eigenvalue)
#[1] 0.9668647 0.9832566 0.9940806
endt=Sys.time()
time_elapse=endt-startt
print(time_elapse)


scaleeps=3
combz3=mfpcaram$scores
dimz2=dim(combz3)[2]
if (dimz2<=2){
  
  minPts2=4
  
}

if (dimz2>2){
  minPts2=dimz2+1
}
dist2=kNNdist(combz3, k = minPts2-1)
#sort distance and use elbow to find optimal epsilon
distdataelbow2=data.frame(sort(dist2))
distdataelbow2$index=1:(dim(combz3)[1])
ipoint2 <- elbow(data = distdataelbow2)
epsoptimal2=ipoint2$sort.dist2._selected*scaleeps
#dbscan results
res2 <- dbscan(combz3, eps =epsoptimal2 , minPts = 25)
#visualize
tdata2=data.frame(combz3)
tdata2$Cluster=as.factor(res2$cluster)
tpd2 <- ggplot(tdata2,aes(X1,X2,colour = Cluster)) + geom_point(aes(shape=Cluster),size=3)+ggtitle(paste0("Users Cluster Results",'\n',"(",dim(tdata2)[1]," Subjects",")")) +
  xlab("MFPCA Score1") + ylab("MFPCA Score2")+ theme(plot.title = element_text(hjust = 0.5))+
  
  theme(text=element_text(size = 20))
tpd2

#plot latent curve by cluster
vec0=c(which(tdata2$Cluster==0))
vec1=c(which(tdata2$Cluster==1))
vec2=c(which(tdata2$Cluster==2))

zlatent=array(c(trueest$Est$EstimateZ_i1,trueest$Est$EstimateZ_i2),dim=c(1000,2000,2))
phat=array(c(trueest$Est$Estimatep_i1,trueest$Est$Estimatep_i2,trueest$Est$Estimatep_i3),
           dim=c(1000,2000,3))


plot_latentfn(zlatent,phat,
              0,1,tclusterdata$Cluster,0,seq(0.1,0.9,length=2000),tclusterdata)

plot_latentfn(zlatent,phat,
              0,1,tclusterdata$Cluster,1,seq(0.1,0.9,length=2000),tclusterdata)
plot_latentfn(zlatent,phat,
              0,1,tclusterdata$Cluster,2,seq(0.1,0.9,length=2000),tclusterdata)





########################################################################################################
#Function to plot Twitter latent curves by cluster
#Input zlatent: 3 d array n by datapoints by Q-1
#phat: 3 d array n by datapoints by Q
#cluster: clustering results
#labelnum: the latent curves of the cluster label
#argval: time interval
#tclusterdata: the clustering results
#Output: Latent curves by cluster
plot_latentfn=function(zlatent,phat,st,et,cluster,labelnum,argval,tclusterdata){
  vectort=c(which(tclusterdata$Cluster==labelnum))
  datapoints=dim(zlatent)[2]
  
  n=length(vectort)
  ########################################################################################
  
  #plot 1: truth and smoothed
  ########################################################################################
  #entire block is  n*t dimension
  #truth Z, P and smoothed Z, P
  numz=dim(zlatent)[3]
  nump=dim(phat)[3]
  zest=zlatent[vectort,,]
  pest=phat[vectort,,]
  
  meanz=colMeans(zest)
  for (i in 1:numz){
    matplot(argval, t(zest[,,i]),
            type='l', lty=1, col="light grey",
            
            main=mtext(bquote("Estimated Latent Cruve ")),
            xlab="Day", ylab="Value",cex.lab = 2,cex.axis = 2,cex.main=2,xaxt="n")
    lines(argval,meanz[,i] ,
          type='l', lty=2, lwd=2, col = "red")
    axis(1,                         # Define x-axis manually
         at = c(0,476/2000,952/2000,1428/2000,1904/2000),cex.lab = 2,cex.axis = 2,cex.main=2,
         labels = c("0","5","10","15","20")) }
  #p_i
  
  meanp=colMeans(pest)   #t*nump
  
  for (i in 1:nump){
    matplot(argval, t(pest[,,i]),
            type='l', lty=1, col="light grey",
            main=mtext(bquote("Estimated Probability Cruve ")),
            xlab="Day", ylab='Value',cex.lab = 2,cex.axis = 2,cex.main=2,xaxt="n")
    
    lines(argval,meanp[,i] ,
          type='l', lty=2, lwd=2, col = "red")
    axis(1,
         at = c(0,476/2000,952/2000,1428/2000,1904/2000),cex.lab = 2,cex.axis = 2,cex.main=2,
         labels = c("0","5","10","15","20"))
    
  }
  
}


#Function to get estimated latent curves from Twitter Data
#Input: functional dummy curves X_i1, X_i2 and X_i3 are n by datapoints matrix with binary values
#ps=1
#k: the number of the eigen functions
#t: time interval
get_estimate_latent=function(X_i1,X_i2,X_i3,ps,k,t){
  
  k=k
  
  #k=3  #number of eigen functions
  q=3  #level of the categorical level
  
  #recover Z_i1 hat using Z_i[1,all j, all n] only related to p1
  Z_i1hat=Z_ihat(X_i1,t)
  #recover Z_i2 hat using X_i[2,all j, all n] only related to p2 
  Z_i2hat=Z_ihat(X_i2,t)
  Z_i3hat=Z_ihat(X_i3,t)
  
  
  if(ps==1){ 
    
    Z_i1hatstar=Z_i1hat+log(1+exp(Z_i3hat))-log(1+exp(Z_i1hat))-Z_i3hat
    Z_i2hatstar=Z_i2hat+log(1+exp(Z_i3hat))-log(1+exp(Z_i2hat))-Z_i3hat
  }
  
  
  
  p_i1hat=p_ihat(Z_i1hatstar,Z_i2hatstar)$p_1hatmatrix
  p_i2hat=p_ihat(Z_i1hatstar,Z_i2hatstar)$p_2hatmatrix
  p_i3hat=p_ihat(Z_i1hatstar,Z_i2hatstar)$p_3hatmatrix
  
  
  # truel=list("Truecatcurve"=X_nt)
  est=list("EstimateZ_i1"=Z_i1hatstar,"EstimateZ_i2"=Z_i2hatstar,"Estimatep_i1"=p_i1hat,"Estimatep_i2"=p_i2hat,"Estimatep_i3"=p_i3hat)
  # return(list("Trueth"=truel,"Est"=est))
  return(list("Est"=est))
}


#Function to find MFPCA scores
#Input: datagenerated is the data generated using the function generate_data_scenario
#M is the initial starting dimension of MFPCA scores
#q=3: number of the total category
#n: number of the subject
#datapoints: the number of the data points per subject
#st: starting time
#et: end time
#t: time interval for the observations

#Output:MFPCA scores, Latent curves, eigen functions, eigen values and mean functions
mfpca_ram=function(datagenerated,M,datapoints,q=3,n,st,et,t){
  #gather information
  t=t
  st=st
  et=et
  Z_i1hatstar=datagenerated$Est$EstimateZ_i1
  Z_i2hatstar=datagenerated$Est$EstimateZ_i2
  n=dim(Z_i1hatstar)[1]
  Zihatstar=array(c(Z_i1hatstar, Z_i2hatstar), c(n,datapoints,2))
  #form functional data object
  fdarange = c(st, et)
  fdabasis = fda::create.bspline.basis(fdarange, 25, 4)
  fdatime = seq(st, et, length = datapoints)
  Zihatstarnew = apply(Zihatstar, c(1, 3), t)
  fdafd = fda::smooth.basis(fdatime, Zihatstarnew, fdabasis)$fd
  nharm = M
  #use ramsay MFPCA 
  fdapcaList = fda::pca.fd(fdafd, nharm)
  #find the dimension when the pve is at least 0.95
  finalM=which(cumsum(fdapcaList$values)/sum(fdapcaList$values)>=0.95)[1]
  fdapcaListfinal = fda::pca.fd(fdafd, finalM)
  scores_z = apply(fdapcaListfinal$scores, c(1, 2), sum)
  values=fdapcaListfinal$values
  meanfn=fdapcaList$meanfd
  eigenfn=fdapcaList$harmonics
  
  #if one dimension explains more than 0.95 pve, collect first two dimensions for visualization purpose
  dims=dim(scores_z)[2]
  if (dims<2){
    finalM=2
    fdapcaListfinal = fda::pca.fd(fdafd, finalM)
    scores_z = apply(fdapcaListfinal$scores, c(1, 2), sum)
    values=fdapcaListfinal$values
    meanfn=fdapcaList$meanfd
    eigenfn=fdapcaList$harmonics
  }
  
  return(list("scores"=scores_z,
              "latentcurve"=Zihatstar,"eigenvalue"=values,"meanfn"=meanfn,
              "eigenfn"=eigenfn))}

#Function to output the smooth probability curves using latent curves Z

#New function to output predicted L-1 latent curves:
#input observed latent curves Z_it and return smoothed p_ihat
#Output: smoothed probability curves
p_ihat=function(Z_i1app,Z_i2app){
  denom=(1+exp(Z_i1app)+exp(Z_i2app))
  p_i1h=exp(Z_i1app)/denom
  p_i2h=exp(Z_i2app)/denom
  p_i3h=1/denom
  return(list("p_1hatmatrix"=p_i1h,"p_2hatmatrix"=p_i2h,"p_3hatmatrix"=p_i3h))
}