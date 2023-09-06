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

#Function to find the L2 distance between two latent curves
#Input: yy- 1D vector (true curve)
#       yy2-1D vector (estimated curve)
#Output: scalar-L2 distance
trapzfnum <- function(yy,yy2)
{
  st=0.0001
  et=1
  x=seq(st,et,length=5000)
  xx=seq(st,et,length=length(yy))
  y1 <- cubicspline(xx, yy,x)
  y2 <- cubicspline(xx, yy2,x)
  out=sqrt(trapz(x, (y1-y2)^2) )
  return(out)
}

#Function to use trapzfnum function and find L2 distance for 2D array, n of them
mse_bw_matrix <- function(truecurve,estcurve)
{
  trapzfnum <- function(yy,yy2)
  {
    st=0.0001
    et=1
    x=seq(st,et,length=5000)
    xx=seq(st,et,length=length(yy))
    y1 <- cubicspline(xx, yy,x)
    y2 <- cubicspline(xx, yy2,x)
    out=sqrt(trapz(x, (y1-y2)^2) )
    return(out)
  }
  
  n=dim(truecurve)[2]
  # datapoints=dim(truecurve)[1]
  # mseall=c(0)
  ######could probably use apply function here it's also subject level
  mseall <- foreach(i = 1:n, .combine = c, .packages = c("pracma")) %dorng% {
    return(rbind(trapzfnum(truecurve[,i], estcurve[,i])))
  }
  
  return(mseall)
}


#Function to find the Hellinger distance between two probability curves
#Input: yy- 1D vector (true curve)
#       yy2-1D vector (estimated curve)
#Output: scalar-Hellinger distance
trapzfnump <- function(yy,yy2)
{
  st=0.0001
  et=1
  x=seq(st,et,length=5000)
  xx=seq(st,et,length=length(yy))
  y1 <- cubicspline(xx, yy,x)
  y1[y1<0]=0
  y2 <- cubicspline(xx, yy2,x)
  y2[y2<0]=0
  out=sqrt(trapz(x, (sqrt(y1)-sqrt(y2))^2) )
  return(out)
}

#Function to use trapzfnump function and find Hellinger distance for 2D array, n of them
mse_bw_matrixp <- function(truecurve,estcurve)
{
  trapzfnump <- function(yy,yy2)
  {
    st=0.0001
    et=1
    x=seq(st,et,length=5000)
    xx=seq(st,et,length=length(yy))
    y1 <- cubicspline(xx, yy,x)
    y1[y1<0]=0
    y2 <- cubicspline(xx, yy2,x)
    y2[y2<0]=0
    out=sqrt(trapz(x, (sqrt(y1)-sqrt(y2))^2) )
    return(out)
  }
  
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
    return(rbind(trapzfnump(truecurve[,i], estcurve[,i])/sqrt_2))
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
                              scenario, num_replicas, est_choice, run_hellinger)
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
    time_elapsed[[row_name]] <<- append(time_elapsed[[row_name]], this_section_time)
    last_time <<- this_time
  }

  for(replica_idx in 1:num_replicas)
  {
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
      tolcat <- table(categ_func_data_list$W[, indv])
      catorder <- order(tolcat, decreasing = TRUE)
      numcat <- length(catorder)
      refcat <- catorder[numcat]
      count_iter <- 0
      while ((min(as.numeric(tolcat)) == 1 && categ_func_data_list$W[, indv][timeseries_length] == refcat && count_iter < 100)
             || (min(as.numeric(tolcat)) == 1 && categ_func_data_list$W[, indv][1] == refcat && count_iter < 100))
      {
        count_iter <- count_iter + 1
        #new_cluster_data <- GenerateClusterData(setting_choice, "A", 3, 5, timeseries_length)

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
    } # end for(indv in 1:num_indvs)

    #Estimation
    cat("EstimateCategFuncData", replica_idx, "\n")
    timestamps01 <- seq(from = 0.0001, to = 1, length=timeseries_length)
    timeKeeperStart("Xiaoxia")
    categFD_est <- EstimateCategFuncData(est_choice, timestamps01, categ_func_data_list$W)
    timeKeeperNext()

    if (run_hellinger)
    {
      # evaluate performance Z and P
      rmse1_temp <- c(by(mse_bw_matrix(Z1,categFD_est$Z1_est) , true_cluster, mean))
      rmse2_temp <- c(by(mse_bw_matrix(Z2,categFD_est$Z2_est), true_cluster, mean))
      rmse[replica_idx, ,] <- rbind(rmse1_temp,rmse2_temp )

      error.p1<- mse_bw_matrixp(p1,categFD_est$p1_est)
      error.p2<- mse_bw_matrixp(p2,categFD_est$p2_est)
      error.p3<- mse_bw_matrixp(p3,categFD_est$p3_est)


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
  } # END of "for(replica_idx in 1:num_replicas)'

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
  
  save(time_elapsed, file=paste("time_elapsed_", num_indvs, "_", timeseries_length, "_",
                                scenario, "_", num_replicas, "_", est_choice, "_", run_hellinger, ".RData", sep=""))
  
  return_vals <- list("cluster_table_true"=cluster_table_true,"cluster_table_est"=cluster_table_est,
       "cluster_table_est_se"=cluster_table_est_se)
  
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
#' @param choice "probit", "binomial", or "
EstimateCategFuncData <- function(choice, timestamps01, W, basis_size=25, method="ML")
{
  if(choice == "probit"){
    X <- GetXFromW(W)
    return(EstimateCategFuncData_probit(timestamps01, X, basis_size=25, method="ML"))
  }else if(choice == "binormial"){
    X <- GetXFromW(W)
    return(EstimateCategFuncData_binorm(timestamps01, X, basis_size=25, method="ML"))
  }else if(choice == "multinormial"){
    return(EstimateCategFuncData_multinormial(timestamps01, W, basis_size=25, method="ML"))
  }
}

#'Function to estimate z and p using wood_multinormial
#'

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
    p1 <- fit_binom$fitted.values[,1]
    p2 <- fit_binom$fitted.values[,2]
    prob[i,,] <- cbind(p1, p2, 1-p1-p2)
    
    z1<- fit_binom$linear.predictors[,1]
    z2<- fit_binom$linear.predictors[,2]
    Z<- cbind(Z, c(z1,z2))
    
  }
  
  return(list(Z1_est=Z[1:timeseries_length ,], Z2_est=Z[1:timeseries_length +timeseries_length ,], 
              p1_est=t(prob[,,1]), p2_est=t(prob[,,2]), p3_est=t(prob[,,3]) ))
}




EstimateCategFuncData_probit <- function(timestamps01, X, basis_size=25, method="ML")
{
  num_indv<- dim(X)[1]
  timeseries_length<- dim(X)[2]
  category_count <- dim(X)[3]

  Z<-NULL
  p<-array(0, c(num_indv, timeseries_length, category_count))
  # for (indv in 1:num_indv){
  #   #i=93
  #   #basis_size=25
  #   x1<- X[indv,,1]
  #   x2<- X[indv,,2]
  #   x3<- X[indv,,3]
  # 
  #   if (timeseries_length<=301 && sum(x1)/timeseries_length<0.004){
  # 
  #     gam_result_1 <- RunGam(timestamps01, x1, "probit", basis_size, method)
  #     p1 <- gam_result_1$prob
  #     p1_linpred <- gam_result_1$linpred
  # 
  #     gam_result_2 <- RunGam(timestamps01, x2, "probit", basis_size, method)
  #     p2 <- gam_result_2$prob
  #     p2_linpred <- gam_result_2$linpred
  # 
  #     gam_result_3 <- RunGam(timestamps01, x3, "probit", basis_size, method)
  #     p3 <- gam_result_3$prob
  #     p3_linpred <- gam_result_3$linpred
  #     denominator_p <- 1+exp(p3_linpred)
  #     z1<- (p1_linpred-p3_linpred)-log( (1+exp(p1_linpred))/(denominator_p))
  #     z2<- (p2_linpred-p3_linpred)-log( (1+exp(p2_linpred))/(denominator_p))
  #   }else{
  #     gam_result_1 <- RunGam(timestamps01, x1, "binomial", basis_size, method)
  #     p1 <- gam_result_1$prob
  #     p1_linpred <- gam_result_1$linpred
  # 
  #     gam_result_2 <- RunGam(timestamps01, x2, "binomial", basis_size, method)
  #     p2 <- gam_result_2$prob
  #     p2_linpred <- gam_result_2$linpred
  # 
  #     gam_result_3 <- RunGam(timestamps01, x3, "binomial", basis_size, method)
  #     p3 <- gam_result_3$prob
  #     p3_linpred <- gam_result_3$linpred
  # 
  #     # estimate the latent curves Z
  #     exp_p3_linepred <- exp(p3_linpred)
  #     z1 <- (p1_linpred-p3_linpred)-log( (1+exp(p1_linpred))/(1+exp_p3_linepred))
  #     z2 <- (p2_linpred-p3_linpred)-log( (1+exp(p2_linpred))/(1+exp_p3_linepred))
  # 
  # 
  #   } # end if special case for probit
  # 
  #   Z <- cbind(Z, c(z1,z2))
  #   psum <- (p1+p2+p3)
  #   p[indv,,] <- cbind(p1/psum, p2/psum, p3/psum)
  # }
  # return(list(Z1_est=Z[1:timeseries_length,],
  #             Z2_est=Z[1:timeseries_length+timeseries_length,],
  #             p1_est=t(p[,,1]),
  #             p2_est=t(p[,,2]),
  #             p3_est=t(p[,,3]) ))
  
  zp <- foreach (indv = 1:num_indv, .combine = cbind, .init = NULL, .packages = c("mgcv")) %dorng%
  {
    source("R/acj/run_gam_function.R")
    
    x1<- X[indv,,1]
    x2<- X[indv,,2]
    x3<- X[indv,,3]
    
    if (timeseries_length<=301 && sum(x1)/timeseries_length<0.004){
      
      gam_result_1 <- RunGam(timestamps01, x1, "probit", basis_size, method)
      p1 <- gam_result_1$prob
      p1_linpred <- gam_result_1$linpred
      
      gam_result_2 <- RunGam(timestamps01, x2, "probit", basis_size, method)
      p2 <- gam_result_2$prob
      p2_linpred <- gam_result_2$linpred
      
      gam_result_3 <- RunGam(timestamps01, x3, "probit", basis_size, method)
      p3 <- gam_result_3$prob
      p3_linpred <- gam_result_3$linpred
      denominator_p <- 1+exp(p3_linpred)
      z1<- (p1_linpred-p3_linpred)-log( (1+exp(p1_linpred))/(denominator_p))
      z2<- (p2_linpred-p3_linpred)-log( (1+exp(p2_linpred))/(denominator_p))
    }else{
      gam_result_1 <- RunGam(timestamps01, x1, "binomial", basis_size, method)
      p1 <- gam_result_1$prob
      p1_linpred <- gam_result_1$linpred
      
      gam_result_2 <- RunGam(timestamps01, x2, "binomial", basis_size, method)
      p2 <- gam_result_2$prob
      p2_linpred <- gam_result_2$linpred
      
      gam_result_3 <- RunGam(timestamps01, x3, "binomial", basis_size, method)
      p3 <- gam_result_3$prob
      p3_linpred <- gam_result_3$linpred
      
      # estimate the latent curves Z
      denominator_p <- 1 + exp(p3_linpred)
      z1 <- (p1_linpred-p3_linpred)-log( (1+exp(p1_linpred))/(denominator_p))
      z2 <- (p2_linpred-p3_linpred)-log( (1+exp(p2_linpred))/(denominator_p))
      
      
    } # end if special case for probit
    
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



EstimateCategFuncData_binorm <- function(timestamps01, X, basis_size=25, method="ML")
{
  num_indv<- dim(X)[1]
  timeseries_length<- dim(X)[2]
  category_count <- dim(X)[3]

  Z<-NULL
  p<-array(0, c(num_indv, timeseries_length, category_count))
  ##########################
  ###########################
  # num_indv is the subject and this step is done by subject level
  # can parallel
  #############################
  #############################
#   for (indv in 1:num_indv)
#   {
#     x1<- X[indv,,1]
#     x2<- X[indv,,2]
#     x3<- X[indv,,3]
#
#     # fit the Binom model
#     ###################################################
#     ###################################################
#     ##updated estimation
#
#     gam_result_1 <- RunGam(timestamps01, x1, "binomial", basis_size, method)
#     p1 <- gam_result_1$prob
#     p1_linpred <- gam_result_1$linpred
#
#     gam_result_2 <- RunGam(timestamps01, x2, "binomial", basis_size, method)
#     p2 <- gam_result_2$prob
#     p2_linpred <- gam_result_2$linpred
#
#     gam_result_3 <- RunGam(timestamps01, x3, "binomial", basis_size, method)
#     p3 <- gam_result_3$prob
#     p3_linpred <- gam_result_3$linpred
#
#     # estimate the latent tranjecotries Z
#     exp_p3_linepred <- exp(p3_linpred)
#     z1 <- (p1_linpred-p3_linpred)-log( (1+exp(p1_linpred))/(1+exp_p3_linepred))
#     z2 <- (p2_linpred-p3_linpred)-log( (1+exp(p2_linpred))/(1+exp_p3_linepred))
#
#     Z <- cbind(Z, c(z1,z2))
#     psum <- (p1+p2+p3)
#     p[indv,,] <- cbind(p1/psum, p2/psum, p3/psum)
#   }
# return(p)
# return(list(Z=Z,p=p))


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

# EXECUTION:

set.seed(123)

# }) # profvis end

# ////////////// A Scenario BEGIN ///////////////////////////

n100t300A=ClusterSimulation(100,300,"A",100,"probit",TRUE)
n100t750A=ClusterSimulation(100,750,"A",100,"probit",TRUE)
n100t2000A=ClusterSimulation(100,2000,"A",100,"probit",TRUE)
n500t300A=ClusterSimulation(500,300,"A",100,"probit",TRUE)

n500t750A=ClusterSimulation(500,750,"A",100,"probit",TRUE)
n500t2000A=ClusterSimulation(500,2000,"A",100,"probit",TRUE)

n1000t300A=ClusterSimulation(1000,300,"A",100,"probit",TRUE)
n1000t750A=ClusterSimulation(1000,750,"A",100,"probit",TRUE)
n1000t2000A=ClusterSimulation(1000,2000,"A",100,"probit",TRUE)

true_tableA=rbind(n100t300A$cluster_table_true,n100t750A$cluster_table_true,n100t2000A$cluster_table_true,
                  n500t300A$cluster_table_true,n500t750A$cluster_table_true,n500t2000A$cluster_table_true,
                  n1000t300A$cluster_table_true,n1000t750A$cluster_table_true,n1000t2000A$cluster_table_true)
rownames(true_tableA)=c("n100t300","n100t750","n100t2000",
                        "n500t300","n500t750","n500t2000",
                        "n1000t300","n1000t750","n1000t2000")

est_tableA=rbind(n100t300A$cluster_table_est,n100t750A$cluster_table_est,n100t2000A$cluster_table_est,
                 n500t300A$cluster_table_est,n500t750A$cluster_table_est,n500t2000A$cluster_table_est,
                 n1000t300A$cluster_table_est,n1000t750A$cluster_table_est,n1000t2000A$cluster_table_est)
rownames(est_tableA)=c("n100t300","n100t750","n100t2000",
                       "n500t300","n500t750","n500t2000",
                       "n1000t300","n1000t750","n1000t2000")


est_tableA_se=rbind(n100t300A$cluster_table_est_se,n100t750A$cluster_table_est_se,n100t2000A$cluster_table_est_se,
                    n500t300A$cluster_table_est_se,n500t750A$cluster_table_est_se,n500t2000A$cluster_table_est_se,
                    n1000t300A$cluster_table_est_se,n1000t750A$cluster_table_est_se,n1000t2000A$cluster_table_est_se)
rownames(est_tableA_se)=c("n100t300","n100t750","n100t2000",
                          "n500t300","n500t750","n500t2000",
                          "n1000t300","n1000t750","n1000t2000")

save(true_tableA,est_tableA,est_tableA_se,file="A_clustering_before_mse_tableA.RData")

mse_tableA=rbind(
  c(n100t300A$mse[1,1],n100t300A$mse[2,1],n100t300A$hellinger[1,1],n100t300A$hellinger[2,1],n100t300A$hellinger[3,1]),
  c(n100t750A$mse[1,1],n100t750A$mse[2,1],n100t750A$hellinger[1,1],n100t750A$hellinger[2,1],n100t750A$hellinger[3,1]),
  c(n100t2000A$mse[1,1],n100t2000A$mse[2,1],n100t2000A$hellinger[1,1],n100t2000A$hellinger[2,1],n100t2000A$hellinger[3,1]),

  c(n100t300A$mse[1,2],n100t300A$mse[2,2],n100t300A$hellinger[1,2],n100t300A$hellinger[2,2],n100t300A$hellinger[3,2]),
  c(n100t750A$mse[1,2],n100t750A$mse[2,2],n100t750A$hellinger[1,2],n100t750A$hellinger[2,2],n100t750A$hellinger[3,2]),
  c(n100t2000A$mse[1,2],n100t2000A$mse[2,2],n100t2000A$hellinger[1,2],n100t2000A$hellinger[2,2],n100t2000A$hellinger[3,2]),

  c(n100t300A$mse[1,3],n100t300A$mse[2,3],n100t300A$hellinger[1,3],n100t300A$hellinger[2,3],n100t300A$hellinger[3,3]),
  c(n100t750A$mse[1,3],n100t750A$mse[2,3],n100t750A$hellinger[1,3],n100t750A$hellinger[2,3],n100t750A$hellinger[3,3]),
  c(n100t2000A$mse[1,3],n100t2000A$mse[2,3],n100t2000A$hellinger[1,3],n100t2000A$hellinger[2,3],n100t2000A$hellinger[3,3]),
)


rownames(mse_tableA)=c("s1n100t300","s1n100t750","s1n100t2000",
                       "s2n100t300","s2n100t750","s2n100t2000",
                       "s3n100t300","s3n100t750","s3n100t2000")
colnames(mse_tableA)=c("z1","z2","p1","p2","p3")
save(true_tableA,est_tableA,est_tableA_se,mse_tableA,file="A_clustering.RData")

# ////////////// A Scenario END ///////////////////////////
#scenario B
#####################################
###scenario B
##########scenarioA
#cluster_simulation=function(n,m,scenario,mc_sims)
n100t300B=ClusterSimulation(100,300,"B",100,"probit",TRUE)
n100t750B=ClusterSimulation(100,750,"B",100,"probit",TRUE)
n100t2000B=ClusterSimulation(100,2000,"B",100,"probit",TRUE)


n500t300B=ClusterSimulation(500,300,"B",100,"probit",TRUE)
n500t750B=ClusterSimulation(500,750,"B",100,"probit",TRUE)
n500t2000B=ClusterSimulation(500,2000,"B",100,"probit",TRUE)


n1000t300B=ClusterSimulation(1000,300,"B",100,"probit",TRUE)
n1000t750B=ClusterSimulation(1000,750,"B",100,"probit",TRUE)
n1000t2000B=ClusterSimulation(1000,2000,"B",100,"probit",TRUE)


true_tableB=rbind(n100t300B$cluster_table_true,n100t750B$cluster_table_true,n100t2000B$cluster_table_true,
                  n500t300B$cluster_table_true,n500t750B$cluster_table_true,n500t2000B$cluster_table_true,
                  n1000t300B$cluster_table_true,n1000t750B$cluster_table_true,n1000t2000B$cluster_table_true)
rownames(true_tableB)=c("n100t300","n100t750","n100t2000",
                        "n500t300","n500t750","n500t2000",
                        "n1000t300","n1000t750","n1000t2000")

est_tableB=rbind(n100t300B$cluster_table_est,n100t750B$cluster_table_est,n100t2000B$cluster_table_est,
                 n500t300B$cluster_table_est,n500t750B$cluster_table_est,n500t2000B$cluster_table_est,
                 n1000t300B$cluster_table_est,n1000t750B$cluster_table_est,n1000t2000B$cluster_table_est)
rownames(est_tableB)=c("n100t300","n100t750","n100t2000",
                       "n500t300","n500t750","n500t2000",
                       "n1000t300","n1000t750","n1000t2000")


est_tableB_se=rbind(n100t300B$cluster_table_est_se,n100t750B$cluster_table_est_se,n100t2000B$cluster_table_est_se,
                    n500t300B$cluster_table_est_se,n500t750B$cluster_table_est_se,n500t2000B$cluster_table_est_se,
                    n1000t300B$cluster_table_est_se,n1000t750B$cluster_table_est_se,n1000t2000B$cluster_table_est_se)
rownames(est_tableB_se)=c("n100t300","n100t750","n100t2000",
                          "n500t300","n500t750","n500t2000",
                          "n1000t300","n1000t750","n1000t2000")


save(true_tableB,est_tableB,est_tableB_se,file="B_clustering.RData")

mse_tableB=rbind(
  c(n100t300B$mse[1,1],n100t300B$mse[2,1],n100t300B$hellinger[1,1],n100t300B$hellinger[2,1],n100t300B$hellinger[3,1]),
  c(n100t750B$mse[1,1],n100t750B$mse[2,1],n100t750B$hellinger[1,1],n100t750B$hellinger[2,1],n100t750B$hellinger[3,1]),
  c(n100t2000B$mse[1,1],n100t2000B$mse[2,1],n100t2000B$hellinger[1,1],n100t2000B$hellinger[2,1],n100t2000B$hellinger[3,1]),

  c(n100t300B$mse[1,2],n100t300B$mse[2,2],n100t300B$hellinger[1,2],n100t300B$hellinger[2,2],n100t300B$hellinger[3,2]),
  c(n100t750B$mse[1,2],n100t750B$mse[2,2],n100t750B$hellinger[1,2],n100t750B$hellinger[2,2],n100t750B$hellinger[3,2]),
  c(n100t2000B$mse[1,2],n100t2000B$mse[2,2],n100t2000B$hellinger[1,2],n100t2000B$hellinger[2,2],n100t2000B$hellinger[3,2]),

  c(n100t300B$mse[1,3],n100t300B$mse[2,3],n100t300B$hellinger[1,3],n100t300B$hellinger[2,3],n100t300B$hellinger[3,3]),
  c(n100t750B$mse[1,3],n100t750B$mse[2,3],n100t750B$hellinger[1,3],n100t750B$hellinger[2,3],n100t750B$hellinger[3,3]),
  c(n100t2000B$mse[1,3],n100t2000B$mse[2,3],n100t2000B$hellinger[1,3],n100t2000B$hellinger[2,3],n100t2000B$hellinger[3,3]),
)


rownames(mse_tableB)=c("s1n100t300","s1n100t750","s1n100t2000",
                       "s2n100t300","s2n100t750","s2n100t2000",
                       "s3n100t300","s3n100t750","s3n100t2000")
colnames(mse_tableB)=c("z1","z2","p1","p2","p3")
save(true_tableB,est_tableB,est_tableB_se,mse_tableB,file="B_clustering.RData")

n100t300C=ClusterSimulation(100,300,"C",100,"probit",TRUE)
n100t750C=ClusterSimulation(100,750,"C",100,"probit",TRUE)
n100t2000C=ClusterSimulation(100,2000,"C",100,"probit",TRUE)


n500t300C=ClusterSimulation(500,300,"C",100,"probit",TRUE)
n500t750C=ClusterSimulation(500,750,"C",100,"probit",TRUE)
n500t2000C=ClusterSimulation(500,2000,"C",100,"probit",TRUE)


n1000t300C=ClusterSimulation(1000,300,"C",100,"probit",TRUE)
n1000t750C=ClusterSimulation(1000,750,"C",100,"probit",TRUE)
n1000t2000C=ClusterSimulation(1000,2000,"C",100,"probit",TRUE)


true_tableC=rbind(n100t300C$cluster_table_true,n100t750C$cluster_table_true,n100t2000C$cluster_table_true,
                  n500t300C$cluster_table_true,n500t750C$cluster_table_true,n500t2000C$cluster_table_true,
                  n1000t300C$cluster_table_true,n1000t750C$cluster_table_true,n1000t2000C$cluster_table_true)
rownames(true_tableC)=c("n100t300","n100t750","n100t2000",
                        "n500t300","n500t750","n500t2000",
                        "n1000t300","n1000t750","n1000t2000")

est_tableC=rbind(n100t300C$cluster_table_est,n100t750C$cluster_table_est,n100t2000C$cluster_table_est,
                 n500t300C$cluster_table_est,n500t750C$cluster_table_est,n500t2000C$cluster_table_est,
                 n1000t300C$cluster_table_est,n1000t750C$cluster_table_est,n1000t2000C$cluster_table_est)
rownames(est_tableC)=c("n100t300","n100t750","n100t2000",
                       "n500t300","n500t750","n500t2000",
                       "n1000t300","n1000t750","n1000t2000")


est_tableC_se=rbind(n100t300C$cluster_table_est_se,n100t750C$cluster_table_est_se,n100t2000C$cluster_table_est_se,
                    n500t300C$cluster_table_est_se,n500t750C$cluster_table_est_se,n500t2000C$cluster_table_est_se,
                    n1000t300C$cluster_table_est_se,n1000t750C$cluster_table_est_se,n1000t2000C$cluster_table_est_se)
rownames(est_tableC_se)=c("n100t300","n100t750","n100t2000",
                          "n500t300","n500t750","n500t2000",
                          "n1000t300","n1000t750","n1000t2000")


save(true_tableC,est_tableC,est_tableC_se,file="C_clustering.RData")

mse_tableC=rbind(
  c(n100t300C$mse[1,1],n100t300C$mse[2,1],n100t300C$hellinger[1,1],n100t300C$hellinger[2,1],n100t300C$hellinger[3,1]),
  c(n100t750C$mse[1,1],n100t750C$mse[2,1],n100t750C$hellinger[1,1],n100t750C$hellinger[2,1],n100t750C$hellinger[3,1]),
  c(n100t2000C$mse[1,1],n100t2000C$mse[2,1],n100t2000C$hellinger[1,1],n100t2000C$hellinger[2,1],n100t2000C$hellinger[3,1]),

  c(n100t300C$mse[1,2],n100t300C$mse[2,2],n100t300C$hellinger[1,2],n100t300C$hellinger[2,2],n100t300C$hellinger[3,2]),
  c(n100t750C$mse[1,2],n100t750C$mse[2,2],n100t750C$hellinger[1,2],n100t750C$hellinger[2,2],n100t750C$hellinger[3,2]),
  c(n100t2000C$mse[1,2],n100t2000C$mse[2,2],n100t2000C$hellinger[1,2],n100t2000C$hellinger[2,2],n100t2000C$hellinger[3,2]),

  c(n100t300C$mse[1,3],n100t300C$mse[2,3],n100t300C$hellinger[1,3],n100t300C$hellinger[2,3],n100t300C$hellinger[3,3]),
  c(n100t750C$mse[1,3],n100t750C$mse[2,3],n100t750C$hellinger[1,3],n100t750C$hellinger[2,3],n100t750C$hellinger[3,3]),
  c(n100t2000C$mse[1,3],n100t2000C$mse[2,3],n100t2000C$hellinger[1,3],n100t2000C$hellinger[2,3],n100t2000C$hellinger[3,3]),
)


rownames(mse_tableC)=c("s1n100t300","s1n100t750","s1n100t2000",
                       "s2n100t300","s2n100t750","s2n100t2000",
                       "s3n100t300","s3n100t750","s3n100t2000")
colnames(mse_tableC)=c("z1","z2","p1","p2","p3")
save(true_tableC,est_tableC,est_tableC_se,mse_tableC,file="B_clustering.RData")
end_time <- Sys.time()
print(end_time - start_time)

if(run_parallel)
{
  parallel::stopCluster(cl = my.cluster)
  initialized_parallel <- FALSE
}
