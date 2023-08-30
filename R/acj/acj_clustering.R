library(profvis)

library(mgcv)

library(foreach)
library(doParallel)
library(doRNG)

n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
cat("Parellel Registered: ", foreach::getDoParRegistered())

# set.seed(123)
# x <- foreach(indv = 1:10, .combine=rbind) %do% {
#   # list("n"=array(rnorm(6), c(2,3)), "a"=sqrt(indv))
#   # c(indv,sqrt(indv),indv*indv)
#   cbind(c(1:10) * indv * 1.0, array(30 + indv,10))
# }
# dim(x)
# typeof(x)
# x
# array(x)
#
#
# set.seed(123)
# zz <- NULL
# pp <- array(0, c(10, 20))
# for(indv in 1:10)
# {
#   pp[indv,] <- cbind(c(1:10) * indv, c(31:40))
#   # zz <- cbind(zz, c(1:10) * indv)
# }
# typeof(pp)
# pp
# array(pp)

# print(x)
#
# x1 <- as.vector(x[,1])
# x2 <- as.vector(x[,2])
# x3 <- as.vector(x[,3])
# print(x1)
# print(x2)
# print(x3)
#
# parallel::stopCluster(cl = my.cluster)

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
    print(cat("count_iter: ", count_iter))
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
  n=dim(truecurve)[2]
  datapoints=dim(truecurve)[1]
  mseall=c(0)
  ######could probably use apply function here it's also subject level
  for (i in 1:n)
  {
    mseall[i]=trapzfnum(truecurve[,i],estcurve[,i])
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
  n=dim(truecurve)[2]
  datapoints= dim(truecurve)[1]
  mseall=c(0)
  ######could probably use apply function here it's also subject level
  for (i in 1:n)
  {
    mseall[i]=trapzfnump(truecurve[,i],estcurve[,i])/sqrt(2)
  }
  return(mseall)
}

ClusterSimulation <- function(num_indvs, timeseries_length,
                              scenario, num_replicas,
                              seed_cluster=1230, seed_cfd=9876)
{
  cat("Cluster Simulation\nNum Indvs:\t", num_indvs,
      "\nTimeseries Len:\t", timeseries_length,
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

  for(replica_idx in 1:num_replicas)
  {
    cat("\nReplica: ", replica_idx)

    # generate clusters
    set.seed(seed_cluster + 100 * replica_idx)
    cat("\nCluster", replica_idx, " --> seed: ", seed_cluster + 100 * replica_idx, "\n")

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
    set.seed(seed_cfd + 100 * replica_idx)
    cat("\nCategFD", replica_idx, " --> seed: ", seed_cfd + 100 * replica_idx)

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
        new_cluster_data <- GenerateClusterData(setting_choice, "A", 3, 5, timeseries_length)

        new_prob_curves <- list(p1 = new_cluster_data$p1, p2 = new_cluster_data$p2, p3 = new_cluster_data$p3)
        new_categ_func_data_list <- GenerateCategFuncData( new_prob_curves )

        # what is this 3 ?? arbitrarily chosen?
        categ_func_data_list$W[, indv] <- new_categ_func_data_list$W[, 3]
        Z1[, indv] <- new_cluster_data$Z1[, 3] # latent curves Z1 and Z2
        Z2[, indv] <- new_cluster_data$Z2[, 3]

        categ_func_data_list$X[indv, , ] <- 0
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
    timestamps01 <- seq(from = 0.0001, to = 1, length=timeseries_length)
    categFD_est <- EstimateCategFuncData(timestamps01, X=categ_func_data_list$X)

    if (num_indvs == 100 && scenario == "A")
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
    mfpca_est <- extract_scores_UNIVFPCA(mZ1=categFD_est$Z1_est, mZ2=categFD_est$Z2_est, tt=timestamps01, PVE=0.95)

    #KMEANS
    true_kmeans_temp <- kmeans_cluster(data=mfpca_true$scores)$label
    est_kmeans_temp <- kmeans_cluster(data=mfpca_est$scores)$label

    #FADP
    true_fadp_temp <- fadp_cluster(mZ1=Z1, mZ2=Z2, tt=timestamps01)$label
    est_fadp_temp <- fadp_cluster(mZ1=categFD_est$Z1_est, mZ2=categFD_est$Z2_est, tt=timestamps01)$label

    #dbscan
    true_dbscan_temp <- dbscan_cluster(data=mfpca_true$scores,1)$label
    est_dbscan_temp <- dbscan_cluster(data=mfpca_est$scores,1)$label

    #record results
    true_dbscan<- cbind(true_dbscan, true_dbscan_temp)
    est_dbscan <- cbind(est_dbscan, est_dbscan_temp)

    true_kmeans<- cbind(true_kmeans, true_kmeans_temp)
    est_kmeans <- cbind(est_kmeans,  est_kmeans_temp)

    true_fadp<- cbind(true_fadp, true_fadp_temp)
    est_fadp <- cbind(est_fadp, est_fadp_temp)
    print(replica_idx)

  } # END of "for(replica_idx in 1:num_replicas)'

  cat("\n replicas done \n")

  if (num_indvs == 100 && scenario == "A")
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



  cluster_table_true=c(mean(results$true_fadp_ri),mean(results$true_fadp_ari), mean(results$true_fadp_cpn),
                       mean(results$true_kmeans_ri),mean(results$true_kmeans_ari),mean(results$true_kmeans_cpn),

                       mean(results$true_dbscan_ri),
                       mean(results$true_dbscan_ari),
                       mean(results$true_dbscan_cpn)
  )
  names(cluster_table_true)=c("fadp RI","fadp ARI","fadp cpn",
                              "kmeans RI","kmeans ARI","kmeans cpn",
                              "dbscan RI","dbscan ARI","dbscan cpn"
  )


  cluster_table_est=c(mean(results$est_fadp_ri),mean(results$est_fadp_ari),

                      mean(results$est_fadp_cpn),


                      mean(results$est_kmeans_ri),mean(results$est_kmeans_ari),
                      mean(results$est_kmeans_cpn),
                      mean(results$est_dbscan_ri),
                      mean(results$est_dbscan_ari),mean(results$est_dbscan_cpn)
  )
  names(cluster_table_est)=c("fadp RI","fadp ARI","fadp cpn",
                             "kmeans RI","kmeans ARI","kmeans cpn",
                             "dbscan RI","dbscan ARI","dbscan cpn"
  )


  cluster_table_est_se=c(sd(results$est_fadp_ri)/sqrt(num_replicas),sd(results$est_fadp_ari)/sqrt(num_replicas),

                         sd(results$est_fadp_cpn)/sqrt(num_replicas),


                         sd(results$est_kmeans_ri)/sqrt(num_replicas),sd(results$est_kmeans_ari)/sqrt(num_replicas),
                         sd(results$est_kmeans_cpn)/sqrt(num_replicas),
                         sd(results$est_dbscan_ri)/sqrt(num_replicas),
                         sd(results$est_dbscan_ari)/sqrt(num_replicas),sd(results$est_dbscan_cpn)/sqrt(num_replicas)
  )
  names(cluster_table_est_se)=c("fadp RI","fadp ARI","fadp cpn",
                                "kmeans RI","kmeans ARI","kmeans cpn",
                                "dbscan RI","dbscan ARI","dbscan cpn"
  )

  if (num_indvs == 100 && scenario == "A"){
    return(list("cluster_table_true"=cluster_table_true,"cluster_table_est"=cluster_table_est,
                "cluster_table_est_se"=cluster_table_est_se,"mse"=mse_sim,"hellinger"=hellinger_sim))

  }

  else{
    return(list("cluster_table_true"=cluster_table_true,"cluster_table_est"=cluster_table_est,
                "cluster_table_est_se"=cluster_table_est_se))

  }
}

#' Runs gam on the given data and returns results
#' Fits a binomial to describe the given in_x
#' @param timestamps01 current time value
#' @param in_x binary series
#' @return fit values and linear predictors both with length of time_series length
RunGam <- function(timestamps01, in_x, basis_size, method)
{
  basis_size_rev_1 <- max(min(round(min(sum(in_x), sum(1-in_x))/2), basis_size ), 5)

  fit_binom_1 <- gam(in_x~s(timestamps01, bs = "cr", m=2, k = basis_size_rev_1),
                     family=binomial(link="probit"), method = method,
                     control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                     optimizer=c("outer","bfgs"))
  p1 <- fit_binom_1$fitted.values
  p1_linpred <- fit_binom_1$linear.predictors

  return(list(prob=p1, linpred=p1_linpred))
}

EstimateCategFuncData <- function(timestamps01, X=NULL, W=NULL, basis_size=25, method="ML", sim=TRUE)
{
  if(is.null(X)) sim<-FALSE

  if(!sim)
  {
    num_indv<- ncol(W)
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
  }

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
#     gam_result_1 <- RunGam(timestamps01, x1, basis_size, method)
#     p1 <- gam_result_1$prob
#     p1_linpred <- gam_result_1$linpred
#
#     gam_result_2 <- RunGam(timestamps01, x2, basis_size, method)
#     p2 <- gam_result_2$prob
#     p2_linpred <- gam_result_2$linpred
#
#     gam_result_3 <- RunGam(timestamps01, x3, basis_size, method)
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

  #///////REDEF FOR PAR
  RunGam <- function(timestamps01, in_x, basis_size, method)
  {
    basis_size_rev_1 <- max(min(round(min(sum(in_x), sum(1-in_x))/2), basis_size ), 5)

    fit_binom_1 <- gam(in_x~s(timestamps01, bs = "cr", m=2, k = basis_size_rev_1),
                       family=binomial(link="probit"), method = method,
                       control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                       optimizer=c("outer","bfgs"))
    p1 <- fit_binom_1$fitted.values
    p1_linpred <- fit_binom_1$linear.predictors

    return(list(prob=p1, linpred=p1_linpred))
  }#///////REDEF FOR PAR

  zp <- foreach (indv = 1:num_indv, .combine = cbind, .init = NULL, .packages = c("mgcv")) %dorng%
  {
    x1<- X[indv,,1]
    x2<- X[indv,,2]
    x3<- X[indv,,3]

    # fit the Binom model
    ###################################################
    ###################################################
    ##updated estimation

    gam_result_1 <- RunGam(timestamps01, x1, basis_size, method)
    p1 <- gam_result_1$prob
    p1_linpred <- gam_result_1$linpred

    gam_result_2 <- RunGam(timestamps01, x2, basis_size, method)
    p2 <- gam_result_2$prob
    p2_linpred <- gam_result_2$linpred

    gam_result_3 <- RunGam(timestamps01, x3, basis_size, method)
    p3 <- gam_result_3$prob
    p3_linpred <- gam_result_3$linpred

    # estimate the latent tranjecotries Z
    exp_p3_linepred <- exp(p3_linpred)
    z1 <- (p1_linpred-p3_linpred)-log( (1+exp(p1_linpred))/(1+exp_p3_linepred))
    z2 <- (p2_linpred-p3_linpred)-log( (1+exp(p2_linpred))/(1+exp_p3_linepred))
    psum <- p1 + p2 + p3

    # Z <- cbind(Z, c(z1,z2))

    # p[indv,,] <- cbind(p1/psum, p2/psum, p3/psum)
    # return(cbind(p1/psum, p2/psum, p3/psum))
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
#' @param k description?? number of psi functions
#' @param mu_1 description??
#' @param mu_2 description??
#' @param score_vals description??
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
  cat("octave", num_indvs, k, num_indvs * k, "\n")
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

result_par <- ClusterSimulation(100, 300, "A", 2)

print(result_par)

# }) # profvis end

parallel::stopCluster(cl = my.cluster)
