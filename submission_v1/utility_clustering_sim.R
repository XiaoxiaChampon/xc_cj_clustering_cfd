############################################################
# Article: Functions to Perform Simulation for Clustering Categorical Functional Data
#           File that contains all the functions necessary to estimate categorical FD with 3 CATEGORIES! 
# Author:  Xiaoxia Champon, Ana-Maria Staicu, Chathura Jayalath
# Date: 09/25/2023
##############################################################


source("utility_estimation.R")
source("utility_generate.R")
source("run_gam_function.R")


##load the libraries
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
}

# profvis({


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
    cat(row_name, "calc time taken:", capture.output(this_section_time), "\n")
    time_elapsed[[row_name]] <<- append(time_elapsed[[row_name]], this_section_time)
    last_time <<- this_time
  }
  
  ######
  Z_true_curve <- Z_est_curve <- list()
  p_true_curve <- p_est_curve <- list()
  W_cfd <- list()
  for(replica_idx in 1:num_replicas)
  {
    cat("\nNum Indvs:", num_indvs,
        "\tTimeseries Len:", timeseries_length,
        "\tScenario:", scenario,
        "\tNum Replicas:", num_replicas)
    cat("\nReplica: ", replica_idx, "\n")
    
    # generate clusters
    cat("Cluster", replica_idx, "\n")
    
    cluster_f1 <- GenerateClusterData(1, scenario, 3, cluster_allocation[1], timeseries_length)
    cluster_f2 <- GenerateClusterData(2, scenario, 3, cluster_allocation[2], timeseries_length)
    cluster_f3 <- GenerateClusterData(3, scenario, 3, cluster_allocation[3], timeseries_length)
    
    # Recover the latent Gaussian process 
    Z1 <- cbind(cluster_f1$Z1, cluster_f2$Z1, cluster_f3$Z1)
    Z2 <- cbind(cluster_f1$Z2, cluster_f2$Z2, cluster_f3$Z2)
    
    
    # Recover the true probability curves 
    p1 <- cbind(cluster_f1$p1, cluster_f2$p1, cluster_f3$p1)
    p2 <- cbind(cluster_f1$p2, cluster_f2$p2, cluster_f3$p2)
    p3 <- cbind(cluster_f1$p3, cluster_f2$p3, cluster_f3$p3)
    prob_curves <- list(p1=p1, p2=p2, p3=p3)
    #########9/11/2023  
    Z_true_curve[[replica_idx]] <- array(c(Z1,Z2),dim=c(timeseries_length,num_indvs,2))
    p_true_curve[[replica_idx]] <- array(c(p1,p2,p3),dim=c(timeseries_length,num_indvs,3))
    ############
    
    
    # generate categFuncData
    cat("CategFD", replica_idx, "\n")
    
    categ_func_data_list <- GenerateCategFuncData(prob_curves)
    
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
      
      
      # 1. check weather one category only appears 1 time and is it in the end of the timeseries
      # 2. OR is it appearing only one time in the begining
      # 3. OR if the category is less thant Q total category
      tolcat <- table(categ_func_data_list$W[, indv])
      catorder <- order(tolcat, decreasing = TRUE)
      numcat <- length(catorder)
      refcat <- catorder[numcat]
      count_iter <- 0
      while (count_iter < 100 && ((min(as.numeric(tolcat)) == 1 && categ_func_data_list$W[, indv][timeseries_length] == refcat)
             || (min(as.numeric(tolcat)) == 1 && categ_func_data_list$W[, indv][1] == refcat)
             || (length(categ_func_data_list$W[, indv]) < numcat)
      ))
      {
        count_iter <- count_iter + 1
        
        new_cluster_data <- GenerateClusterData(setting_choice, scenario, 3, 5, timeseries_length)
        
        new_prob_curves <- list(p1 = new_cluster_data$p1, p2 = new_cluster_data$p2, p3 = new_cluster_data$p3)
        new_categ_func_data_list <- GenerateCategFuncData( new_prob_curves )
        
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
    #####################9/11/2023
    Z_est_curve[[replica_idx]]=array(c(categFD_est$Z1_est,categFD_est$Z2_est),dim=c(timeseries_length,num_indvs,2))
    p_est_curve[[replica_idx]]=array(c(categFD_est$p1_est,categFD_est$p2_est,categFD_est$p3_est),dim=c(timeseries_length,num_indvs,3))
    W_cfd[[replica_idx]]=categ_func_data_list$W
    #####################9/11/2023
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
  
  return_vals <- list("cluster_table_true"=cluster_table_true,
                      "cluster_table_est"=cluster_table_est,
                      "cluster_table_est_se"=cluster_table_est_se,
                      "Z_true_curves"=Z_true_curve,
                      "Z_est_curves"=Z_est_curve,
                      "p_true_curves"=p_true_curve,
                      "p_est_curves"=p_est_curve,
                      "W_cfd"=W_cfd)
  
  if (run_hellinger)
  {
    return_vals$mse <- mse_sim
    return_vals$hellinger <- hellinger_sim
  }
  return(return_vals)
}

#############
#Function to systematic sampling
getsys = function(N,n){
  k = floor(N/n)
  r = sample(1:k, 1)
  seq(r, r + k*(n-1), k)
}

#systsample<- data[getsys(nrow(data), 100), ]

#################
ClusterSimulation_Sub <- function(num_indvs, timeseries_length,
                              scenario, num_replicas, est_choice, run_hellinger,
                              nt2000data)
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
  for(replica_idx in 1:num_replicas)
  {
    cat("\nNum Indvs:", num_indvs,
        "\tTimeseries Len:", timeseries_length,
        "\tScenario:", scenario,
        "\tNum Replicas:", num_replicas)
    cat("\nReplica: ", replica_idx, "\n")
    
    # generate clusters
    cat("Cluster", replica_idx, "\n")
    
    
    # generate categFuncData
    cat("CategFD", replica_idx, "\n")
    sample_cluster1 <- sample(1:cluster_allocation[1],num_indvs*occur_fraction[1],replace = FALSE)
    sample_cluster2 <- sample(1:cluster_allocation[2],num_indvs*occur_fraction[2],replace = FALSE)
    sample_cluster3 <- sample(1:cluster_allocation[3],num_indvs*occur_fraction[3],replace = FALSE)
    nsample=c(sample_cluster1,cluster_allocation[1]+sample_cluster2,
              cluster_allocation[1]+cluster_allocation[2]+sample_cluster3 )
    tsample=getsys(2000, timeseries_length)
    
    list_number <- 9
    
    Z1=nt2000data$Z_true_curves[[list_number]][[num_replicas]][tsample,nsample,1]
    Z2=nt2000data$Z_true_curves[[list_number]][[num_replicas]][tsample,nsample,2]
    
    
    # Recover the true probability curves 
    p1 <- nt2000data$p_true_curves[[list_number]][[num_replicas]][tsample,nsample,1]
    p2 <- nt2000data$p_true_curves[[list_number]][[num_replicas]][tsample,nsample,2]
    p3 <- nt2000data$p_true_curves[[list_number]][[num_replicas]][tsample,nsample,3]
    prob_curves <- list(p1=p1, p2=p2, p3=p3)
    #########9/11/2023  
    Z_true_curve[[replica_idx]]=array(c(Z1,Z2),dim=c(timeseries_length,num_indvs,2))
    p_true_curve[[replica_idx]]=array(c(p1,p2,p3),dim=c(timeseries_length,num_indvs,3))
    ############
    
    # 
    # Z1_est=nt2000data$Z_est_curves[[1]][[num_replicas]][tsample,nsample,1]
    # Z2_est=nt2000data$Z_est_curves[[1]][[num_replicas]][tsample,nsample,2]
  
  
    #Estimation
    cat("EstimateCategFuncData", replica_idx, "\n")
    #timestamps_2000 <- seq(from = 0.0001, to = 1, length=timeseries_length)
    timestamps_2000 <- seq(from = 0.0001, to = 1, length=2000)
    timestamps01 <- timestamps_2000[tsample]
    timeKeeperStart("Xiaoxia")
    categ_func_data_list_W =nt2000data$W_cfd[[list_number]][[num_replicas]][tsample,nsample]
    categFD_est <- EstimateCategFuncData(est_choice, timestamps01,categ_func_data_list_W  )
    #####################9/11/2023
    Z_est_curve[[replica_idx]]=array(c(categFD_est$Z1_est,categFD_est$Z2_est),dim=c(timeseries_length,num_indvs,2))
    p_est_curve[[replica_idx]]=array(c(categFD_est$p1_est,categFD_est$p2_est,categFD_est$p3_est),dim=c(timeseries_length,num_indvs,3))
    W_cfd[[replica_idx]]=categ_func_data_list_W
    #####################9/11/2023
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
    cfd_scores <- cfda_score_function(categ_func_data_list_W, nCores,timestamps01,basis_num )
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
  
  return_vals <- list("cluster_table_true"=cluster_table_true,
                      "cluster_table_est"=cluster_table_est,
                      "cluster_table_est_se"=cluster_table_est_se,
                      "Z_true_curves"=Z_true_curve,
                      "Z_est_curves"=Z_est_curve,
                      "p_true_curves"=p_true_curve,
                      "p_est_curves"=p_est_curve,
                      "W_cfd"=W_cfd)
  
  if (run_hellinger)
  {
    return_vals$mse <- mse_sim
    return_vals$hellinger <- hellinger_sim
  }
  return(return_vals)
}
#######################################################


RunExperiment <- function(scenario, num_replicas, est_choice)
{
  # scenario="C"
  # num_replicas=3
  # est_choice="probit"
  n100t300C <- ClusterSimulation(100,300,scenario,num_replicas,est_choice,TRUE)
  n100t750C <- ClusterSimulation(100,750,scenario,num_replicas,est_choice,TRUE)
  n100t2000C <- ClusterSimulation(100,2000,scenario,num_replicas,est_choice,TRUE)
  
  
  n500t300C <- ClusterSimulation(500,300,scenario,num_replicas,est_choice,TRUE)
  n500t750C <- ClusterSimulation(500,750,scenario,num_replicas,est_choice,TRUE)
  n500t2000C <- ClusterSimulation(500,2000,scenario,num_replicas,est_choice,TRUE)
  
  
  n1000t300C <- ClusterSimulation(1000,300,scenario,num_replicas,est_choice,TRUE)
  n1000t750C <- ClusterSimulation(1000,750,scenario,num_replicas,est_choice,TRUE)
  n1000t2000C <- ClusterSimulation(1000,2000,scenario,num_replicas,est_choice,TRUE)
  
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
  
  
  save(true_tableC,est_tableC,est_tableC_se,file=paste(scenario,num_replicas,est_choice,"_beforeMSE_clustering.RData",sep="_"))
  
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
  
  save(true_est_w_data,file=paste(scenario,num_replicas,est_choice,"true_est_w_data_clustering.RData",sep="_"))
}


RunExperimentSub <- function(scenario, num_replicas, est_choice,n100t2000data)
{
  # scenario="C"
  # num_replicas=3
  # est_choice="probit"
  n100t300C <- ClusterSimulation_Sub(100,300,scenario,num_replicas,est_choice,TRUE,n1000t2000data)
  n100t750C <- ClusterSimulation_Sub(100,750,scenario,num_replicas,est_choice,TRUE,n1000t2000data)
  n100t2000C <- n1000t2000data
  
  
  n500t300C <- ClusterSimulation_Sub(500,300,scenario,num_replicas,est_choice,TRUE,n1000t2000data)
  n500t750C <- ClusterSimulation_Sub(500,750,scenario,num_replicas,est_choice,TRUE,n1000t2000data)
  n500t2000C <- n1000t2000data
  
  
  n1000t300C <- ClusterSimulation_Sub(1000,300,scenario,num_replicas,est_choice,TRUE,n1000t2000data)
  n1000t750C <- ClusterSimulation_Sub(1000,750,scenario,num_replicas,est_choice,TRUE,n1000t2000data)
  n1000t2000C <- n1000t2000data
  
  
  true_tableC <- rbind(n100t300C$cluster_table_true,n100t750C$cluster_table_true,n100t2000C$cluster_table_true[3,],
                       n500t300C$cluster_table_true,n500t750C$cluster_table_true,n500t2000C$cluster_table_true[6,],
                       n1000t300C$cluster_table_true,n1000t750C$cluster_table_true,n1000t2000C$cluster_table_true[9,])
  rownames(true_tableC) <- c("n100t300","n100t750","n100t2000",
                             "n500t300","n500t750","n500t2000",
                             "n1000t300","n1000t750","n1000t2000")
  
  est_tableC <- rbind(n100t300C$cluster_table_est,n100t750C$cluster_table_est,n100t2000C$cluster_table_est[3,],
                      n500t300C$cluster_table_est,n500t750C$cluster_table_est,n500t2000C$cluster_table_est[6,],
                      n1000t300C$cluster_table_est,n1000t750C$cluster_table_est,n1000t2000C$cluster_table_est[9,])
  rownames(est_tableC) <- c("n100t300","n100t750","n100t2000",
                            "n500t300","n500t750","n500t2000",
                            "n1000t300","n1000t750","n1000t2000")
  
  
  est_tableC_se <- rbind(n100t300C$cluster_table_est_se,n100t750C$cluster_table_est_se,n100t2000C$cluster_table_est_se[3,],
                         n500t300C$cluster_table_est_se,n500t750C$cluster_table_est_se,n500t2000C$cluster_table_est_se[6,],
                         n1000t300C$cluster_table_est_se,n1000t750C$cluster_table_est_se,n1000t2000C$cluster_table_est_se[9,])
  rownames(est_tableC_se) <- c("n100t300","n100t750","n100t2000",
                               "n500t300","n500t750","n500t2000",
                               "n1000t300","n1000t750","n1000t2000")
  
  
  save(true_tableC,est_tableC,est_tableC_se,file=paste(scenario,num_replicas,est_choice,"_beforeMSE_clustering.RData",sep="_"))
  
  mse_tableC100 <- rbind(
    c(n100t300C$mse[1,1],n100t300C$mse[2,1],n100t300C$hellinger[1,1],n100t300C$hellinger[2,1],n100t300C$hellinger[3,1]),
    c(n100t750C$mse[1,1],n100t750C$mse[2,1],n100t750C$hellinger[1,1],n100t750C$hellinger[2,1],n100t750C$hellinger[3,1]),
    c(n100t2000C$mse_tableC100[3,1],n100t2000C$mse_tableC100[3,2],n100t2000C$mse_tableC100[3,3],n100t2000C$mse_tableC100[3,4],n100t2000C$mse_tableC100[3,5]),
    
    c(n100t300C$mse[1,2],n100t300C$mse[2,2],n100t300C$hellinger[1,2],n100t300C$hellinger[2,2],n100t300C$hellinger[3,2]),
    c(n100t750C$mse[1,2],n100t750C$mse[2,2],n100t750C$hellinger[1,2],n100t750C$hellinger[2,2],n100t750C$hellinger[3,2]),
    c(n100t2000C$mse_tableC100[6,1],n100t2000C$mse_tableC100[6,2],n100t2000C$mse_tableC100[6,3],n100t2000C$mse_tableC100[6,4],n100t2000C$mse_tableC100[6,5]),
    
    c(n100t300C$mse[1,3],n100t300C$mse[2,3],n100t300C$hellinger[1,3],n100t300C$hellinger[2,3],n100t300C$hellinger[3,3]),
    c(n100t750C$mse[1,3],n100t750C$mse[2,3],n100t750C$hellinger[1,3],n100t750C$hellinger[2,3],n100t750C$hellinger[3,3]),
    c(n100t2000C$mse_tableC100[9,1],n100t2000C$mse_tableC100[9,3],n100t2000C$mse_tableC100[9,3],n100t2000C$mse_tableC100[9,4],n100t2000C$mse_tableC100[9,5])
  )
  
  
  rownames(mse_tableC100) <- c("s1n100t300","s1n100t750","s1n100t2000",
                               "s2n100t300","s2n100t750","s2n100t2000",
                               "s3n100t300","s3n100t750","s3n100t2000")
  colnames(mse_tableC100) <- c("z1","z2","p1","p2","p3")
  
  mse_tableC500 <- rbind(
    c(n500t300C$mse[1,1],n500t300C$mse[2,1],n500t300C$hellinger[1,1],n500t300C$hellinger[2,1],n500t300C$hellinger[3,1]),
    c(n500t750C$mse[1,1],n500t750C$mse[2,1],n500t750C$hellinger[1,1],n500t750C$hellinger[2,1],n500t750C$hellinger[3,1]),
    c(n500t2000C$mse_tableC500[3,1],n500t2000C$mse_tableC500[3,2],n500t2000C$mse_tableC500[3,3],n500t2000C$mse_tableC500[3,4],n500t2000C$mse_tableC500[3,5]),
    
    c(n500t300C$mse[1,2],n500t300C$mse[2,2],n500t300C$hellinger[1,2],n500t300C$hellinger[2,2],n500t300C$hellinger[3,2]),
    c(n500t750C$mse[1,2],n500t750C$mse[2,2],n500t750C$hellinger[1,2],n500t750C$hellinger[2,2],n500t750C$hellinger[3,2]),
    c(n500t2000C$mse_tableC500[6,1],n500t2000C$mse_tableC500[6,2],n500t2000C$mse_tableC500[6,3],n500t2000C$mse_tableC500[6,4],n500t2000C$mse_tableC500[6,5]),
    
    c(n500t300C$mse[1,3],n500t300C$mse[2,3],n500t300C$hellinger[1,3],n500t300C$hellinger[2,3],n500t300C$hellinger[3,3]),
    c(n500t750C$mse[1,3],n500t750C$mse[2,3],n500t750C$hellinger[1,3],n500t750C$hellinger[2,3],n500t750C$hellinger[3,3]),
    c(n500t2000C$mse_tableC500[9,1],n500t2000C$mse_tableC500[9,2],n500t2000C$mse_tableC500[9,3],n500t2000C$mse_tableC500[9,4],n500t2000C$mse_tableC500[9,5])
  )
  
  
  rownames(mse_tableC500) <- c("s1n500t300","s1n500t750","s1n500t2000",
                               "s2n500t300","s2n500t750","s2n500t2000",
                               "s3n500t300","s3n500t750","s3n500t2000")
  colnames(mse_tableC500) <- c("z1","z2","p1","p2","p3")
  
  
  mse_tableC1000 <- rbind(
    c(n1000t300C$mse[1,1],n1000t300C$mse[2,1],n1000t300C$hellinger[1,1],n1000t300C$hellinger[2,1],n1000t300C$hellinger[3,1]),
    c(n1000t750C$mse[1,1],n1000t750C$mse[2,1],n1000t750C$hellinger[1,1],n1000t750C$hellinger[2,1],n1000t750C$hellinger[3,1]),
    c(n1000t2000C$mse_tableC1000[3,1],n1000t2000C$mse_tableC1000[3,2],n1000t2000C$mse_tableC1000[3,3],n1000t2000C$mse_tableC1000[3,4],n1000t2000C$mse_tableC1000[3,5]),
    
    c(n1000t300C$mse[1,2],n1000t300C$mse[2,2],n1000t300C$hellinger[1,2],n1000t300C$hellinger[2,2],n1000t300C$hellinger[3,2]),
    c(n1000t750C$mse[1,2],n1000t750C$mse[2,2],n1000t750C$hellinger[1,2],n1000t750C$hellinger[2,2],n1000t750C$hellinger[3,2]),
    c(n1000t2000C$mse_tableC1000[6,1],n1000t2000C$mse_tableC1000[6,2],n1000t2000C$mse_tableC1000[6,3],n1000t2000C$mse_tableC1000[6,4],n1000t2000C$mse_tableC1000[6,5]),
    
    c(n1000t300C$mse[1,3],n1000t300C$mse[2,3],n1000t300C$hellinger[1,3],n1000t300C$hellinger[2,3],n1000t300C$hellinger[3,3]),
    c(n1000t750C$mse[1,3],n1000t750C$mse[2,3],n1000t750C$hellinger[1,3],n1000t750C$hellinger[2,3],n1000t750C$hellinger[3,3]),
    c(n1000t2000C$mse_tableC1000[9,1],n1000t2000C$mse_tableC1000[9,2],n1000t2000C$mse_tableC1000[9,3],n1000t2000C$mse_tableC1000[9,4],n1000t2000C$mse_tableC1000[9,5])
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
                                                n100t2000C$Z_true_curves[[3]],
                                                n500t300C$Z_true_curves,
                                                n500t750C$Z_true_curves,
                                                n500t2000C$Z_true_curves[[6]],
                                                n1000t300C$Z_true_curves,
                                                n1000t750C$Z_true_curves,
                                                n1000t2000C$Z_true_curves[[9]]),
                           "Z_est_curves"=list(n100t300C$Z_est_curves,
                                               n100t750C$Z_est_curves,
                                               n100t2000C$Z_est_curves[[3]],
                                               n500t300C$Z_est_curves,
                                               n500t750C$Z_est_curves,
                                               n500t2000C$Z_est_curves[[6]],
                                               n1000t300C$Z_est_curves,
                                               n1000t750C$Z_est_curves,
                                               n1000t2000C$Z_est_curves[[9]]),
                           "p_true_curves"=list(n100t300C$p_true_curves,
                                                n100t750C$p_true_curves,
                                                n100t2000C$p_true_curves[[3]],
                                                n500t300C$p_true_curves,
                                                n500t750C$p_true_curves,
                                                n500t2000C$p_true_curves[[6]],
                                                n1000t300C$p_true_curves,
                                                n1000t750C$p_true_curves,
                                                n1000t2000C$p_true_curves[[9]]),
                           "p_est_curves"=list(n100t300C$p_est_curves,
                                               n100t750C$p_est_curves,
                                               n100t2000C$p_est_curves[[3]],
                                               n500t300C$p_est_curves,
                                               n500t750C$p_est_curves,
                                               n500t2000C$p_est_curves[[6]],
                                               n1000t300C$p_est_curves,
                                               n1000t750C$p_est_curves,
                                               n1000t2000C$p_est_curves[[9]]),
                           "W_cfd"=list(n100t300C$W_cfd,
                                        n100t750C$W_cfd,
                                        n100t2000C$W_cfd[[3]],
                                        n500t300C$W_cfd,
                                        n500t750C$W_cfd,
                                        n500t2000C$W_cfd[[6]],
                                        n1000t300C$W_cfd,
                                        n1000t750C$W_cfd,
                                        n1000t2000C$W_cfd[[9]])
  ))
  save(true_est_w_data,file=paste(scenario,num_replicas,est_choice,"true_est_w_data_clustering.RData",sep="_"))
}
set.seed(123)
C_2_probit <- RunExperiment("C",2,"probit")


if(run_parallel)
{
  parallel::stopCluster(cl = my.cluster)
  initialized_parallel <- FALSE
}
