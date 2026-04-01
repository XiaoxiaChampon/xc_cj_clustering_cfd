############################################################
# Copyright 2024 Xiaoxia Champon

# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the “Software”), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:

# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.
############################################################
# Purpose: Adding weekend effects in modeling the probability and simulation from real data
# Author:  Xiaoxia Champon
# Date: 11/20/2024
##############################################################
#for Moore-Penrose inverse
library(MASS)

# For: profiling and visualization of profiling
#library(profvis)

# For: gam 
library(mgcv)

# For: cubicspline
library(pracma)

# For: fpca.face
library(refund)

# For: FADPclust
#library(FADPclust)

#For smooth basis
library(fda)

# For: kNNdist
library(dbscan)

# For: elbow
 # devtools::install_github("ahasverus/elbow")
 # library(elbow)

# For: rand.index
library(fossil)

# For: cfda method
#library(cfda)

# For: gather method
library(tidyverse)

#FOR best number of clusters
library(NbClust)

# ---- For: parallelization ----
# For: foreach loop
library(foreach)

###########
library(optparse)

# Define options
option_list <- list(
  make_option(c("-j", "--jobid"), type="integer", default=123,
              help="Job Index", metavar="JOBID"),
  make_option(c("-n", "--numcpus"), type="integer", default=32,
              help="Num CPUs", metavar="NUMCPUS"),
  #make_option(c("-s", "--subjects"), type="integer", default=100,
       #       help="Num Subs", metavar="NUMSUB"),
  
   make_option(c("-r", "--replicas"), type="integer", default=100,
               help="Num Replicas", metavar="NUMREPLICAS")
)

#####need for hazel
# Create parser and parse options
parser <- OptionParser(option_list=option_list)
options <- parse_args(parser)

# options_jobid <- options$jobid
# options_numcpus <- options$numcpus
# options_replicas <- options$replicas
#num_subjects <- options$subjects


options_jobid <- 1
options_numcpus <- 9
options_replicas <- 2
# options_subjects <- 100
###########
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
  # n.cores <- parallel::detectCores()
  n.cores <- options_numcpus
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  doParallel::registerDoParallel(cl = my.cluster)
  cat("Parellel Registered: ", foreach::getDoParRegistered(), " (num cores=", n.cores, ")\n")
  initialized_parallel <- TRUE
  
  # registerDoRNG(123) # ///<<<< THIS CREATES THE ERROR FOR FADPClust !!!
}



#' Create directories
if (!dir.exists("outputs")){
  dir.create("outputs")
}
if (!dir.exists("outputs/Twitter_RMSE")){
  dir.create("outputs/Twitter_RMSE")
}



#estimate z and p
# source("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/acj/gam_weekends.R")
# source("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/acj/time_track_function.R")
source("gam_weekends.R")
source("time_track_function.R")
source("hazel_function.R")


ClusterSimulation <- function(num_subjects, timeseries_length_new, t, eigen_score_final,Z_after_final,
                              num_replicas, est_choice, run_hellinger, temp_folder,options_numcpus)
{
  # timeseries_length_new <- 100
  # num_replicas <- 2
  # est_choice <- "binomial"
  # run_hellinger <- TRUE
  # num_subjects <- 100
  cat("Cluster Simulation\nNum Indvs:\t", num_subjects,
      "\nTimeseries Len:\t", timeseries_length_new,
      "\nNum Replicas:\t", num_replicas)
  
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
    cat("\nNum Indvs:", num_subjects,
        "\tTimeseries Len:", timeseries_length_new,
        "\tNum Replicas:", num_replicas)
    cat("\nReplica: ", replica_idx, "\n")
    
    # generate clusters
    # set.seed(seed_cluster + 100 * replica_idx)
    # cat("\nCluster", replica_idx, " --> seed: ", seed_cluster + 100 * replica_idx, "\n")
    cat("Cluster", replica_idx, "\n")
    
    ZP_data <- GenerateZPData (t, eigen_score_final,Z_after_final)
    Z1 <- ZP_data$Z1
    Z2 <- ZP_data$Z2
    
    # Recover the true probability curves --> could there be more than 3 ???
    p1 <- ZP_data$p1
    p2 <- ZP_data$p2
    p3 <- ZP_data$p3
    
    
    true_cluster <- kmeans_cluster (data=ZP_data$eigen_score_new)$ label
    cat("\nTrue Cluster:\n", true_cluster)
    
    true_cluster_db <- dbscan_cluster (data=ZP_data$eigen_score_new, 3.5)$ label
    cat("\nTrue Cluster DB:\n", true_cluster_db)
    
    #(1) True cluster
    ##find the cluster
    # true_cluster_kmeans <- kmeans_cluster (data=ZP_data$eigen_score_new)
    # cat("\nTrue Cluster Kmeans:\n", true_cluster_kmeans)
    # 
    # true_cluster_db <- dbscan_cluster (data=ZP_data$eigen_score_new, 3.5)
    # cat("\nTrue Cluster DB:\n", true_cluster_db)
    # 
    # true_cluster_fadp <- fadp_cluster (Z1,
    #                                    Z2,
    #                                    t)
    # cat("\nTrue Cluster FADP:\n", true_cluster_fadp)
    # 
    # cfda_scores <- cfda_score_function (W_matrix, options_numcpus,t, basis_size )
    # 
    # true_cluster_cfda <- dbscan_cluster (cfda_scores, 3.5)
    # cat("\nTrue Cluster CFDA:\n", true_cluster_cfda)
    
    
    #############################
    prob_curves <- list(p1=p1, p2=p2, p3=p3)
    num_indv_all <- dim(eigen_score_final$scores)[1]
    timeseries_length <- length(t)
    
    
    categ_func_data_list <- GenerateCategFuncData(prob_curves)
    
    # # what is Q vals ? better name???
    # Q_vals <- unique(c(categ_func_data_list$W))
    # if(is.numeric(Q_vals))
    # {
    #   Q_vals <- sort(Q_vals)
    # }
    
    
    #####################################
    #subset the user
    #672 1344 2016
    individual_sample <- sample(1:num_indv_all, num_subjects)
    num_indvs <- num_subjects
    timestamps01 <- t
    
    st=timestamps01[1]
    et=tail(timestamps01,n=1)
    timestamps01_new=seq(st,et,length=timeseries_length_new)
    
    vec <- matrix(1:num_indvs, nrow=num_indvs, ncol=1)
    
    
    Z1_est <- apply(vec, 1, function(x) {cubicspline(timestamps01, 
                                                     Z1[,individual_sample],
                                                     timestamps01_new)})
    
    Z2_est <- apply(vec, 1, function(x) {cubicspline(timestamps01, 
                                                     Z2[,individual_sample],
                                                     timestamps01_new)})
    p1_est <- apply(vec, 1, function(x) {cubicspline(timestamps01, 
                                                     p1[,individual_sample],
                                                     timestamps01_new)})
    
    p2_est <- apply(vec, 1, function(x) {cubicspline(timestamps01, 
                                                     p2[,individual_sample],
                                                     timestamps01_new)})
    p3_est <- apply(vec, 1, function(x) {cubicspline(timestamps01, 
                                                     p3[,individual_sample],
                                                     timestamps01_new)})
    prob_curves <- list(p1=p1_est, p2=p2_est, p3=p3_est)
    #########9/11/2023  
     # Z_true_curve=array(c(Z1,Z2),dim=c(timeseries_length_new,num_indv_all,2))
     # p_true_curve=array(c(p1,p2,p3),dim=c(timeseries_length,num_indv_all,3))
    Z_true_curve[[replica_idx]]=array(c(Z1_est,Z2_est),dim=c(timeseries_length_new,num_indvs,2))
    p_true_curve[[replica_idx]]=array(c(p1_est,p2_est,p3_est),dim=c(timeseries_length_new,num_indvs,3))
    ############
    
    
    # generate categFuncData
    # set.seed(seed_cfd + 100 * replica_idx)
    # cat("\nCategFD", replica_idx, " --> seed: ", seed_cfd + 100 * replica_idx, "\n")
    cat("CategFD", replica_idx, "\n")
    
    categ_func_data_list <- GenerateCategFuncData(prob_curves)
    
    #what is Q vals ? better name???
    Q_vals <- unique(c(categ_func_data_list$W))
    if(is.numeric(Q_vals))
    {
      Q_vals <- sort(Q_vals)
    }
    for(indv in 1:num_indvs){
    # better names for the following variables ???
    # 1. check weather one category only appears 1 time and is it in the end of the timeseries
    # 2. OR is it appearing only one time in the begining
    # 3. OR if the category is less thant Q total category
    # In general, one category only occurs 2 times
    # If timepoints=300, one category only occurs less than 4 times 3/300=0.01
    #If timepoints=750, one category only occurs less than 6 times 5/750=0.0067
    tolcat <- table(categ_func_data_list$W[,indv])
    #tolcat <- table(categ_func_data_list$W[,indv])
    catorder <- order(tolcat, decreasing = TRUE)
    numcat <- length(catorder)
    refcat <- catorder[numcat]
    count_iter <- 0
    while (count_iter < 100 && 
           ( (numcat < length(Q_vals))
             # ||(timeseries_length==300  && min(as.numeric(tolcat)) < 4)
             # ||(timeseries_length==750  && min(as.numeric(tolcat)) < 10)
             
             ||(timeseries_length<=672  && min(as.numeric(tolcat)) < 4)
             ||(timeseries_length<=1344  && min(as.numeric(tolcat)) < 10)
           )
    )
    {
      count_iter <- count_iter + 1
      
      #new_cluster_data <- GenerateClusterData(setting_choice, scenario, 3, 5, timeseries_length)
      new_cluster_data <- GenerateZPData (t, eigen_score_final,Z_after_final)
      
      new_prob_curves <- list(p1 = new_cluster_data$p1, p2 = new_cluster_data$p2, p3 = new_cluster_data$p3)
      new_categ_func_data_list <- GenerateCategFuncData( new_prob_curves )
      
      # what is this 3 ?? arbitrarily chosen?
      categ_func_data_list$W[, indv] <- new_categ_func_data_list$W[, 3]
      Z1_est[, indv] <- new_cluster_data$Z1[, 3] # latent curves Z1 and Z2
      categ_func_data_list$X[indv, , ] <- 0
      Z2_est[, indv] <- new_cluster_data$Z2[, 3]
      
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
    timestamps01 <- timestamps01_new
    timeKeeperStart("Xiaoxia")
    categFD_est <- EstimateCategFuncData(est_choice, timestamps01_new, categ_func_data_list$W)
    #####################9/11/2023
    Z_est_curve[[replica_idx]]=array(c(categFD_est$Z1_est,categFD_est$Z2_est),dim=c(timeseries_length_new,num_indvs,2))
    p_est_curve[[replica_idx]]=array(c(categFD_est$p1_est,categFD_est$p2_est,categFD_est$p3_est),dim=c(timeseries_length_new,num_indvs,3))
    W_cfd[[replica_idx]]=categ_func_data_list$W
    #####################9/11/2023
    timeKeeperNext()
    
    if (run_hellinger)
    {
      # evaluate performance Z and P
      rmse1_temp <- c(by(mse_bw_matrix(Z1_est, categFD_est$Z1_est, timestamps01_new) , true_cluster[individual_sample], mean))
      rmse2_temp <- c(by(mse_bw_matrix(Z2_est, categFD_est$Z2_est, timestamps01_new), true_cluster[individual_sample], mean))
      rmse[replica_idx, ,] <- rbind(rmse1_temp, rmse2_temp )
      
      error.p1 <- mse_bw_matrixp(p1_est, categFD_est$p1_est, timestamps01_new)
      error.p2 <- mse_bw_matrixp(p2_est, categFD_est$p2_est, timestamps01_new)
      error.p3 <- mse_bw_matrixp(p3_est, categFD_est$p3_est, timestamps01_new)
      
      
      hellinger[replica_idx, ,] <-  rbind( c(by(error.p1, true_cluster[individual_sample], mean)),
                                           c(by(error.p2, true_cluster[individual_sample], mean)),
                                           c(by(error.p3, true_cluster[individual_sample], mean)))
      
    }
    
    # Record clustering performance. USE only Z
    #**** extract_scores_UNIVFPCA is the FASTEST
    
    mfpca_true <-extract_scores_UNIVFPCA(mZ1=Z1_est, mZ2=Z2_est,  tt=timestamps01_new , PVE=0.95)
    #plot(  scores_true$scores[, 1:2])
    timeKeeperStart("univfpca")
    mfpca_est <- extract_scores_UNIVFPCA(mZ1=categFD_est$Z1_est, mZ2=categFD_est$Z2_est, tt=timestamps01_new, PVE=0.95)
    timeKeeperNext()
    
    #KMEANS
    true_kmeans_temp <- kmeans_cluster(data=mfpca_true$scores)$label
    timeKeeperStart("kmeans")
    est_kmeans_temp <- kmeans_cluster(data=mfpca_est$scores)$label
    timeKeeperNext()
    
    #FADP
    true_fadp_temp <- fadp_cluster(mZ1=Z1_est, mZ2=Z2_est, tt=timestamps01)$label
    timeKeeperStart("fadp")
    est_fadp_temp <- fadp_cluster(mZ1=categFD_est$Z1_est, mZ2=categFD_est$Z2_est, tt=timestamps01)$label
    timeKeeperNext()
    #dbscan
    true_dbscan_temp <- dbscan_cluster(data=mfpca_true$scores,3.5)$label
    timeKeeperStart("dbscan")
    est_dbscan_temp <- dbscan_cluster(data=mfpca_est$scores,3.5)$label
    timeKeeperNext()
    ##dbscan cfda
    ##########no cfda
    basis_num <- 10
    #nCores= n.cores
    true_dbscan_temp_cfda <- true_cluster_db
    parallel::stopCluster(cl = my.cluster)
    timeKeeperStart("cfd")
    cfd_scores <- cfda_score_function(categ_func_data_list$W, options_numcpus,timestamps01_new,basis_num )
    my.cluster <- parallel::makeCluster(options_numcpus, type = "PSOCK")
    doParallel::registerDoParallel(cl = my.cluster)
    cat("Parellel Registered: ", foreach::getDoParRegistered(), "\n")
    est_dbscan_temp_cfda <- dbscan_cluster(data=cfd_scores,3.5)$label
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
                                                                      num_replicas, "_", est_choice, "_", run_hellinger, ".RData", sep="")))
  
  return_vals <- list("cluster_table_true"=cluster_table_true,
                      "cluster_table_est"=cluster_table_est,
                      "cluster_table_est_se"=cluster_table_est_se,
                      "Z_true_curves"=Z_true_curve,
                      "Z_est_curves"=Z_est_curve,
                      "p_true_curves"=p_true_curve,
                      "p_est_curves"=p_est_curve,
                      "W_cfd"=W_cfd)
  
  save(return_vals, file=file.path(temp_folder, paste("TwitterSim_", num_indvs, "_", timeseries_length_new, "_",
                                                     num_replicas, "_", est_choice, "_", run_hellinger, ".RData", sep="")))
  
  if (run_hellinger)
  {
    return_vals$mse <- mse_sim
    return_vals$hellinger <- hellinger_sim
  }
  return(return_vals)
}


RunExperiment <- function( num_replicas, est_choice, eigen_score_final,Z_after_final, 
                           options_numcpus,t,some_identifier="noid")
{
  temp_folder <- file.path("outputs", "Twitter_RMSE", paste( num_replicas, "_", est_choice, "_", some_identifier, sep=""))
  # Empty the directory if it exists
  if(dir.exists(temp_folder)){
    unlink(temp_folder, recursive = TRUE)
  }
  dir.create(temp_folder)
  print(temp_folder)
  
  #ClusterSimulation <- function(num_subjects, timeseries_length_new, t, eigen_score_final,Z_after_final,
  #num_replicas, est_choice, run_hellinger, temp_folder,options_numcpus)
  # n100t300C <- ClusterSimulation(100,672,t,eigen_score_final,Z_after_final,num_replicas,est_choice,TRUE,temp_folder,options_numcpus)
  # n100t750C <- ClusterSimulation(100,1344,scenario,num_replicas,est_choice,TRUE,temp_folder,options_numcpus)
  # n100t2000C <- ClusterSimulation(100,2016,scenario,num_replicas,est_choice,TRUE,temp_folder,options_numcpus)
  # 
  # 
  # n500t300C <- ClusterSimulation(500,672,t,eigen_score_final,Z_after_final,num_replicas,est_choice,TRUE,temp_folder,options_numcpus)
  # n500t750C <- ClusterSimulation(500,1344,t,eigen_score_final,Z_after_final,num_replicas,est_choice,TRUE,temp_folder,options_numcpus)
  # n500t2000C <- ClusterSimulation(500,2016,t,eigen_score_final,Z_after_final,num_replicas,est_choice,TRUE,temp_folder,options_numcpus)
  # 
  # 
  # n1000t300C <- ClusterSimulation(1000,672,t,eigen_score_final,Z_after_final,num_replicas,est_choice,TRUE,temp_folder,options_numcpus)
  # n1000t750C <- ClusterSimulation(1000,1344,t,eigen_score_final,Z_after_final,num_replicas,est_choice,TRUE,temp_folder,options_numcpus)
  # n1000t2000C <- ClusterSimulation(1000,2016,st,eigen_score_final,Z_after_final,num_replicas,est_choice,TRUE,temp_folder,options_numcpus)
  # 
  
  n100t300C <- ClusterSimulation(100,600,t,eigen_score_final,Z_after_final,num_replicas,est_choice,TRUE,temp_folder,options_numcpus)
  n100t750C <- ClusterSimulation(100,700,t,eigen_score_final,Z_after_final,num_replicas,est_choice,TRUE,temp_folder,options_numcpus)
  n100t2000C <- ClusterSimulation(100,800,t,eigen_score_final,Z_after_final,num_replicas,est_choice,TRUE,temp_folder,options_numcpus)


  n500t300C <- ClusterSimulation(300,600,t,eigen_score_final,Z_after_final,num_replicas,est_choice,TRUE,temp_folder,options_numcpus)
  n500t750C <- ClusterSimulation(300,700,t,eigen_score_final,Z_after_final,num_replicas,est_choice,TRUE,temp_folder,options_numcpus)
  n500t2000C <- ClusterSimulation(300,800,t,eigen_score_final,Z_after_final,num_replicas,est_choice,TRUE,temp_folder,options_numcpus)


  n1000t300C <- ClusterSimulation(500,600,t,eigen_score_final,Z_after_final,num_replicas,est_choice,TRUE,temp_folder,options_numcpus)
  n1000t750C <- ClusterSimulation(500,700,t,eigen_score_final,Z_after_final,num_replicas,est_choice,TRUE,temp_folder,options_numcpus)
  n1000t2000C <- ClusterSimulation(500,800,t,eigen_score_final,Z_after_final,num_replicas,est_choice,TRUE,temp_folder,options_numcpus)

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
  
  save(true_est_w_data,file=file.path("outputs", paste(num_replicas,est_choice,some_identifier,"true_est_w_data_clustering.RData",sep="_")))
  
  return(true_est_w_data)
}


timeKeeperStart("Xiaoxia")

#RunExperiment <- function( num_replicas, est_choice, eigen_score_final,Z_after_final, 
#                           options_numcpus,t,some_identifier="noid")
load("Twitter_RMSE.RData")

#load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/Twitter_RMSE.RData")
  set.seed(123)
  A_2_binomial <- RunExperiment(2,"binomial",eigen_score_final,Z_after_final,
                                options_numcpus,t,"test")
  save(A_2_binomial, file ="A_2_binomial.RData" )
  timeKeeperNext()
  # num_replicas <- 2
  # est_choice <- "binomial"
  # some_identifier= "test"
   