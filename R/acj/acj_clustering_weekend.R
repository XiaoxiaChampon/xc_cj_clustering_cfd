
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
# Purpose: Adding weekend effects in Twitter Simulation
# Author:  Xiaoxia Champon
# Date: 11/06/2024
##############################################################

#load("Twitter_figure_label_mul_weekend.RData")


#length(vec1) #3743, 43, 28
#length(vec2)
#length(vec0)
#3830, 5
#head(tclusterdata)
#ksi1       ksi2 weekend_vector1 V4 Cluster

# vec1_score1_mean <- mean(tclusterdata[vec1,1])
# vec1_score2_mean <- mean(tclusterdata[vec1,2])
# 
# vec1_score1_sd <- sd(tclusterdata[vec1,1])
# vec1_score2_sd <- sd(tclusterdata[vec1,2])
# 
# vec1_weekend1_mean <- mean(tclusterdata[vec1,3])
# vec1_weekend2_mean <- mean(tclusterdata[vec1,4])
# 
# vec1_weekend1_sd <- sd(tclusterdata[vec1,3])
# vec1_weekend2_sd <- sd(tclusterdata[vec1,4])
# ####################################
# vec2_score1_mean <- mean(tclusterdata[vec2,1])
# vec2_score2_mean <- mean(tclusterdata[vec2,2])
# 
# vec2_score1_sd <- sd(tclusterdata[vec2,1])
# vec2_score2_sd <- sd(tclusterdata[vec2,2])
# 
# vec2_weekend1_mean <- mean(tclusterdata[vec2,3])
# vec2_weekend2_mean <- mean(tclusterdata[vec2,4])
# 
# vec2_weekend1_sd <- sd(tclusterdata[vec2,3])
# vec2_weekend2_sd <- sd(tclusterdata[vec2,4])
# ####################################
# vec0_score1_mean <- mean(tclusterdata[vec0,1])
# vec0_score2_mean <- mean(tclusterdata[vec0,2])
# 
# vec0_score1_sd <- sd(tclusterdata[vec0,1])
# vec0_score2_sd <- sd(tclusterdata[vec0,2])
# 
# vec0_weekend1_mean <- mean(tclusterdata[vec0,3])
# vec0_weekend2_mean <- mean(tclusterdata[vec0,4])
# 
# vec0_weekend1_sd <- sd(tclusterdata[vec0,3])
# vec0_weekend2_sd <- sd(tclusterdata[vec0,4])
# ##################################
# eigenf_func <- eigen_score$Phi
# weekend_vector <- as.factor(c(rep(c(rep(0,5*24*60/20),rep(1,2*24*60/20)),4)))
# 
# save(vec1_score1_mean, vec1_score2_mean, vec1_score1_sd, vec1_score2_sd,
#      vec1_weekend1_mean, vec1_weekend2_mean, vec1_weekend1_sd, vec1_weekend2_sd,
#      vec2_score1_mean, vec2_score2_mean, vec2_score1_sd, vec2_score2_sd,
#      vec2_weekend1_mean, vec2_weekend2_mean, vec2_weekend1_sd, vec2_weekend2_sd,
#      vec0_score1_mean, vec0_score2_mean, vec0_score1_sd, vec0_score2_sd,
#      vec0_weekend1_mean, vec0_weekend2_mean, vec0_weekend1_sd, vec0_weekend2_sd,
#      eigenf_func,   file = "Twiiter_Sim_Input_New.RData" )
#load("Twitter_eigen_weekend.RData")
#Z_after <-  t(sqrt(length(timestamps01))*eigen_score $scores%*%ginv(eigen_score$Phi))
###########################################
load("Twiiter_Sim_Input_New.RData")
source("hazel_function.R")
source("time_track_function.R")
############################################
# For: profiling and visualization of profiling
# library(profvis)
# profiling_result <- profvis({
#   
#for eval.fd
library(fda)

#for semimetric.basis
library(fda.usc)

#for cluster.stats
library(fpc)

# For: gam 
library(mgcv)

# For: cubicspline
library(pracma)

# For: fpca.face
library(refund)

# For: FADPclust
#library(FADPclust)

# For: kNNdist
library(dbscan)

# For: elbow
#devtools::install_github("ahasverus/elbow")
#library(elbow)

# For: rand.index
library(fossil)

# For: cfda method
#library(cfda)

# For: gather method
library(tidyverse)

###########
library(optparse)

# Define options
option_list <- list(
  make_option(c("-j", "--jobid"), type="integer", default=123,
              help="Job Index", metavar="JOBID"),
  make_option(c("-n", "--numcpus"), type="integer", default=32,
              help="Num CPUs", metavar="NUMCPUS"),
  make_option(c("-r", "--replicas"), type="integer", default=100,
              help="Num Replicas", metavar="NUMREPLICAS")
  # make_option(c("-s", "--subjects"), type="integer", default=100,
  #             help="Num Subjects/Individuals", metavar="NUMSUBJECTS"),
  # make_option(c("-b", "--boots"), type="integer", default=100,
  #             help="Num Bootstraps", metavar="NUMBOOTS"),
  # make_option(c("-l", "--timelength"), type="integer", default=90,
  #             help="Time Length", metavar="TIMELENGTH")
)

#####need for hazel
# Create parser and parse options
parser <- OptionParser(option_list=option_list)
options <- parse_args(parser)

options_jobid <- options$jobid
options_numcpus <- options$numcpus
options_replicas <- options$replicas
#options_subjects <- options$subjects
# options_boots <- options$boots
# options_timelength <- options$timelength
#####################
# 
options_jobid <- 1
options_numcpus <- 9
options_replicas <- 1



# options_subjects <- 100
# options_boots <- 100
# options_timelength <- 90
# Use the options
cat("Job Idx:", options_jobid, "\n")
cat("Num CPUs:", options_numcpus, "\n")
cat("Num Replicas:", options_replicas, "\n")
#cat("Num Subjects:", options_subjects, "\n")


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
  #n.cores <- parallel::detectCores() - 1
  n.cores <- options_numcpus
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
  
  out_dbscan <- dbscan::dbscan(data, eps =epsoptimal , minPts = minPts)
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
    source("trapzfnum_function.R")
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
    source("trapzfnum_function.R")
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
EstimateCategFuncData <- function(choice, timestamps01, W, basis_size=25, method="ML", indvs_to_eval)
{
  if(choice == "probit"){
    X <- GetXFromW(W)
    return(EstimateCategFuncData_probit(timestamps01, X, basis_size, method, 1/150))
  }else if(choice == "binomial"){
    X <- GetXFromW(W)
    return(EstimateCategFuncData_binorm(timestamps01, X, basis_size, method))
  }else if(choice == "multinomial"){
    #return(EstimateCategFuncData_multinormial(timestamps01, W, basis_size, method))
    return(EstimateCategFuncData_multinormial_weekend_parallel(timestamps01, W, basis_size, method, indvs_to_eval))
  }
}

#' Fit multinomial from W categorical data to extra smooth latent curves
#' Fits a binomial to describe the given in_x
#' @param timestamps01 current time value
#' @param W : 2D categorical data matrix, t * n dimension, n is the number of subjects, t is the time
#' @return list: fit values and linear predictors both with length of time_series length, Z is t*n, p: t*n

EstimateCategFuncData_multinormial_weekend_parallel <- function(timestamps01, W, basis_size=25, method="ML", indvs_to_eval)
{
  # num_indv <- ncol(W)
  timeseries_length <-  nrow(W)
  category_count <- length(unique(c(W)))
  #weekend_vector <- as.factor(c(rep(c(rep(0,480),rep(1,192)),4))[1:timeseries_length])
  weekend_vector <- as.factor(c(rep(c(rep(0,5*24*60/20),rep(1,2*24*60/20)),4)))[1:timeseries_length]
  Z <- NULL
  # prob<-array(0, c(num_indv, timeseries_length , category_count))
  # weekend_vector_coef <- matrix(0, num_indv, category_count-1)
  
  
  # Run parallel only if its more than the number of CPUs
  if(options_numcpus > length(indvs_to_eval)){
    # browser()
    Z_P_WeekendCoef <- NULL
    for (idx in indvs_to_eval) {
      # browser()
      print(idx)
      try_result <- tryCatch({
        return_val <- numeric(timeseries_length*(category_count+category_count-1)+category_count-1)
        if ( length(table(W[,idx])) == category_count) {
          
          fit_binom<-gam(list(W[,idx]-1~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector,
                              ~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector
          ),
          family=multinom(K=category_count-1), method = method,
          control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
          optimizer=c("outer","bfgs")) 
          #####
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
          if (names(table(W[,idx]))[2]=="3"){
            W[,idx][W[,idx]==3] <- 2
            basis_size_rev <- max(min(round(min(unname(table(W[,idx])[2]), sum(1-unname(table(W[,idx])[2])))/2), basis_size ), 5)
            fit_binom <- gam(W[,idx]-1~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                             family = "binomial", method = method,
                             control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                             optimizer=c("outer","bfgs"))
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
            ##########################
            p3 <- predict(fit_binom,type="response")
            p1 <- 1-p3
            p2 <- rep(0,timeseries_length)
            
          }else {
            basis_size_rev <- max(min(round(min(unname(table(W[,idx])[2]), sum(1-unname(table(W[,idx])[2])))/2), basis_size ), 5)
            fit_binom <- gam(W[,idx]-1~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                             family = "binomial", method = method,
                             control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                             optimizer=c("outer","bfgs"))
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
            p2 <- predict(fit_binom,type="response")
            p1 <- 1-p2
            p3 <- rep(0,timeseries_length)
          } # sub if else
        } # big if else
        
        #2t*n matrix
        # Z<- cbind(Z, c(z1,z2))
        # ##find probability
        # Z_cbind=cbind(z1,z2)
        # exp_z=exp(Z_cbind)
        # denominator_p=1+exp_z[,1]+exp_z[,2]
        # p1 <- exp_z[,1]/denominator_p
        # p2 <- exp_z[,2]/denominator_p
        # p3=1/denominator_p
        #3D matrix t*n*category 
        #prob[idx,,] <- cbind(p1, p2, p3)
        # 5*t +2 length
        return_val <- c(z1,z2,p1,p2,p3, weekend_vector_coef)
        print("Good")
        ret <- list(idx, TRUE, return_val)
        ret
      },
      error = function(e){
        print(paste("Error in iteration",idx, ":", e$message))
        ret <- list(idx, FALSE, e$message)
        ret
      })
      
      # return(try_result)
      Z_P_WeekendCoef[[idx]] <- try_result
      
    } # end of for loop
  
  }else{
  
  # Z_P_WeekendCoef <- NULL
  # for (idx in indvs_to_eval) {
    
  Z_P_WeekendCoef <- foreach(idx = indvs_to_eval , .packages = c("mgcv")) %dorng% {
    
    #T_rep <- foreach(this_row = 1:5) %dorng%
    # { #source("./source_code/R/data_generator.R")
    
      # print(idx)
    try_result <- tryCatch({
      return_val <- numeric(timeseries_length*(category_count+category_count-1)+category_count-1)
      if ( length(table(W[,idx])) == category_count) {
        
        fit_binom<-gam(list(W[,idx]-1~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector,
                            ~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector
        ),
        family=multinom(K=category_count-1), method = method,
        control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
        optimizer=c("outer","bfgs")) 
        #####
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
        if (names(table(W[,idx]))[2]=="3"){
          W[,idx][W[,idx]==3] <- 2
          basis_size_rev <- max(min(round(min(unname(table(W[,idx])[2]), sum(1-unname(table(W[,idx])[2])))/2), basis_size ), 5)
          fit_binom <- gam(W[,idx]-1~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                           family = "binomial", method = method,
                           control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                           optimizer=c("outer","bfgs"))
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
          ##########################
          p3 <- predict(fit_binom,type="response")
          p1 <- 1-p3
          p2 <- rep(0,timeseries_length)
          
        }else {
          basis_size_rev <- max(min(round(min(unname(table(W[,idx])[2]), sum(1-unname(table(W[,idx])[2])))/2), basis_size ), 5)
          fit_binom <- gam(W[,idx]-1~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                           family = "binomial", method = method,
                           control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                           optimizer=c("outer","bfgs"))
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
          p2 <- predict(fit_binom,type="response")
          p1 <- 1-p2
          p3 <- rep(0,timeseries_length)
        } # sub if else
      } # big if else
      
      #2t*n matrix
      # Z<- cbind(Z, c(z1,z2))
      # ##find probability
      # Z_cbind=cbind(z1,z2)
      # exp_z=exp(Z_cbind)
      # denominator_p=1+exp_z[,1]+exp_z[,2]
      # p1 <- exp_z[,1]/denominator_p
      # p2 <- exp_z[,2]/denominator_p
      # p3=1/denominator_p
      #3D matrix t*n*category 
      #prob[idx,,] <- cbind(p1, p2, p3)
      # 5*t +2 length
      return_val <- c(z1,z2,p1,p2,p3, weekend_vector_coef)
      print("Good")
      ret <- list(idx, TRUE, return_val)
      ret
    },
    error = function(e){
      print(paste("Error in iteration",idx, ":", e$message))
      ret <- list(idx, FALSE, e$message)
      ret
    })
      
      return(try_result)
      # Z_P_WeekendCoef[[idx]] <- try_result
      
    } # end of for loop
  }
  
  # browser()
  non_null_ZPW <- Filter(Negate(is.null), Z_P_WeekendCoef)
  bad_local_indexes <- which(sapply(non_null_ZPW, function(x) x[[2]] == FALSE))
  bad_indexes <- c()
  if(0 < length(bad_local_indexes)){
    bad_indexes <- sapply(bad_local_indexes, function(x) non_null_ZPW[[x]][[1]])
  }
  
  
  return(list(bad_indexes=bad_indexes, partial_outcome=Z_P_WeekendCoef, indvs_to_eval=indvs_to_eval))
  
  # Z_P_WeekendCoef <- do.call(rbind, Z_P_WeekendCoef)
  # #t*n
  # partial_outcome <- list(Z1_est=t(Z_P_WeekendCoef[,1:timeseries_length]), 
  #      Z2_est=t(Z_P_WeekendCoef[,(1+timeseries_length):(2*timeseries_length)]),
  #      p1_est=t(Z_P_WeekendCoef[,(1+timeseries_length*2):(3*timeseries_length)]), 
  #      p2_est=t(Z_P_WeekendCoef[,(1+timeseries_length*3):(4*timeseries_length)]), 
  #      p3_est=t(Z_P_WeekendCoef[,(1+timeseries_length*4):(5*timeseries_length)]) ,
  #      weekend_vector_coef = Z_P_WeekendCoef[,(1+timeseries_length*(category_count+category_count-1)):dim(Z_P_WeekendCoef)[2]])
  # 
  # partial_return <- list(partial_outcome=Z_P_WeekendCoef,
  #             bad_indexes=bad_indexes,
  #             indvs_evaled=indvs_to_eval)
  # 
  # return(partial_return)
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
#k=2
#setting = 1
GenerateClusterData <- function(setting, scenario, k, num_indvs, timeseries_length, weekend_vector, eigenf_func)
{
  setting_object <- GetMuAndScore(setting, scenario, k, num_indvs, weekend_vector)
  cluster_f <- GenerateClusterDataScenario(num_indvs,
                                           timeseries_length,
                                           k,
                                           mu_1 = setting_object$mu_1,
                                           mu_2 = setting_object$mu_2,
                                           score_vals = setting_object$score_vals,
                                           weekend_columns = setting_object$weekend_columns,
                                           eigenf_func = eigenf_func)
  return (cluster_f )
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
#'
# cluster_allocation <- c(10,20,30)
 # abc2 <- GetMuAndScore (2,"B" , 2, 30, weekend_vector)
 # abc1 <- GetMuAndScore (1,"B" , 2, 50, weekend_vector)
 # abc3 <- GetMuAndScore (3,"B" , 2, 20, weekend_vector)
 # abc1$score_vals
 # abc2$score_vals
 # abc3$score_vals
# setting <- 2
# scenario <- "B"
# k=2
# weekend_vector <- as.factor(c(rep(c(rep(0,5*24*60/20),rep(1,2*24*60/20)),4)))[1:timeseries_length]
GetMuAndScore <- function(setting, scenario, k,  num_indvs, weekend_vector)
{

  #all_score_values = rep(0, 4)
  
  if(1 == setting)
  { 
    weekend_1 <- rnorm(num_indvs, vec1_weekend1_mean*0.1, vec1_weekend1_sd*0.1)
    mu_1 <-  outer(as.numeric(weekend_vector), weekend_1, "*") 
    weekend_2 <- rnorm(num_indvs, vec1_weekend2_mean, vec1_weekend2_sd*0.1)
    mu_2 <- outer(as.numeric(weekend_vector), weekend_2, "*") 
    
    score_front <- c(vec1_score1_mean*0.1, vec1_score1_sd*0.1,vec1_score2_mean, vec1_score2_sd*0.1)
  } else if(2 == setting)
  {
    weekend_1 <- rnorm(num_indvs, vec2_weekend1_mean*0.025, vec2_weekend1_sd*0.05)
    mu_1 <-  outer(as.numeric(weekend_vector), weekend_1, "*") 
    weekend_2 <- rnorm(num_indvs, vec2_weekend2_mean*0.1, vec2_weekend2_sd)
    mu_2 <- outer(as.numeric(weekend_vector), weekend_2, "*") 
    
    score_front <- c(vec2_score1_mean*0.025, vec2_score1_sd*0.05,vec2_score2_mean*0.1, vec2_score2_sd)
  } else if(3 == setting)
  {
    weekend_1 <- rnorm(num_indvs, vec0_weekend1_mean*0.015, vec0_weekend1_sd*0.025)
    mu_1 <-  outer(as.numeric(weekend_vector), weekend_1, "*") 
    weekend_2 <- rnorm(num_indvs, vec0_weekend2_mean*0.05, vec0_weekend2_sd*0.025)
    mu_2 <- outer(as.numeric(weekend_vector), weekend_2, "*") 
    
    score_front <- c(vec0_score1_mean*0.015, vec0_score1_sd*0.025,vec0_score2_mean*0.05, vec0_score2_sd*0.025)
  }
  
  # for(idx in 1:length(score_front))
  # {
  #   all_score_values[idx] <- score_front[idx]
  # }
  # 
  weekend_columns <- cbind(weekend_1, weekend_2)
  return(list("mu_1" = mu_1, "mu_2" = mu_2, "score_vals" = score_front,
              "weekend_columns" = weekend_columns))
}


GetMuAndScore <- function(setting, scenario, k,  num_indvs, weekend_vector)
{

  #all_score_values = rep(0, 4)
  
  if(1 == setting)
  { 
    weekend_1 <- rnorm(num_indvs, vec1_weekend1_mean*0.1, vec1_weekend1_sd*0.1)
    mu_1 <-  outer(as.numeric(weekend_vector), weekend_1, "*") 
    weekend_2 <- rnorm(num_indvs, vec1_weekend2_mean, vec1_weekend2_sd*0.1)
    mu_2 <- outer(as.numeric(weekend_vector), weekend_2, "*") 
    
    score_front <- c(vec1_score1_mean*0.1, vec1_score1_sd*0.1,vec1_score2_mean, vec1_score2_sd*0.1)
  } else if(2 == setting)
  {
    weekend_1 <- rnorm(num_indvs, vec2_weekend1_mean*0.025, vec2_weekend1_sd*0.05)
    mu_1 <-  outer(as.numeric(weekend_vector), weekend_1, "*") 
    weekend_2 <- rnorm(num_indvs, vec2_weekend2_mean*0.1, vec2_weekend2_sd)
    mu_2 <- outer(as.numeric(weekend_vector), weekend_2, "*") 
    
    score_front <- c(vec2_score1_mean*0.025, vec2_score1_sd*0.05,vec2_score2_mean*0.1, vec2_score2_sd)
  } else if(3 == setting)
  {
    weekend_1 <- rnorm(num_indvs, vec0_weekend1_mean*0.015, vec0_weekend1_sd*0.025)
    mu_1 <-  outer(as.numeric(weekend_vector), weekend_1, "*") 
    weekend_2 <- rnorm(num_indvs, vec0_weekend2_mean*0.05, vec0_weekend2_sd*0.025)
    mu_2 <- outer(as.numeric(weekend_vector), weekend_2, "*") 
    
    score_front <- c(vec0_score1_mean*0.015, vec0_score1_sd*0.025,vec0_score2_mean*0.05, vec0_score2_sd*0.025)
  }
  
  # for(idx in 1:length(score_front))
  # {
  #   all_score_values[idx] <- score_front[idx]
  # }
  # 
  weekend_columns <- cbind(weekend_1, weekend_2)
  return(list("mu_1" = mu_1, "mu_2" = mu_2, "score_vals" = score_front,
              "weekend_columns" = weekend_columns))
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
                                        k = 2,
                                        mu_1,
                                        mu_2,
                                        score_vals,
                                        weekend_columns,
                                        eigenf_func)
{
  #timestamps01 <- seq(from = 0.0001, to = 1, length=timeseries_length)
  
  # noise octaves
  # cat("octave", num_indvs, k, num_indvs * k, "\n")
  scores_standard <- matrix(rnorm(num_indvs * k), ncol = k)
  scores <- scores_standard %*% diag(sqrt(score_vals[c(2,4)]))
  scores <- sweep(scores, 2, score_vals[c(1,3)], "+")
  #
  BIG_mu <- rbind(mu_1, mu_2)
  #BIG_phi <- PsiFunc(k, timestamps01)
  BIG_phi <- eigenf_func
  
  #Z_after <-  t(sqrt(length(timestamps01))*eigen_score $scores%*%ginv(eigen_score$Phi))

  Z <- t(sqrt(timeseries_length)*scores%*%ginv(eigenf_func) )+ BIG_mu
  #Z <- BIG_phi %*% t(scores) + BIG_mu
  Z1 <- Z[1:timeseries_length, ]
  Z2 <- Z[1:timeseries_length + timeseries_length, ]
  
  ###add 1.9,1.5
  # Z1 <- Z1 + (max(Z1)-min(Z1))*1.85
  # Z2 <- Z2 + (max(Z2)-min(Z2))*1.65
  
  
  expZ1 <- exp(Z1)
  expZ2 <- exp(Z2)
  denom <- 1 + expZ1 + expZ2
  p1 <- expZ1 / denom
  p2 <- expZ2 / denom
  p3 <- 1 / denom
  
  # vectorize for future work!!!
  return(list(Z1 = Z1, Z2 = Z2,
              p1 = p1, p2 = p2, p3 = p3,
              MEAN = BIG_mu, PHI = BIG_phi, MFPC = scores , weekend_columns = weekend_columns))
}

#' Psi function
#'
PsiFunc <- function(klen, timestamps01)
{
  psi_k1 <- sapply(c(1:klen), function(i) sin((2 * i + 1) * pi * timestamps01))
  psi_k2 <- sapply(c(1:klen), function(i) cos(2 * i * pi * timestamps01))
  return(rbind(psi_k1, psi_k2))
}


#cluster_number <- c (3743, 43, 16, 28)
#proportion_cluster <- cluster_number/sum(cluster_number)
#0.977284595 0.011227154 0.004177546 0.007310705


# num_indvs = 100
# timeseries_length = 504

# num_indvs = 100
# timeseries_length = 1008

# num_indvs = 100
#timeseries_length = 1512
# timeseries_length = 2016

# num_indvs = 500
# timeseries_length = 504

# num_indvs = 500
# timeseries_length = 1008

# num_indvs = 500
# timeseries_length = 1512
# timeseries_length = 2016

# num_indvs = 1000
# timeseries_length = 504

# num_indvs = 1000
# timeseries_length = 1008

# num_indvs = 1000
#timeseries_length = 1512
# timeseries_length = 2016

scenario = "B"
num_replicas = 1
est_choice = "multinomial"
run_hellinger = TRUE
some_identifier = "test"
temp_folder = temp_folder <- file.path("outputs", "clustersims", paste(scenario, "_", num_replicas, "_", est_choice, "_", some_identifier, sep=""))

RegenIndv <- function(regen_anyway, indv, cluster_allocation, num_indvs, out_categ_func_data_list, out_weekend_columns,
                      out_Z1, out_Z2, out_total_regens,
                      scenario, timeseries_length, 
                      weekend_vector, eigenf_func, Q_vals){
  # if(indv == 75){ browser() }
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
  
  tolcat <- table(out_categ_func_data_list$W[,indv])
  catorder <- order(tolcat, decreasing = TRUE)
  numcat <- length(catorder)
  refcat <- catorder[numcat]
  count_iter <- 0
  while (regen_anyway || 
         (count_iter < 100 && 
         ( 
           (length(tolcat) == 1)
           || (("1" %in% names(table(out_categ_func_data_list$W[,indv])) == FALSE)==TRUE)
         ))
  )
  {
    regen_anyway <- FALSE
    count_iter <- count_iter + 1
    
    #new_cluster_data <- GenerateClusterData(setting_choice, scenario, 3, 5, timeseries_length)
    new_cluster_data <- GenerateClusterData(setting_choice, scenario, 2, 5, timeseries_length,  weekend_vector, eigenf_func)
    out_weekend_columns[indv,] <-  new_cluster_data$weekend_columns[3,]
    
    
    new_prob_curves <- list(p1 = new_cluster_data$p1, p2 = new_cluster_data$p2, p3 = new_cluster_data$p3)
    new_categ_func_data_list <- GenerateCategFuncData( new_prob_curves )
    
    # what is this 3 ?? arbitrarily chosen?
    out_categ_func_data_list$W[, indv] <- new_categ_func_data_list$W[, 3]
    out_Z1[, indv] <- new_cluster_data$Z1[, 3] # latent curves Z1 and Z2
    out_categ_func_data_list$X[indv, , ] <- 0
    out_Z2[, indv] <- new_cluster_data$Z2[, 3]
    
    for (this_time in 1:timeseries_length)
    {
      out_categ_func_data_list$X[indv, this_time, which(Q_vals == out_categ_func_data_list$W[, indv][this_time])] <- 1
    }
    
    tolcat <- table(out_categ_func_data_list$W[, indv])
    catorder <- order(tolcat, decreasing = TRUE)
    numcat <- length(catorder)
    refcat <- catorder[numcat]
  } # end while
  out_total_regens <- out_total_regens + count_iter
  cat("Indivdual", indv, "with regeneration", out_total_regens, "\n")
  
  return_values <- list(categ_func_data_list=out_categ_func_data_list, 
                        weekend_columns=out_weekend_columns,
                        Z1=out_Z1, Z2=out_Z2, total_regens=out_total_regens)
}

ClusterSimulation <- function(num_indvs, timeseries_length,
                              scenario, num_replicas, est_choice, run_hellinger, temp_folder,
                              eigenf_func_input = eigenf_func)
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
  
  timestamps01 <- seq(0.001,0.99,length=(dim(eigenf_func_input)[1]/2))
 
  weekend_vector <- as.factor(c(rep(c(rep(0,5*24*60/20),rep(1,2*24*60/20)),4)))[1:timeseries_length]
  timestamps01_new <- timestamps01
  timeseries_length_new <- length(timestamps01_new)
  timestamps01 <- seq(0.001,0.99,length=timeseries_length)
  eigenf_func_new <- eigenf_func_input
  eigenf_func_sim <- matrix(0,nrow = 2*timeseries_length, ncol=2)
  eigenf_func_sim[1:timeseries_length,1] <- cubicspline(timestamps01_new, eigenf_func_new[1:timeseries_length_new,1],timestamps01 )
  eigenf_func_sim[1:timeseries_length,2] <- cubicspline(timestamps01_new, eigenf_func_new[1:timeseries_length_new,2],timestamps01 )
  eigenf_func_sim[(1+timeseries_length):(timeseries_length*2),1] <- cubicspline(timestamps01_new, 
                                                                            eigenf_func_new[(1+timeseries_length_new):(timeseries_length_new*2),1],
                                                                           timestamps01 )
  eigenf_func_sim[(1+timeseries_length):(timeseries_length*2),2] <- cubicspline(timestamps01_new, 
                                                                            eigenf_func_new[(1+timeseries_length_new):(timeseries_length_new*2),2],
                                                                           timestamps01 )
  
  eigenf_func <- eigenf_func_sim
  
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
    
    # cluster_f1 <- GenerateClusterData(1, scenario, 3, cluster_allocation[1], timeseries_length)
    # cluster_f2 <- GenerateClusterData(2, scenario, 3, cluster_allocation[2], timeseries_length)
    # cluster_f3 <- GenerateClusterData(3, scenario, 3, cluster_allocation[3], timeseries_length)
    
    cluster_f1 <- GenerateClusterData(1, scenario, 2, cluster_allocation[1], timeseries_length, weekend_vector,eigenf_func)
    cluster_f2 <- GenerateClusterData(2, scenario, 2, cluster_allocation[2], timeseries_length,  weekend_vector, eigenf_func)
    cluster_f3 <- GenerateClusterData(3, scenario, 2, cluster_allocation[3], timeseries_length,  weekend_vector, eigenf_func)
    weekend_columns <- rbind(cluster_f1$weekend_columns,
                             cluster_f2$weekend_columns,
                             cluster_f3$weekend_columns)
    
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
    # browser()
    # what is Q vals ? better name???
    Q_vals <- unique(c(categ_func_data_list$W))
    if(is.numeric(Q_vals))
    {
      Q_vals <- sort(Q_vals)
    }
    
    # I need to know what this loop is meant to do !??? maybe there is a better way
    for(indv in 1:num_indvs)
    {
      # if(indv %in% 1:cluster_allocation[1])
      # {
      #   setting_choice <- 1
      # }
      # if(indv %in% (cluster_allocation[1] + 1):(cluster_allocation[1] + cluster_allocation[2]))
      # {
      #   setting_choice <- 2
      # }
      # if(indv %in% (cluster_allocation[1] + cluster_allocation[2] + 1):num_indvs)
      # {
      #   setting_choice <- 3
      # }
      
      # better names for the following variables ???
      # 1. check weather one category only appears 1 time and is it in the end of the timeseries
      # 2. OR is it appearing only one time in the begining
      # 3. OR if the category is less thant Q total category
      # In general, one category only occurs 2 times
      # If timepoints=300, one category only occurs less than 4 times 3/300=0.01
      #If timepoints=750, one category only occurs less than 6 times 5/750=0.0067
      
      # tolcat <- table(categ_func_data_list$W[,indv])
      # catorder <- order(tolcat, decreasing = TRUE)
      # numcat <- length(catorder)
      # refcat <- catorder[numcat]
      # count_iter <- 0
      # while (count_iter < 100 && 
      #        ( 
      #          (length(tolcat) == 1)
      #          || (("1" %in% names(table(categ_func_data_list$W[,indv])) == FALSE)==TRUE)
      #        )
      # )
      # {
      #   count_iter <- count_iter + 1
      #   
      #   #new_cluster_data <- GenerateClusterData(setting_choice, scenario, 3, 5, timeseries_length)
      #   new_cluster_data <- GenerateClusterData(setting_choice, scenario, 2, 5, timeseries_length,  weekend_vector, eigenf_func)
      #   weekend_columns[indv,] <-  new_cluster_data$weekend_columns[3,]
      #   
      #   
      #   new_prob_curves <- list(p1 = new_cluster_data$p1, p2 = new_cluster_data$p2, p3 = new_cluster_data$p3)
      #   new_categ_func_data_list <- GenerateCategFuncData( new_prob_curves )
      #   
      #   # what is this 3 ?? arbitrarily chosen?
      #   categ_func_data_list$W[, indv] <- new_categ_func_data_list$W[, 3]
      #   Z1[, indv] <- new_cluster_data$Z1[, 3] # latent curves Z1 and Z2
      #   categ_func_data_list$X[indv, , ] <- 0
      #   Z2[, indv] <- new_cluster_data$Z2[, 3]
      #   
      #   for (this_time in 1:timeseries_length)
      #   {
      #     categ_func_data_list$X[indv, this_time, which(Q_vals == categ_func_data_list$W[, indv][this_time])] <- 1
      #   }
      #   
      #   tolcat <- table(categ_func_data_list$W[, indv])
      #   catorder <- order(tolcat, decreasing = TRUE)
      #   numcat <- length(catorder)
      #   refcat <- catorder[numcat]
      # } # end while
      # total_regens <- total_regens + count_iter
      # cat("Indivdual", indv, "with regeneration",total_regens, "\n")
      
      regen_indv_vals <- RegenIndv(FALSE, indv, cluster_allocation, num_indvs, categ_func_data_list, weekend_columns,
                            Z1, Z2, total_regens,
                            scenario, timeseries_length, 
                            weekend_vector, eigenf_func, Q_vals)
      categ_func_data_list <- regen_indv_vals$categ_func_data_list
      weekend_columns <- regen_indv_vals$weekend_columns
      Z1 <- regen_indv_vals$Z1
      Z2 <- regen_indv_vals$Z2
      total_regens <- regen_indv_vals$total_regens
      
    } # end for(indv in 1:num_indvs)
    
    #Estimation
    cat("EstimateCategFuncData", replica_idx, "\n")
    #timestamps01 <- seq(from = 0.001, to = 0.99, length=timeseries_length)
    cat("length of time series", length(timestamps01), "\n")
    cat("Estimate Start with W shape", dim(categ_func_data_list$W), "\n")
    timeKeeperStart("X method")
    print(which(categ_func_data_list$W[,52] == 2))
    
    # Build estimate while regenerating
    est_outcome <- EstimateCategFuncData(est_choice, timestamps01, categ_func_data_list$W, indvs_to_eval = 1:num_indvs)
    current_outcome <- est_outcome$partial_outcome
    count_est_iter <- 0
    while(0 < length(est_outcome$bad_indexes) && count_est_iter < 1000){
      # browser()
      # regen all bad indexes
      for(indv in est_outcome$bad_indexes){
        count_est_iter <- count_est_iter + 1
        regen_indv_vals <- RegenIndv(TRUE, indv, cluster_allocation, num_indvs, categ_func_data_list, weekend_columns,
                                     Z1, Z2, total_regens,
                                     scenario, timeseries_length, 
                                     weekend_vector, eigenf_func, Q_vals)
        categ_func_data_list <- regen_indv_vals$categ_func_data_list
        weekend_columns <- regen_indv_vals$weekend_columns
        Z1 <- regen_indv_vals$Z1
        Z2 <- regen_indv_vals$Z2
        total_regens <- regen_indv_vals$total_regens
      }
      print(which(categ_func_data_list$W[,52] == 2))
      
      est_outcome <- EstimateCategFuncData(est_choice, timestamps01, categ_func_data_list$W, indvs_to_eval = est_outcome$bad_indexes)
      
      for (idx in est_outcome$indvs_to_eval) {
        current_outcome[[idx]] <- est_outcome$partial_outcome[[idx]]
      }
      # current_outcome <- est_outcome$partial_outcome
    }
    cat("Count of estimation iterations:", count_est_iter, "\n")
    # browser()
    Z_P_WeekendCoef <- lapply(current_outcome, function(x) x[[3]])
    Z_P_WeekendCoef <- do.call(rbind, Z_P_WeekendCoef)
    #t*n
    category_count <- length(unique(c(categ_func_data_list$W)))
    categFD_est <- list(Z1_est=t(Z_P_WeekendCoef[,1:timeseries_length]),
                            Z2_est=t(Z_P_WeekendCoef[,(1+timeseries_length):(2*timeseries_length)]),
                            p1_est=t(Z_P_WeekendCoef[,(1+timeseries_length*2):(3*timeseries_length)]),
                            p2_est=t(Z_P_WeekendCoef[,(1+timeseries_length*3):(4*timeseries_length)]),
                            p3_est=t(Z_P_WeekendCoef[,(1+timeseries_length*4):(5*timeseries_length)]) ,
                            weekend_vector_coef = Z_P_WeekendCoef[,(1+timeseries_length*(category_count+category_count-1)):dim(Z_P_WeekendCoef)[2]])
    
    
    cat("Estimation finish", dim(categ_func_data_list$W), "\n")
    # browser()
    #####################9/11/2023
    Z_est_curve[[replica_idx]]=array(c(categFD_est$Z1_est,categFD_est$Z2_est),dim=c(timeseries_length,num_indvs,2))
    p_est_curve[[replica_idx]]=array(c(categFD_est$p1_est,categFD_est$p2_est,categFD_est$p3_est),dim=c(timeseries_length,num_indvs,3))
    W_cfd[[replica_idx]]=categ_func_data_list$W
    #####################9/11/2023
    timeKeeperNext()
    #save(categFD_est, file = "ZPW_500_672.RData" )
    #load("ZPW_500_672.RData")
    
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
      # [,1]      [,2]      [,3]
      # [1,] 0.2435779 0.1695561 0.1365889
      # [2,] 0.4602496 0.3658547 0.2503850
      # [3,] 0.2711508 0.2180779 0.2187912
      
    }
    
    # Record clustering performance. USE only Z
    #**** extract_scores_UNIVFPCA is the FASTEST
    
    mfpca_true <-extract_scores_UNIVFPCA(mZ1=Z1, mZ2=Z2,  tt=timestamps01 , PVE=0.95)
    #plot(  scores_true$scores[, 1:2])
    timeKeeperStart("univfpca")
    mfpca_est <- extract_scores_UNIVFPCA(mZ1=categFD_est$Z1_est, mZ2=categFD_est$Z2_est, tt=timestamps01, PVE=0.95)
    timeKeeperNext()
    
    mfpca_true$scores <- cbind(mfpca_true$scores, weekend_columns)
    mfpca_est$scores <-  cbind(mfpca_est$scores, categFD_est$weekend_vector_coef)
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
    # true_dbscan_temp <- dbscan_cluster(data=mfpca_true$scores,1)$label
    # timeKeeperStart("dbscan")
    # est_dbscan_temp <- dbscan_cluster(data=mfpca_est$scores,1)$label
    
    true_dbscan_temp <- dbscan_cluster(data=mfpca_true$scores,4)$label
    timeKeeperStart("dbscan")
    est_dbscan_temp <- dbscan_cluster(data=mfpca_est$scores,4)$label
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
    #est_dbscan_temp_cfda <- dbscan_cluster(data=cfd_scores,1)$label
    est_dbscan_temp_cfda <- dbscan_cluster(data=cfd_scores,4)$label
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


RunExperiment <- function(scenario, num_replicas, est_choice, some_identifier="noid", eigenf_func_input = eigenf_func)
{
  temp_folder <- file.path("outputs", "clustersims", paste(scenario, "_", num_replicas, "_", est_choice, "_", some_identifier, sep=""))
  # Empty the directory if it exists
  if(dir.exists(temp_folder)){
    unlink(temp_folder, recursive = TRUE)
  }
  dir.create(temp_folder)
  print(temp_folder)
  
  
  # n100t300C <- ClusterSimulation(100,672,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
  # n100t750C <- ClusterSimulation(100,1344,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
  # n100t2000C <- ClusterSimulation(100,2016,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
  # 
  # 
  # n500t300C <- ClusterSimulation(500,672,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
  # n500t750C <- ClusterSimulation(500,1344,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
  # n500t2000C <- ClusterSimulation(500,2016,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
  # 
  # 
  # n1000t300C <- ClusterSimulation(1000,672,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
  # n1000t750C <- ClusterSimulation(1000,1344,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
  # n1000t2000C <- ClusterSimulation(1000,2016,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
  # 
  
  # n100t300C <- ClusterSimulation(100,504,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
  # n100t750C <- ClusterSimulation(100,1008,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
  # n100t2000C <- ClusterSimulation(100,1512,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
  
  
  # n500t300C <- ClusterSimulation(500,504,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
  # n500t750C <- ClusterSimulation(500,1008,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
  # n500t2000C <- ClusterSimulation(500,1512,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
  
  
  n1000t300C <- ClusterSimulation(1000,504,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
  n1000t750C <- ClusterSimulation(1000,1008,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
  n1000t2000C <- ClusterSimulation(1000,1512,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
  
  
  
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

# EXECUTION:

# }) # profvis end



# set.seed(123)
# A_2_probit <- RunExperiment("A",2,"probit","test")
# 
# set.seed(123)
# A_2_binomial <- RunExperiment("A",2,"binomial","test")


# set.seed(123)
# B_2_multinomial <- RunExperiment("B", options_replicas,"multinomial","paper1")
# save(B_2_multinomial, file = "Hazel_mul_B2.RData")
# save(B_2_multinomial,file=file.path("outputs", paste(options_jobid, options_replicas,"Hazel_mul_B2.RData",sep="_")))

set.seed(123)
n1000t300C <- ClusterSimulation(1000,2016,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
save(n1000t300C , file = "n1000t300C.RData")
n1000t750C <- ClusterSimulation(2000,2016,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
save(n1000t750C , file = "n1000t750C.RData")
n1000t2000C <- ClusterSimulation(3000,2016,scenario,num_replicas,est_choice,TRUE,temp_folder, eigenf_func_input = eigenf_func)
save(n1000t2000C , file = "n1000t2000C.RData")

#})

#htmlwidgets::saveWidget(profiling_result, "profiling_result.html")

#save(C_2_probit,file="C_2_probit.RData")
# set.seed(123)
# A_100_probit <- RunExperiment("A",100,"probit")
# 
# 



if(run_parallel)
{
  parallel::stopCluster(cl = my.cluster)
  initialized_parallel <- FALSE
}

#n=100 t=672
# Xiaoxia calc time taken: Time difference of 14.96065 secs 
# count_iter:  1 
# count_iter:  2 
# count_iter:  3 
# count_iter:  4 
# count_iter:  5 
# univfpca calc time taken: Time difference of 0.01541018 secs 
# kmeans calc time taken: Time difference of 0.01435995 secs 
# fadp calc time taken: Time difference of 0.5094869 secs 
# dbscan calc time taken: Time difference of 0.001412153 secs 
# Parellel Registered:  TRUE 
# cfd calc time taken: Time difference of 3.078088 mins 
