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
#library(pracma)

# For: fpca.face
library(refund)

# For: FADPclust
#library(FADPclust)

# For: kNNdist
#library(dbscan)

# For: elbow
# devtools::install_github("ahasverus/elbow")
# library(elbow)

# For: rand.index
#library(fossil)

# For: cfda method
#library(cfda)

# For: gather method
#library(tidyverse)

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
               help="Num CPUs", metavar="NUMCPUS")
  # make_option(c("-r", "--replicas"), type="integer", default=100,
  #             help="Num Replicas", metavar="NUMREPLICAS")
)

#####need for hazel
# Create parser and parse options
parser <- OptionParser(option_list=option_list)
options <- parse_args(parser)

options_jobid <- options$jobid
options_numcpus <- options$numcpus
# options_replicas <- options$replicas
# #options_subjects <- options$subjects

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



#' #' Create directories 
#' if (!dir.exists("outputs")){
#'   dir.create("outputs")
#' }
#' if (!dir.exists("outputs/clustersims")){
#'   dir.create("outputs/clustersims")
#' }


#February 14, 2017 to March 16, 2017, Tuesday, Thursday
#first read the three values and make W
# X_i1_full <- load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/submission_code/X_i1full.RData")
# X_i2_full <- load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/submission_code/X_i2full.RData")
# X_i3_full <- load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/submission_code/X_i3full.RData")

X_i1_full <- load("X_i1full.RData")
X_i2_full <- load("X_i2full.RData")
X_i3_full <- load("X_i3full.RData")

no_tweet_data <- X_i1f
tweet_no_mention_data <- X_i2f
tweet_mention_data <- X_i3f

rm(X_i1f)
rm(X_i2f)
rm(X_i3f)

#estimate z and p
# source("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/acj/gam_weekends.R")
# source("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/acj/time_track_function.R")
source("gam_weekends.R")
source("time_track_function.R")
num_indv <- dim(no_tweet_data)[1]
timestamps01 <- as.numeric(noquote(colnames(no_tweet_data)))
basis_size <- 25 
timeseries_length <- length(timestamps01)
method <- "REML"
category_count <- 3
threshold_probability <- 0.004

timeKeeperStart("Xiaoxia")
zp <- foreach (indv = 1:num_indv, .combine = cbind, .init = NULL, .packages = c("mgcv", "refund")) %dorng%
  {
    #source("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/acj/gam_weekends.R")
    source("gam_weekends.R")
    
    
    # x1<- X[indv,,1]
    # x2<- X[indv,,2]
    # x3<- X[indv,,3]
    
    # x1 <- no_tweet_data[indv,]
    # x2 <- tweet_no_mention_data[indv,]
    # x3 <- tweet_mention_data[indv,]
    #15 minutes, one day 24 hours is 96 data points, 5 days is 480 data points, two weekends days are 192 data points
    # one week is 480+192 = 672 data points
    #two week is 1344 data points
    # three week is 2016 data points
    #four week is 2688 data points
    #maximum 2000 data points is 286 weeks
    #weekend_vector <- c(rep(c(rep(0,480),rep(1,192)),4))[1:length(x1)]
    
    weekend_vector <- c(rep(c(rep(0,480),rep(1,192)),4))[96:(96+length(no_tweet_data)-1)]
    
   
    probit_binom <- function(x_binary){
      if (sum(x_binary)/timeseries_length < threshold_probability){
        gam_result_binary <- RunGam_Day(timestamps01,weekend_vector, x_binary, "probit", basis_size, method)
        p_binary <- gam_result_binary$prob
        p_binary_linpred <- gam_result_binary$linpred
      }else{
        gam_result_binary <- RunGam_Day(timestamps01, weekend_vector, x_binary, "binomial", basis_size, method)
        p_binary <- gam_result_binary$prob
        p_binary_linpred <- gam_result_binary$linpred
      }
      return(list("p_binary"=p_binary,"p_binary_linpred"=p_binary_linpred))
    }
    
    #r_1 <- probit_binom(x1)
    r_1 <- probit_binom(unlist(no_tweet_data[indv,]))
    p1 <- r_1$p_binary
    p1_linpred <- r_1$p_binary_linpred
    
    r_2 <- probit_binom(unlist(tweet_no_mention_data[indv,]))
    p2 <- r_2$p_binary
    p2_linpred <- r_2$p_binary_linpred
    
    r_3 <- probit_binom(unlist(tweet_mention_data[indv,]))
    p3 <- r_3$p_binary
    p3_linpred <- r_3$p_binary_linpred
    
    # gam_result_1 <- RunGam_Day(timestamps01,  weekend_vector, unlist(no_tweet_data[indv,]), "binomial", basis_size, method)
    # p1 <- gam_result_1$prob
    # p1_linpred <- gam_result_1$linpred
    # 
    # gam_result_2 <- RunGam_Day(timestamps01,  weekend_vector, unlist(tweet_no_mention_data[indv,]), "binomial", basis_size, method)
    # p2 <- gam_result_2$prob
    # p2_linpred <- gam_result_2$linpred
    # 
    # gam_result_3 <- RunGam_Day(timestamps01,  weekend_vector, unlist(tweet_mention_data[indv,]), "binomial", basis_size, method)
    # p3 <- gam_result_3$prob
    # p3_linpred <- gam_result_3$linpred
    
    # estimate the latent tranjecotries Z
    denominator_p <- 1 + exp(p3_linpred)
    z1 <- (p1_linpred-p3_linpred)-log( (1+exp(p1_linpred))/(denominator_p))
    z2 <- (p2_linpred-p3_linpred)-log( (1+exp(p2_linpred))/(denominator_p))
    
    psum <- p1 + p2 + p3
    return(c(c(z1,z2), cbind(p1/psum, p2/psum, p3/psum)))
    #return(c(c(z1,z2), cbind(p1/(p1 + p2 + p3), p2/(p1 + p2 + p3), p3/(p1 + p2 + p3))))
  }
# Unravel the two variables from zp
z_rows_count <- timeseries_length * 2
Z <- array(zp[1:z_rows_count, ], c(z_rows_count, num_indv))
p <- array(t(matrix(zp[(z_rows_count + 1):dim(zp)[1], ], ncol=num_indv)), c(num_indv, timeseries_length, category_count))
timeKeeperNext()  
Z1_est <- Z[1:timeseries_length,]
Z2_est <- Z[1:timeseries_length+timeseries_length,]

p1_est <- t(p[,,1])
p2_est <- t(p[,,2])
p3_est <- t(p[,,3])

#save(Z1_est, Z2_est, p1_est, p2_est, p3_est, timestamps01 , index_remove, file = "Twitter_Est.RData")
save(Z1_est, Z2_est, p1_est, p2_est, p3_est, timestamps01 ,  file = "Twitter_Est.RData")
load("Twitter_Est.RData")

save(Z1_est, Z2_est, p1_est, p2_est, p3_est, timestamps01 ,  file = "Twitter_Est_probit_binom.RData")

load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/acj/Twitter_Est.RData")

#find columns that has na vlues
cols_with_na <- apply(Z1_est, 2, function(x) any(is.na(x)))
which_columns <- which(cols_with_na)  # Get column indices
print(which_columns)
#probit 1942 2309 3033 3274 3469

cols_with_inf <- apply(Z1_est, 2, function(x) any((x==-Inf)))
which_columns_inf <- which(cols_with_inf)  # Get column indices
print(which_columns_inf)
# 1942 2309 2365 2721 3033 3274 3469
#probit 2309 2365 2721

cols_with_inf1 <- apply(Z1_est, 2, function(x) any((x==Inf)))
which_columns1 <- which(cols_with_inf1)  # Get column indices
print(which_columns1)
#probit 3511 of the individuals have inf

cols_with_inf2 <- apply(Z2_est, 2, function(x) any((x==-Inf)))
which_columns2 <- which(cols_with_inf2)  # Get column indices
print(which_columns2)
# probit [1]   63  622  637 2282 2316 2389 2496 2572 2588 2608 2743 2771 3007 3177 3436 3595 3643 3649 3656 3670 3680
#[22] 3698 3702 3704 3722 3723 3736 3765 3766 3773 3782 3789 3797 3801 3804 3805 3806 3807 3815 3820 3838 3839
#[43] 3845 3846 3849 3850 3860 3861 3862 3867 3871

cols_with_inf22 <- apply(Z2_est, 2, function(x) any((x==Inf)))
which_columns22 <- which(cols_with_inf22)  # Get column indices
print(which_columns22)
#probit 3481 individuals have Inf


par(mfrow = c(1,2))
matplot(timestamps01,Z1_est[,-c(which_columns)],type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "Z1")
matplot(timestamps01,Z2_est[,-c(which_columns)],type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "Z2")

cols_with_large <- apply(Z1_est, 2, function(x) any((x>200)))
which_columns_large <- which(cols_with_large)  # Get column indices
print(which_columns_large)
#380

cols_with_largerz2 <- apply(Z2_est, 2, function(x) any((x>250 )))
which_columnsz2 <- which(cols_with_largerz2)  # Get column indices
print(which_columnsz2)
#380

cols_with_largerz22 <- apply(Z2_est, 2, function(x) any((x< -400 )))
which_columnsz22 <- which(cols_with_largerz22)  # Get column indices
print(which_columnsz22)

save(timestamps01, Z1_est, Z2_est,which_columns, which_columns_large, which_columnsz22, file ="twitter_est_extreme.RData")
par(mfrow = c(1,2))
matplot(timestamps01,Z1_est[,-c(which_columns, which_columns_large, which_columnsz22)],type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "Z1")
matplot(timestamps01,Z2_est[,-c(which_columns, which_columns_large, which_columnsz22)],type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "Z2")

par(mfrow = c(1,2))
matplot(t,Z1_est[,-c(which_columns, which_columns_large, which_columnsz22)],type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "Z1")
matplot(t,Z2_est[,-c(which_columns, which_columns_large, which_columnsz22)],type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "Z2")

par(mfrow = c(1,3))
matplot(t,p1_est[,-c(which_columns, which_columns_large, which_columnsz22)],type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "p1",ylim=c(0,1))
matplot(t,p2_est[,-c(which_columns, which_columns_large, which_columnsz22)],type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "p2",ylim=c(0,1))
matplot(t,p3_est[,-c(which_columns, which_columns_large, which_columnsz22)],type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "p3",ylim=c(0,1))



par(mfrow = c(1,3))
matplot(t,logit_p(p1_est[,-c(which_columns, which_columns_large, which_columnsz22)]),type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "logit p1")
matplot(t,logit_p(p2_est[,-c(which_columns, which_columns_large, which_columnsz22)]),type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "logit p2")
matplot(t,logit_p(p3_est[,-c(which_columns, which_columns_large, which_columnsz22)]),type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "logit p3")

save(timestamps01, Z1_est, Z2_est,which_columns, which_columns_large, which_columnsz22, file ="twitter_est_extreme_t_final.RData")
#################################
#generate new Z
extract_scores_UNIVFPCA <- function (mZ1,mZ2, tt , PVE=0.95)
{
  # mZ1 <- Z1_est[,-c(index_remove)]
  # mZ2 <- Z2_est[,-c(index_remove)]
  
  # mZ1 <- Z1_est
  # mZ2 <- Z2_est
  # tt <- timestamps01
  
  # mZ1 <- Z1_est[,-c(index_remove)][,-c(1,2718)][,-c(540,2717)][,-2715]
  # mZ2 <- Z2_est[,-c(index_remove)][,-c(1,2718)][,-c(540,2717)][,-2715]
  
  # mZ1 <- Z1_est[,-c(which_columns)]
  # mZ2 <-  Z2_est[,-c(which_columns)]
  m<- nrow(mZ1)
  n<-ncol(mZ1)
  
  out1 <- refund::fpca.face(Y=t(mZ1), argvals =tt, pve = 0.99)
  out2 <- refund::fpca.face(Y=t(mZ2), argvals =tt, pve = 0.99)
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
#simulation
load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/time_stamp.RData")
timestamps01  <- t
#3865, 5 score dimention
#eigen_score <- extract_scores_UNIVFPCA (Z1_est[,-c(which_columns)],Z2_est[,-c(which_columns)], timestamps01 , PVE=0.95)
eigen_score <- extract_scores_UNIVFPCA (Z1_est[,-c(which_columns, which_columns_large, which_columnsz22)],
                                        Z2_est[,-c(which_columns, which_columns_large, which_columnsz22)], 
                                        t , PVE=0.95)

save(eigen_score, t, file = "Twitter_eigen.RData")

cumsum(c(sd(eigen_score $scores[,1]),sd(eigen_score $scores[,2]),sd(eigen_score $scores[,3]),
         sd(eigen_score $scores[,4]),sd(eigen_score $scores[,5])))/sum(c(sd(eigen_score $scores[,1]),sd(eigen_score $scores[,2]),sd(eigen_score $scores[,3]),
                                                                       sd(eigen_score $scores[,4]),sd(eigen_score $scores[,5])))
#[1] 0.3888913 0.6441416 0.7973338 0.9027171 1.0000000
#n * 672*2=1344, n*t
library(MASS)
Z_after <-  t(sqrt(length(t))*eigen_score $scores%*%ginv(eigen_score$Phi))


# load("/Users/xzhao17/Documents/project1mac/Twitterhapf.RData")
# 
par(mfrow = c(1,2))
matplot(timestamps01,Z1_est[,-c(which_columns, which_columns_large, which_columnsz22)]-Z_after[1:2000,],type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "Z1")
#recover

matplot(timestamps01,Z2_est[,-c(which_columns, which_columns_large, which_columnsz22)]-Z_after[2001:4000,],type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "Z2")


cols_with_larger_diff <- apply(Z1_est[,-c(which_columns, which_columns_large, which_columnsz22)]-Z_after[1:2000,], 
                               2, function(x) any((x< -90 )))
which_columns_diff <- which(cols_with_larger_diff)  # Get column indices
print(which_columns_diff)
#541
cols_with_larger_diff2 <- apply(Z2_est[,-c(which_columns, which_columns_large, which_columnsz22)]-Z_after[2001:4000,], 
                               2, function(x) any((x< -100 )))
which_columns_diff2 <- which(cols_with_larger_diff2)  # Get column indices
print(which_columns_diff2)
#2819

##########################
par(mfrow = c(1,2))
matplot(timestamps01,Z1_est[,c(which_columns_diff,which_columns_diff2)],type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "Z1")
#recover

matplot(timestamps01,Z2_est[,c(which_columns_diff,which_columns_diff2)],type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "Z2")

###########################



#Z_after: 2t*n
save(eigen_score, Z_after, which_columns_diff, which_columns_diff2, Z1_est,Z2_est,
     t, which_columns, which_columns_large,which_columnsz22, file = "Twitter_simulation_data.RData" )

eigen_score_final <- extract_scores_UNIVFPCA (Z1_est[,-c(which_columns, which_columns_large, which_columnsz22, 
                                                   which_columns_diff, which_columns_diff2)],
                                        Z2_est[,-c(which_columns, which_columns_large, which_columnsz22,
                                                   which_columns_diff, which_columns_diff2)], 
                                        t , PVE=0.95)
Z_after_final <-  t(sqrt(length(t))*eigen_score_final$scores%*%ginv(eigen_score_final$Phi))
save(eigen_score_final, Z_after_final, which_columns_diff, which_columns_diff2, Z1_est,Z2_est,
     t, which_columns, which_columns_large,which_columnsz22,which_columns_diff, which_columns_diff2,
     file = "Twitter_simulation_data_final.RData" )


par(mfrow = c(1,2))
matplot(t,Z1_est[,-c(which_columns, which_columns_large, which_columnsz22,which_columns_diff, which_columns_diff2)]-Z_after_final[1:2000,],type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "Z1")
#recover

matplot(t,Z2_est[,-c(which_columns, which_columns_large, which_columnsz22,which_columns_diff, which_columns_diff2)]-Z_after_final[2001:4000,],type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "Z2")
####################################

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


library(ggplot2)
library(elbow)
library(dbscan)
tclusterdata=data.frame(eigen_score_final$scores)
est_dbscan <- dbscan_cluster(data=eigen_score_final $scores,3.8)$label
tclusterdata$Cluster=as.factor(est_dbscan)
colnames(tclusterdata)[1:2] =c("ksi1","ksi2")
tps2 <- ggplot(tclusterdata,aes(ksi1,ksi2,colour = Cluster)) + geom_point(aes(shape=Cluster),size=3)+ggtitle(paste0("Cluster Results",'\n',"(",dim(tclusterdata)[1]," Subjects",")")) +
  xlab(expression('Score '* widehat(xi[i1]))) + ylab(expression('Score '* widehat(xi[i2])))+ theme(plot.title = element_text(hjust = 0.5))+#theme(text = element_text(size = 20))+
  #theme(legend.text = element_text(size = 30))  
  #theme(axis.text=element_text(size = 20),axis.title = element_text(size =30))
  theme(text=element_text(size = 20))
tps2
table(est_dbscan)
est_dbscan
# 0    1    2    3 
# 43 3792   11   17 
dim(tclusterdata)[1]
#3863 



combinedscore_z <- eigen_score_final$scores
dimz=dim(combinedscore_z)[2] #find how many features data has
if (dimz<=2){minPts=4} #this is what user defined. if you could automatically find a good one.that might solve all my problem. nOW I use default value from literature minPts=4 for 2 features and 2*features+1 for data has more than 2 features
if (dimz>2){minPts=2*dimz+1}
dist=kNNdist(combinedscore_z, k =minPts-1) #step 1 sort the distance. input 2D array you want to cluster, and then minPts-1
#ninty5p=quantile(dist, probs = pct)
scaleep=3.5  #this set to be 1, that is the scale you multiple the optimal epsilon, I set it to 12.84 just because of my twitter data is very noisy. this is the how much you want to multiple the optimal radius. increase the circle basically. instead of the optimal epsilon circle. I increase the circle so my cluster is meaingful otherwise if smalla circle, too many clusters are created
#########change to max increase
distdataelbow=data.frame(sort(dist)) #sort the distance
distdataelbow$index=1:(dim(combinedscore_z)[1])
ipoint <- elbow(data = distdataelbow)
epsoptimal=(ipoint$sort.dist._selected)*scaleep #these three lines find the elbow 

resnospike <- dbscan(eigen_score_final$scores, eps =epsoptimal , minPts = 20)

tclusterdatanospke=data.frame(combinedscore_z[,])
tclusterdatanospke$Cluster=as.factor(resnospike$cluster)
library(ggplot2)
tpsnospike <- ggplot(tclusterdatanospke,aes(X1,X2,colour = Cluster)) + geom_point(aes(shape=Cluster),size=3)+ggtitle(paste0("Twitter Users Cluster Results",'\n',"(",dim(tclusterdatanospke)[1]," Subjects",")")) +
  #xlab(expression('Score '* widehat(xi[i1]))) + ylab(expression('Score '* widehat(xi[i2])))
  xlab('MFPCA Score1 ') + ylab('MFPCA Score2 ')+ theme(plot.title = element_text(hjust = 0.5))+#theme(text = element_text(size = 20))+
  #theme(legend.text = element_text(size = 30))  
  #theme(axis.text=element_text(size = 20),axis.title = element_text(size =30))
  theme(text=element_text(size = 20))+theme(legend.position = c(0.92,0.8))
tpsnospike

table(tclusterdatanospke$Cluster)
# 0    1    2    3 
# 96 3533  211   23 
####################################
#plot
vec0=c(which(tclusterdatanospke$Cluster==0))
vec1=c(which(tclusterdatanospke$Cluster==1))
vec2=c(which(tclusterdatanospke$Cluster==2))
vec3=c(which(tclusterdatanospke$Cluster==3))

#2000*3683
Z1_figure <- Z1_est[,-c(which_columns, which_columns_large, which_columnsz22, 
                        which_columns_diff, which_columns_diff2)]

Z2_figure <- Z2_est[,-c(which_columns, which_columns_large, which_columnsz22, 
                        which_columns_diff, which_columns_diff2)]

p1_figure <- p1_est[,-c(which_columns, which_columns_large, which_columnsz22, 
                        which_columns_diff, which_columns_diff2)]
p2_figure <- p2_est[,-c(which_columns, which_columns_large, which_columnsz22, 
                        which_columns_diff, which_columns_diff2)]
p3_figure <- p3_est[,-c(which_columns, which_columns_large, which_columnsz22, 
                        which_columns_diff, which_columns_diff2)]
argval <- t

save(Z1_figure, Z2_figure,
     p1_figure, p2_figure,
     p3_figure, argval, t,
     vec0, vec1, vec2, vec3,
     tclusterdatanospke,z_min,z_max,W_matrix,
     file = "Twitter_new_figure_label.RData")
# min(Z1_figure[,vec3])
# [1] 2.237043
# > max(Z1_figure[,vec3])
# [1] 11.3504
load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/Twitter_new_figure_label.RData")
###########################################################################################
#graph W by cluster

#W is n * t
W_matrix_after <- W_matrix[-c(which_columns, which_columns_large, which_columnsz22, 
                               which_columns_diff, which_columns_diff2),]
#sample one user per cluster
user_index <- c(sample(vec1,1),sample(vec2,1),sample(vec3,1))

unname(W_matrix_after[user_index[1],start_number:end_number])
unname(W_matrix_after[user_index[2],start_number:end_number])
unname(W_matrix_after[user_index[3],start_number:end_number])

save(user_index, W_matrix_after, start_number, end_number, t,
     file = "User_Cluster_W.RData")

#################graph for each user
graph_user_category <- function(index_number){
  user_data <- data.frame(
    Time = t[start_number:end_number],  # 10 days
    Category = as.factor(W_matrix_after[user_index[index_number],start_number:end_number] ) # Random categories
  )
  library(dplyr)
  # Rename categories using recode
  user_data <- user_data %>%
    mutate(Category = recode(Category, "1" = "No Tweet",
                             "2" = "Tweet Non Ref Brand",
                             "3" = "Tweet on Ref Brand"))
  user_fig <- ggplot(user_data, aes(x = Time, y = Category)) +
    geom_tile(aes(fill = Category), color = "white") +
    scale_fill_manual(values = c("No Tweet" = "red", 
                                 "Tweet Non Ref Brand" = "green",
                                 "Tweet on Ref Brand" = "blue")) +
    labs(title = "", x = "", y = "") +
    theme_minimal()+
    scale_x_continuous(
      breaks = c(min(t[start_number:end_number]), max(t[start_number:end_number])),          # Specify positions for the labels
      labels = c("March 11, 2017","March 12, 2017")  # Specify corresponding labels
    ) +
    theme(legend.position = "none",
          text=element_text(size = 15))
  print( user_fig)
}

graph_user_category(1)
graph_user_category(2)
graph_user_category(3)


user_data <- data.frame(
  Time = t[start_number:end_number],  # 10 days
  Category = as.factor(W_matrix_after[user_index[3],start_number:end_number] ) # Random categories
)
library(dplyr)
# Rename categories using recode
user_data <- user_data %>%
  mutate(Category = recode(Category, "1" = "No Tweet",
                           "2" = "Tweet Non Ref Brand",
                           "3" = "Tweet on Ref Brand"))
ggplot(user_data, aes(x = Time, y = Category)) +
  geom_tile(aes(fill = Category), color = "white") +
  scale_fill_manual(values = c("No Tweet" = "red", 
                               "Tweet Non Ref Brand" = "green",
                               "Tweet on Ref Brand" = "blue")) +
  labs(title = "", x = "", y = "") +
  theme_minimal()+
  scale_x_continuous(
    breaks = c(min(t[start_number:end_number]), max(t[start_number:end_number])),          # Specify positions for the labels
    labels = c("March 11, 2017","March 12, 2017")  # Specify corresponding labels
  ) +
 # theme(legend.position = "none",
        #text=element_text(size = 15))
theme(legend.position = "bottom",
      text=element_text(size = 20))+
guides(color = guide_legend(ncol = 2))

#one day 72 data points
#pch=15, square, pch=16 filled circle, pch=17, triangle, pch=18 is diamond
# Assign symbols based on y values

 
 # pch_values <- ifelse(y == 1, 18,  # Square for y <= 3
 #                      +  ifelse(y ==2 , 4,  # Filled Circle for 3 < y <= 6
 #                                +   8))  # Triangle for y > 6
# # Plot using matplot
# par(mfrow=c(1,3))
# # start_number <- 10
# # end_number <- start_number + 10
# plot(t[start_number:end_number], t(W_matrix_after[user_index[1],start_number:end_number]), pch = pch_values, type = "p", col = 1:3, lty = 1,
#         main = "", xlab = "", ylab = "",xaxt = "n",ylim=c(0.75,1.25),
#      yaxt = "n")
# # Add custom y-axis with only 3 values
# # axis(2, at = c(1, 2, 3),  #cex.lab = 1.5,cex.axis = 2,cex.main=2,
# #      labels = c("No Tweet", "Tweet Non Ref Brand", "Tweet on Ref Brand"))
# axis(2, at = c(1),  #cex.lab = 1.5,cex.axis = 2,cex.main=2,
#      labels = c(""))
# 
# axis(1,                         # Define x-axis manually
#      #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
#      at = c(min(t[start_number:end_number]), max(t[start_number:end_number])),
#      #cex.lab = 1.5,cex.axis = 2,cex.main=2,
#      labels = c("March 11, 2017","March 12, 2017"))
# 
# # Add a legend
# #legend("topright", legend = c("No Tweet", "Tweet Non Ref Brand", "Tweet on Ref Brand"),
# #       pch = pch_values, col = c(1,1,1))
# plot(t[start_number:end_number], t(W_matrix_after[user_index[2],start_number:end_number]), pch = pch_values, type = "p", col = 1:3, lty = 1,
#      main = "", xlab = "", ylab = "",xaxt = "n",ylim=c(0.75,1.25),
#      yaxt = "n")
# # axis(2, at = c(1, 2, 3),  #cex.lab = 1.5,cex.axis = 2,cex.main=2,
# #      labels = c("No Tweet", "Tweet Non Ref Brand", "Tweet on Ref Brand"))
# axis(2, at = c(1),  #cex.lab = 1.5,cex.axis = 2,cex.main=2,
#      labels = c(""))
# 
# axis(1,                         # Define x-axis manually
#      #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
#      at = c(min(t[start_number:end_number]), max(t[start_number:end_number])),
#      #cex.lab = 1.5,cex.axis = 2,cex.main=2,
#      labels = c("March 11, 2017","March 12, 2017"))
# 
# plot(t[start_number:end_number], t(W_matrix_after[user_index[3],start_number:end_number]), pch = pch_values, type = "p", col = 1:3, lty = 1,
#      main = "", xlab = "", ylab = "",xaxt = "n",ylim=c(0.75,1.25),
#      yaxt = "n")
# axis(2, at = c(1),  #cex.lab = 1.5,cex.axis = 2,cex.main=2,
#      labels = c(""))
# 
# axis(1,                         # Define x-axis manually
#      #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
#      at = c(min(t[start_number:end_number]), max(t[start_number:end_number])),
#      #cex.lab = 1.5,cex.axis = 2,cex.main=2,
#      labels = c("March 11, 2017","March 12, 2017"))

#############
#how to graph category time series
# Example categorical time series data
# data <- data.frame(
#   Time = as.Date("2024-01-01") + 0:9,  # 10 days
#   Category = sample(c("A", "B", "C"), 10, replace = TRUE)  # Random categories
# )
# print(data)
# library(ggplot2)
# 
# ggplot(data, aes(x = Time, y = Category, group = 1)) +
#   geom_line(aes(color = Category)) +
#   geom_point(aes(shape = Category), size = 3) +
#   scale_y_discrete(limits = c("A", "B", "C")) +  # Ensure consistent order
#   labs(title = "Categorical Time Series", x = "Time", y = "Category") +
#   theme_minimal()
# ggplot(data, aes(x = Time, y = Category)) +
#   geom_tile(aes(fill = Category), color = "white") +
#   scale_fill_manual(values = c("A" = "red", "B" = "blue", "C" = "green")) +
#   labs(title = "Categorical Heatmap Over Time", x = "Time", y = "Category") +
#   theme_minimal()
# #only show xaxis at certain point with labels
# # Plot with specific x-axis labels
# ggplot(data, aes(x = x, y = y)) +
#   geom_line() +
#   scale_x_continuous(
#     breaks = c(2, 5, 8),          # Specify positions for the labels
#     labels = c("Two", "Five", "Eight")  # Specify corresponding labels
#   ) +
#   labs(title = "Custom X-axis Labels", x = "Custom Points", y = "Value") +
#   theme_minimal()
###########################################################################################


load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/Twitter_new_figure_label.RData")
# max(Z1_figure)
# [1] 180.3087
# > min(Z1_figure)
# [1] -0.3315023
# > max(Z2_figure)
# [1] 174.6473
# > min(Z2_figure)
# [1] -261.4091
logit_p <- function(p){log(p / (1 - p))}

p1_figure <- logit_p(p1_est[,-c(which_columns, which_columns_large, which_columnsz22, 
                        which_columns_diff, which_columns_diff2)])
p2_figure <- logit_p(p2_est[,-c(which_columns, which_columns_large, which_columnsz22, 
                        which_columns_diff, which_columns_diff2)])
p3_figure <- logit_p(p3_est[,-c(which_columns, which_columns_large, which_columnsz22, 
                        which_columns_diff, which_columns_diff2)])

z_min <- min(c(Z1_figure,Z2_figure))
z_max <- max(c(Z1_figure,Z2_figure))
p_min <- min (c (p1_figure, p2_figure, p3_figure))
p_max <- max (c (p1_figure, p2_figure, p3_figure))

save(t,argval,Z1_figure, Z2_figure,
     p1_figure, p2_figure,
     p3_figure, argval, t,
     vec0, vec1, vec2, vec3,
     tclusterdatanospke,z_min,
     z_max,p_min,p_max,
     file = "Twiiter_figure_logit.RData")

setwd("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd")
load("Twiiter_figure_logit.RData")
#par(mfrow=c(4,5))
par(mfrow=c(1,5))
##vec1
meanz=colMeans(t(Z1_figure)[vec1,])
matplot(t, t(t(Z1_figure)[vec1,]),
        type='l', lty=1, col="light grey",
       
        main=expression(hat("Z")^1*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min,z_max))
lines(argval,meanz ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 

meanz2=colMeans(t(Z2_figure)[vec1,])
matplot(t, t(t(Z2_figure)[vec1,]),
        type='l', lty=1, col="light grey",
        
        main=expression(hat("Z")^2*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min,z_max))
lines(argval,meanz2 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 

#######################################

meanp=colMeans(t(p1_figure)[vec1,])
matplot(t, t(t(p1_figure)[vec1,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        main=expression(hat("p")^1*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min,p_max))
lines(argval,meanp ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 

meanp2=colMeans(t(p2_figure)[vec1,])
matplot(t, t(t(p2_figure)[vec1,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        main=expression(hat("p")^2*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min,p_max))
lines(argval,meanp2 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 


meanp3=colMeans(t(p3_figure)[vec1,])
matplot(t,t(t(p3_figure)[vec1,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        main=expression(hat("p")^3*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min,p_max))
lines(argval,meanp3 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 
##vec2
argval=t
meanz=colMeans(t(Z1_figure)[vec2,])
matplot(t, t(t(Z1_figure)[vec2,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        #main=expression(hat("Z")^1*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min,z_max))
lines(argval,meanz ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 

meanz2=colMeans(t(Z2_figure)[vec2,])
matplot(t, t(t(Z2_figure)[vec2,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        #main=expression(hat("Z")^2*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min,z_max))
lines(argval,meanz2 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 

#######################################

meanp=colMeans(t(p1_figure)[vec2,])
matplot(t, t(t(p1_figure)[vec2,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        #main=expression(hat("p")^1*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min,p_max))
lines(argval,meanp ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 

meanp2=colMeans(t(p2_figure)[vec2,])
matplot(t, t(t(p2_figure)[vec2,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        #main=expression(hat("p")^2*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min,p_max))
lines(argval,meanp2 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 


meanp3=colMeans(t(p3_figure)[vec2,])
matplot(t, t(t(p3_figure)[vec2,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        #main=expression(hat("p")^3*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min,p_max))
lines(argval,meanp3 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 
#vec3
meanz=colMeans(t(Z1_figure)[vec3,])
matplot(t, t(t(Z1_figure)[vec3,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        #main=expression(hat("Z")^1*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min,z_max))
lines(argval,meanz ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 

meanz2=colMeans(t(Z2_figure)[vec3,])
matplot(t, t(t(Z2_figure)[vec3,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
       # main=expression(hat("Z")^2*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min,z_max))
lines(argval,meanz2 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 

#######################################
meanp=colMeans(t(p1_figure)[vec3,])
matplot(t, t(t(p1_figure)[vec3,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        #main=expression(hat("p")^1*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min,p_max))
lines(argval,meanp ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 

meanp2=colMeans(t(p2_figure)[vec3,])
matplot(t, t(t(p2_figure)[vec3,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        #main=expression(hat("p")^2*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min,p_max))
lines(argval,meanp2 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 


meanp3=colMeans(t(p3_figure)[vec3,])
matplot(t, t(t(p3_figure)[vec3,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        #main=expression(hat("p")^3*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min,p_max))
lines(argval,meanp3 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 2,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 
#vec0
meanz=colMeans(t(Z1_figure)[vec0,])
matplot(t, t(t(Z1_figure)[vec0,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        #main=expression(hat("Z")^1*(t)),
        xlab="Day", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min,z_max))
lines(argval,meanz ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 

meanz2=colMeans(t(Z2_figure)[vec0,])
matplot(t, t(t(Z2_figure)[vec0,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        # cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(-12,3),main=expression(hat("Z")^2*(t)),cex.main=2          main=mtext(bquote(" ")),
        xlab="Day",
        #main=expression(hat("Z")^2*(t)), 
        ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min,z_max))
lines(argval,meanz2 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 

#######################################
meanp=colMeans(t(p1_figure)[vec0,])
matplot(t, t(t(p1_figure)[vec0,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        #main=expression(hat("p")^1*(t)),
        xlab="Day", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min,p_max))
lines(argval,meanp ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 

meanp2=colMeans(t(p2_figure)[vec0,])
matplot(t, t(t(p2_figure)[vec0,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        #main=expression(hat("p")^2*(t)),
        xlab="Day", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min,p_max))
lines(argval,meanp2 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 


meanp3=colMeans(t(p3_figure)[vec0,])
matplot(t, t(t(p3_figure)[vec0,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        #main=expression(hat("p")^3*(t)),
        xlab="Day", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min,p_max))
lines(argval,meanp3 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 

###############
#plot just p3 for all clusters
par(mfrow=c(1,4))
#cluster 1, cluster 2, cluster 3, noise 
#vec1
meanp3=colMeans(t(p3_figure)[vec1,])
matplot(t, t(t(p3_figure)[vec1,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        
        xlab="Day", ylab="",cex.lab = 2,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(min(t(p3_figure)[vec1,]),max(t(p3_figure)[vec1,])),main="Cluster One")
lines(argval,meanp3 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 
#vec2
meanp3=colMeans(t(p3_figure)[vec2,])
matplot(t, t(t(p3_figure)[vec2,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        
        xlab="Day", ylab="",cex.lab = 2,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(min(t(p3_figure)[vec2,]),max(t(p3_figure)[vec2,])),main="Cluster Two")
lines(argval,meanp3 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 
#vec3
meanp3=colMeans(t(p3_figure)[vec3,])
matplot(t, t(t(p3_figure)[vec3,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        
        xlab="Day", ylab="",cex.lab = 2,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(min(t(p3_figure)[vec3,]),max(t(p3_figure)[vec3,])),main="Cluster Three")
lines(argval,meanp3 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 2,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))
#vec0
meanp3=colMeans(t(p3_figure)[vec0,])
matplot(t, t(t(p3_figure)[vec0,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
        
        main="Noise Group",
        xlab="Day", ylab="",cex.lab = 2,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(min(t(p3_figure)[vec0,]),max(t(p3_figure)[vec0,])))
lines(argval,meanp3 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     cex.lab = 2,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30")) 
########################################################
#######kmeans
reskmeansall = NbClust::NbClust(data =  combinedscore_z, diss = NULL,
                                distance = "euclidean", min.nc = 2, max.nc = 5,
                                method = "kmeans",index="silhouette")
tclusterdataall=data.frame(combinedscore_z)
tclusterdataall$Cluster=as.factor(reskmeansall$Best.partition)

library(ggplot2)
tpsall <- ggplot(tclusterdataall,aes(X1,X2,colour = Cluster)) + geom_point(aes(shape=Cluster),size=3)+ggtitle(paste0("K-means Cluster Results",'\n',"(",dim(tclusterdataall)[1]," Subjects",")")) +
  #xlab(expression('Score '* widehat(xi[i1]))) + ylab(expression('Score '* widehat(xi[i2])))
  xlab('MFPCA Score1 ') + ylab('MFPCA Score2 ')+ theme(plot.title = element_text(hjust = 0.5))+#theme(text = element_text(size = 20))+
  #theme(legend.text = element_text(size = 30))  
  #theme(axis.text=element_text(size = 20),axis.title = element_text(size =30))
  theme(text=element_text(size = 20))
tpsall


# 
# 
# par(mfrow = c(1,2))
# matplot(timestamps01,Z1_est,type='l', lty=1, col="light grey",
#         xlab = "Time", ylab = "Value", main = "Z1")
# # matlines(timestamps01,Z_after[1:672,],type='l', lty=1, col="light blue")
# # 
# # 
# matplot(timestamps01,Z2_est,type='l', lty=1, col="light grey",
#         xlab = "Time", ylab = "Value", main = "Z2")
# #recover
# matlines(timestamps01,Z_after[673:1344,],type='l', lty=1, col="light blue")


# par(mfrow = c(1,2))
# matplot(timestamps01,Z1_est[,-c(index_remove)][,2718],type='l', lty=1, col="light grey",
#         xlab = "Time", ylab = "Value", main = "Z1")
# matplot(timestamps01,Z2_est[,-c(index_remove)][,2718],type='l', lty=1, col="light grey",
#         xlab = "Time", ylab = "Value", main = "Z2")



# par(mfrow = c(1,2))
# matplot(timestamps01,Z1_est[,-c(index_remove)][,1],type='l', lty=1, col="light grey",
#         xlab = "Time", ylab = "Value", main = "Z1")
# matplot(timestamps01,Z2_est[,-c(index_remove)][,1],type='l', lty=1, col="light grey",
#         xlab = "Time", ylab = "Value", main = "Z2")
# 
# 
# par(mfrow = c(1,2))
# matplot(timestamps01,Z1_est[,-c(index_remove)][,-c(1,2718)][,-c(540,2717)][,-2715],type='l', lty=1, col="light grey",
#         xlab = "Time", ylab = "Value", main = "Z1")
# matplot(timestamps01,Z2_est[,-c(index_remove)][,-c(1,2718)][,-c(540,2717)][,-2715],type='l', lty=1, col="light grey",
#         xlab = "Time", ylab = "Value", main = "Z2")




# 
# options(max.print = 1000000)
# #user 3274
# index_z2 <- c(0)
# for (i in 1:dim(Z2_est)[2]){
#   if (length(which(Z2_est[,i]< -250000))>=1){
#     print(paste0("user", i,"index", which(Z2_est[,i]< -250000)[1:5]))
#     index_z2[i] <- i
#   }
#   
# }
# 
# #1
# index_z2 <- na.omit( index_z2)[-1]
# index_z2
# 
# #3469, 3274, 3033
# index_less <- c(0)
# for (i in 1:dim(Z1_est)[2]){
#   if (length(which(Z1_est[,i]< -200))>=1){
#     print(paste0("user", i,"index", which(Z1_est[,i]< -200)[1:5]))
#     index_less[i] <- i
#   }
#   
# }
# 
# #7
# index_less <- na.omit( index_less)[-1]
# index_less
# 
# index_more <- c(0)
# for (i in 1:dim(Z1_est)[2]){
#   if (length(which(Z1_est[,i]> 200))>=1){
#     print(paste0("user", i,"index", which(Z1_est[,i]> 200)))
#     index_more[i] <- i
#   }
# }
# #1
# index_more <- na.omit( index_more)[-1]
# index_more
# 
# 
# 
# #############
# index_remove <- c(index_z2, index_less, index_more)[-1]
# index_remove
# 
# eigen_score <- extract_scores_UNIVFPCA (Z1_est[,-c(index_remove)][,-c(1,2718)][,-c(540,2717)][,-2715],
#                                         Z2_est[,-c(index_remove)][,-c(1,2718)][,-c(540,2717)][,-2715], timestamps01 , PVE=0.95)
# 
# max(eigen_score$scores[,1])
# min(eigen_score$scores[,1])
# which(eigen_score$scores[,1]==max(eigen_score$scores[,1]))
# which(eigen_score$scores[,1]==min(eigen_score$scores[,1]))
# 
# which(eigen_score$scores[,2]==max(eigen_score$scores[,2]))
# 
# plot(eigen_score$scores[,1],eigen_score$scores[,2])
# 
# plot(eigen_score$scores[703,1],eigen_score$scores[703,2])
# 
# 
# save(eigen_score, file = "Twitter_eigen.RData")


###latent curves
#n*2t
#Z_after <-  t(sqrt(length(timestamps01))*eigen_score$scores%*%ginv(eigen_score$Phi))
