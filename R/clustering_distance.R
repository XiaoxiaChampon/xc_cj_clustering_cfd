#Document for Chathura 8/14/2023
##########################################################
############################################################
# Article: Functions Needed for Clustering Categorical Functional Data
#           File that contains all the functions necessary to generate data 
# Author:  Ana-Maria Staicu
# Date: 07/01/2023
##############################################################

#Function to generate the LATENT PROCESS and EVENT PROBAB
#INPUT 
#k - number of eigen functions
#n - number of subjects, 
#m - number of timepoints per curve. Regular grid (same for every)
#Q - numebr of categories. Fixed to 3
#setting  - 1,2,3 settings considered in the ms 
#scenario - A or B used for clustering
#
#OUTPUT
# list with Z1,Z2, p1,p3, p3, mu, Phi, scores

#m - number of time points


generate_data_scenario=function(k=3,n,m, setting=3,  scenario="A", Q=3){
  # k=3;n=100;m=250; setting=3;   scenario="A"; Q=3
  tt <- seq(from = 0.0001, to = 1, length=m)
  k=k
  q=Q  #number of categories
  
  score.var <- rep(0, k) # assumes k=3
  
  # Define MEAN and SCORE VARIANCES for the 3 settings 
  if (setting==1){
    if (scenario=="A" ){ 
      mu_1=function(t)   -1+2*t+2*t^2  
      mu_2=function(t)   -2.5+exp(t*2) 
      #
      score.var[1]<- 1
      score.var[2]<- 1/2
      score.var[3]<- 1/4  }
    
    if (scenario=="B" ){
      mu_1=function(t) -1+2*t+2*t^2
      mu_2=function(t)  -0.5+ exp(t*2) 
      score.var[1]<- 1
      score.var[2]<- 1/2
      score.var[3]<- 1/4  }
    if (scenario=="C"){
      mu_1=function(t)   -1+2*t+2*t^2
      mu_2=function(t)   -2.5+exp(t*2)

      score.var[1]<- 50
      score.var[2]<- 25
      score.var[3]<- 5  }
    
  }
  
  if(setting==2){
    mu_1=function(t) 4*t^2 -1.2        
    mu_2=function(t) 4*t^2 -3.5          
    score.var[1]<- 1
    score.var[2]<- 1/2
    score.var[3]<- 1/4}
  
  if(setting==3){
    mu_1=function(t) -2.2 + 4*t^2 
    mu_2=function(t) -7+6*t^2 
    
   # mu_2=function(t) -7+6*t^2-2
     score.var[1]<- 1
     score.var[2]<- 1/4
     score.var[3]<- 1/16

    }
   
  
  # Define BASIS which is assume the same for all settings
  #*** changed basis to avoid meeting halfway in the same point
  psi_fn=function(k,t){
    psi_k1=sapply(c(1:k), function(i) sin((2*i+1)*pi*t ))
    psi_k2=sapply(c(1:k), function(i) cos(2*i*pi*t ))  
    rbind(psi_k1, psi_k2)
  }
  
  #*** define the matrix of MFPC scores 
  
  scores_standard <- matrix(rnorm(n*k), ncol=k) 
  scores<- scores_standard %*% diag(sqrt(score.var))
  
  #*** construct Z 
  
  BIG_mu <- c(mu_1(tt), mu_2(tt))
  BIG_phi <- psi_fn(k=k, t=tt) 
  
  Z <- BIG_phi %*% t(scores) +BIG_mu 
  Z1<- Z[1:m,]
  Z2<- Z[1:m+m,]
  
 
  
  # RECOVER probabilities for each category
  expZ1 <- exp(Z1)
  expZ2 <- exp(Z2)
  denom <- 1+ expZ1+  expZ2
  p1 <- expZ1/ denom
  p2 <- expZ2/ denom
  #p3 <- 1-p1-p2
  p3 <- 1/denom
  
  #OUTPUT  
  return(list(Z1=Z1, Z2=Z2, p1=p1, p2=p2, p3=p3, MEAN=BIG_mu, PHI=BIG_phi,  MFPC =scores))
}




#Function to generate Categ Functional Data, given the latent event probabilities
# INPUT
# p list of probability matrices that are m by n with probabilities for n subjects, 
#        each observed m times; elements called p1, p2, p3
#
# OUTPUT
# categ FD - m by n with Q=3  categories 
generate_CategFD_scenario=function(p, seed=NULL){
  Q=length(p) # assume Q=3
  p1=p$p1; p2=p$p2; p3=p$p3
  n <- ncol(p1)
  m<- nrow(p1)
  
  if(!is.null(seed)) set.seed(seed)
  
  W <- matrix(0, ncol=n, nrow=m) # identify new object - matrix of dimensions m by n
  X_array <- array(0, c(n,m, Q))
  for (i in c(1:n)){
    X <-sapply(c(1:m) , function(j) rmultinom(n=1, size=1, prob = c(p1[j,i],p2[j,i], p3[j,i]) ) )
    W[,i] <- apply(X, 2, which.max)
    X_array[i,,] <- t(X)
  }
  
  return(list(X=X_array, W=W))
}



############################################################
# Article: Functions Needed for Clustering Categorical Functional Data
#           File that contains all the functions necessary to estimate categorical FD with 3 CATEGORIES! 
# Author:  Ana-Maria Staicu
# Date: 07/01/2023
##############################################################


# estimate_categFD - function to recover the latent trajectories
# 
# INPUT
# X the n by m by Q array with values 0 and 1 obtained from the categ fn DATA (this is needed if SIM=true)
# W the n by m matrix of categories - REQUIRED if SIM=FALSE
# tt - the grid points at which the curves are observed. Regular + Balanced sampling desig
# sim = indicator with values TRUE (if we perform simulation) and FALSE (if we apply it to categ FD)
#       default TRUE
#
#OUTPUT
# list with Z1_est and Z2_est; p1_est, p2_est, p3_est

estimate_categFD <- function(X=NULL, W=NULL, tt, basis_size=25, method="ML", sim=TRUE){

  if(is.null(X)) sim<-FALSE
  
  if(!sim){
    n<- ncol(W)
    m <-nrow(W)
    Q<- length(unique(c(W)))
    Q_vals <- unique(c(W))
    if(is.numeric(Q_vals)) Q_vals<- sort(Q_vals)
    
    X<- array(0, c(n,m,Q))
    for(i in 1:n){
      w1<-W[,i]
      for(j in 1:m) X[i,j, which(Q_vals==w1[j])] <- 1
    }
  }
  
  n<- dim(X)[1]
  m<- dim(X)[2]
  Q <- dim(X)[3]
  
  Z<-NULL
  #change from 
  p<-array(0, c(n, m, Q))
  ##########################
  ###########################
  #n is the subject and this step is done by subject level
  #can parallel
  #############################
  #############################
  for (i in 1:n){
    x1<- X[i,,1]
    x2<- X[i,,2]
    x3<- X[i,,3]
    
    # fit the Binom model
    
    ###################################################
    ###################################################
    ##updated estimation
   
      basis_size_rev_1<-max(min(round(  min(sum(x1), sum(1-x1))/2), basis_size ), 5)

      fit_binom_1<- gam(x1~s(tt, bs = "cr", m=2, k = basis_size_rev_1),
                        #family="binomial", method = method,
                        family=binomial(link="probit"), method = method,
                        control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                        optimizer=c("outer","bfgs"))
    p1 <- fit_binom_1$fitted.values
    p1_linpred <- fit_binom_1$linear.predictors
    
  
      basis_size_rev_2<-max(min(round(  min(sum(x2), sum(1-x2))/2), basis_size ), 5)

      fit_binom_2<- gam(x2~s(tt, bs = "cr", m=2, k = basis_size_rev_2),
                        family=binomial(link="probit"), method = method,
                        control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                        optimizer=c("outer","bfgs"))
  

    p2 <- fit_binom_2$fitted.values
    p2_linpred <- fit_binom_2$linear.predictors
 
    basis_size_rev_3<-max(min(round(  min(sum(x3), sum(1-x3))/2), basis_size ), 5)

    fit_binom_3<- gam(x3~s(tt, bs = "cr", m=2, k = basis_size_rev_3),
                        family=binomial(link="probit"), method = method,
                        control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                        optimizer=c("outer","bfgs"))
    
    p3 <- fit_binom_3$fitted.values
    p3_linpred <- fit_binom_3$linear.predictors
   
    
    
    # estimate the latent tranjecotries Z
    z1<- (p1_linpred-p3_linpred)-log( (1+exp(p1_linpred))/(1+exp(p3_linpred)))
    z2<- (p2_linpred-p3_linpred)-log( (1+exp(p2_linpred))/(1+exp(p3_linpred)))
    
    Z<- cbind(Z, c(z1,z2))
    p[i,,] <- cbind(p1/(p1+p2+p3), p2/(p1+p2+p3), p3/(p1+p2+p3))
    
  }
  
  return(list(Z1_est=Z[1:m,], Z2_est=Z[1:m+m,], 
              p1_est=t(p[,,1]), p2_est=t(p[,,2]), p3_est=t(p[,,3]) ))
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


extract_scores_UNIVFPCA <- function (mZ1,mZ2, tt=tt , PVE=0.95){
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

#Function to find the L2 distance between two latent curves
#Input: yy- 1D vector (true curve)
#       yy2-1D vector (estimated curve)
#Output: scalar-L2 distance
trapzfnum <- function(yy,yy2){
  
  st=0.0001
  et=1
  x=seq(st,et,length=5000)
  xx=seq(st,et,length=length(yy))
  y1 <- cubicspline(xx, yy,x)
  y2 <- cubicspline(xx, yy2,x)
  out=sqrt(trapz(x, (y1-y2)^2) )
  out
  
}

#Function to find the Hellinger distance between two probability curves
#Input: yy- 1D vector (true curve)
#       yy2-1D vector (estimated curve)
#Output: scalar-Hellinger distance
trapzfnump <- function(yy,yy2){
  st=0.0001
  et=1
  x=seq(st,et,length=5000)
  xx=seq(st,et,length=length(yy))
  y1 <- cubicspline(xx, yy,x)
  y1[y1<0]=0
  y2 <- cubicspline(xx, yy2,x)
  y2[y2<0]=0
  out=sqrt(trapz(x, (sqrt(y1)-sqrt(y2))^2) )
  out
  
}
#Function to use trapzfnum function and find L2 distance for 2D array, n of them
mse_bw_matrix=function(truecurve,estcurve){
  n=dim(truecurve)[2]
  datapoints=dim(truecurve)[1]
  mseall=c(0)
  ######could probably use apply function here it's also subject level
  for (i in 1:n){
    mseall[i]=trapzfnum(truecurve[,i],estcurve[,i])
    
  }
  mseall
  
}

#Function to use trapzfnump function and find Hellinger distance for 2D array, n of them
mse_bw_matrixp=function(truecurve,estcurve){
  n=dim(truecurve)[2]
  datapoints= dim(truecurve)[1]
  mseall=c(0)
  ######could probably use apply function here it's also subject level
  for (i in 1:n){
    mseall[i]=trapzfnump(truecurve[,i],estcurve[,i])/sqrt(2)
   
    
  }
  mseall
}
############################################################
# Article: Functions Needed for Clustering Categorical Functional Data
#           Main file to generate data and obtain simulation results
# Author:  Ana-Maria Staicu
# Date: 07/01/2023
##############################################################


#Load the libraries
#############################################################
library(xtable)
library(fda)
library(refund)
library(Matrix)
library(MASS)
library(arm)
library(mgcv)
library(readr)
library(funData)
library(MFPCA)
library(purrr)
library(tidyverse)
library(gridExtra)
library(dbscan)
library(tidyr)
library(dplyr)
library(fossil)
#devtools::install_github("ahasverus/elbow")
library(elbow)
library(pdfCluster)
library(FADPclust)
library(pracma)


###################################################################################
#This can be distributed to different  simulation
# Function to Cluster and Calculate the estimation accuracy and clustering
#Input n: scalar-number of subjects/observations
#      m: scalr-number of time points per curve
#      scenario: options are "A", "B" or "C", which uses different mean function and scores
#      mc_sims: the number of the monte carlo simulations

#Output: list-if n==100, also calculate the accuracy mse and Hellinger distance
#             else, list(cluster_table_true-clustering results based on true scores,
#                        "cluster_table_est"-clustering results based on estimated scores,
#                         "cluster_table_est_se"=standared errors)
cluster_simulation=function(n,m,scenario,mc_sims){
  n=1000
  # m=2000
  # scenario="C"
  # mc_sims=5
  if (scenario=="A"){
    p1=0.75
    p2=0.22
    p3=0.03
  }
  
  if (scenario=="B"){
    p1=0.5
    p2=0.3
    p3=0.2
  }
  
  if (scenario=="C"){
 
    
    # p1=0.4
    # p2=0.5
    # p3=0.1
    p1=0.1
    p2=0.6
    p3=0.3
  }
  
  n1<-n*p1; 
  n2<- n*p2; 
  n3<-n*p3
  true_cluster <- c(rep(1, n1), rep(2, n2), rep(3, n3)) #this is what we know whenw
  
#when we generate the data, but if it's real data. We don't know the true cluster.
  #how can we evaluate the dbscan at each step to see which parameter is the best?
  #that's my question
  #how to create teh true_cluster? how do you create the labesl for each subjects
  #in the real data?
  true_cluster_db <- c(rep(1, n1), rep(2, n2), rep(0, n3))
  tt <- seq(from = 0.0001, to = 1, length=m) # gridpoint t
  
  # seed for latent process
  seed1<-1230
  # seed for categFD, given the latent process
  seed0<- 9876

  
   rmse<- array(0, c(mc_sims, 2, 3))       #2 components and 3 settings
   hellinger <- array(0, c(mc_sims, 3, 3)) #3 events probab and 3 settings
  
  true_kmeans <- est_kmeans  <- NULL # records clustering membership on the TRUE SCORES/ESTIMATED SCORES
  true_fadp <- est_fadp <-  NULL # records clustering membership on the TRUE Z/ESTIMATED Z
  true_dbscan <- est_dbscan <-  NULL

  for(ii in 1:mc_sims){
    #ii=2
    set.seed(seed1+100*ii)

    # generate clusters
    CLSTF1 <- generate_data_scenario(k=3,n=n1, m=m, setting=1, scenario, Q=3)
    CLSTF2 <- generate_data_scenario(k=3,n=n2, m=m, setting=2, scenario,  Q=3)
    CLSTF3 <- generate_data_scenario(k=3,n=n3, m=m, setting=3, scenario,  Q=3)
    
    
    # Recover the latent Gaussian process
    Z1<-cbind(CLSTF1$Z1,CLSTF2$Z1 , CLSTF3$Z1  )
    Z2<-cbind(CLSTF1$Z2,CLSTF2$Z2 , CLSTF3$Z2  )
    
    # Recover the true probab curves 
    p1<-cbind(CLSTF1$p1,CLSTF2$p1 , CLSTF3$p1  )
    p2<-cbind(CLSTF1$p2,CLSTF2$p2 , CLSTF3$p2  )
    p3<-cbind(CLSTF1$p3,CLSTF2$p3 , CLSTF3$p3  )
    p<-list(p1=p1 ,p2=p2, p3=p3)
    
    #generate categFD
    set.seed(seed0+ii*100)
    out <- generate_CategFD_scenario(p=p)
    
    X<- out$X  
    W<- out$W
    
    
    #######add 7/25/2023 while loop to update cluster 3
    ##################################################
    Q_vals <- unique(c(W))
    
    if(is.numeric(Q_vals)) 
      Q_vals<- sort(Q_vals)
    ###############add to exclude the rare event at the end
    ################
    #W m*n  n=n1 72 +n2 22+n3. n1+1 : (n1+n2)
    
    for (clstr_2_idx in 1:n) {
      if (clstr_2_idx %in% 1:n1 ){
        setting_choice=1
      }
      
      if (clstr_2_idx %in% (n1+1):(n1+n2) ){
        setting_choice=2
      }
      
      if (clstr_2_idx %in% (n1+n2+1):(n) ){
        setting_choice=3
      }
      
      
      tolcat = table(W[, clstr_2_idx])
      catorder = order(tolcat, decreasing = TRUE)
      numcat = length(catorder)
      refcat = catorder[numcat]
      Wnew = matrix(0, nrow = m, ncol = 5)
      count_iter = 0
      while ((min(as.numeric(tolcat)) == 1 && W[, clstr_2_idx][m] == refcat && count_iter < 100)
             || (min(as.numeric(tolcat)) == 1 && W[, clstr_2_idx][1] == refcat && count_iter < 100))
      {
        count_iter = count_iter + 1
        print(cat("count_iter: ", count_iter))
        CLSTF22 = generate_data_scenario(
          k = 3,
          n = 5,
          m = m,
          setting = setting_choice,
          scenario = "A",
          Q = 3
        )
        
        # 
        # Recover the true probab curves
        p11 <- CLSTF22$p1
        p21 <- CLSTF22$p2
        p31 <- CLSTF22$p3
        
        p111 <- list(p1 = p11 , p2 = p21, p3 = p31)
        out1 <- generate_CategFD_scenario(p = p111)
        
        #X <- out$X #n,m, q
        Wnew <- out1$W #m,n.  m*5
        W[, clstr_2_idx] = Wnew[, 3]
        ###
        # Recover the latent Gaussian process
        Z1[,clstr_2_idx] <- CLSTF22$Z1[,3]
        Z2[,clstr_2_idx] <- CLSTF22$Z2[,3]
        ###########
        X[clstr_2_idx,,]=0
        w1 = W[, clstr_2_idx]
        for (j in 1:m)
          X[clstr_2_idx, j, which(Q_vals == w1[j])] <- 1
        
        tolcat = table(W[, clstr_2_idx])
        catorder = order(tolcat, decreasing = TRUE)
        numcat = length(catorder)
        refcat = catorder[numcat]
        Wnew = matrix(0, nrow = m, ncol = 5)
        
      }
      X = X
    }
    
    
    ###########
    
    #ESTIMATION
    categFD_est <- estimate_categFD(X=X, tt=tt)
    
    # recover latent process
    Z1_est=categFD_est$Z1_est
    Z2_est=categFD_est$Z2_est
    ###
    par(mfrow=c(2,2))
    matplot(seq(0.0001,1,length=m),Z2[,1:n1],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,xlim=c(0,1),ylim=c(-40,15))
    matlines(seq(0.0001,1,length=m),Z2[,(n1+1):(n1+n2)],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
    matlines(seq(0.0001,1,length=m),Z2[,(n1+n2+1):n],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
    
    
    # #plot first truelatent curves p_2 by clusters
    matplot(seq(0.0001,1,length=m),Z2_est[,1:n1],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(-40,15))
    matlines(seq(0.0001,1,length=m),Z2_est[,(n1+1):(n1+n2)],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
    matlines(seq(0.0001,1,length=m),Z2_est[,(n1+n2+1):n],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
    
    
    
    
    #recover =probabilities for each category
    p1_est<-categFD_est$p1_est
    p2_est<-categFD_est$p2_est
    p3_est<-categFD_est$p3_est
    
  if ( n==100 && scenario=="A"){
    # evaluate performance Z and P
    rmse1_temp <- c(by(mse_bw_matrix(Z1,Z1_est) , true_cluster, mean))
    rmse2_temp <- c(by(mse_bw_matrix(Z2,Z2_est), true_cluster, mean))
    rmse[ii, ,] <- rbind(rmse1_temp,rmse2_temp )
    
    error.p1<- mse_bw_matrixp(p1,p1_est)
    error.p2<- mse_bw_matrixp(p2,p2_est)
    error.p3<- mse_bw_matrixp(p3,p3_est)
    
    
    hellinger[ii, ,] <-  rbind( c(by(error.p1, true_cluster, mean)),
                                c(by(error.p2, true_cluster, mean)),
                                c(by(error.p3, true_cluster, mean)))
    
  }
   
  
    
    # Record clustering performance. USE only Z
    #**** extract_scores_UNIVFPCA is the FASTEST
    mfpca_true <-extract_scores_UNIVFPCA(mZ1=Z1, mZ2=Z2,  tt=tt , PVE=0.95)
    #plot(  scores_true$scores[, 1:2])
    mfpca_est <- extract_scores_UNIVFPCA(mZ1=Z1_est, mZ2=Z2_est, tt=tt , PVE=0.95)
    
    
    plot(mfpca_true$scores[,1],mfpca_true$scores[,2])
    plot(mfpca_est$scores[,1],mfpca_est$scores[,2])
    
    #KMEANS
    true_kmeans_temp <- kmeans_cluster(data=mfpca_true$scores)$label
    est_kmeans_temp <- kmeans_cluster(data=mfpca_est$scores)$label
    
    #FADP
    true_fadp_temp <- fadp_cluster(mZ1=Z1, mZ2=Z2, tt=tt)$label
    est_fadp_temp <- fadp_cluster(mZ1=Z1_est, mZ2=Z2_est, tt=tt)$label
    
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
    print(ii)
  }
  
  
  if (n==100 && scenario=="A"){
    # Assess accuracy in the simulation study
    mse_sim= apply(rmse, c(2,3), mean)
    hellinger_sim= apply(hellinger, c(2,3), mean)
    
  }

  
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
  
  
  cluster_table_est_se=c(sd(results$est_fadp_ri)/sqrt(mc_sims),sd(results$est_fadp_ari)/sqrt(mc_sims),
                         
                         sd(results$est_fadp_cpn)/sqrt(mc_sims),
                         
                         
                         sd(results$est_kmeans_ri)/sqrt(mc_sims),sd(results$est_kmeans_ari)/sqrt(mc_sims), 
                         sd(results$est_kmeans_cpn)/sqrt(mc_sims),
                         sd(results$est_dbscan_ri)/sqrt(mc_sims),
                         sd(results$est_dbscan_ari)/sqrt(mc_sims),sd(results$est_dbscan_cpn)/sqrt(mc_sims)
  )
  names(cluster_table_est_se)=c("fadp RI","fadp ARI","fadp cpn",
                                "kmeans RI","kmeans ARI","kmeans cpn",
                                "dbscan RI","dbscan ARI","dbscan cpn"
  )
  
  if (n==100 && scenario=="A"){
    return(list("cluster_table_true"=cluster_table_true,"cluster_table_est"=cluster_table_est,
                "cluster_table_est_se"=cluster_table_est_se,"mse"=mse_sim,"hellinger"=hellinger_sim))
    
  }
  
 else{
   return(list("cluster_table_true"=cluster_table_true,"cluster_table_est"=cluster_table_est,
               "cluster_table_est_se"=cluster_table_est_se))
   
 }
}


##########scenarioA
#cluster_simulation=function(n,m,scenario,mc_sims)
#n: number of subjects, m: number of time points, scenario: "A", "B" or "C", mc_sims: number of simulations
n100t300A=cluster_simulation(100,300,"A",50)
n100t750A=cluster_simulation(100,750,"A",50)
n100t2000A=cluster_simulation(100,2000,"A",50)


n500t300A=cluster_simulation(500,300,"A",50)
n500t750A=cluster_simulation(500,750,"A",50)
n500t2000A=cluster_simulation(500,2000,"A",50)


n1000t300A=cluster_simulation(1000,300,"A",50)
n1000t750A=cluster_simulation(1000,750,"A",50)
n1000t2000A=cluster_simulation(1000,2000,"A",50)


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


mse_tableA=rbind(
  cbind(n100t300A$mse[1,1],n100t300A$mse[2,1],n100t300A$hellinger[1,1],n100t300A$hellinger[2,1],n100t300A$hellinger[3,1]),
  cbind(n100t750A$mse[1,1],n100t750A$mse[2,1],n100t750A$hellinger[1,1],n100t750A$hellinger[2,1],n100t750A$hellinger[3,1]),
  cbind(n100t2000A$mse[1,1],n100t2000A$mse[2,1],n100t2000A$hellinger[1,1],n100t2000A$hellinger[2,1],n100t2000A$hellinger[3,1]),
  
  cbind(n100t300A$mse[1,2],n100t300A$mse[2,2],n100t300A$hellinger[1,2],n100t300A$hellinger[2,2],n100t300A$hellinger[3,2]),
  cbind(n100t750A$mse[1,2],n100t750A$mse[2,2],n100t750A$hellinger[1,2],n100t750A$hellinger[2,2],n100t750A$hellinger[3,2]),
  cbind(n100t2000A$mse[1,2],n100t2000A$mse[2,2],n100t2000A$hellinger[1,2],n100t2000A$hellinger[2,2],n100t2000A$hellinger[3,2]),
  
  cbind(n100t300A$mse[1,3],n100t300A$mse[2,3],n100t300A$hellinger[1,3],n100t300A$hellinger[2,3],n100t300A$hellinger[3,3]),
  cbind(n100t750A$mse[1,3],n100t750A$mse[2,3],n100t750A$hellinger[1,3],n100t750A$hellinger[2,3],n100t750A$hellinger[3,3]),
  cbind(n100t2000A$mse[1,3],n100t2000A$mse[2,3],n100t2000A$hellinger[1,3],n100t2000A$hellinger[2,3],n100t2000A$hellinger[3,3]),
  )


rownames(mse_tableA)=c("s1n100t300","s1n100t750","s1n100t2000",
                          "s2n100t300","s2n100t750","s2n100t2000",
                          "s3n100t300","s3n100t750","s3n100t2000") 
colnames(mse_tableA)=c("z1","z2","p1","p2","p3")
save(true_tableA,est_tableA,est_tableA_se,mse_tableA,file="A_clustering.RData")

# load("staicu_Anew20_final.RData")
# library(xtable)
# xtable(est_tableA)


#####################################
###scenario B
##########scenarioA
#cluster_simulation=function(n,m,scenario,mc_sims)
n100t300B=cluster_simulation(100,300,"B",50)
n100t750B=cluster_simulation(100,750,"B",50)
n100t2000B=cluster_simulation(100,2000,"B",50)


n500t300B=cluster_simulation(500,300,"B",50)
n500t750B=cluster_simulation(500,750,"B",50)
n500t2000B=cluster_simulation(500,2000,"B",50)


n1000t300B=cluster_simulation(1000,300,"B",50)
n1000t750B=cluster_simulation(1000,750,"B",50)
n1000t2000B=cluster_simulation(1000,2000,"B",50)


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



#cluster_simulation=function(n,m,scenario,mc_sims)
n100t300C=cluster_simulation(100,300,"C",50)
n100t750C=cluster_simulation(100,750,"C",50)
n100t2000C=cluster_simulation(100,2000,"C",50)


n500t300C=cluster_simulation(500,300,"C",50)
n500t750C=cluster_simulation(500,750,"C",50)
n500t2000C=cluster_simulation(500,2000,"C",50)


n1000t300C=cluster_simulation(1000,300,"C",50)
n1000t750C=cluster_simulation(1000,750,"C",50)
n1000t2000C=cluster_simulation(1000,2000,"C",50)


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

