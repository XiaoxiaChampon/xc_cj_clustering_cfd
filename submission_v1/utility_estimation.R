############################################################
# Article: Functions Needed for Clustering Categorical Functional Data
#           File that contains all the functions necessary to estimate categorical FD with 3 CATEGORIES! 
# Author:  Xiaoxia Champon, Ana-Maria Staicu, Chathura Jayalath
# Date: 09/25/2023
##############################################################



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
#' @param  yy- 1D vector (true curve)
#' @param        yy2-1D vector (estimated curve)
#' @return : scalar-L2 distance
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
#' @param  truecurve- 2D array (true curve)
#' @param   estcurve-2D array (estimated curve)
#' @return : scalar-L2 distance

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
#output: scores matrix for that specific Nth simulation   n*M     M is the column number of scores
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
    return(EstimateCategFuncData_probit(timestamps01, X, basis_size=25, method="ML"))
  }else if(choice == "binormial"){
    X <- GetXFromW(W)
    return(EstimateCategFuncData_binorm(timestamps01, X, basis_size=25, method="ML"))
  }else if(choice == "multinormial"){
    return(EstimateCategFuncData_multinormial(timestamps01, W, basis_size=25, method="ML"))
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


#' @param timestamps01, 1D array, time interval that cfd is observed
#' @param  X: 3D array, t*n*Q, t: the number of time points, n: the number of individuals, Q: the number of categories
#' @param  basis_size=25, the number of basis function used 
#' @param  method="ML"
#' @return list of 2D array: True Z curves , Est Z curves, True p curves, Est p curves
#'                           all have dimension t*n

EstimateCategFuncData_probit <- function(timestamps01, X, basis_size=25, method="ML")
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



