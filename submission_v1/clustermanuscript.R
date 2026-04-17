#Manuscript: Clustering Categorical functional data with application to social media
#Purpose: code to reproduce the results from the manuscript and demonstrate the catfda R package
#Date : 4/6/2023
#Part I: list required libraries
#Part II: define all the functions
#Part III: sample code to reproduce Figures, Table 1 RMSE and Table 2 clustering results for the simulated data in the manuscript
#Part IV: Twitter Data pre processing
#Part V: Clustering results for Twiiter Data from the manuscript
#Part VI: catfdcluster function from catfda R package to cluster categorical functional data
#Part VII: sample code to use function catfdcluster from the catfda package 
#Note:
#-for demonstration purpose we only include one combination with pre defined n and datapoints in Part III 
# to reproduce Figure 2, Figure 3 and its clustering results
#-users can adjust n and datapoints to reproduce the results of all the combinations in the manuscript
#-users can use the acc_table_sim function to simulate each one of the 9 combinations N times by setting the seedlength=N

#part I: load required libraries
##################################################################################################
##################################################################################################
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
library(cfda)
library(tidyr)
library(dplyr)
library(factoextra)
library(NbClust)
library(fossil)
library(FADPclust) 
library(ggplot2)
devtools::install_github("ahasverus/elbow")
library(elbow)
library(pdfCluster)
library(scales)
####
#Twitter data collection and processing libraries
#Set up collection
library(twitteR)
library(devtools)
library(readr)


#part II: define functions
#######################################################################################################
#######################################################################################################
#Function to return the logit
logit <- function(x){
  log(x)-log(1-x)
}

######################################################################################################
#Function: wrapper function for gam() which outputs the fitted values
#
#Inputs: 
# z : index z = 1,...,N 
# Curves : N x D matrix of observed binary series
# tt : grid of timepoints going from 0 to 1 with D observations
# k : number of basis functions
# method: method used to evaluate the gam
#
#Output: 
# Fitted values from the gam function for subject z 
#
#####
#get smoothed curves
regression_g = function(z, Curves, tt, k=25, method="ML"){   #changed from 10 to 25, changed to 100
  z1 = Curves[z,]
  gam1 <- gam(z1~s(tt, bs = "cr", m=2, k = k),
              family="binomial", method = method,control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),optimizer=c("outer","bfgs"))
  return(gam1$fitted.values)
}


######################################################################################################################
#Function to output the smooth curves using link function

#New function to output predicted L-1 latent curves:
#input observed binary curves X_it and return smoothed Z_ihat
#tt: vector time interval
#Output: smoothed latent curves
Z_ihat=function(Curves,tt){
  N=dim(Curves)[1]
  vec = matrix(1:(N), ncol = 1)
  smoothed_x = logit(t(apply(vec, 1, function(x) regression_g(x, Curves, tt))))
  smoothed_x
}
########################################################################################################################


#Function to generate scores
#input mu:scalar- mean of the scores, sd: scalar -standard deviation
#output: score
score=function(mu,sd){
  rnorm(1,mean=mu,sd=sd)
}

#############################################################################################################
#Function to generate categorical curves
#Input:
#n number of subjects
#datapoints
#sparse=5 and sparse=7 yes   sparse=0 no
#scorevar=1 bigger var , scorevar=3 smaller var
#ps=1 equation (4) in the manuscript to generate latent curves Z_il from V_iq
#k  #number of eigen functions
#q  #level of the categorical level
#st: starting time
#et: end time

#Output
#True: True categorical curves, true latent curves, true eigen functions, true probability curves
#Est:  Estimated Latent curves, estimated probability curves

generate_data_scenario=function(k,n,datapoints,sparse,scorevar,ps,seed=123,st,et){
  k=k
  seed=seed
  st=st
  et=et
  scorevar=scorevar
  #k=3  #number of eigen functions
  q=3  #level of the categorical level
  
  if (sparse==0){
    mu_1=function(t){
      #mean function mu1
      1+4*t
    }
    mu_2=function(t){
      #mean function mu2
      0.97+6*t^2
    }
    
    p_ihat=function(Z_i1app,Z_i2app){
      denom=(1+exp(Z_i1app)+exp(Z_i2app))
      p_i1h=exp(Z_i1app)/denom
      p_i2h=exp(Z_i2app)/denom
      p_i3h=1- p_i1h- p_i2h
      return(list("p_1hatmatrix"=p_i1h,"p_2hatmatrix"=p_i2h,"p_3hatmatrix"=p_i3h))
    }
    
  }
  
  
  if(sparse==5){
    mu_1=function(t){
      3.8+4*t^2-5  
      
    }
    mu_2=function(t){
      1.5+4*t^2-5   
      
    }
    
    p_ihat=function(Z_i1app,Z_i2app){
      denom=(1+exp(Z_i1app)+exp(Z_i2app))
      p_i1h=exp(Z_i1app)/denom
      p_i2h=exp(Z_i2app)/denom
      p_i3h=1/denom
      # p_i3h=1-p_i1h- p_i2h
      return(list("p_1hatmatrix"=p_i1h,"p_2hatmatrix"=p_i2h,"p_3hatmatrix"=p_i3h))
    }
  }
  
  if(sparse==7){
    mu_1=function(t){
      3.8+4*t^2-6  
      
    }
    mu_2=function(t){
      0.97+6*t^2-8
    }
    
    
    p_ihat=function(Z_i1app,Z_i2app){
      denom=(1+exp(Z_i1app)+exp(Z_i2app))
      p_i1h=exp(Z_i1app)/denom
      p_i2h=exp(Z_i2app)/denom
      p_i3h=1/denom
      # p_i3h=1-p_i1h- p_i2h
      return(list("p_1hatmatrix"=p_i1h,"p_2hatmatrix"=p_i2h,"p_3hatmatrix"=p_i3h))
    }
  }
  
  mu_vec=rep(0,k)
  #eigen function
  psi_fn=function(k){
    
    psi_k1=matrix(rep(1,length(t)*k),ncol=k)  
    psi_k2=matrix(rep(1,length(t)*k),ncol=k) 
    for (i in 1:k) {
      psi_k1[,i]=sin(2*i*pi*t )
      psi_k2[,i]=cos(2*i*pi*t )
    }
    list("psi_k1"=psi_k1,"psi_k2"=psi_k2)
  }
  
  
  t=seq(from = st,to = et, length=datapoints)
  
  X_i=array(0,dim=c(q,datapoints,n))  #multinormial results: row is level q, column is time points, n is the number of subjects, each column only has one row of 1 and every other rows are 0
  X_nt=matrix(rep(1,n*length(t)),nrow=n,ncol=length(t))  #true observations of categorical-valued outcome, each row represent one subject, columns represent time points
  score_matrix=matrix(rep(1,n*k),nrow=n,ncol=k)  #row is number of subjects, column is the number of eigen functions
  psi_score_matrix_1=matrix(rep(1,n*length(t)),ncol=n)  #dim: length(t)*nsubjects
  psi_score_matrix_2=matrix(rep(1,n*length(t)),ncol=n)
  Z_i1=matrix(rep(1,n*length(t)),nrow=n)  #True latent curves1:row is n subjects, col is t time points
  Z_i2=matrix(rep(1,n*length(t)),nrow=n) #True latent curve 2
  p_i1=matrix(rep(0,n*length(t)),nrow=n)  #True p_i1
  p_i2=matrix(rep(0,n*length(t)),nrow=n)  #True p_i2
  p_i3=matrix(rep(0,n*length(t)),nrow=n)  #True p_i3
  for (i in 1:n){
    set.seed(seed+i)
    if (k==3){
      if (scorevar==1){
        #score varies based on i
        score_1=score(0,sqrt(1))
        score_2=score(0,sqrt(1/2))
        score_3=score(0,1/2)
      }
      
      if (scorevar==3){
        #score varies based on i
        score_1=score(0,1)
        score_2=score(0,sqrt(1/4))
        score_3=score(0,1/4)
      }
      
      
      score_vector=cbind(score_1,score_2,score_3)
    }
    
    psi_k1=psi_fn(k)$psi_k1
    psi_k2=psi_fn(k)$psi_k2
    
    #Z varies based on i
    #psi t*k, score: t*k,  psi%*%t(score)
    psi_score_matrix_1[,i]=psi_k1%*%t(score_vector)
    Z_i1[i,]=mu_1(t)+psi_score_matrix_1[,i]
    
    psi_score_matrix_2[,i]=psi_k2%*%t(score_vector)
    Z_i2[i,]=mu_2(t)+psi_score_matrix_2[,i]
    
    
    #p varies based on i
    denominator=(1+exp(as.vector(Z_i1[i,]))+exp(as.vector(Z_i2[i,])))
    p_i1[i,]=(exp(as.vector(Z_i1[i,])))/denominator
    p_i2[i,]=(exp(as.vector(Z_i2[i,])))/denominator
    p_i3[i,]=1-p_i1[i,]-p_i2[i,]
    
    
    #X_i varies based on i,j
    
    for (j in 1:length(t)){
      X_i[,j,i]=rmultinom(n=1, size=1, prob=c(p_i1[i,j],p_i2[i,j],p_i3[i,j]))
    }
    
    #X_it varies based on i
    X_it=c(1)
    for (j in 1:length(t)){
      X_it[j]=as.vector(which(X_i[,j,i] == 1))
    }
    X_nt[i,]=X_it
    
    #collect score matrix
    score_matrix[i,]=score_vector
  }
  
  #collect value and graph
  #collect first two rows of observed binary curves
  X_i1=t(X_i[1,,])  #all n row subjects , t columns values related to p1
  X_i2=t(X_i[2,,]) #all n row subjects , t columns values related to p2
  X_i3=t(X_i[3,,]) #all n row subjects , t columns values related to p3
  
  #recover Z_i1 hat using X_i[1,all j, all n] only related to p1
  Z_i1hat=Z_ihat(X_i1,t)
  #recover Z_i2 hat using X_i[2,all j, all n] only related to p2 
  Z_i2hat=Z_ihat(X_i2,t)
  #recover Z_i3 hat using X_i[2,all j, all n] only related to p2 
  Z_i3hat=Z_ihat(X_i3,t)
  
  
  #recover Q-1  Latent curves
  if(ps==1){ 
    Z_i1hatstar=Z_i1hat+log(1+exp(Z_i3hat))-log(1+exp(Z_i1hat))-Z_i3hat
    Z_i2hatstar=Z_i2hat+log(1+exp(Z_i3hat))-log(1+exp(Z_i2hat))-Z_i3hat
  }
  

  p_i1hat=p_ihat(Z_i1hatstar,Z_i2hatstar)$p_1hatmatrix
  p_i2hat=p_ihat(Z_i1hatstar,Z_i2hatstar)$p_2hatmatrix
  p_i3hat=p_ihat(Z_i1hatstar,Z_i2hatstar)$p_3hatmatrix
  

  truel=list("TrueZ_i1"=Z_i1,"TrueZ_i2"=Z_i2,"Truep_i1"=p_i1,"Truep_i2"=p_i2,"Truep_i3"=p_i3,"Truescore_matrix"=score_matrix,"Truecatcurve"=X_nt,"comp1"=psi_k1,"comp2"=psi_k2)
  est=list("EstimateZ_i1"=Z_i1hatstar,"EstimateZ_i2"=Z_i2hatstar,"Estimatep_i1"=p_i1hat,"Estimatep_i2"=p_i2hat,"Estimatep_i3"=p_i3hat)
  return(list("Trueth"=truel,"Est"=est))
}

#######################################################################################################
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
#########################################################################################################
####
#cfda method to extract scores
####
#write a function to produce the scores
#input one  X_nt matrix n*t: 1 of  Nth simulation, n*t* N
#output scores matrix for that specific Nth simulation   n*M     M is the column number of scores
cfda_score=function(cfda_data,datapoints, M,st,et){
  t=seq(from = st,to = et, length=datapoints)
  times=seq(1,datapoints,by=1)
  colnames(cfda_data)=paste0("time_",times)
  time=rep(t,times=dim(cfda_data)[1])
  subjects=seq(1,dim(cfda_data)[1],by=1)
  cfda_new=data.frame(cbind(cfda_data,subjects))
  #reshape the data into the long format
  data_long <- gather(cfda_new, time_t, state, paste0("time_",times[1]):paste0("time_",times[datapoints]),factor_key=TRUE)
  cfda_d=data_long%>%arrange(subjects)%>%dplyr::rename(id=subjects)%>%dplyr::select(id,state)
  cfda_final=data.frame(cbind(cfda_d,time))
  Tmax=max(cfda_final$time)
  cfda_cut=cut_data(cfda_final, Tmax = Tmax)
  m <- 10
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
  fmca <- compute_optimal_encoding(cfda_cut, b, nCores = 2, verbose = FALSE)
  #choose the number of scores that explians at least 0.95 pve
  nPc95 <-which(cumsum(prop.table(fmca$eigenvalues)) > 0.95)[1]
  cfda_score=fmca$pc[, 1:nPc95]
}

#####################################################################################################
#function to output average score matrix for N simulations
#input the data set n*t*N, n number of subjects, t number of data points,N number of simulations
#M: number of the column of scores matrix
#output average score matrix n*M
cfda_score_apply=function(data,datapoints, M ,st, et){
  data_used=data$AverageTrue$Truecurves
  n=dim(data_used)[1]
  datapoints=dim(data_used)[2]
  N=dim(data_used)[3]
  M=M
  st=st
  et=et
  cfda_score_matrix=apply(data_used,3, cfda_score ,datapoints, M, st, et)
}

#######################################################################################################
#Function to get the kmeans and dbscan clustering results with one specific combination using pre defined n and datapoints 
#q: number of category 3
#k: number of true eigen functions 3
#knnum: the number of the neighbours when calculating the knn distance
#spline1D: the number of splines needed and we set it
#as 25 in all the simulations, 
#Mcfda is the column of initial mfpca scores, 
#pct is the percent of the distance,
#pctcfda is the percent of the cfda distance,
#minPts is the initial minPts, 
#min.nc is the minimum number of clusters,
#max.nc is the maximum number of clusters
acc_table_graph=function(seed,k,n,datapoints,q,sparse1,sparse2,sparse3,scorevar1,scorevar2,ps,st,et,
                         spline1D,Mcfda,knnum,pct, pctcfda,minPts,min.nc,max.nc){
  seed=seed
  n=n
  q=q
  k=k
  knnum=knnum
  minPts=minPts
  
  onep=0.5
  twop=0.3
  threep=0.2
  # onep=0.75
  # twop=0.22
  # threep=0.03
  clustern50t250=generate_data_scenario(k=k,n=round(onep*n),datapoints,sparse1,scorevar1,ps,seed=seed,st,et)
  clustern50t250p2=generate_data_scenario(k=k,n=round(twop*n),datapoints,sparse2,scorevar2,ps,seed=seed,st,et)
  clustern50t250p3=generate_data_scenario(k=k,n=round(n*threep),datapoints,sparse3,scorevar2,ps,seed=seed,st,et)
  
  
  truen100t250=mapply(rbind,clustern50t250$Trueth,clustern50t250p2$Trueth,clustern50t250p3$Trueth,SIMPLIFY=FALSE,USE.NAMES = TRUE)
  estn100t250=mapply(rbind,clustern50t250$Est,clustern50t250p2$Est,clustern50t250p3$Est,SIMPLIFY=FALSE,USE.NAMES = TRUE)
  combn100t250=list("Trueth"=truen100t250, "Est"=estn100t250)
  
  sp=Sys.time()
  combn100t250mfpca=mfpca_ram(combn100t250,2,datapoints,q=3,
                                       n,st,et,seq(st,et,length=datapoints))
  
  ept=Sys.time()
  tp=ept-sp
  
  datalist=list("1sthalf"=clustern50t250,"2ndhalf"=clustern50t250p2,"combmfpca"=combn100t250mfpca)
  
  ######
  dbst=Sys.time()
  combinedscore_z <- combn100t250mfpca$scores
  min.nc=2
  max.nc=5
  dimz=dim(combinedscore_z)[2]
  ###########
  if (dimz<=2){
    
    minPts=4
    
  }
  
  
  if (dimz>2){

    minPts=dimz+1
  }
  #dist=kNNdist(combinedscore_z, k = knnum)
  dist=kNNdist(combinedscore_z, k = minPts-1)
  #########change to max increase
  distdataelbow=data.frame(sort(dist))
  distdataelbow$index=1:(dim(combinedscore_z)[1])
  ipoint <- elbow(data = distdataelbow)
  epsoptimal=ipoint$sort.dist._selected*3
  #epsoptimal=sortdist[which.max(diff(sortdist))]
  
  ##############################
  
  #res <- dbscan(combinedscore_z, eps =ninty5p , minPts = minPts)
  res <- dbscan(combinedscore_z, eps =epsoptimal , minPts = minPts)
  true_label <- c(rep(1,round(onep*n)),rep(2,round(twop*n)),rep(3,round(threep*n)))
  correctp=rand.index(true_label, res$cluster)
  correctpadj=adj.rand.index(true_label, res$cluster)
  #############add accuracy for noise group
  c3=which(true_label==3)
  cc1=which(res$cluster==1)
  cc2=which(res$cluster==2)
  cc3=which(res$cluster==3)
  noisep=max(length(which(c3 %in% cc1)),length(which(c3 %in% cc2)),length(which(c3 %in% cc3)))/length(c3)
  ##############
  dbet=Sys.time()
  dbt=tp+dbet-dbst
  ###################
  # k-means
  kst=Sys.time()
  reskmeans = NbClust::NbClust(data =  combinedscore_z, diss = NULL,
                               distance = "euclidean", min.nc = min.nc, max.nc = max.nc,
                               method = "kmeans",index="silhouette")
  correctk=rand.index(true_label, reskmeans$Best.partition)
  correctkadj=adj.rand.index(true_label, reskmeans$Best.partition)
  ########
  cc1k=which(reskmeans$Best.partition==1)
  cc2k=which(reskmeans$Best.partition==2)
  cc3k=which(reskmeans$Best.partition==3)
  noisepk=max(length(which(c3 %in% cc1k)),length(which(c3 %in% cc2k)),length(which(c3 %in% cc3k)))/length(c3)
  #######
  ket=Sys.time()
  kt=tp+ket-kst
  
  acctable=cbind(correctp,correctk,correctpadj,correctkadj,noisep,noisepk)
  timedbk=list("dbt"=dbt,"kt"=kt)
  return(list("time"=timedbk,"acctable"=acctable))
}
########################################################################################################
#Function to simulate the specific n and datapoints combination N times
#input seedlength: the number of simulations N, 
#spline1D: the number of splines needed and we set it as 25 in all the simulations, 
#Mcfda is the column of initial mfpca scores, 
#pct is the percent of the distance,
#pctcfda is the percent of the cfda distance, 
#minPts is the initial minPts,
#min.nc is the minimum number of clusters, 
#max.nc is the maximum number of clusters
##Output : sample results for one simulation rexp, 
#time is the running time, acc is the average rand index as well as the adjusted rand index, 
#timese is the standard error of running time and 
#accse is the standard error for rand index and adjusted rand index
acc_table_sim=function(seedlength,k,n,datapoints,q,sparse1,sparse2,sparse3,scorevar1,scorevar2,ps,
                       st,et,spline1D,Mcfda,knnum,pct, pctcfda,minPts,min.nc,max.nc){
  
  onep=0.5
  twop=0.3
  threep=0.2
  # onep=0.75
  # twop=0.22
  # threep=0.03

  seedvec = matrix((1:(seedlength))+123, ncol = 1)
  results = apply(seedvec, 1, function(x) acc_table_graph(x, k,n,datapoints,q,sparse1,sparse2,sparse3,scorevar1,scorevar2,ps,
                                                          st,et,spline1D,Mcfda,knnum,pct, pctcfda,minPts,min.nc,max.nc))
  resultsexp=results[[1]]
  tcn100t250=sapply(results,"[[",1)
  accn100t250=sapply(results,"[[",2)
  timeave=apply(matrix(unlist(tcn100t250),nrow=2,ncol=seedlength),1,mean)
  timese=apply(matrix(unlist(tcn100t250),nrow=2,ncol=seedlength),1,sd)/sqrt(round(onep*n)+round(twop*n)+round(threep*n))
  accmean=apply(accn100t250,1,mean)
  accse=apply(accn100t250,1,sd)/sqrt(round(onep*n)+round(twop*n)+round(threep*n))
  return(list("rexp"=resultsexp,"time"=timeave,
              "acc"=accmean,"timese"=timese,"accse"=accse))
}
#########################################################################################################
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
#########################################################################################################
#Function to input the Latent curves Z and output the probability curves in catfda package
#Input: 3 D array, n * datapoints * (Q-1) where Q is the number of the category 
#Output: 3 D array, n * datapoints * Q probability curves, one per category
phatf=function (Zlatent) {
  numcat = dim(Zlatent)[3]+1
  n = dim(Zlatent)[1]
  nt = dim(Zlatent)[2]
  phatarray = array(0, dim = c(n, nt, numcat))
  denomarray = array(0, dim = c(n, nt, (numcat - 1)))
  for (i in 1:(numcat - 1)) {
    denomarray[, , i] = exp(Zlatent[, , i])
    
  }
  demsum = 1 + apply(denomarray, c(1, 2), sum)
  for (i in 1:(numcat - 1)) {
    phatarray[, , i] = exp(Zlatent[, , i])/demsum
  }
  sump = apply(phatarray, c(1, 2), sum)
  phatarray[, , numcat] = 1 - sump
  return(phat = phatarray)
}
########################################################################################################
#Function to make a functional data object
#input ufdata is observed curves, t is the time interval
#output is functional data object
mfundata=function(ufdata,t){
  mvdata=funData::funData(argvals = list(t), X = ufdata)
  mvdata
}
##########################################################################################
#Function to reproduce Results from Table 1 in the manuscript
#Input: true curve, estimated curve in n times datapoints dimension
#Output: mse and relative root MSE
mse_bw_matrix=function(truecurve,estcurve,datapoints,n){
  mse=sum(apply((truecurve-estcurve)^2,1, function(x) {sum(x)/datapoints}))/n
  gv=(truecurve-estcurve)^2
  gvn=gv/apply(gv,1,function(x){sqrt(sum(x^2))})
  relmse=sum(apply(gvn,1, function(x) {sum(x)/datapoints}))/n
  return(list("mse"=mse,"relmse"=relmse))
}
############################################################################################
#Function to get the relative mse between the curves
acc_table_mse=function(seed,k,n,datapoints,q,sparse1,sparse2,sparse3,scorevar1,scorevar2,ps,st,et,
                         spline1D,Mcfda,knnum,pct, pctcfda,minPts,min.nc,max.nc){
  seed=seed
  clustern50t250=generate_data_scenario(k=k,n=n,datapoints,sparse1,scorevar1,ps,seed=seed,st,et)
  clustern50t250p2=generate_data_scenario(k=k,n=n,datapoints,sparse2,scorevar1,ps,seed=seed,st,et)
  clustern50t250p3=generate_data_scenario(k=k,n=n,datapoints,sparse3,scorevar2,ps,seed=seed,st,et)
  
  #######################################
  Z_i1hatstar=clustern50t250$Est$EstimateZ_i1
  Z_i2hatstar=clustern50t250$Est$EstimateZ_i2
  
  Zihatone=clustern50t250$Trueth$TrueZ_i1
  Zihattwo=clustern50t250$Trueth$TrueZ_i2
  
  
  simp_i1=clustern50t250$Trueth$Truep_i1
  simp_i2=clustern50t250$Trueth$Truep_i2
  simp_i3=clustern50t250$Trueth$Truep_i3
  simp_i1_est=clustern50t250$Est$Estimatep_i1
  simp_i2_est=clustern50t250$Est$Estimatep_i2
  simp_i3_est=clustern50t250$Est$Estimatep_i3
  
  
  z11=mse_bw_matrix(Z_i1hatstar,Zihatone,datapoints,n)$relmse
  z12=mse_bw_matrix(Z_i2hatstar,Zihattwo,datapoints,n)$relmse
  p11=mse_bw_matrix(simp_i1,simp_i1_est,datapoints,n)$relmse
  p12=mse_bw_matrix(simp_i2,simp_i2_est,datapoints,n)$relmse
  p13=mse_bw_matrix(simp_i3,simp_i3_est,datapoints,n)$relmse
  #########################################
  Z_i1hatstar2=clustern50t250p2$Est$EstimateZ_i1
  Z_i2hatstar2=clustern50t250p2$Est$EstimateZ_i2
  
  Zihatone2=clustern50t250p2$Trueth$TrueZ_i1
  Zihattwo2=clustern50t250p2$Trueth$TrueZ_i2
  
  
  simp_i12=clustern50t250p2$Trueth$Truep_i1
  simp_i22=clustern50t250p2$Trueth$Truep_i2
  simp_i32=clustern50t250p2$Trueth$Truep_i3
  simp_i1_est2=clustern50t250p2$Est$Estimatep_i1
  simp_i2_est2=clustern50t250p2$Est$Estimatep_i2
  simp_i3_est2=clustern50t250p2$Est$Estimatep_i3
  
  
  z21=mse_bw_matrix(Z_i1hatstar2,Zihatone2,datapoints,n)$relmse
  z22=mse_bw_matrix(Z_i2hatstar2,Zihattwo2,datapoints,n)$relmse
  p21=mse_bw_matrix(simp_i12,simp_i1_est2,datapoints,n)$relmse
  p22=mse_bw_matrix(simp_i22,simp_i2_est2,datapoints,n)$relmse
  p23=mse_bw_matrix(simp_i32,simp_i3_est2,datapoints,n)$relmse
  
  ############################################
  Z_i1hatstar3=clustern50t250p3$Est$EstimateZ_i1
  Z_i2hatstar3=clustern50t250p3$Est$EstimateZ_i2
  
  Zihatone3=clustern50t250p3$Trueth$TrueZ_i1
  Zihattwo3=clustern50t250p3$Trueth$TrueZ_i2
  
  
  
  simp_i13=clustern50t250p3$Trueth$Truep_i1
  simp_i23=clustern50t250p3$Trueth$Truep_i2
  simp_i33=clustern50t250p3$Trueth$Truep_i3
  simp_i1_est3=clustern50t250p3$Est$Estimatep_i1
  simp_i2_est3=clustern50t250p3$Est$Estimatep_i2
  simp_i3_est3=clustern50t250p3$Est$Estimatep_i3
  
  
  z31=mse_bw_matrix(Z_i1hatstar3,Zihatone3,datapoints,n)$relmse
  z32=mse_bw_matrix(Z_i2hatstar3,Zihattwo3,datapoints,n)$relmse
  p31=mse_bw_matrix(simp_i13,simp_i1_est3,datapoints,n)$relmse
  p32=mse_bw_matrix(simp_i23,simp_i2_est3,datapoints,n)$relmse
  p33=mse_bw_matrix(simp_i33,simp_i3_est3,datapoints,n)$relmse
  #################################################
  acctable=cbind(z11,z12,p11,p12,p13)
  acctablep=cbind(z21,z22,p21,p22,p23)
  acctable3=cbind(z31,z32,p31,p32,p33)
  #############################################################################
  
  return(list("acctable"=acctable,"acctable2"=acctablep,"acctable3"=acctable3))
}
################################################################################
#Function to simulate the mse results for the N times where N is the seedlength
acc_table_simmse=function(seedlength,k,n,datapoints,q,sparse1,sparse2,sparse3,scorevar1,scorevar2,ps,
                       st,et,spline1D,Mcfda,knnum,pct, pctcfda,minPts,min.nc,max.nc){
  seedvec = matrix((1:(seedlength))+123, ncol = 1)
  results = apply(seedvec, 1, function(x) acc_table_graph(x, k,n,datapoints,q,sparse1,sparse2,sparse3,scorevar1,scorevar2,ps,
                                                          st,et,spline1D,Mcfda,knnum,pct, pctcfda,minPts,min.nc,max.nc))
  #resultsexp=results[[1]]
  accn100t2501=sapply(results,"[[",1)
  accn100t2502=sapply(results,"[[",2)
  accn100t2503=sapply(results,"[[",3)
  onezp=apply(accn100t2501,1,mean)
  twozp=apply(accn100t2502,1,mean)
  threezp=apply(accn100t2503,1,mean)
  
  onesd=apply(accn100t2501,1,sd)/sqrt(n)
  twosd=apply(accn100t2502,1,sd)/sqrt(n)
  threesd=apply(accn100t2503,1,sd)/sqrt(n)
  return(list("onezp"=onezp,"twozp"=twozp,"threezp"=threezp,"onesd"=onesd,"twosd"=twosd,"threesd"=threesd))
}
##########################################################################################################
##########################################################################################################
#Part III: sample code to reproduce the clusering results for Table 2 from the manuscript
#for demonstration purpose, this is the only one choice of n and datapoints combination, one simulation
#for Figure 2 and Figure 3 in the manuscript, n=1000, t=2000
#The results in Table 2 from manuscript include 9 different combination where n=100,500, 1000, datapoints=250, 500, 2000
#And each combination is repeated 100 times for all methods except only 5 times for cfdadbscan method

#generate data and plot the true, estimate curves from manuscript setting one, two and three
###################################
#this can be changed based on one of the 9 different combinations from the manuscript
#keep everything else fixed as it is to reproduce the results from the manuscript
n=1000 #total number of curves
datapoints=2000 #number of time points
scaleeps=2.6 #this is the scale for optimal epsilon, it's set to 2.6 for n=1000, datapoints=2000; 
#it can vary based on the domain knowledge 
###################################
p1=0.5 #proportion of the curves from cluster one
p2=0.3 #proportion of the curves from cluster two
p3=0.2 #proportion of the curves from cluster three

###########
seed=123 #seed
sparse1=0
scorevar1=1
sparse2=5
sparse3=7
scorevar3=3
st=0.01
et=0.99
k=3
ps=1
cluster1=generate_data_scenario(k=k,n=round(n*p1),datapoints,sparse1,
                                   scorevar1,ps,seed=seed,st,et)


cluster2=generate_data_scenario(k=k,n=round(n*p2),datapoints,sparse2,
                                   scorevar1,ps,seed=seed,st,et)
cluster3=generate_data_scenario(k=k,n=round(n*p3),datapoints,sparse3,
                                   scorevar3,ps,seed=seed,st,et)
truet=mapply(rbind,cluster1$Trueth,cluster2$Trueth,cluster3$Trueth,SIMPLIFY=FALSE,USE.NAMES = TRUE)
estt=mapply(rbind,cluster1$Est,cluster2$Est,cluster3$Est,SIMPLIFY=FALSE,USE.NAMES = TRUE)
combt=list("Trueth"=truet, "Est"=estt)
catfddata=combt$Trueth$Truecatcurve
save(catfddata,file="catfddata.RData")


######
#plot first true latent curves Z_1 by clusters
par(mfrow=c(2,2))
matplot(seq(st,et,length=datapoints),t(combt$Trueth$TrueZ_i1)[,1:round(n*p1)],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(-5.5,10))
matlines(seq(st,et,length=datapoints),t(combt$Trueth$TrueZ_i1)[,(round(n*p1)+1):(round(n*p1)+round(n*p2))],col="green",type = "l",xlab="",ylab="",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
matlines(seq(st,et,length=datapoints),t(combt$Trueth$TrueZ_i1)[,(round(n*p1)+round(n*p2)+1):n],col="blue",type = "l",xlab="",ylab="",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)

#plot second true latent curves Z_2. by clusters
matplot(seq(0,1,length=datapoints),t(combt$Trueth$TrueZ_i2)[,1:round(n*p1)],col="red",type = "l",xlab="",ylab="",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(-10,10))
matlines(seq(0,1,length=datapoints),t(combt$Trueth$TrueZ_i2)[,(round(n*p1)+1):(round(n*p1)+round(n*p2))],col="green",type = "l",xlab="",ylab="",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
matlines(seq(0,1,length=datapoints),t(combt$Trueth$TrueZ_i2)[,(round(n*p1)+round(n*p2)+1):n],col="blue",type = "l",xlab="",ylab="",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
######

par(mfrow=c(2,2))
######
#plot first estimated latent curves Z_1 by clusters
matplot(seq(st,et,length=datapoints),t(combt$Est$EstimateZ_i1)[,1:round(n*p1)],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(-10,10))
matlines(seq(st,et,length=datapoints),t(combt$Est$EstimateZ_i1)[,(round(n*p1)+1):(round(n*p1)+round(n*p2))],col="green",type = "l",xlab="",ylab="",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
matlines(seq(0,1,length=datapoints),t(combt$Est$EstimateZ_i1)[,(round(n*p1)+round(n*p2)+1):n],col="blue",type = "l",xlab="",ylab="",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)


#plot first estimated latent curves Z_2 by clusters
matplot(seq(0,1,length=datapoints),t(combt$Est$EstimateZ_i2)[,1:round(n*p1)],col="red",type = "l",xlab="",ylab="",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(-10,10))
matlines(seq(0,1,length=datapoints),t(combt$Est$EstimateZ_i2)[,(round(n*p1)+1):(round(n*p1)+round(n*p2))],col="green",type = "l",xlab="",ylab="",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
matlines(seq(0,1,length=datapoints),t(combt$Est$EstimateZ_i2)[,(round(n*p1)+round(n*p2)+1):n],col="blue",type = "l",xlab="",ylab="",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)

#####
#plot first truelatent curves p_1 by clusters
par(mfrow=c(2,3))
matplot(seq(st,et,length=datapoints),t(combt$Trueth$Truep_i1)[,1:20],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(0,1),xlim=c(0,1))
matlines(seq(st,et,length=datapoints),t(combt$Trueth$Truep_i1)[,(500+1):(520)],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
matlines(seq(st,et,length=datapoints),t(combt$Trueth$Truep_i1)[,(800+1):820],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)

#plot first truelatent curves p_2 by clusters
matplot(seq(st,et,length=datapoints),t(combt$Trueth$Truep_i2)[,1:20],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(0,1))
matlines(seq(st,et,length=datapoints),t(combt$Trueth$Truep_i2)[,(500+1):(520)],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
matlines(seq(st,et,length=datapoints),t(combt$Trueth$Truep_i2)[,(800+1):820],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)

##plot first truelatent curves p_3 by clusters
matplot(seq(st,et,length=datapoints),t(combt$Trueth$Truep_i3)[,1:20],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(0,1))
matlines(seq(st,et,length=datapoints),t(combt$Trueth$Truep_i3)[,(500+1):(520)],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
matlines(seq(st,et,length=datapoints),t(combt$Trueth$Truep_i3)[,(800+1):820],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)


#selected 20 subjects
#plot first estimated latent curves p_1 by clusters
matplot(seq(st,et,length=datapoints),t(combt$Est$Estimatep_i1)[,1:20],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(0,1))
matlines(seq(st,et,length=datapoints),t(combt$Est$Estimatep_i1)[,(500+1):(520)],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
matlines(seq(0,1,length=datapoints),t(combt$Est$Estimatep_i1)[,(800+1):820],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
######

#plot first estimated latent curves p_2 by clusters
matplot(seq(st,et,length=datapoints),t(combt$Est$Estimatep_i2)[,1:20],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(0,1))
matlines(seq(st,et,length=datapoints),t(combt$Est$Estimatep_i2)[,(500+1):(520)],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
matlines(seq(0,1,length=datapoints),t(combt$Est$Estimatep_i2)[,(800+1):820],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
######

#plot first estimated latent curves p_3 by clusters
matplot(seq(st,et,length=datapoints),t(combt$Est$Estimatep_i3)[,1:20],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(0,1))
matlines(seq(st,et,length=datapoints),t(combt$Est$Estimatep_i3)[,(500+1):(520)],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
matlines(seq(0,1,length=datapoints),t(combt$Est$Estimatep_i3)[,(800+1):820],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
######
#All p curves
#plot first estimated latent curves p_1 by clusters
matplot(seq(st,et,length=datapoints),t(combt$Est$Estimatep_i1)[,1:round(n*p1)],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(0,1))
matlines(seq(st,et,length=datapoints),t(combt$Est$Estimatep_i1)[,(round(n*p1)+1):(round(n*p1)+round(n*p2))],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
matlines(seq(0,1,length=datapoints),t(combt$Est$Estimatep_i1)[,(round(n*p1)+round(n*p2)+1):n],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
######
#plot first truelatent curves p_2 by clusters
matplot(seq(st,et,length=datapoints),t(combt$Trueth$Truep_i2)[,1:round(n*p1)],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(0,1))
matlines(seq(st,et,length=datapoints),t(combt$Trueth$Truep_i2)[,(round(n*p1)+1):(round(n*p1)+round(n*p2))],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
matlines(seq(st,et,length=datapoints),t(combt$Trueth$Truep_i2)[,(round(n*p1)+round(n*p2)+1):n],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
#plot first estimated latent curves p_2 by clusters
matplot(seq(st,et,length=datapoints),t(combt$Est$Estimatep_i2)[,1:round(n*p1)],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(0,1))
matlines(seq(st,et,length=datapoints),t(combt$Est$Estimatep_i2)[,(round(n*p1)+1):(round(n*p1)+round(n*p2))],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
matlines(seq(0,1,length=datapoints),t(combt$Est$Estimatep_i2)[,(round(n*p1)+round(n*p2)+1):n],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
######
##plot first truelatent curves p_3 by clusters
matplot(seq(st,et,length=datapoints),t(combt$Trueth$Truep_i3)[,1:round(n*p1)],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(0,1))
matlines(seq(st,et,length=datapoints),t(combt$Trueth$Truep_i3)[,(round(n*p1)+1):(round(n*p1)+round(n*p2))],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
matlines(seq(st,et,length=datapoints),t(combt$Trueth$Truep_i3)[,(round(n*p1)+round(n*p2)+1):n],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
#plot first estimated latent curves p_3 by clusters
matplot(seq(st,et,length=datapoints),t(combt$Est$Estimatep_i3)[,1:round(n*p1)],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(0,1))
matlines(seq(st,et,length=datapoints),t(combt$Est$Estimatep_i3)[,(round(n*p1)+1):(round(n*p1)+round(n*p2))],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
matlines(seq(0,1,length=datapoints),t(combt$Est$Estimatep_i3)[,(round(n*p1)+round(n*p2)+1):n],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
######
########
#plot the categories by cluster
par(mfrow=c(3,1))
hist((combt$Trueth$Truecatcurve[1:500,]),main="Categories for Cluster One",xlab="",cex.lab=1.5, cex.axis=2)
hist((combt$Trueth$Truecatcurve[501:800,]),main="Categories for Cluster Two",xlab="",cex.lab=1.5, cex.axis=2)
hist((combt$Trueth$Truecatcurve[801:1000,]),main="Categories for Cluster Three",xlab="",cex.lab=1.5, cex.axis=2)

#plot category pie chart
###########
#Cluster One
pie1=as.data.frame(table((combt$Trueth$Truecatcurve[1:500,])))
pie1$pct=round(pie1$Freq/(sum(pie1$Freq)),2)
pie1$labels= scales::percent(pie1$pct)
colnames(pie1)[colnames(pie1) == "Var1"] <- "Category"

ggplot(pie1, aes(x = "", y = pct, fill = Category)) +
  geom_col() +
  geom_text(aes(label = labels),size=10,
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+ggtitle("Cluster One")+theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size = 20))
############
#Cluster Two
pie2=as.data.frame(table((combt$Trueth$Truecatcurve[501:800,])))
pie2$pct=round(pie2$Freq/(sum(pie2$Freq)),2)
pie2$labels= scales::percent(pie2$pct)
colnames(pie2)[colnames(pie2) == "Var1"] <- "Category"

ggplot(pie2, aes(x = "", y = pct, fill = Category)) +
  geom_col() +
  geom_text(aes(label = labels),size=10,
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+ggtitle("Cluster Two")+theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size = 20))
##########
pie3=as.data.frame(table((combt$Trueth$Truecatcurve[801:1000,])))
pie3$pct=round(pie3$Freq/(sum(pie3$Freq)),2)
pie3$labels= scales::percent(pie3$pct)
colnames(pie3)[colnames(pie3) == "Var1"] <- "Category"

ggplot(pie3, aes(x = "", y = pct, fill = Category)) +
  geom_col() +
  geom_text(aes(label = labels),size=10,
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+ggtitle("Cluster Three")+theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size = 20))

########
#Table 1 RMSE results
#collect true an estimated curves from cluster one
Z_i1hatstar=cluster1$Est$EstimateZ_i1
Z_i2hatstar=cluster1$Est$EstimateZ_i2

Zihatone=cluster1$Trueth$TrueZ_i1
Zihattwo=cluster1$Trueth$TrueZ_i2


simp_i1=cluster1$Trueth$Truep_i1
simp_i2=cluster1$Trueth$Truep_i2
simp_i3=cluster1$Trueth$Truep_i3
simp_i1_est=cluster1$Est$Estimatep_i1
simp_i2_est=cluster1$Est$Estimatep_i2
simp_i3_est=cluster1$Est$Estimatep_i3

#Find RMSE for Z and p curves
z11=mse_bw_matrix(Z_i1hatstar,Zihatone,datapoints,n)$mse
z12=mse_bw_matrix(Z_i2hatstar,Zihattwo,datapoints,n)$mse
p11=mse_bw_matrix(simp_i1,simp_i1_est,datapoints,n)$mse
p12=mse_bw_matrix(simp_i2,simp_i2_est,datapoints,n)$mse
p13=mse_bw_matrix(simp_i3,simp_i3_est,datapoints,n)$mse
#########################################
#collect true an estimated curves from cluster two
Z_i1hatstar2=cluster2$Est$EstimateZ_i1
Z_i2hatstar2=cluster2$Est$EstimateZ_i2

Zihatone2=cluster2$Trueth$TrueZ_i1
Zihattwo2=cluster2$Trueth$TrueZ_i2


simp_i12=cluster2$Trueth$Truep_i1
simp_i22=cluster2$Trueth$Truep_i2
simp_i32=cluster2$Trueth$Truep_i3
simp_i1_est2=cluster2$Est$Estimatep_i1
simp_i2_est2=cluster2$Est$Estimatep_i2
simp_i3_est2=cluster2$Est$Estimatep_i3

#Find RMSE for Z and p curves
z21=mse_bw_matrix(Z_i1hatstar2,Zihatone2,datapoints,n)$mse
z22=mse_bw_matrix(Z_i2hatstar2,Zihattwo2,datapoints,n)$mse
p21=mse_bw_matrix(simp_i12,simp_i1_est2,datapoints,n)$mse
p22=mse_bw_matrix(simp_i22,simp_i2_est2,datapoints,n)$mse
p23=mse_bw_matrix(simp_i32,simp_i3_est2,datapoints,n)$mse
############################################
#collect true an estimated curves from cluster three
Z_i1hatstar3=cluster3$Est$EstimateZ_i1
Z_i2hatstar3=cluster3$Est$EstimateZ_i2

Zihatone3=cluster3$Trueth$TrueZ_i1
Zihattwo3=cluster3$Trueth$TrueZ_i2



simp_i13=cluster3$Trueth$Truep_i1
simp_i23=cluster3$Trueth$Truep_i2
simp_i33=cluster3$Trueth$Truep_i3
simp_i1_est3=cluster3$Est$Estimatep_i1
simp_i2_est3=cluster3$Est$Estimatep_i2
simp_i3_est3=cluster3$Est$Estimatep_i3


# #Find RMSE for Z and p curves
z31=mse_bw_matrix(Z_i1hatstar3,Zihatone3,datapoints,n)$mse
z32=mse_bw_matrix(Z_i2hatstar3,Zihattwo3,datapoints,n)$mse
p31=mse_bw_matrix(simp_i13,simp_i1_est3,datapoints,n)$mse
p32=mse_bw_matrix(simp_i23,simp_i2_est3,datapoints,n)$mse
p33=mse_bw_matrix(simp_i33,simp_i3_est3,datapoints,n)$mse
####################
#One combination: n=1000, datapoints=2000 Table 1 result
acctable=cbind(z11,z12,p11,p12,p13)
acctablep=cbind(z21,z22,p21,p22,p23)
acctable3=cbind(z31,z32,p31,p32,p33)
#setting one
round(acctable,2)
#setting two
round(acctablep,2)
#setting three
round(acctable3,2)
#############################################################################

########
#use MFPCA to find scores and cluster using k-means
combz3=mfpca_ram(combt,2,datapoints,q=3,
                         n,st,et,seq(st,et,length=datapoints))$scores
##k-means results
tdatak2=data.frame(combz3)
reskmeans2 = NbClust::NbClust(data =  combz3, diss = NULL,
                              distance = "euclidean", min.nc = 2, max.nc = 5,
                              method = "kmeans",index="silhouette")
true_label <- c(rep(1,round(n*p1)),rep(2,round(n*p2)),rep(3,round(n*p3)))
correctk2=rand.index(true_label, reskmeans2$Best.partition)
correctk22=adj.rand.index(true_label, reskmeans2$Best.partition)
tdatak2$Cluster=as.factor(reskmeans2$Best.partition)
tpdk2 <- ggplot(tdatak2,aes(X1,X2,colour = Cluster)) + geom_point(aes(shape=Cluster),size=3)+ggtitle(paste0(" Users Cluster Results",'\n',"(",dim(tdatak2)[1]," Subjects",")")) +
  xlab("MFPCA Score1") + ylab("MFPCA Score2")+ theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size = 20))
fviz_nbclust(combz3,method = c("silhouette"),kmeans)+theme(text=element_text(size = 20))+
labs(subtitle = "Silhouette method")

########
c3=which(true_label==3)
cc3k=which(reskmeans2$Best.partition==3)
noisepk=length(cc3k)/length(c3)
#######


######
tpdk2
###dbscan results
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
res2 <- dbscan(combz3, eps =epsoptimal2 , minPts = minPts2)
#visualize
correctp2=rand.index(true_label, res2$cluster)
correctp22=adj.rand.index(true_label, res2$cluster)
tdata2=data.frame(combz3)
tdata2$Cluster=as.factor(res2$cluster)
tpd2 <- ggplot(tdata2,aes(X1,X2,colour = Cluster)) + geom_point(aes(shape=Cluster),size=3)+ggtitle(paste0("Users Cluster Results",'\n',"(",dim(tdata2)[1]," Subjects",")")) +
  xlab("MFPCA Score1") + ylab("MFPCA Score2")+ theme(plot.title = element_text(hjust = 0.5))+
  
  theme(text=element_text(size = 20))
tpd2
#####
#############add accuracy for noise group

cc1=which(res2$cluster==1)
cc2=which(res2$cluster==2)
cc3=which(res2$cluster==3)
noisep=max(length(which(c3 %in% cc1)),length(which(c3 %in% cc2)),length(which(c3 %in% cc3)))/length(c3)
##############


#collect rand index and adjusted rand index for both methods
datapd=c(correctk2,correctp2)
datapd1=c(correctk22,correctp22)
datapd2=c(noisepk,noisep)
df_tabled=rbind(datapd,datapd1,datapd2)
accdata=data.frame(df_tabled)
colnames(accdata)=c("k-means","dbscan")
rownames(accdata)=c("rand index","adjusted rand index","correct noise percentage")
#results for just kmeans and dbscan
accdata

#########
#use FADP to cluster 
#collect estimated curves
Zihatone=combt$Est$EstimateZ_i1
Zihattwo=combt$Est$EstimateZ_i2
fdarange = c(st, et)
fdabasis = fda::create.bspline.basis(fdarange, 25, 4)
fdatime = seq(st, et, length = datapoints)
fdafd1 = fda::smooth.basis(fdatime, t(Zihatone), fdabasis)$fd
fdafd2 = fda::smooth.basis(fdatime, t(Zihattwo), fdabasis)$fd
FADlist=list(fdafd1,fdafd2)
FADP <- FADPclust(fdata = FADlist, cluster = 2:5, method = "FADP1",
                  proportion = seq(0.02, 0.2, 0.02), f.cut = 0.15,
                  pve = 0.9, stats = "Avg.silhouette")
accz1z2=rand.index(true_label, FADP$clust)
accz1z2adj=adj.rand.index(true_label, FADP$clust)
#rand index
accz1z2
#adjusted rand index
accz1z2adj
#########
#add results for FADP
accdata$FADP=c(accz1z2,accz1z2adj)
accdata
#cfda results
Mcfda=2
combz4=cfda_score(combt$Trueth$Truecatcurve,datapoints=datapoints,M=Mcfda,st,et)
dimz4=dim(combz4)[2]
if (dimz4<=2){
  
  minPts4=4
  
}

if (dimz4>2){
  minPts4=dimz4+1
}
dist4=kNNdist(combz4, k = minPts4-1)
#using elbow to find optimal epsilon
distdataelbow4=data.frame(sort(dist4))
distdataelbow4$index=1:(dim(combz4)[1])
ipoint4 <- elbow(data = distdataelbow4)
epsoptimal4=ipoint4$sort.dist2._selected*scaleeps
res4 <- dbscan(combz4, eps =epsoptimal4 , minPts = minPts4)
correctp3=rand.index(true_label, res4$cluster)
correctp33=adj.rand.index(true_label, res4$cluster)
tdata4=data.frame(combz4)
tdata4$Cluster=as.factor(res4$cluster)
tpd4 <- ggplot(tdata4,aes(X1,X2,colour = Cluster)) + geom_point(aes(shape=Cluster),size=3)+ggtitle(paste0("Users Cluster Results",'\n',"(",dim(tdata4)[1]," Subjects",")")) +
  xlab("MFPCA Score1") + ylab("MFPCA Score2")+ theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size = 20))
tpd4
###summarize all four methods with the rand index and adjusted rand index results into a table
datapd3=c(correctp3,correctp33)
accdata$cfdadbscan=datapd3
#One combination n=1000, datapoints=2000 Table 2 result
accdata

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
# import pandas as pd
# import matplotlib.pyplot as plt
# import datetime
# import os
# import re
# import glob
# import os
# 
# #%%
# 
# absolute_path = os.path.dirname(__file__)
# #%%
# profile_df = pd.read_csv(absolute_path + '\\twitter_data\\Users_profile_data.csv')
# profile_df = profile_df.dropna()
# 
# #%%
# file = absolute_path + '\\twitter_data\\Tweets_of_4776users_CSVs\\'
# csv_files = glob.glob(os.path.join(file, "*.csv"))
# df_tweets = pd.DataFrame()
# 
# for f in csv_files:
#   temp = pd.read_csv(f, encoding='latin-1')
# df_tweets = df_tweets.append(temp, ignore_index=True)
# 
# df_tweets = df_tweets.add_suffix('_tweeted_data')
# print(df_tweets)
# #%%
# df_tweets['UserID_tweeted_data'] = df_tweets['UserID_tweeted_data'].str.strip('@')
# #%%
# def parser(fn):
#   txt = open(fn).read()
# preparse = re.findall('"\d+",\d+,".+?","[^\t\n\r\f\v]+?","[^\t\n\r\f\v]+?"', txt, re.DOTALL)
# parsedRows = []
# for line in preparse:
#   columns = line.split(',')
# output = {}
# 
# output['Index'] = columns.pop(0)
# output['Tweet_ID'] = columns.pop(0)
# output['UserID'] = columns.pop()
# output['Company'] = columns.pop()
# output['Text'] = ','.join(columns)
# 
# parsedRows.append(output)
# 
# return parsedRows
# 
# #parsed = [t.split(',') for t in preparse]
# #print(parsed)
# 
# #%%
# 
# data = parser(absolute_path + r'\twitter_data\EDITED - Reference_tweet_data.txt')
# reference_df = pd.DataFrame(data)
# 
# reference_df['Index'] = reference_df['Index'].str.strip('"')
# reference_df['Index'] = reference_df['Index'].astype(int)
# reference_df['UserID'] = reference_df['UserID'].str.strip('"')
# reference_df =reference_df.add_suffix('_ref_tweets')
# print ("\nUnique values :  \n",reference_df.nunique())
# 
# duplicated = reference_df[reference_df.duplicated(subset='UserID_ref_tweets', keep=False)]
# reference_without_duplicates = reference_df.drop_duplicates(
#   subset = ['UserID_ref_tweets', 'Company_ref_tweets'],
#   keep = 'last').reset_index(drop=True)
# #%%
# df = pd.merge(reference_without_duplicates, profile_df, left_on='UserID_ref_tweets', right_on='user_id')
# print ("\nFeatures : \n" ,df.columns.tolist())
# print ("\nMissing values :  ", df.isnull().sum().values.sum())
# print ("\nUnique values :  \n",df.nunique())
# 
# #%%
# finalDF = pd.merge(df_tweets,df, left_on='UserID_tweeted_data', right_on='UserID_ref_tweets')
# print ("\nFeatures : \n" ,finalDF.columns.tolist())
# print ("\nMissing values :  ", finalDF.isnull().sum().values.sum())
# print ("\nUnique values :  \n",finalDF.nunique())
# 
# #%%
# finalDF['text_tweeted_data'] = finalDF['text_tweeted_data'].astype(str)
# finalDF['Company_ref_tweets'] = finalDF['Company_ref_tweets'].astype(str)
# finalDF['Company_ref_tweets'] = finalDF['Company_ref_tweets'].str.strip('"')
# finalDF['mention'] = finalDF.apply(lambda x: x.Company_ref_tweets in x.text_tweeted_data, axis=1)
# print(finalDF['mention'].value_counts())
# 
# #%%
# finalDF['mention_type'] = 2
# finalDF.loc[finalDF['mention']==False,'mention_type'] = 1
# print(finalDF['mention_type'].value_counts())
# 
# #%%
# finalDF['DateTime2_tweeted_data'] = pd.to_datetime(finalDF['DateTime2_tweeted_data'])
# finalDF['DateTime_tweeted_data'] = finalDF.apply(lambda x: pd.Timestamp(x.DateTime2_tweeted_data), axis=1)
# finalDF['DateTime_tweeted_data'] = finalDF.apply(lambda x: x.DateTime_tweeted_data.timestamp(), axis=1)
# #%%
# """WRITE FORMATTED DATA TO CSV"""
# finalDF.to_csv(absolute_path + r'\final_output.csv')
# 
# #%%
# """IMPORT FORMATTED DATA"""
# finalDF = pd.read_csv(r'C:\Users\Rob\OneDrive\NCSU PhD\twitter\final_output.csv')
# 
# print(finalDF.columns)
# 
# #%%
# """REMOVE EXTRA COLUMNS"""
# final_column_list = ['DateTime2_tweeted_data', 'Index_ref_tweets', 'mention_type']
# finalDF = finalDF[finalDF.columns.intersection(final_column_list)]
# finalDF['mention_type'] = finalDF['mention_type'].astype(int)
# 
# #%%
# """lAST MONTH OF ACTIVITY"""
# 
# finalDF['DateTime2_tweeted_data'] = pd.to_datetime(finalDF['DateTime2_tweeted_data'])
# month_ago = finalDF['DateTime2_tweeted_data'].max() - datetime.timedelta(days=30)
# 
# last_30_days_DF = finalDF.query('DateTime2_tweeted_data > @month_ago')
# print(last_30_days_DF.head())
# 
# #%%
# 
# def round_time(dt: datetime, unit=20):
#   #seconds = dt - dt.date()
#   #unit_seconds = unit.total_seconds()
#   #rounded_seconds = seconds - (seconds % unit_seconds)
#   #return dt.date() + rounded_seconds
#   return dt.replace(second=0, microsecond=0, minute=dt.minute-(dt.minute%unit), hour=dt.hour)
# 
# """Bin Each User"""
# each_user = last_30_days_DF.groupby('Index_ref_tweets')
# output = pd.DataFrame()
# counter = 0
# start_time=round_time(dt=last_30_days_DF['DateTime2_tweeted_data'].min())
# end_time=round_time(dt=last_30_days_DF['DateTime2_tweeted_data'].max())
# print(start_time)
# print(end_time)
# 
# t_index = pd.DatetimeIndex(pd.date_range(start=start_time, end=end_time, freq="20min"))
# 
# for user_name, user_data in each_user:
#   print(f"Read {counter} out of {len(each_user)}")
# #df_bin=df_subset.set_index('DateTime2_tweeted_data',drop=False)
# binuser=user_data.groupby(pd.Grouper(key='DateTime2_tweeted_data', freq='20min')).max()
# binuser = binuser.reindex(t_index)
# binuser['mention_type'] = binuser['mention_type'].fillna(5)
# binuser['Index_ref_tweets'] = binuser['Index_ref_tweets'].fillna(user_name)
# binuser = binuser.transpose()
# print(binuser)
# output = output.append(binuser, ignore_index=True)
# counter += 1
# 
# #%%
# """SPLIT INTO THREE DATAFRAMES"""
# """Mentions"""
# mentions = output.replace(to_replace=[1, 5, 2], value=[0, 0, 1])
# mentions.to_csv(absolute_path + r'\mentions_dataset.csv')
# 
# no_mentions = output.replace(to_replace=[1, 5, 2], value=[1, 0, 0])
# no_mentions.to_csv(absolute_path + r'\no_mentions_dataset.csv')
# 
# no_tweets = output.replace(to_replace=[1, 5, 2], value=[0, 1, 0])
# no_tweets.to_csv(absolute_path + r'\no_tweet_dataset.csv')
#Output the three csv files: One file for each category
#######################################################################################################
#R code to pre process three csv files in the form of the categorical functional data
#Please download all the files from processedTwitterData and make sure the path is correct
X_i1r=read.csv(file='no_tweet_dataset.csv')[1:1000,-1]   #3875, 2162  to 1000,2161  1st col originally the id X
X_i2r=read.csv(file='no_mentions_dataset.csv')[1:1000,-1]  
X_i3r=read.csv(file='mentions_dataset.csv')[1:1000,-1]
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

#load("Twitterram.RData")
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
res2 <- dbscan(combz3, eps =epsoptimal2 , minPts = minPts2)
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


#############################################################################################################
############################################################################################################
#Part VI: catfdcluster function from catfda R package to cluster categorical functional data
#######################################################################
# Purpose: clustering categorical functional data
# Author:blank
# Date: 1/21/2023
#####################################################################

#\name{catfdcluster}
# \alias{catfdcluster}
# \title{cluster categorical functional data}
# 
# \description{
#   This function clusters categorical functional data, where each series are categorical values 
#observed on a dense interval. It uses kmeans or DBSCAN to cluster on the multivariate functional 
#principal component scores extracted from the multivariate latent cruves, which induce the categorical functional data. Ramsay method concatenates the multivariate latent curves and the final scores are the linear combinations of the scores from each component. Happ method uses univariate expansion and can consider the weight of each component or different domains.
#   
# }
# \usage{
#   catfdcluster(catfd,argval,splines1D,M,knnum,pct,minPts,max.nc,min.nc,method,numdim,scaleeps)
# }
# 
# 
# \arguments{
#   \itemize{
#     \item{catfd}{: matrix n rows is the number of individuals, T columns is the number of time total time points. Each row represents a categorical series observed on a dense interval at each point t. }
#     
#     \item {argval} {: vector length T. The observational time interval.}
#     \item {splines1D} {: scalar the number of univariate basis needed. Usually 25.}
#     \item{M}{: scalar The number of principa components. }
#     \item{knnum}{: scalar The number of neighbours to calculate the distance in DBSCAN}
#     \item{pct}{: scalar between 0 and 1: The percentage of the distance used to determin the epsilon in DBSCAN. }
#     \item{minPts}{: scalar The minimum number of points needed to define a core point in DBSCAN}
#     \item{max.nc}{: scalar maximum number of clusters wanted in Kmeans}
#     \item{min.nc}{: scalar minimum number of clusters wanted in Kmeans}
#     \item{method}{: "ramsay" or "happ" Ramsay concatenates the multivariate latent curves and the final scores are the linear combinations of the scores from each component. Happ can account the weight among different components and handle different domains.}
#     \item{numdim}{: "all" or "two". "all" uses all the dimensions of the mfpca scores that explain at least 95 percent of the variance, "two" uses only the first two dimensions}
#     \item{scaleeps}{: scalar The scale of the optimal epsilon, 1 means optimal epsilon from elbow method}
#   }
#   
# }
# 
# 
# \details{
#   See Clustering of categorical valued functional data with application to social media
# }
# \value{
#   \itemize{
#     \item{scores}{: matrix n rows is the number of individuals, k columns is the number of principal components}
#     \item{dbcluster}{: vector length n is the clustering results for each individual using dbscan}
#     \item{dbtable}{: table summary of the dbscan clustering results}
#     \item{kcluster}{: vector length n is the clustering results for each individual using kmeans}
#     \item{kmeantable}{: table summary of the kmeans clustering results}
#     \item{latentcurve}{: 3 D array n rows is the number of individuals, T columns is the total number of observational points, l: Total category Q minus 1 latent curves}
#     \item{meanfn}{: lists Q minus 1 vectors each has length T. The mean function for each one of the latent curves}
#     \item{eigenvalue}{: vector length k is each one of the eigen values from multivariate functional principal component analysis}
#     \item{eigenfn}{: matrix T rows is the number of total observational points, k columns is the number of principal components}
#     \item{probcurve}{: 3 D array n rows is the number of individual, T columns is the number of observational times, Q: the number of categories}
#   }
#   
#   
# }
# \references{
#   Ramsay, J. O., Silverman, B. W. (2005). Functional Data Analysis. Springer. ISBN: 9780387400808
#   
#   Happ, Clara, and Sonja Greven. Multivariate functional principal component analysis for data observed on different (dimensional) domains. Journal of the American Statistical Association 113.522 (2018): 649-659.
# }
# \author{For this manuscript : blank
# }
# 
# 
# \seealso{
#   
#   Clustering of categorical valued functional data with application to social media
# }
# \examples{
#   catclust=catfdcluster(matrix(sample(c(0,1,2),100*250,replace=TRUE),nrow=100,ncol=250),seq(0,1,length=250),25,3,3,0.9,4,5,2,"happ","two",scaleeps=1)
# }


catfdcluster=function(catfd,argval,splines1D,M,knnum,pct,minPts,max.nc,min.nc,method,numdim,scaleeps){
  st=min(argval)
  et=max(argval)
  datapoints=dim(catfd)[2]
  tolcat=table(catfd)
  catorder=order(tolcat,decreasing = TRUE)
  numcat=length(catorder)
  refcat=catorder[numcat]
  refmat=catfd
  refmat[refmat!=refcat]=0
  refmat[refmat==refcat]=1
  nsub=dim(catfd)[1]
  ntime=dim(catfd)[2]
  subdata=array(data=NA,c(nsub,ntime,numcat))
  for (i in 1:numcat){
    datacopy=catfd
    datacopy[datacopy!=i]=0
    datacopy[datacopy==i]=1
    subdata[,,i]=datacopy
  }
  t=seq(st,et,length=datapoints)
  #input observed X_i1 or X_i2 binary curves and return smoothed Z_i1hat, Z_i2hat
  Zihat=array(data=NA,c(nsub,ntime,numcat))
  for (i in 1:numcat){
    datacopy=subdata
    Zihat[,,i]=Z_ihat(datacopy[,,i],t)
  }
  
  Zihatstar=array(data=NA,c(nsub,ntime,numcat-1))
  for (i in 1:(numcat-1)){
    datacopy=Zihat
    Zihatstar[,,i]=Zihat[,,i]+log(1+exp(Zihat[,,numcat]))-log(1+exp(Zihat[,,i]))-Zihat[,,numcat]
  }
  
  phatmat=phatf(Zihatstar)
  
  
  
  #############################################
  if (method=="ramsay"){
    fdarange = c(st, et)
    fdabasis = fda::create.bspline.basis(fdarange, splines1D, 4)
    fdatime = seq(st, et, length = datapoints)
    Zihatstarnew = apply(Zihatstar, c(1, 3), t)
    fdafd = fda::smooth.basis(fdatime, Zihatstarnew, fdabasis)$fd
    nharm = M
    fdapcaList = fda::pca.fd(fdafd, nharm)
    
    finalM=which(cumsum(fdapcaList$values)/sum(fdapcaList$values)>=0.95)[1]
    fdapcaListfinal = fda::pca.fd(fdafd, finalM)
    scores_z = apply(fdapcaListfinal$scores, c(1, 2), sum)
    values=fdapcaListfinal$values
    meanfn=fdapcaList$meanfd
    eigenfn=fdapcaList$harmonics
    
    dims=dim(scores_z)[2]
    if (dims<2){
      #finalM=which(cumsum(fdapcaList$values)/sum(fdapcaList$values)>=0.98)[1]
      finalM=2
      fdapcaListfinal = fda::pca.fd(fdafd, finalM)
      scores_z = apply(fdapcaListfinal$scores, c(1, 2), sum)
      values=fdapcaListfinal$values
      meanfn=fdapcaList$meanfd
      eigenfn=fdapcaList$harmonics
    }
    
  }
  
  if (method=="happ"){
    vecapply=matrix(1:(dim(Zihatstar)[3]),ncol=1)
    mfdataused=apply(vecapply,1,function(x) {mfundata(Zihatstar[,,x],t)})
    mvdata=funData::multiFunData(mfdataused)
    
    
    uniexpan=list()
    # MFPCA based on univariate FPCA Z_ihat
    for (i in 1:(numcat-1)){
      uniexpan[[i]]=list(type = "splines1D", k = splines1D)
    }
    
    # MFPCA based on univariate FPCA Z_ihat
    uFPCA <- MFPCA::MFPCA(mvdata, M = M, uniExpansions = uniexpan)
    
    finalM=which(cumsum(uFPCA$values)/sum(uFPCA$values)>=0.95)[1]
    uFPCA <- MFPCA::MFPCA(mvdata, M = finalM, uniExpansions = uniexpan)
    scores_z=uFPCA$scores
    values=uFPCA$values
    eigenfn=uFPCA$functions
    meanfn=uFPCA$meanFunction}
  
  dimz=dim(scores_z)[2]
  if (dimz<=2){minPts=4}
  if (dimz>2){minPts=dimz+1}
  dist=dbscan::kNNdist(scores_z, k = minPts-1)
  ninty5p=quantile(dist, probs = pct)
  
  distdataelbow=data.frame(sort(dist))
  distdataelbow$index=1:(dim(scores_z)[1])
  ipoint <- elbow(data = distdataelbow)
  epsoptimal=ipoint$sort.dist._selected*scaleeps
  ##############################
  
  distdata=data.frame(sort(dist))
  distdata$index=1:dim(distdata)[1]
  
  ##############################################
  
  if (numdim=="all"){
    #dbscan
    res <- dbscan::dbscan(scores_z, eps =epsoptimal , minPts = minPts)
    
    clustertable=table(res$cluster)
    
    #########Kmeans
    reskmeans=NbClust::NbClust(data = scores_z, diss = NULL, distance = "euclidean",
                               min.nc = min.nc, max.nc = max.nc, method = "kmeans",index="silhouette")
    
    clustertablek=table(reskmeans$Best.partition)
  }
  if (numdim=="two"){
    #dbscan
    res <- dbscan::dbscan(scores_z[,1:2], eps =epsoptimal , minPts = minPts)
    
    clustertable=table(res$cluster)
    
    #########Kmeans
    reskmeans=NbClust::NbClust(data = scores_z[,1:2], diss = NULL, distance = "euclidean",
                               min.nc = min.nc, max.nc = max.nc, method = "kmeans",index="silhouette")
    clustertablek=table(reskmeans$Best.partition)
  }
  return(list("scores"=scores_z,"dbcluster"=res$cluster,"dbtable"=clustertable,
              "kcluster"=reskmeans$Best.partition,"kmeantable"=clustertablek,
              "latentcurve"=Zihatstar,"eigenvalue"=values,"meanfn"=meanfn,
              "eigenfn"=eigenfn,"probcurve"=phatmat))}
################################################################################################
###########################################################################################
#Part VII sample code to cluster categorical funcitonal data using the catfdcluster function from catffd package
#here n=1000, datapoints=2000, Q=3 categories and scaleeps=2.6
#Please download the catfddata.Rdata file
load("catfdadata.RData")
st=Sys.time()
catfdclust=catfdcluster(catfddata,seq(0,1,length=2000),25,3,3,0.9,4,5,2,"ramsay","all",scaleeps=2.6)

et=Sys.time()
timeeplase=et-st
timeeplase

##plot the results
clusterdata=data.frame(catfdclust$scores[,c(1,2)])
clusterdata$dbCluster=as.factor(catfdclust$dbcluster)
clusterdata$kCluster=as.factor(catfdclust$kcluster)
clustervisual <- ggplot(clusterdata,aes(X1,X2,colour = dbCluster)) + geom_point(aes(shape=dbCluster),size=3)+ggtitle(paste0("Users Cluster Results",'\n',"(",dim(clusterdata)[1]," Subjects",")")) +
  xlab("MFPCA Score1") + ylab("MFPCA Score2")+ theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size = 20))
clustervisual

clustervisualk <- ggplot(clusterdata,aes(X1,X2,colour = kCluster)) + geom_point(aes(shape=kCluster),size=3)+ggtitle(paste0("Users Cluster Results",'\n',"(",dim(clusterdata)[1]," Subjects",")")) +
  xlab("MFPCA Score1") + ylab("MFPCA Score2")+ theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size = 20))
clustervisualk






