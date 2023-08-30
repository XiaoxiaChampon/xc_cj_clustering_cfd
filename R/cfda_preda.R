#2/22/2023

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
# library(knitr)
# library(kableExtra)
library(dbscan)
library(cfda)
library(tidyr)
library(dplyr)
# library( fdm2id)
library(fossil)


#Function to return the logit
logit <- function(x){
  return(log(x/(1-x)))
}

#####
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
# Fitted values from the game function for subject z 
#
#####
#get smoothed curves
regression_g = function(z, Curves, tt, k=10, method="ML"){   #changed from 10 to 25
  z1 = Curves[z,]
  gam1 <- gam(z1~s(tt, bs = "cr", m=2, k = k),
              family="binomial", method = method,control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),optimizer=c("outer","bfgs"))
  return(gam1$fitted.values)
}




#####
#Function: wrapper function for bayesglm() which outputs expected coefficients for the gaussian responses
#
#Inputs: 
# z : index z = 1,...,N 
# dta : data frame contain the subject id, mean function, eigenfunctions and latent continuous curves
# lm_structure: formula displaying the bayesglm function
# eigen_vals1: eigen_values for the ksi coefficients
#
#Output: 
# Fitted values from the game function for subject z 
#get bayesian scores

regression_bf2 = function(z, dta, lm_structure, eigen_vals1){
  bayesglm1 = bayesglm(lm_structure,
                       # family = binomial(link = "logit"),
                       family=gaussian,
                       data = subset(dta, dta$id==z),
                       prior.scale = sqrt(eigen_vals1), #set scales
                       prior.df = Inf, #normal priors
                       scaled = F ) #Do not adjust the scales
  return(bayesglm1$coefficients)
}

######################################################################################################################
#original before 4/6/2022
##
#Step 1 of the proposed method in Anthony paper one: smoothing using link function
##

#New function to output predicted L-1 latent curves:
#input observed X_i1 or X_i2 binary curves and return smoothed Z_i1hat, Z_i2hat
Z_ihat=function(Curves_train,tt){
  N_train=dim(Curves_train)[1]
  vec = matrix(1:(N_train), ncol = 1)
  smoothed_x = logit(t(apply(vec, 1, function(x) regression_g(x, Curves_train, tt))))
  smoothed_x
}
########################################################################################################################

#####################################################################################################
# 4//6/ 2022
##
#Step 1 of the proposed method in Anthony paper one: smoothing using link function
##

#New function to output predicted L latent probability curves:
#input observed X_i1 or X_i2 binary curves and return smoothed p_i1hat, p_i2hat


#########################################################################################################




score=function(mu,sd){
  rnorm(1,mean=mu,sd=sd)
}

#n number of subjects
#datapoints
#sparse=1 yes   sparse=0 no
#scorevar=2 bigger var , scorevar=1 smaller var
#ps=1 find z1,z2,z3, find p1,p2,p3, logp1-logp3
#ps=2 find z1,z2,z3, find p1,p2,p3=1-p1-p2 logp1-logp3
#ps=3 find z1,z2 staicu find z1hat z2hat
#ps=4 find z1,z2 staicu find z1hat z2hat but only use numerator
#k  #number of eigen functions
#q  #level of the categorical level

generate_data_scenario=function(k,n,datapoints,sparse,scorevar,ps,seed=123,st,et){
  k=k
  seed=seed
  st=st
  et=et
  scorevar=scorevar
  #k=3  #number of eigen functions
  q=3  #level of the categorical level
  
  if(sparse==1){
    mu_1=function(t){
      3.8+4*t  #-0.64+4*t    #4.5+4*t   #sparse 3.8+4*t
      
    }
    mu_2=function(t){
      1.5+4*t^2    #0.97+6*t^2
      
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
  
  
  if (sparse==0){
    mu_1=function(t){
      #3.8+4*t  
      -0.64+4*t    #4.5+4*t   #sparse 3.8+4*t
      
    }
    mu_2=function(t){
      #1.5+4*t^2    
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
  
  
  if (sparse==3){
    mu_1=function(t){
      3.8+4*t  #-0.64+4*t    #4.5+4*t   #sparse 3.8+4*t
      
    }
    mu_2=function(t){
      #1.5+4*t^2    
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
  
  
  
  
  mu_vec=rep(0,k)
  
  
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
        score_1=score(0,1)
        score_2=score(0,sqrt(1/2))
        score_3=score(0,1/2)
        
      }
      
      if (scorevar==2){
        #score varies based on i
        score_1=score(0,1)
        score_2=score(0,sqrt(1/3))
        score_3=score(0,1/3)
      }
      
      if (scorevar==3){
        #score varies based on i
        score_1=score(0,1)
        score_2=score(0,sqrt(1/4))
        score_3=score(0,1/4)
      }
      
      
      if (scorevar==4){
        #score varies based on i
        score_1=score(-0.5,1)
        score_2=score(1,sqrt(1/2))
        score_3=score(0.25,1/2)
      }
      
      score_vector=cbind(score_1,score_2,score_3)
    }
    
    
    
    if (k==4){
      
      if (scorevar==1){
        score_1=score(0,1)
        score_2=score(0,sqrt(1/2))
        score_3=score(0,1/2)
        score_4=score(0,sqrt(1/8))
        # cpve=cumsum(c(1,sqrt(1/2),1/2,sqrt(1/8)))/sum(c(1,sqrt(1/2),1/2,sqrt(1/8)))
        # cvar=c(1,sqrt(1/2),1/2,sqrt(1/8))
      }
      
      
      
      if (scorevar==2){
        score_1=score(0,1)
        score_2=score(0,sqrt(1/3))
        score_3=score(0,1/3)
        score_4=score(0,sqrt(1/27))
        # cpve=cumsum(c(1,sqrt(1/3),1/3,sqrt(1/27)))/sum(c(1,sqrt(1/3),1/3,sqrt(1/27)))
        # cvar=c(1,sqrt(1/3),1/3,sqrt(1/27))
      }
      
      
      
      if (scorevar==3){
        score_1=score(0,1)
        score_2=score(0,sqrt(1/4))
        score_3=score(0,1/4)
        score_4=score(0,sqrt(1/64))
        # cpve=cumsum(c(1,sqrt(1/4),1/4,sqrt(1/64)))/sum(c(1,sqrt(1/4),1/4,sqrt(1/64)))
        # cvar=c(1,sqrt(1/4),1/4,sqrt(1/64))
      }
      score_vector=cbind(score_1,score_2,score_3,score_4)
      
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
    
    
    #X_i varies based on i
    #X_i=matrix(rep(1,k*length(t)),nrow=k,ncol=length(t))
    
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
  Z_i3hat=Z_ihat(X_i3,t)
  
  
  
  if(ps==1){ 
    Z_i1hatstar=Z_i1hat+log(1+exp(Z_i3hat))-log(1+exp(Z_i1hat))-Z_i3hat
    Z_i2hatstar=Z_i2hat+log(1+exp(Z_i3hat))-log(1+exp(Z_i2hat))-Z_i3hat
  }
  
  
  if(ps==2){
    smooth_p1=(1+exp(-Z_i1hat))^(-1)
    smooth_p2=(1+exp(-Z_i2hat))^(-1)
    smooth_p3=1-smooth_p1-smooth_p2
    Z_i1hatstar=log(smooth_p1)-log(smooth_p3)
    Z_i2hatstar=log(smooth_p2)-log(smooth_p3)
  }
  
  if(ps==3){
    common_check=1-exp(Z_i1hat+Z_i2hat)-exp(Z_i2hat)
    common_check[common_check<=0]=0.001
    Z_i1hatstar=Z_i1hat+log(1+exp(Z_i2hat))-log(common_check)
    Z_i2hatstar=Z_i2hat+log(1+exp(Z_i1hat))-log(common_check)
  }
  
  
  
  if (ps==4){
    common_check=1-exp(Z_i1hat+Z_i2hat)-exp(Z_i2hat)
    common_check[common_check<=0]=0.001
    Z_i1hatstar=Z_i1hat-log(common_check)
    Z_i2hatstar=Z_i2hat+log(1+exp(Z_i1hat))-log(common_check)
  }
  
  
  
  p_i1hat=p_ihat(Z_i1hatstar,Z_i2hatstar)$p_1hatmatrix
  p_i2hat=p_ihat(Z_i1hatstar,Z_i2hatstar)$p_2hatmatrix
  p_i3hat=p_ihat(Z_i1hatstar,Z_i2hatstar)$p_3hatmatrix
  
  
  # truel=list("TrueX1"=X_i1,"TrueX2"=X_i2,"TrueX3"=X_i3,"TrueZ_i1"=Z_i1,"TrueZ_i2"=Z_i2,"Truep_i1"=p_i1,"Truep_i2"=p_i2,"Truep_i3"=p_i3,"Truescore_matrix"=score_matrix,"Truecatcurve"=X_nt)
  truel=list("TrueZ_i1"=Z_i1,"TrueZ_i2"=Z_i2,"Truep_i1"=p_i1,"Truep_i2"=p_i2,"Truep_i3"=p_i3,"Truescore_matrix"=score_matrix,"Truecatcurve"=X_nt,"comp1"=psi_k1,"comp2"=psi_k2)
  est=list("EstimateZ_i1"=Z_i1hatstar,"EstimateZ_i2"=Z_i2hatstar,"Estimatep_i1"=p_i1hat,"Estimatep_i2"=p_i2hat,"Estimatep_i3"=p_i3hat)
  return(list("Trueth"=truel,"Est"=est))
}



# MFPCA_result
baysglm_result=function(datagenerated,datapoints,q=3,n,st,et,k,splines1D,M){
  datapoints=datapoints
  t=seq(from = st,to = et, length=datapoints)
  #collect True Z_is, p_is and estimated latent curves Z_ihat, p_ihat
  #return(list(truevalues, estvalues))
  
  Z_i1hatstar=datagenerated$Est$EstimateZ_i1
  Z_i2hatstar=datagenerated$Est$EstimateZ_i2
  
  Zihatstar=array(c(Z_i1hatstar,Z_i2hatstar),c(round(0.7*n)+round(0.3*n)+round(0.04*n),datapoints,q-1))
  
  #####get z scores
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
  
  dims=dim(scores_z)[2]
  if (dims<2){
    finalM=which(cumsum(fdapcaList$values)/sum(fdapcaList$values)>=0.99)[1]
    fdapcaListfinal = fda::pca.fd(fdafd, finalM)
    scores_z = apply(fdapcaListfinal$scores, c(1, 2), sum)
  }
  #####################
  
  # fdarange = c(0, 2300)
  # fdabasis = fda::create.bspline.basis(fdarange, 105, 6)
  # fdatime = seq(0, 2300, length = 1401)
  # 
  # fdafd = fda::smooth.basis(fdatime, handwrit, fdabasis)$fd
  # nharm = 3
  # fdapcaList = fda::pca.fd(fdafd, nharm)
  # scores_z = apply(fdapcaList$scores, c(1, 2), sum)
  #######################
  
  
  
  return(list("Bayesian_scores"=scores_z,"values"=fdapcaListfinal$values))
  
}


####
#cfda
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
  data_long <- gather(cfda_new, time_t, state, paste0("time_",times[1]):paste0("time_",times[datapoints]),factor_key=TRUE)
  cfda_d=data_long%>%arrange(subjects)%>%rename(id=subjects)%>%dplyr::select(id,state)
  cfda_final=data.frame(cbind(cfda_d,time))
  Tmax=max(cfda_final$time)
  cfda_cut=cut_data(cfda_final, Tmax = Tmax)
  m <- 10
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
  fmca <- compute_optimal_encoding(cfda_cut, b, nCores = 16, verbose = FALSE)
  nPc90 <-which(cumsum(prop.table(fmca$eigenvalues)) > 0.95)[1]
  cfda_score=fmca$pc[, 1:nPc90]
  cfda_score[,1:M]
}


#function to output average score matrix for N simulations
#input the data set n*t*N, number of data points, M: number of the column of scores matrix
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
  # cfda_score_apply_final=apply(cfda_score_matrix,c(1,2),mean)
  # cfda_score_apply_final
}


#cfda_score=function(cfda_data,datapoints, M,st,et)


# ggsavepath="/ggsaves"
#q: number of category 3
#k: number of true eigen functions 3
acc_table_graph=function(seed,k,n,datapoints,q,sparse1,sparse2,sparse3,scorevar1,scorevar2,ps,st,et,
                         spline1D,Mcfda,knnum,pct, pctcfda,minPts,min.nc,max.nc){
  seed=seed
  n=n
  q=q
  k=k
  knnum=knnum
  minPts=minPts
  
  #check
  ######
  # k=3
  # n=100
  # datapoints=250
  # sparse1=0
  # sparse2=1
  # sparse3=3
  # scorevar1=1
  # scorevar2=2
  # ps=1
  # seed=123
  # st=0.01
  # et=0.99
  ##
  
  clustern50t250=generate_data_scenario(k=k,n=round(0.7*n),datapoints,sparse1,scorevar1,ps,seed=seed,st,et)
  clustern50t250p2=generate_data_scenario(k=k,n=round(0.3*n),datapoints,sparse2,scorevar2,ps,seed=seed,st,et)
  clustern50t250p3=generate_data_scenario(k=k,n=round(n*0.04),datapoints,sparse3,scorevar1,ps,seed=seed,st,et)
  
  
  truen100t250=mapply(rbind,clustern50t250$Trueth,clustern50t250p2$Trueth,clustern50t250p3$Trueth,SIMPLIFY=FALSE,USE.NAMES = TRUE)
  estn100t250=mapply(rbind,clustern50t250$Est,clustern50t250p2$Est,clustern50t250p3$Est,SIMPLIFY=FALSE,USE.NAMES = TRUE)
  combn100t250=list("Trueth"=truen100t250, "Est"=estn100t250)
  
  sp=Sys.time()
  #q=3
  #spline1D=25
  #Mcfda=3
  #combn100t250mfpca=baysglm_result(datagenerated=combn100t250,datapoints,q=q,n=n,st,et,k,spline1D,Mcfda)
  
  combn100t250cfda=cfda_score(truen100t250$Truecatcurve,datapoints=datapoints,M=Mcfda,st,et)
  ept=Sys.time()
  tp=ept-sp
  
  #datalist=list("1sthalf"=clustern50t250,"2ndhalf"=clustern50t250p2,"combmfpca"=combn100t250mfpca)
  
  ######
  dbst=Sys.time()
  combinedscore_z <- combn100t250cfda
  #knnum=3
  #pct=0.9
  #minPts=3
  min.nc=2
  max.nc=5
  dimz=dim(combinedscore_z)[2]
  if (dimz<=2){minPts=4}
  if (dimz>2){minPts=2*dimz+1}
  #dist=kNNdist(combinedscore_z, k = knnum)
  dist=kNNdist(combinedscore_z, k = minPts)
  ninty5p=quantile(dist, probs = pct)
  
  #########change to max increase
  sortdist=sort(dist)
  epsoptimal=sortdist[which.max(diff(sortdist))]
  
  ##############################
  
  #res <- dbscan(combinedscore_z, eps =ninty5p , minPts = minPts)
  res <- dbscan(combinedscore_z, eps =epsoptimal , minPts = minPts)
  true_label <- c(rep(1,round(0.7*n)),rep(2,round(0.3*n)),rep(3,round(0.04*n)))
  correctp=rand.index(true_label, res$cluster)
  dbet=Sys.time()
  dbt=tp+dbet-dbst
  
  kst=Sys.time()
  
  
  reskmeans = NbClust::NbClust(data =  combinedscore_z, diss = NULL, 
                               distance = "euclidean", min.nc = min.nc, max.nc = max.nc, 
                               method = "kmeans",index="silhouette")
  correctk=rand.index(true_label, reskmeans$Best.partition)
  ket=Sys.time()
  kt=tp+ket-kst
  
  acctable=cbind(correctp,correctk)
  
  #return(list("time"=timecomp,"acctable"=acctable,"clustertable"=clustertable,"clustertablecfda"=clustertablecfda,"clist"=clist,"clistcfda"=clistcfda,"datalist"=datalist))
  timedbk=list("dbt"=dbt,"kt"=kt)
  return(list("time"=timedbk,"acctable"=acctable))
}


acc_table_sim=function(seedlength,k,n,datapoints,q,sparse1,sparse2,sparse3,scorevar1,scorevar2,ps,
                       st,et,spline1D,Mcfda,knnum,pct, pctcfda,minPts,min.nc,max.nc){
  seedvec = matrix((1:(seedlength))+123, ncol = 1)
  results = apply(seedvec, 1, function(x) acc_table_graph(x, k,n,datapoints,q,sparse1,sparse2,sparse3,scorevar1,scorevar2,ps,
                                                          st,et,spline1D,Mcfda,knnum,pct, pctcfda,minPts,min.nc,max.nc))
  resultsexp=results[[1]]
  tcn100t250=sapply(results,"[[",1)
  accn100t250=sapply(results,"[[",2)
  timeave=apply(matrix(unlist(tcn100t250),nrow=2,ncol=seedlength),1,mean)
  timese=apply(matrix(unlist(tcn100t250),nrow=2,ncol=seedlength),1,sd)/sqrt(round(0.7*n)+round(0.3*n)+round(0.04*n))
  accmean=apply(accn100t250,1,mean)
  accse=timese=apply(accn100t250,1,sd)/sqrt(round(0.7*n)+round(0.3*n)+round(0.04*n))
  return(list("rexp"=resultsexp,"time"=timeave,
              "acc"=accmean,"timese"=timese,"accse"=accse))
}


#a=Sys.time()
options(warn=-1)
#n=100
#acc_table_sim=function(seedlength,k,n,datapoints,q,sparse1,sparse2,sparse3,scorevar1,scorevar2,ps,
#st,et,spline1D,Mcfda,knnum,pct, pctcfda,minPts,min.nc,max.nc)
###################################
nst=Sys.time()
n100t250=acc_table_sim(5,3,100,250,3,0,1,3,1,2,1,0.01,0.99,
                       25,3,2,0.9, 0.93,4,2,5)

net=Sys.time()
t1=net-nst

# 
nst=Sys.time()
n100t500=acc_table_sim(5,3,100,500,3,0,1,3,1,2,1,0.01,0.99,
                       25,3,2,0.96, 0.93,4,2,5)
net=Sys.time()
t2= net-nst

# 
# 
nst=Sys.time()
n100t2000=acc_table_sim(5,3,100,2000,3,0,1,3,1,2,1,0.01,0.99,
                        25,3,2,0.9, 0.93,4,2,5)
net=Sys.time()
t3=net-nst

#####################################
#n=500
nst=Sys.time()
n500t250=acc_table_sim(5,3,500,250,3,0,1,3,1,2,1,0.01,0.99,
                       25,3,2,0.9, 0.93,4,2,5)
net=Sys.time()
t4=net-nst

# 
nst=Sys.time()
n500t500=acc_table_sim(5,3,500,500,3,0,1,3,1,2,1,0.01,0.99,25,3,2,0.9, 0.93,4,2,5)
net=Sys.time()
t5=net-nst

# 
nst=Sys.time()
n500t2000=acc_table_sim(5,3,500,2000,3,0,1,3,1,2,1,0.01,0.99,
                        25,3,2,0.95, 0.93,4,2,5)
net=Sys.time()
t6=net-nst

#n=1000
#####################################
nst=Sys.time()
n1000t250=acc_table_sim(5,3,1000,250,3,0,1,3,1,2,1,0.01,0.99,
                        25,3,2,0.9, 0.93,4,2,5)
net=Sys.time()
t7=net-nst

nst=Sys.time()
n1000t500=acc_table_sim(5,3,1000,500,3,0,1,3,1,2,1,0.01,0.99,
                        25,3,2,0.9, 0.93,4,2,5)
net=Sys.time()
t8=net-nst

nst=Sys.time()
n1000t2000=acc_table_sim(5,3,1000,2000,3,0,1,3,1,2,1,0.01,0.99,
                         25,3,2,0.95, 0.93,4,2,5)
net=Sys.time()
t9=net-nst

#sresultn100t250=n100t250sim100[[1]]
# sampleacc=n100t250sim100[1][[1]]$acctable
# sampletime=n100t250sim100[1][[1]]$time
#tcn100t250=sapply(n100t250sim100,"[[",1)
#accn100t250=sapply(n100t250sim100,"[[",2)
#b=Sys.time()
#c=b-a
#c 

save(n100t250,n100t500,n100t2000,
     n500t250,n500t500,n500t2000,
     n1000t250,n1000t500,n1000t2000,t1,t2,t3,t4,t5,t6,t7,t8,t9,
     file="Allscorescfda.RData")
#load("Allscores.RData")

####
#test for true and est eigen values
#generate_data_scenario=function(k,n,datapoints,sparse,scorevar,
#ps,seed=123,st,et)


##

# checkte=function(k,n,datapoints,sparse1,sparse2,sparse3,scorevar1,
#                  scorevar2,ps,seed,st,et){
#   clustern50t250=generate_data_scenario(k=k,n=round(0.7*n),datapoints,sparse1,scorevar1,ps,seed=seed,st,et)
#   clustern50t250p2=generate_data_scenario(k=k,n=round(0.3*n),datapoints,sparse2,scorevar2,ps,seed=seed,st,et)
#   clustern50t250p3=generate_data_scenario(k=k,n=round(n*0.04),datapoints,sparse3,scorevar1,ps,seed=seed,st,et)
#   
#   
#   truen100t250=mapply(rbind,clustern50t250$Trueth,clustern50t250p2$Trueth,clustern50t250p3$Trueth,SIMPLIFY=FALSE,USE.NAMES = TRUE)
#   estn100t250=mapply(rbind,clustern50t250$Est,clustern50t250p2$Est,clustern50t250p3$Est,SIMPLIFY=FALSE,USE.NAMES = TRUE)
#   combn100t250=list("Trueth"=truen100t250, "Est"=estn100t250)
#   
#   sp=Sys.time()
#   #q=3
#   #spline1D=25
#   #Mcfda=3
#   combn100t250mfpca=baysglm_result(datagenerated=combn100t250,datapoints,q=q,n=n,st,et,k,spline1D,Mcfda)
#   ept=Sys.time()
#   tp=ept-sp
#   
#   
#   Z_i1hatstar=combn100t250$Trueth$TrueZ_i1
#   Z_i2hatstar=combn100t250$Trueth$TrueZ_i2
#   
#   Zihatstar=array(c(Z_i1hatstar,Z_i2hatstar),c(round(0.7*n)+round(0.3*n)+round(0.04*n),datapoints,q-1))
#   
#   #####get z scores
#   fdarange = c(st, et)
#   fdabasis = fda::create.bspline.basis(fdarange, splines1D, 4)
#   fdatime = seq(st, et, length = datapoints)
#   Zihatstarnew = apply(Zihatstar, c(1, 3), t)
#   fdafd = fda::smooth.basis(fdatime, Zihatstarnew, fdabasis)$fd
#   nharm = M
#   fdapcaList = fda::pca.fd(fdafd, nharm)
#   
#   finalM=which(cumsum(fdapcaList$values)/sum(fdapcaList$values)>=0.98)[1]
#   fdapcaListfinal = fda::pca.fd(fdafd, finalM)
#   scores_z = apply(fdapcaListfinal$scores, c(1, 2), sum)
#   
#   
#   
#   return(list("sdim"=dim(combn100t250mfpca$Bayesian_scores),
#               "cums"=cumsum(combn100t250mfpca$values)/sum(combn100t250mfpca$values),
#               "tsdim"=dim(scores_z),
#               "tcums"=cumsum(fdapcaListfinal$values)/sum(fdapcaListfinal$values)
#               
#               ))
# }

#checkte=function(k,n,datapoints,sparse1,sparse2,sparse3,scorevar1,
#scorevar2,ps,seed,st,et)
# k=3
# #n=100
# #datapoints=250
# sparse1=0
# sparse2=1
# sparse3=3
# scorevar1=1
# scorevar2=2
# ps=1
# seed=123
# st=0.01
# et=0.99
# splines1D=25
# M=2
# checkte(k,100,250,sparse1,sparse2,sparse3,scorevar1,
#                  scorevar2,ps,seed,st,et)
# checkte(k,100,500,sparse1,sparse2,sparse3,scorevar1,
#         scorevar2,ps,seed,st,et)
# checkte(k,100,2000,sparse1,sparse2,sparse3,scorevar1,
#         scorevar2,ps,seed,st,et)