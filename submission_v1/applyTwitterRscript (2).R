#apply to Twitter

load("X_i1.RData")
load("X_i2.RData")
load("X_i3.RData")

# X_i1=X_i1[1:50,1:1000]
# X_i2=X_i2[1:50,1:1000]
# X_i3=X_i3[1:50,1:1000]
#############################

#01_functions
###############################

# 
# library(xtable)
library(fda)
library(refund)
library(Matrix)
library(MASS)
library(arm)
library(mgcv)
# library(readr)
library(funData)
library(MFPCA)
# library(purrr)
library(tidyverse)
library(gridExtra)
# library(knitr)
# library(kableExtra)
 library(dbscan)
# library(cfda)
library(tidyr)
library(dplyr)

############################################################################
#function to to generate and simulate data
##########################################################################




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

#adjust to ML 4/6/2022
# regression_g = function(z, tt, k=25, method="ML"){
#   gam1 <- gam(z~s(tt, bs = "cr", m=2, k = k), 
#               family="binomial", method = method,control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),optimizer=c("outer","bfgs"))
#   return(gam1$fitted.values)
# }


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
  N_train=dim(as.matrix(Curves_train))[1]
  vec = matrix(1:(N_train), ncol = 1)
  smoothed_x = logit(t(apply(vec, 1, function(x) regression_g(x, as.matrix(Curves_train), tt))))
  smoothed_x
}

p_ihat=function(Z_i1app,Z_i2app){
  denom=(1+exp(Z_i1app)+exp(Z_i2app))
  p_i1h=exp(Z_i1app)/denom
  p_i2h=exp(Z_i2app)/denom
  p_i3h=1/denom
  # p_i3h=1-p_i1h- p_i2h
  return(list("p_1hatmatrix"=p_i1h,"p_2hatmatrix"=p_i2h,"p_3hatmatrix"=p_i3h))
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
# datapoints=dim(X_i1)[2]
# n=dim(X_i1)[1]

datapoints=dim(X_i1)[2]
n=dim(X_i1)[1]
# st=as.numeric(colnames(X_i1))[1]
st=as.numeric(colnames(X_i1))[1]
# et=tail(as.numeric(colnames(X_i1)),n=1)
et=tail(as.numeric(colnames(X_i1)),n=1)
t=as.numeric(colnames(X_i1))
get_estimate_latent=function(X_i1,X_i2,X_i3,ps,k,t){
  
  k=k
  
  #k=3  #number of eigen functions
  q=3  #level of the categorical level
  # n=10 #number of subjects
  
  
  #collect value and graph
  #collect first two rows of observed binary curves
  # X_i1=t(X_i[1,,])  #all n row subjects , t columns values related to p1
  # X_i2=t(X_i[2,,]) #all n row subjects , t columns values related to p2
  # X_i3=t(X_i[3,,]) #all n row subjects , t columns values related to p3
  
  #recover Z_i1 hat using Z_i[1,all j, all n] only related to p1
  Z_i1hat=Z_ihat(X_i1,t)
  #recover Z_i2 hat using X_i[2,all j, all n] only related to p2 
  Z_i2hat=Z_ihat(X_i2,t)
  Z_i3hat=Z_ihat(X_i3,t)
  
  
  # Z_i1hatstar=log(exp(Z_i1hat))+log(1+exp(Z_i3hat))-log(1+exp(Z_i1hat))-log(exp(Z_i3hat))
  # Z_i2hatstar=log(exp(Z_i2hat))+log(1+exp(Z_i3hat))-log(1+exp(Z_i2hat))-log(exp(Z_i3hat))
  
  
  # common_check=1-exp(Z_i1hat+Z_i2hat)-exp(Z_i2hat)
  # common_check[common_check<=0]=0.001
  # 
  # 
  # 
  # Z_i1hatstar=Z_i1hat-log(common_check)
  # Z_i2hatstar=Z_i2hat+log(1+exp(Z_i1hat))-log(common_check)
  if(ps==1){ 
    # smooth_p1=(1+exp(-Z_i1hat))^(-1)
    # smooth_p2=(1+exp(-Z_i2hat))^(-1)
    # smooth_p3=(1+exp(-Z_i3hat))^(-1)
    # Z_i1hatstar=log(smooth_p1)-log(smooth_p3)
    # Z_i2hatstar=log(smooth_p2)-log(smooth_p3)
    
    
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
  
  
  # 
  # #plot and compare
  # par(mfrow=c(1,3))
  # matplot(t(p_i1[1:ns,]),col="green")  #true Z_i1 for first subject
  # matlines(t(p_i1hat[1:ns,]),col="red")
  # 
  # matplot(t(p_i2[1:ns,]),col="green")  #true Z_i1 for first subject
  # matlines(t(p_i2hat[1:ns,]),col="red")
  # 
  # matplot(t(p_i3[1:ns,]),col="green")  #true Z_i1 for first subject
  # matlines(t(p_i3hat[1:ns,]),col="red")
  # truel=list("TrueX1"=X_i1,"TrueX2"=X_i2,"TrueX3"=X_i3,"TrueZ_i1"=Z_i1,"TrueZ_i2"=Z_i2,"Truep_i1"=p_i1,"Truep_i2"=p_i2,"Truep_i3"=p_i3,"Truescore_matrix"=score_matrix,"Truecatcurve"=X_nt)
  
  
  ####################
  #need to add it
  ####################
  
  # truel=list("Truecatcurve"=X_nt)
  est=list("EstimateZ_i1"=Z_i1hatstar,"EstimateZ_i2"=Z_i2hatstar,"Estimatep_i1"=p_i1hat,"Estimatep_i2"=p_i2hat,"Estimatep_i3"=p_i3hat)
  # return(list("Trueth"=truel,"Est"=est))
  return(list("Est"=est))
}


#MFPCA and get scores
#must need function can't use baysglm vector size can't be allocated 6.7G
bayesian_scores=function(muhat,eigenvalues,eigenfunctions,Curves_train, MFPCAdim,datapoints,n,st,et,t){
  # t=seq(from = st,to = et, length=datapoints)
  t=t
  st=st
  et=et
  ##
  #STEP 3:
  # Set up and apply Bayesglm framework
  ##
  fit = list(mu = muhat,
             evalues = eigenvalues,
             efunctions = eigenfunctions)
  
  mu_t_hat = fit$mu
  eigen_vals1 = fit$evalues
  eigen_funcs1 = fit$efunctions
  
  #data frame used in bayesglm
  tt=t
  N_train=n  #number of subjects
  D=datapoints  #number of time points
  dta = data.frame(index = rep(tt, N_train*MFPCAdim),
                   value = c(t(Curves_train)),
                   id = rep(1:N_train, each = D*MFPCAdim))
  
  npc = length(eigen_vals1)
  if(npc>1){
    for (z in 1:npc) {
      dta <- cbind(dta, rep(eigen_funcs1[,z], N_train*MFPCAdim))
    }
  }else{
    dta = cbind(dta, matrix(eigen_funcs1, ncol =1))
  }
  
  #assign names to data frame
  names(dta)[4:(4 + npc - 1)] <- c(paste0("ksi", 1:npc))
  #repeat mean function in data frame once per user
  dta$mu = rep(mu_t_hat , N_train)
  
  #get formula for glm
  glm_structure = paste(paste0("ksi", 1:npc), collapse = "+")
  glm_structure = paste("value ~ -1 + offset(mu) +" , glm_structure , sep="")
  #set scale for the glm
  prior_scales_test = eigen_vals1
  
  #Estimate the Scores fodtar the training set
  vec = matrix(1:N_train, ncol = 1)
  scores_train = t(apply(vec, 1, function(x) regression_bf2(x, dta, glm_structure, prior_scales_test)))
  return(scores_train)
  
}


# MFPCA_result
baysglm_result=function(datagenerated,M,datapoints,q=3,n,st,et,t){
  datapoints=datapoints
  # t=seq(from = st,to = et, length=datapoints)
  t=t
  st=st
  et=et
  
  # Z_i1=datagenerated$Trueth$TrueZ_i1
  # Z_i2=datagenerated$Trueth$TrueZ_i2
  Z_i1hatstar=datagenerated$Est$EstimateZ_i1
  Z_i2hatstar=datagenerated$Est$EstimateZ_i2
  
  
  # p_i1=datagenerated$Trueth$Truep_i1
  # p_i2=datagenerated$Trueth$Truep_i2
  # p_i3=datagenerated$Trueth$Truep_i3
  p_i1hat=datagenerated$Est$Estimatep_i1
  p_i2hat=datagenerated$Est$Estimatep_i2
  p_i3hat=datagenerated$Est$Estimatep_i3
  
  
  
  # categorical_curve=truevalues$Truecatcurve  #dim:n=200 subjects, col=t 500 points
  #step 2: MFPCA
  # Z_ihat muldata
  m1=funData(argvals = list(t), X = Z_i1hatstar) #20*500
  m2=funData(argvals = list(t), X = Z_i2hatstar) 
  
  # p_ihat muldata
  m1p=funData(argvals = list(t), X = p_i1hat) #20*500
  m2p=funData(argvals = list(t), X = p_i2hat)
  m3p=funData(argvals = list(t), X = p_i3hat)
  # 
  mvdata=multiFunData(list(m1,m2))
  mvdatap=multiFunData(list(m1p,m2p,m3p))
  M=M
  
  # MFPCA based on univariate FPCA Z_ihat
  # uFPCA <- MFPCA(mvdata, M = M, uniExpansions = list(list(type = "uFPCA"),
  #                                                    list(type = "uFPCA")))
  
  # uFPCA <- MFPCA(mvdata, M = M, uniExpansions = list(list(type = "splines1D", k = 10),
  #                                                           list(type = "splines1D", k = 10)),
  #                  fit = TRUE)
  
  uFPCA <- MFPCA(mvdata, M = M, uniExpansions = list(list(type = "splines1D", k = 10),
                                                     list(type = "splines1D", k = 10))
                 )
  # MFPCA based on univariate FPCA p_ihat
  uFPCAp <- MFPCA(mvdatap, M = M, uniExpansions = list(list(type = "splines1D", k = 10),
                                                       list(type = "splines1D", k = 10),list(type = "splines1D", k = 10)))
  
  # 
   MFPCA_results_Z=list("meanfnZ"=uFPCA$meanFunction,"eigenfnZ"=uFPCA$functions,"eigenvalueZ"=uFPCA$values,"scoresZ"=uFPCA$scores)
  
  
  
  # MFPCA_results_p=list("meanfnp"=uFPCAp$meanFunction,"eigenfnp"=uFPCAp$functions,"eigenvaluep"=uFPCAp$values,"scoresp"=uFPCAp$scores)
  # no need for scores
  #MFPCA_results_Z=list("meanfnZ"=uFPCA$meanFunction,"eigenvalueZ"=uFPCA$values)
  # 
  MFPCA_results_p=list("meanfnp"=uFPCAp$meanFunction,"eigenvaluep"=uFPCAp$values)
  # 
  
  
  
  
  
  ##use bayesian to find scores can't allocate vector size 6.7Gb
  #bayesian scores for pi curves
  ###################################################################################################################
  #first need to combine multivariate data to one big matrix
  #MU1+mu2
  muhatp=c(c(uFPCAp$meanFunction[[1]]@X),c(uFPCAp$meanFunction[[2]]@X),c(uFPCAp$meanFunction[[3]]@X))
  # #eigenvalues=uFPCA$values/500
  eigenvaluesp=uFPCAp$values
  # #eigenfunctions=(rbind(t(uFPCA$functions[[1]]@X),t(uFPCA$functions[[2]]@X),t(uFPCA$functions[[3]]@X)))*sqrt(500)
  eigenfunctionsp=(rbind(t(uFPCAp$functions[[1]]@X),t(uFPCAp$functions[[2]]@X),t(uFPCAp$functions[[3]]@X)))
  # #smoothed probability curves instead of binary curves
  Curves_trainp=cbind(p_i1hat,p_i2hat,p_i3hat)
  q=q
  MFPCAdimp=q  #how many probabilities or how many categories
  
  # 
  # #bayesian scores for Zi curves
  # #first need to combine multivariate data to one big matrix
  # #MU1+mu2
  muhatz=c(c(uFPCA$meanFunction[[1]]@X),c(uFPCA$meanFunction[[2]]@X))
  #eigenvalues=uFPCA$values/500
  eigenvaluesz=uFPCA$values
  #eigenfunctions=(rbind(t(uFPCA$functions[[1]]@X),t(uFPCA$functions[[2]]@X),t(uFPCA$functions[[3]]@X)))*sqrt(500)
  eigenfunctionsz=(rbind(t(uFPCA$functions[[1]]@X),t(uFPCA$functions[[2]]@X)))
  #smoothed Z curves instead of binary curves
  # Curves_trainz=cbind(log(smooth_p1)-log(smooth_p3),log(smooth_p2)-log(smooth_p3))
  
  Curves_trainz=cbind( Z_i1hatstar, Z_i2hatstar)
  MFPCAdimz=q-1  #how many latent Z_i curves, or category-1
  n=n
  st=st
  et=et
  t=t
  # scores_p=bayesian_scores(muhatp,eigenvaluesp,eigenfunctionsp,Curves_trainp, MFPCAdimp)
  #scores_z=bayesian_scores(muhatz,eigenvaluesz,eigenfunctionsz,Curves_trainz, MFPCAdimz,datapoints=datapoints,n=n,st=st,et=et,t)
  # bayesianscores=list("baypscores"=scores_p,"bayzscores"=scores_z)
  #scores_p=bayesian_scores(muhatp,eigenvaluesp,eigenfunctionsp,Curves_trainp, MFPCAdimp,datapoints=datapoints,n=n,st=st,et=et,t)
  # # 
  #pcaZ_i1hat=apply(t(uFPCA$functions[[1]]@X)%*%t(scores_z),2,function(x) c(uFPCA$meanFunction[[1]]@X)+x)  #500*20
  #pcaZ_i2hat=apply(t(uFPCA$functions[[2]]@X)%*%t(scores_z),2,function(x) c(uFPCA$meanFunction[[2]]@X)+x)
  # 
  # pcap_i1hatscore=apply(t(uFPCAp$functions[[1]]@X)%*%t(scores_p),2,function(x) c(uFPCAp$meanFunction[[1]]@X)+x)  #500*20
  # pcap_i2hatscore=apply(t(uFPCAp$functions[[2]]@X)%*%t(scores_p),2,function(x) c(uFPCAp$meanFunction[[2]]@X)+x)
  # pcap_i3hatscore=apply(t(uFPCAp$functions[[3]]@X)%*%t(scores_p),2,function(x) c(uFPCAp$meanFunction[[3]]@X)+x)
  # # 
  # #Bayesian estimate
  #BayEstimate=list("BaypcaZ_i1"=pcaZ_i1hat,"BaypcaZ_i2"=pcaZ_i2hat)
  # 
  
  # return(list("Truth"=truevalues,"Smoothed"=estvalues,"MFPCAZ"=MFPCA_results_Z,"MFPCAp"=MFPCA_results_p,"Bayesian_scores"=bayesianscores,"MFPCAestimate"=MFPCAest,"BayesianEst"=BayEstimate))
  #if no need for MFPCA scores and estimate
  # return(list("MFPCAZ"=MFPCA_results_Z,"MFPCAp"=MFPCA_results_p,"Bayesian_scores"=bayesianscores,"BayesianEst"=BayEstimate))
  #can't allocate vector size 6.7G whe using baysglm
  
   return(list("MFPCAZ"=MFPCA_results_Z,"ufpca"=uFPCA))
  #return(list("MFPCAZ"=MFPCA_results_Z,"MFPCAp"=MFPCA_results_p,"Bayesian_scores"=scores_z,
              #"Bayesian_scoresp"=scores_p,"BayesianEst"=BayEstimate,"eigenvalues"=eigenvaluesz))

}


##ramsay
baysglm_resultram=function(datagenerated,M,datapoints,q=3,n,st,et,t){
  datapoints=datapoints
  # t=seq(from = st,to = et, length=datapoints)
  t=t
  st=st
  et=et
  Z_i1hatstar=datagenerated$Est$EstimateZ_i1
  Z_i2hatstar=datagenerated$Est$EstimateZ_i2
  n=dim(Z_i1hatstar)[1]
  Zihatstar=array(c(Z_i1hatstar, Z_i2hatstar), c(n,datapoints,2))

  fdarange = c(st, et)
  fdabasis = fda::create.bspline.basis(fdarange, 25, 4)
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
  
  return(list("scores"=scores_z,
              "latentcurve"=Zihatstar,"eigenvalue"=values,"meanfn"=meanfn,
              "eigenfn"=eigenfn))
}




startt=Sys.time()
#step1
#get_estimate_latent=function(X_i1,X_i2,X_i3,ps,k,t)
trueest=get_estimate_latent(as.matrix(X_i1),as.matrix(X_i2),as.matrix(X_i3),1,3,t)
#MFPCA
#baysglm_result=function(datagenerated,M,datapoints,q=3,n,st,et,t)
#MFPCAr=baysglm_result(trueest,3,datapoints,q=3,n,st,et,t)
library(MFPCA)
#MFPCArcheck=baysglm_result(trueest,3,2000,q=3,1000,0.1,0.9,seq(0.1,0.9,length=2000))
#MFPCAr4=baysglm_result(trueest,4,datapoints,q=3,n,st,et,t)
twittersim=baysglm_result(trueest,3,2000,q=3,1000,0.1,0.9,seq(0.1,0.9,length=2000))
plot(twittersim$MFPCAZ$scoresZ[,1],twittersim$MFPCAZ$scoresZ[,2])
combinedscore_z=twittersim$MFPCAZ$scoresZ[,1:2]
library(dbscan)
library(ggplot2)
dist=kNNdist(combinedscore_z, k = 3)
distdata=data.frame(sort(dist))
distdata$index=1:dim(distdata)[1]
#pct=0.98.  #happ
pct=0.98
ninty5p=quantile(dist, probs = pct)
ninty5p
#########################5/29/2023
load("Twitterram.RData")
combinedscore_z=mfpcaram$scores
dimz=dim(combinedscore_z)[2]
if (dimz<=2){minPts=4}
if (dimz>2){minPts=2*dimz+1}
library(dbscan)
library(fossil)
devtools::install_github("ahasverus/elbow")
library(elbow)
library(pdfCluster)
dist=kNNdist(combinedscore_z, k = minPts-1)
#ninty5p=quantile(dist, probs = pct)
scaleep=3
#########change to max increase
distdataelbow=data.frame(sort(dist))
distdataelbow$index=1:(dim(combinedscore_z)[1])
ipoint <- elbow(data = distdataelbow)
epsoptimal=(ipoint$sort.dist._selected)*scaleep
#epsoptimal=sortdist[which.max(diff(sortdist))]

##############################

#res <- dbscan(combinedscore_z, eps =ninty5p , minPts = minPts)
res <- dbscan(combinedscore_z, eps =epsoptimal , minPts = 5)

#res <- dbscan(combinedscore_z, eps =ninty5p, minPts = 5)
table(res$cluster)
tclusterdata=data.frame(combinedscore_z)
tclusterdata$Cluster=as.factor(res$cluster)
library(ggplot2)
tps <- ggplot(tclusterdata,aes(X1,X2,colour = Cluster)) + geom_point(aes(shape=Cluster),size=3)+ggtitle(paste0("Twitter Users Cluster Results",'\n',"(",dim(tclusterdata)[1]," Subjects",")")) +
  xlab(expression('Score '* widehat(xi[i1]))) + ylab(expression('Score '* widehat(xi[i2])))+ theme(plot.title = element_text(hjust = 0.5))+#theme(text = element_text(size = 20))+
  #theme(legend.text = element_text(size = 30))  
  #theme(axis.text=element_text(size = 20),axis.title = element_text(size =30))
  theme(text=element_text(size = 20))
tps
vec0=c(which(tclusterdata$Cluster==0))
vec1=c(which(tclusterdata$Cluster==1))
vec2=c(which(tclusterdata$Cluster==2))

knn=3
minPts=5
ninty5p=2.389757 
save(twittersim,tclusterdata,vec0,vec1,vec2,knn,minPts,ninty5p,file="Twittersim.RData")



save(trueest,MFPCAr,MFPCAr4,file="Twitter.RData")
# save(trueest,MFPCAr,file="Twittert.RData")
endt=Sys.time()
time_elapse=endt-startt
print(time_elapse)


#######################
#ramsay
startt=Sys.time()
load("Twitter.RData")
#MFPCArcheck=baysglm_result(trueest,3,2000,q=3,1000,0.1,0.9,seq(0.1,0.9,length=2000))
mfpcaram=baysglm_resultram(trueest,2,2000,q=3,1000,0.1,0.9,seq(0.1,0.9,length=2000))
save(trueest,mfpcaram,file="Twitterram.RData")
##cumsum(mfpcaram$eigenvalue)/sum(mfpcaram$eigenvalue)
#[1] 0.9668647 0.9832566 0.9940806
# save(trueest,MFPCAr,file="Twittert.RData")
endt=Sys.time()
time_elapse=endt-startt
print(time_elapse)
###########################################
load("Twitterram.RData")
#graph k-means results
reskmeanstwitters = NbClust::NbClust(data =  mfpcaram$scores, diss = NULL, 
                             distance = "euclidean", min.nc = 2, max.nc = 5, 
                             method = "kmeans",index="silhouette")
reskmeanstwittergap = NbClust::NbClust(data =  mfpcaram$scores, diss = NULL, 
                                    distance = "euclidean", min.nc = 2, max.nc = 5, 
                                    method = "kmeans",index="gap")
reskmeanstwitterma = NbClust::NbClust(data =  mfpcaram$scores, diss = NULL, 
                                       distance = "euclidean", min.nc = 2, max.nc = 5, 
                                       method = "kmeans",index="marriot")
library(factoextra)
set.seed(123)
#km.res <- kmeans(mfpcaram$scores, 4, nstart = 25)
fviz_nbclust(mfpcaram$scores,method = c("silhouette"),kmeans)
fviz_nbclust(mfpcaram$scores,method = c( "gap_stat"),kmeans)
fviz_nbclust(mfpcaram$scores,method = c( "wss"),kmeans)


tclusterdatas=data.frame(mfpcaram$scores)
tclusterdatas$Cluster=as.factor(reskmeanstwitters$Best.partition)

tclusterdatagap=data.frame(mfpcaram$scores)
tclusterdatagap$Cluster=as.factor(reskmeanstwittergap$Best.partition)

tclusterdatama=data.frame(mfpcaram$scores)
tclusterdatama$Cluster=as.factor(reskmeanstwitterma$Best.partition)


tpsks <- ggplot(tclusterdatas,aes(X1,X2,colour = Cluster)) + geom_point(aes(shape=Cluster),size=3)+ggtitle(paste0("Twitter Users Cluster Results",'\n',"(",dim(tclusterdatas)[1]," Subjects",")")) +
  xlab(expression('Score '* widehat(xi[i1]))) + ylab(expression('Score '* widehat(xi[i2])))+ theme(plot.title = element_text(hjust = 0.5))+#theme(text = element_text(size = 20))+
  #theme(legend.text = element_text(size = 30))  
  #theme(axis.text=element_text(size = 20),axis.title = element_text(size =30))
  theme(text=element_text(size = 20))
tpsks

tpsgap <- ggplot(tclusterdatagap,aes(X1,X2,colour = Cluster)) + geom_point(aes(shape=Cluster),size=3)+ggtitle(paste0("Twitter Users Cluster Results",'\n',"(",dim(tclusterdatagap)[1]," Subjects",")")) +
  xlab(expression('Score '* widehat(xi[i1]))) + ylab(expression('Score '* widehat(xi[i2])))+ theme(plot.title = element_text(hjust = 0.5))+#theme(text = element_text(size = 20))+
  #theme(legend.text = element_text(size = 30))  
  #theme(axis.text=element_text(size = 20),axis.title = element_text(size =30))
  theme(text=element_text(size = 20))
tpsgap

tpsma <- ggplot(tclusterdatama,aes(X1,X2,colour = Cluster)) + geom_point(aes(shape=Cluster),size=3)+ggtitle(paste0("Twitter Users Cluster Results",'\n',"(",dim(tclusterdatama)[1]," Subjects",")")) +
  xlab(expression('Score '* widehat(xi[i1]))) + ylab(expression('Score '* widehat(xi[i2])))+ theme(plot.title = element_text(hjust = 0.5))+#theme(text = element_text(size = 20))+
  #theme(legend.text = element_text(size = 30))  
  #theme(axis.text=element_text(size = 20),axis.title = element_text(size =30))
  theme(text=element_text(size = 20))
tpsma
############################################
#graph dbscan results



##############################
load("Twitter.RData")
#graph
graph_fn_estpca(trueest,MFPCAr,2000,1000,-5,15,-5,10,-0.01,1.1,-0.01,1.1,-0.1,0.5,-10,10,-5,5,st,et,t)
#################################################################################################################################
#need to get rid of the truth
#################################################################################################################################
graph_fn_estpca=function(datagenerated,baysglmresult,datapoints,n,z1low,z1up,z2low,z2up,p1low,p1up,p2low,p2up,p3low,p3up,xslow,xsup,yslow,ysup,st,et,t){
  datapoints=datapoints
  t=t
  st=st
  et=et
  n=n
  ########################################################################################

  #plot 1: truth and smoothed
  ########################################################################################
  #entire block is  n*t dimension
  #truth Z, P and smoothed Z, P
  # Z_i1=dataset$AverageTrue$TruemeanZ1
  # Z_i2=dataset$AverageTrue$TruemeanZ2
  # Z_i1_est=dataset$AverageSmooth$sZ1mean
  # Z_i2_est=dataset$AverageSmooth$sZ2mean
  # p_i1=dataset$AverageTrue$Truemeanp1
  # p_i2=dataset$AverageTrue$Truemeanp2
  # p_i3=dataset$AverageTrue$Truemeanp3
  # p_i1_est=dataset$AverageSmooth$sp1mean
  # p_i2_est=dataset$AverageSmooth$sp2mean
  # p_i3_est=dataset$AverageSmooth$sp3mean


  # Z_i1=datagenerated$Trueth$TrueZ_i1
  # Z_i2=datagenerated$Trueth$TrueZ_i2
  Z_i1_est=datagenerated$Est$EstimateZ_i1
  Z_i2_est=datagenerated$Est$EstimateZ_i2


  # p_i1=datagenerated$Trueth$Truep_i1
  # p_i2=datagenerated$Trueth$Truep_i2
  # p_i3=datagenerated$Trueth$Truep_i3
  p_i1_est=datagenerated$Est$Estimatep_i1
  p_i2_est=datagenerated$Est$Estimatep_i2
  p_i3_est=datagenerated$Est$Estimatep_i3



  # fmeanzi1hat=dataset$AverageMFPCA$pcaz1mean
  # fmeanzi2hat=dataset$AverageMFPCA$pcaz2mean

  fmeanzi1hat=baysglmresult$MFPCAZ$meanfnZ[[1]]@X
  fmeanzi2hat=baysglmresult$MFPCAZ$meanfnZ[[2]]@X

  # fmeanpi1hat=dataset$AverageMFPCA$pcap1mean
  # fmeanpi2hat=dataset$AverageMFPCA$pcap2mean
  # fmeanpi3hat=dataset$AverageMFPCA$pcap3mean


  #bayesian mfpca z and bayesian mfpca p are all t*n
  pcaZ_i1hat=baysglmresult$BayesianEst$BaypcaZ_i1
  pcaZ_i2hat=baysglmresult$BayesianEst$BaypcaZ_i2
  # pcap_i1hat=dataset$AverageBMFPCA$bfpcap1
  # pcap_i2hat=dataset$AverageBMFPCA$bfpcap2
  # pcap_i3hat=dataset$AverageBMFPCA$bfpcap3

  #scores
  # score_matrix=datagenerated$Trueth$Truescore_matrix
  bmfpca_zscores=baysglmresult$Bayesian_scores



  # Z1 2 10 Z2 -1   8  p10.6 1    0  0.35       0   0.04
  #truth and smoothed

  par(mfrow=c(3,2))

  #Z_i
  # meanzi1=apply(t(Z_i1),1,mean)
  # matplot( t(Z_i1),
  #          type='l', lty=1, col="light grey",
  #          main = mtext(bquote("True Latent Cruve "*'Z'[i1])),
  #          xlab="Number of Datapoints", ylab=expression('Z'[i1]),ylim=c(z1low,z1up))  ##sparse z1 y (3,18)
  # lines(meanzi1 ,
  #       type='l', lty=1, lwd=2, col = "red")
  # 
  # 
  # 
  # 
  # meanzi2=apply(t(Z_i2),1,mean)
  # matplot(t(Z_i2),
  #         type='l', lty=1, col="light grey",
  #         main = mtext(bquote("True Latent Cruve "*'Z'[i2])),
  #         xlab="Number of Datapoints", ylab=expression('Z'[i2]),ylim=c(z2low,z2up))   ##sparse z2 y (1,14)
  # lines(meanzi2,
  #       type='l', lty=1, lwd=2, col = "red")
  # 
  # legend(x = "topleft",  horiz = TRUE,        # Position
  #        legend = c("Individual", "Mean"),  # Legend texts
  #        #lty = c(1, 2),           # Line types
  #        col = c("grey", "red"),           # Line colors
  #        lwd = 2)                 # Line width


  meanzi1hat=apply(t(Z_i1_est),1,mean)
  matplot( t(Z_i1_est),
           type='l', lty=1, col="light grey",
           main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i1])),
           xlab="Number of Datapoints", ylab=expression(widehat('Z')[i1]),ylim=c(z1low,z1up))
  lines(meanzi1hat ,
        type='l', lty=1, lwd=2, col = "red")



  meanzi2hat=apply(t(Z_i2_est),1,mean)
  matplot( t(Z_i2_est),
           type='l', lty=1, col="light grey",
           main = mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i2])),
           xlab="Number of Datapoints", ylab=expression(widehat('Z')[i2]),ylim=c(z2low,z2up))
  lines(meanzi2hat ,
        type='l', lty=1, lwd=2, col = "red")

  legend(x = "topleft",  horiz = TRUE,        # Position
         legend = c("Individual", "Mean"),  # Legend texts
         #lty = c(1, 2),           # Line types
         col = c("grey", "red"),           # Line colors
         lwd = 2)
  ##########################################################
  #March 14
  #plot truth mean and esitmated mean on the same scale
  #sparse mean need lengend to be on bottomright for Z: ylim=c(-1.4,2) xlim=c(-4,4)
  #non sparse mean need lengend to be on topleft for Z: ylim=c(-0.6,0.6) xlim=c(-1.5,1.5)
  #########################################################
  #Z1
  # matplot( t(Z_i1),
  #          type='l', lty=1, col="light blue",
  #          main = mtext(bquote("True and Estimated "*'Z'[i1])),
  #          xlab="Number of Datapoints", ylab=expression('Z'[i1]),ylim=c(z1low,z1up))  ##sparse z1 y (3,18)

  matplot(t(Z_i1_est),
           type='l', lty=1, col="light coral",
           # main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i1])),
           xlab="Number of Datapoints", ylab=expression(widehat('Z')[i1]))

  lines(meanzi1hat,type='l', lwd=6,lty=1, col="dark blue",
        # main=mtext(bquote(" Mean Latent Cruves "*'Z'[i1])),
        xlab="Number of Datapoints", ylab=expression("Average"*'Z'[i1]))
  # lines(meanzi1 ,
  #       type='l', lty=1, lwd=6, col = "dark red")

  ##########################################################
  #Z2
  # matplot(t(Z_i2),
  #         type='l', lty=1, col="light blue",
  #         main = mtext(bquote("True and Estimated "*'Z'[i2])),
  #         xlab="Number of Datapoints", ylab=expression('Z'[i2]),ylim=c(z2low,z2up))  ##sparse z2 y (1,14)

  matplot( t(Z_i2_est),
            type='l', lty=1, col="light coral",
            # main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i2])),
            xlab="Number of Datapoints", ylab=expression(widehat('Z')[i2]))

  lines(meanzi2hat,type='l', lty=1,lwd=6, col="dark blue",
        # main=mtext(bquote(" Mean Latent Cruves "*'Z'[i2])),
        xlab="Number of Datapoints", ylab=expression("Average"*'Z'[i2]))
  # lines(meanzi2,
  #       type='l', lty=1, lwd=6, col = "dark red")
  # # legend(x = "topleft",          # Position
  # #        legend = c("True Mean", "Estimated Mean"),  # Legend texts
  # #        #lty = c(1, 2),           # Line types
  # #        col = c("grey", "red"),           # Line colors
  # #        lwd = 2)                 # Line width
  # 
  # legend("topleft",  horiz = TRUE,        # Position
  #        legend = c("True Mean", "Estimated Mean"),  # Legend texts
  #        #lty = c(1, 2),           # Line types
  #        col = c("dark blue", "dark red"),           # Line colors
  #        lwd = 2)                 # Line width


  ####################################################################################
  par(mfrow=c(3,3))
  #p_i

  # pi1=apply(t(p_i1),1,mean)
  # matplot( t(p_i1),
  #          type='l', lty=1, col="light grey",
  #          main = mtext(bquote("True Latent Probability Cruve "*'p'[i1])),
  #          xlab="Number of Datapoints", ylab=expression('p'[i1]),ylim = c(p1low,p1up))     #sparse p1 y 0.75-0.95
  # lines(pi1 ,
  #       type='l', lty=1, lwd=2, col = "red")
  # 
  # 
  # pi2=apply(t(p_i2),1,mean)
  # matplot( t(p_i2),
  #          type='l', lty=1, col="light grey",
  #          main = mtext(bquote("True Latent Probability Cruve "*'p'[i2])),
  #          xlab="Number of Datapoints", ylab=expression('p'[i2]),ylim=c(p2low,p2up))    #sparse p1 y 0045-0.2
  # lines(pi2 ,
  #       type='l', lty=1, lwd=2, col = "red")
  # 
  # pi3=apply(t(p_i3),1,mean)
  # matplot( t(p_i3),
  #          type='l', lty=1, col="light grey",
  #          main = mtext(bquote("True Latent Probability Cruve "*'p'[i3])),
  #          xlab="Number of Datapoints", ylab=expression('p'[i3]),ylim=c(p3low,p3up))
  # lines(pi3 ,
  #       type='l', lty=1, lwd=2, col = "red")


  pi1hat=apply(t(p_i1_est),1,mean)
  matplot( t(p_i1_est),
           type='l', lty=1, col="light grey",
           main=mtext(bquote("Estimated Latent Cruve "*widehat('p')[i1])),
           xlab="Number of Datapoints", ylab=expression(widehat('p')[i1]),ylim = c(p1low,p1up))
  lines(pi1hat ,
        type='l', lty=1, lwd=2, col = "red")


  pi2hat=apply(t(p_i2_est),1,mean)
  matplot(t(p_i2_est),
          type='l', lty=1, col="light grey",
          main=mtext(bquote("Estimated Latent Cruve "*widehat('p')[i2])),
          xlab="Number of Datapoints", ylab=expression(widehat('p')[i2]),ylim=c(p2low,p2up))
  lines(pi2hat ,
        type='l', lty=1, lwd=2, col = "red")

  pi3hat=apply(t(p_i3_est),1,mean)
  matplot(t(p_i3_est),
          type='l', lty=1, col="light grey",
          main=mtext(bquote("Estimated Latent Cruve "*widehat('p')[i3])),
          xlab="Number of Datapoints", ylab=expression(widehat('p')[i3]),ylim=c(p3low,p3up))
  lines(pi3hat ,
        type='l', lty=1, lwd=2, col = "red")
  legend(x = "topright",          # Position
         legend = c("Individual", "Mean"),  # Legend texts
         #lty = c(1, 2),           # Line types
         col = c("grey", "red"),           # Line colors
         lwd = 2)                 # Line width


  #####################################################################################################
  #plot mean and estimated mean on the same scale
  #########################################################################################
  #p1
  # matplot(t(p_i1),
  #         type='l', lty=1, col="light blue",
  #         main = mtext(bquote("True and Estimated  "*'p'[i1])),
  #         xlab="Number of Datapoints", ylab=expression('p'[i1]),ylim = c(p1low,p1up))     #sparse p1 y 0.75-0.95
  matplot( t(p_i1_est),
            type='l', lty=1, col="light coral",
            # main=mtext(bquote("Estimated Latent Cruve "*widehat('p')[i1])),
            xlab="Number of Datapoints", ylab=expression(widehat('p')[i1]))

  # lines(pi1,type='l', lty=1,lwd=6, col="dark blue",
  #       # main=mtext(bquote("Mean Latent Cruve "*'p'[i1])),
  #       xlab="Number of Datapoints", ylab=expression("Average"*'p'[i1]))
  lines(pi1hat ,
        type='l', lty=1, lwd=6, col = "dark red")


  #p2
  # matplot( t(p_i2),
  #          type='l', lty=1, col="light blue",
  #          main = mtext(bquote("True True and Estimated "*'p'[i2])),
  #          xlab="Number of Datapoints", ylab=expression('p'[i2]),ylim=c(p2low,p2up))    #sparse p1 y 0045-0.2

  matplot( t(p_i2_est),
            type='l', lty=1, col="light coral",
            # main=mtext(bquote("Estimated Latent Cruve "*widehat('p')[i2])),
            xlab="Number of Datapoints", ylab=expression(widehat('p')[i2]))

  # lines(pi2,type='l', lty=1,lwd=6,  col="dark blue",
  #       # main=mtext(bquote("Mean Latent Cruve "*'p'[i2])),
  #       xlab="Number of Datapoints", ylab=expression("Average"*'p'[i2]))
  lines(pi2hat ,
        type='l', lty=1, lwd=6, col = "dark red")

  #p3
  # matplot( t(p_i3),
  #          type='l', lty=1, col="light blue",
  #          main = mtext(bquote("True and Estimated "*'p'[i3])),
  #          xlab="Number of Datapoints", ylab=expression('p'[i3]),ylim=c(p3low,p3up))
  matplot( t(p_i3_est),
            type='l', lty=1, col="light coral",
            # main=mtext(bquote("Estimated Latent Cruve "*widehat('p')[i3])),
            xlab="Number of Datapoints", ylab=expression(widehat('p')[i3]))


  # lines(pi3,type='l', lty=1, lwd=6, col="dark blue",
  #       # main=mtext(bquote("Mean Latent Cruve "*'p'[i3])),
  #       xlab="Number of Datapoints", ylab=expression("Average"*'p'[i3]))
  lines(pi3hat,
        type='l', lty=1, lwd=6, col = "dark red")
  # legend(x = "topleft",          # Position
  #        legend = c("True Mean", "Estimated Mean"),  # Legend texts
  #        #lty = c(1, 2),           # Line types
  #        col = c("grey", "red"),           # Line colors
  #        lwd = 2)                 # Line width
  #


  legend(x = "topright",          # Position
         legend = c("True Mean", "Estimated Mean"),  # Legend texts
         #lty = c(1, 2),           # Line types
         col = c("dark blue", "dark red"),           # Line colors
         lwd = 2)                 # Line width


  ###############################################################################################
  #plot2 truth and bayesian mfpca
  #4/8/2022 truth and mfpca
  ###############################################################################################

  #truth and bayesian estimate
  par(mfrow=c(3,2))
  #Z_i
  #meanzi1150=apply(t(Z_i1_150),1,mean)
  # matplot( t(Z_i1),
  #          type='l', lty=1, col="light grey",
  #          main = mtext(bquote('True Latent Curve '*' Z'[i1])),
  #          xlab="Number of Datapoints", ylab=expression('Z'[i1]),ylim=c(z1low,z2up))     #sparse z1
  # lines(meanzi1 ,
  #       type='l', lty=1, lwd=2, col = "red")
  # 
  # legend(x = "topleft",  horiz = TRUE,        # Position
  #        legend = c("Individual", "Mean"),  # Legend texts
  #        #lty = c(1, 2),           # Line types
  #        col = c("grey", "red"),           # Line colors
  #        lwd = 2)                 # Line width
  # 
  # 
  # 
  # 
  # #meanzi2150=apply(t(Z_i2_150),1,mean)
  # matplot( t(Z_i2),
  #          type='l', lty=1, col="light grey",
  #          main = mtext(bquote('True Latent Curve '*' Z'[i2])),
  #          xlab="Number of Datapoints", ylab=expression('Z'[i2]),ylim=c(z2low,z2up))   #sparse z2
  # lines(meanzi2 ,
  #       type='l', lty=1, lwd=2, col = "red")





  matplot( pcaZ_i1hat,
           type='l', lty=1, col="light grey",
           main = mtext(bquote('Recovered Latent Curve  '* widehat('Z')[i1])),
           xlab="Number of Datapoints", ylab=expression(widehat("Z")[i1]),ylim=c(z1low,z1up))
  lines(fmeanzi1hat ,
        type='l', lty=1, lwd=2, col = "red")


  matplot( pcaZ_i2hat,
           type='l', lty=1, col="light grey",
           main = mtext(bquote('Recovered Latent Curve  '* widehat('Z')[i2])),
           xlab="Number of Datapoints", ylab=expression(widehat("Z")[i2]),ylim=c(z2low,z2up))
  lines(fmeanzi2hat ,
        type='l', lty=1, lwd=2, col = "red")




  ###############################################################################
  #plot truth and estimated mean on the same scale
  ##############################################################################
  #z1
  # matplot( t(Z_i1),
  #          type='l', lty=1, col="light blue",
  #          main = mtext(bquote('True and Recovered '*' Z'[i1])),
  #          xlab="Number of Datapoints", ylab=expression('Z'[i1]),ylim=c(z1low,z1up))     #sparse z1
  # matlines( pcaZ_i1hat,
  #           type='l', lty=1, col="light coral",
  #           # main = mtext(bquote('Recovered Latent Curve  '* widehat('Z')[i1])),
  #           xlab="Number of Datapoints", ylab=expression(widehat("Z")[i1]))
  # 
  # lines( meanzi1,
  #        type='l', lty=1, lwd=6,col="dark blue",
  #        # main = mtext(bquote(' Mean Latent Curve '*' Z'[i1])),
  #        xlab="Number of Datapoints", ylab=expression('Average'*'Z'[i1]))
  # lines(fmeanzi1hat ,
  #       type='l', lty=1, lwd=6, col = "dark red")


  #########################################################################
  #z2
  # matplot(t(Z_i2),
  #         type='l', lty=1, col="light blue",
  #         main = mtext(bquote('True and Recovered '*' Z'[i2])),
  #         xlab="Number of Datapoints", ylab=expression('Z'[i2]),ylim=c(z2low,z2up))   #sparse z2
  # 
  # matlines( pcaZ_i2hat,
  #           type='l', lty=1, col="light coral",
  #           # main = mtext(bquote('Recovered Latent Curve  '* widehat('Z')[i2])),
  #           xlab="Number of Datapoints", ylab=expression(widehat("Z")[i2]))
  # 
  # lines(meanzi2 ,
  #       type='l', lty=1, lwd=6, col = "dark blue")
  # 
  # lines( fmeanzi2hat ,
  #        type='l', lty=1, lwd=6, col="dark red",
  #        # main = mtext(bquote(' Mean Latent Curve '*' Z'[i2])),
  #        xlab="Number of Datapoints", ylab=expression('Average'*'Z'[i2]))
  # 
  # legend("topleft",   horiz = TRUE,       # Position
  #        legend = c("True Mean","MFPCA Mean"),  # Legend texts
  #        #lty = c(1, 2),           # Line types
  #        col = c("dark blue", "dark red"),           # Line colors
  #        lwd = 2)                 # Line width

  # legend(x = "topleft",          # Position
  #        legend = c("True Mean", "MFPCA Mean"),  # Legend texts
  #        #lty = c(1, 2),           # Line types
  #        col = c("grey", "red"),           # Line colors
  #        lwd = 2)                 # Line width


  ###############################################################
  #plot scores
  ##############################################################
  #all score matrix are n*k  n is subjet, k is number of eigen functions
  #   par(mfrow=c(2,2))
  #   n=n
  #   #true Z scores
  #   score_matrix=dataset$AverageTrue$Truemeanscore
  #   plot(score_matrix[,1:2],xlab=expression('True Score '* xi[i]*' [ ,1]'),ylab=expression('True Score '* xi[i]*'[ ,2]'),col=4,
  #        cex.lab=1.5,cex.axis=1.5,cex=1, main="Z score")
  #   text(score_matrix[,1],score_matrix[,2],labels=seq(1,n,length=n),cex=1)
  #
  #
  #   # #MFPCA p scores
  #   # plot(out.50.150.10$AverageMFPCA$pcapscore[,1:2],xlab='MFPCA PC p Score 1',ylab='MFPCA PC p Score 2',col=4,
  #   # 	cex.lab=1.5,cex.axis=1.5,cex=1)
  #   # text(out.50.150.10$AverageMFPCA$pcapscore[,1],out.50.150.10$AverageMFPCA$pcapscore[,2],labels=seq(1,n,length=n),cex=1)
  #
  #
  #   #Bayesian MFPCA p scores
  #   plot(dataset$AverageBMFPCA$bfpcapscore[,1:2],xlab=expression('Recovered Score '* widehat(xi[i])*' [ ,1]'),ylab=expression('Bayesian Score '*widehat(xi[i])*'[ ,2]'),col=4,
  #        cex.lab=1.5,cex.axis=1.5,cex=1, main="Recovered Probability Score")
  #   text(dataset$AverageBMFPCA$bfpcapscore[,1],dataset$AverageBMFPCA$bfpcapscore[,2],labels=seq(1,n,length=n),cex=1)
  #
  #
  #   #MFPCA Z scores
  #   # plot(out.50.150.10$AverageMFPCA$pcazscore[,1:2],xlab='MFPCA PC Z Score 1',ylab='MFPCA PC Z Score 2',col=4,
  #   # 	cex.lab=1.5,cex.axis=1.5,cex=1)
  #   # text(out.50.150.10$AverageMFPCA$pcazscore[,1],out.50.150.10$AverageMFPCA$pcazscore[,2],labels=seq(1,n,length=n),cex=1)
  #
  #
  #   #Bayesian MFPCA Z scores
  #   plot(dataset$AverageBMFPCA$bfpcazscore[,1:2],xlab=expression('Recovered Score '* widehat(xi[i])*' [ ,1]'),ylab=expression('Bayesian Score '*widehat(xi[i])*'[ ,2]'),col=4,
  #        cex.lab=1.5,cex.axis=1.5,cex=1, main="Recovered Z Score")
  #
  #
  # text(dataset$AverageBMFPCA$bfpcazscore[,1],dataset$AverageBMFPCA$bfpcazscore[,2],labels=seq(1,n,length=n),cex=1)



  ################################################################################################################
  #plot true scores and recovered(MFPCA+Baysglm) scores on the same scale

  ################################################################################################################

  # n=n
  #true Z scores
  par(mfrow=c(1,1))
  # score_matrix=datagenerated$Trueth$Truescore_matrix
  # bmfpca_zscores=baysglmresult$Bayesian_scores
  # plot(score_matrix[,1:2],xlab=expression( xi[i]*' [ ,1]'),ylab=expression(xi[i]*'[ ,2]'),col="grey",
  #      xlim=c(xslow,xsup), ylim=c(yslow,ysup),main=" Scores")
  bmfpca_zscores=baysglmresult$Bayesian_scores
  plot(bmfpca_zscores[,1:2],xlab=expression( xi[i]*' [ ,1]'),ylab=expression(xi[i]*'[ ,2]'),col="grey",
       xlim=c(xslow,xsup), ylim=c(yslow,ysup),main=" Scores")

  # plot(score_matrix[,1:2],xlab=expression( xi[i]*' [ ,1]'),ylab=expression(xi[i]*'[ ,2]'),col="grey",
  #      cex.lab=1.5,cex.axis=1.5,cex=1, main=" Scores")
  # text(score_matrix[,1],score_matrix[,2],labels=seq(1,n,length=n),cex=1)


  # #MFPCA p scores
  # plot(out.50.150.10$AverageMFPCA$pcapscore[,1:2],xlab='MFPCA PC p Score 1',ylab='MFPCA PC p Score 2',col=4,
  # 	cex.lab=1.5,cex.axis=1.5,cex=1)
  # text(out.50.150.10$AverageMFPCA$pcapscore[,1],out.50.150.10$AverageMFPCA$pcapscore[,2],labels=seq(1,n,length=n),cex=1)


  #Bayesian MFPCA Z scores
  # bmfpca_zscores=baysglmresult$Bayesian_scores
  # points(bmfpca_zscores[,1:2],col="blue")
  # 
  # legend(x = "topleft",          # Position
  #        legend = c("True Z score", "Recovered Z score"),  # Legend texts
  #        #lty = c(1, 2),           # Line types
  #        col = c("grey", "blue"),           # Line colors
  #        # lwd = 2
  #        pch = c(1,1),
  #        bty = "n")                 # Line width

}

graph_scores=function(datagenerated,baysglmresult,datapoints,n,z1low,z1up,z2low,z2up,p1low,p1up,p2low,p2up,p3low,p3up,xslow,xsup,yslow,ysup,st,et,t){
  par(mfrow=c(1,1))
  # score_matrix=datagenerated$Trueth$Truescore_matrix
  # bmfpca_zscores=baysglmresult$Bayesian_scores
  # plot(score_matrix[,1:2],xlab=expression( xi[i]*' [ ,1]'),ylab=expression(xi[i]*'[ ,2]'),col="grey",
  #      xlim=c(xslow,xsup), ylim=c(yslow,ysup),main=" Scores")
  bmfpca_zscores=baysglmresult$Bayesian_scores
  plot(bmfpca_zscores[,1:2],xlab=expression( xi[i]*' [ ,1]'),ylab=expression(xi[i]*'[ ,2]'),col="blue",
       xlim=c(xslow,xsup), ylim=c(yslow,ysup),main=" Scores")  
}

graph_scores(trueest,MFPCAr,2000,1000,-5,15,-5,10,-0.01,1.1,-0.01,1.1,-0.1,0.5,-20,40,-5,25,st,et,t)

graph_scores(trueest,MFPCAr4,2000,1000,-5,15,-5,10,-0.01,1.1,-0.01,1.1,-0.1,0.5,-20,40,-5,25,st,et,t)

vectort=c(which(rest$cluster==0))
vectort=c(which(res$cluster==0))
graph_zp(trueest,MFPCAr,vectort,2000,
                  z1low,z1up,z2low,z2up,p1low,p1up,
         p2low,p2up,p3low,p3up,xslow,xsup,yslow,ysup,st,et,t)
  
graph_zp=function(datagenerated,baysglmresult,vectort,datapoints,
                  z1low,z1up,z2low,z2up,p1low,p1up,p2low,p2up,p3low,p3up,xslow,xsup,yslow,ysup,st,et,t){
  
    datapoints=datapoints
    t=t
    st=st
    et=et
    n=length(vectort)
    ########################################################################################
    
    #plot 1: truth and smoothed
    ########################################################################################
    #entire block is  n*t dimension
    #truth Z, P and smoothed Z, P
   
    Z_i1_est=datagenerated$Est$EstimateZ_i1[vectort,]
    Z_i2_est=datagenerated$Est$EstimateZ_i2[vectort,]
    
    p_i1_est=datagenerated$Est$Estimatep_i1[vectort,]
    p_i2_est=datagenerated$Est$Estimatep_i2[vectort,]
    p_i3_est=datagenerated$Est$Estimatep_i3[vectort,]

    
    fmeanzi1hat=baysglmresult$MFPCAZ$meanfnZ[[1]]@X
    fmeanzi2hat=baysglmresult$MFPCAZ$meanfnZ[[2]]@X
    
   
    
    
    #bayesian mfpca z and bayesian mfpca p are all t*n
    pcaZ_i1hat=baysglmresult$BayesianEst$BaypcaZ_i1[vectort,]
    pcaZ_i2hat=baysglmresult$BayesianEst$BaypcaZ_i2[vectort,]
    # pcap_i1hat=dataset$AverageBMFPCA$bfpcap1
    # pcap_i2hat=dataset$AverageBMFPCA$bfpcap2
    # pcap_i3hat=dataset$AverageBMFPCA$bfpcap3
    
    #scores
    bmfpca_zscores=baysglmresult$Bayesian_scores[vectort,]
    
    
    
    # Z1 2 10 Z2 -1   8  p10.6 1    0  0.35       0   0.04
    #truth and smoothed
    
    par(mfrow=c(1,2))
    
    #Z_i
    meanzi1hat=apply(t(Z_i1_est),1,mean)
    matplot( t(Z_i1_est),
             type='l', lty=1, col="light grey",
             main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i1])),
             xlab="Number of Datapoints", ylab=expression(widehat('Z')[i1]),ylim=c(z1low,z1up))
    lines(meanzi1hat ,
          type='l', lty=1, lwd=2, col = "red")
    
    
    
    meanzi2hat=apply(t(Z_i2_est),1,mean)
    matplot( t(Z_i2_est),
             type='l', lty=1, col="light grey",
             main = mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i2])),
             xlab="Number of Datapoints", ylab=expression(widehat('Z')[i2]),ylim=c(z2low,z2up))
    lines(meanzi2hat ,
          type='l', lty=1, lwd=2, col = "red")
    
    legend(x = "topleft",  horiz = TRUE,        # Position
           legend = c("Individual", "Mean"),  # Legend texts
           #lty = c(1, 2),           # Line types
           col = c("grey", "red"),           # Line colors
           lwd = 2)
  
    
    ####################################################################################
    par(mfrow=c(1,3))
    #p_i
    
    pi1hat=apply(t(p_i1_est),1,mean)
    matplot( t(p_i1_est),
             type='l', lty=1, col="light grey",
             main=mtext(bquote("Estimated Latent Cruve "*widehat('p')[i1])),
             xlab="Number of Datapoints", ylab=expression(widehat('p')[i1]),ylim = c(p1low,p1up))
    lines(pi1hat ,
          type='l', lty=1, lwd=2, col = "red")
    
    
    pi2hat=apply(t(p_i2_est),1,mean)
    matplot(t(p_i2_est),
            type='l', lty=1, col="light grey",
            main=mtext(bquote("Estimated Latent Cruve "*widehat('p')[i2])),
            xlab="Number of Datapoints", ylab=expression(widehat('p')[i2]),ylim=c(p2low,p2up))
    lines(pi2hat ,
          type='l', lty=1, lwd=2, col = "red")
    
    pi3hat=apply(t(p_i3_est),1,mean)
    matplot(t(p_i3_est),
            type='l', lty=1, col="light grey",
            main=mtext(bquote("Estimated Latent Cruve "*widehat('p')[i3])),
            xlab="Number of Datapoints", ylab=expression(widehat('p')[i3]),ylim=c(p3low,p3up))
    lines(pi3hat ,
          type='l', lty=1, lwd=2, col = "red")
    legend(x = "topright",          # Position
           legend = c("Individual", "Mean"),  # Legend texts
           #lty = c(1, 2),           # Line types
           col = c("grey", "red"),           # Line colors
           lwd = 2)                 # Line width
    
    
    ###############################################################################################
    #plot2 truth and bayesian mfpca
    #4/8/2022 truth and mfpca
    ###############################################################################################
    
    #truth and bayesian estimate
    par(mfrow=c(1,2))
    #Z_i
    
    
    matplot( pcaZ_i1hat,
             type='l', lty=1, col="light grey",
             main = mtext(bquote('Recovered Latent Curve  '* widehat('Z')[i1])),
             xlab="Number of Datapoints", ylab=expression(widehat("Z")[i1]),ylim=c(z1low,z1up))
    lines(fmeanzi1hat ,
          type='l', lty=1, lwd=2, col = "red")
    
    
    matplot( pcaZ_i2hat,
             type='l', lty=1, col="light grey",
             main = mtext(bquote('Recovered Latent Curve  '* widehat('Z')[i2])),
             xlab="Number of Datapoints", ylab=expression(widehat("Z")[i2]),ylim=c(z2low,z2up))
    lines(fmeanzi2hat ,
          type='l', lty=1, lwd=2, col = "red")
    
    
    
    
    ###############################################################
    #plot scores
    ##############################################################
    # n=n
    #true Z scores
    par(mfrow=c(1,1))
    bmfpca_zscores=baysglmresult$Bayesian_scores
    plot(bmfpca_zscores[,1:2],xlab=expression( xi[i]*' [ ,1]'),ylab=expression(xi[i]*'[ ,2]'),col="grey",
         xlim=c(xslow,xsup), ylim=c(yslow,ysup),main=" Scores")
  
}



dist=kNNdist(combinedscore_z, k = knnum)
ninty5p=quantile(dist, probs = pct)



#distdata=data.frame(sort(dist))
#distdata$index=1:dim(distdata)[1]
# dp <- ggplot(distdata,aes(index,sort.dist.)) + geom_line()+ggtitle(paste0(knnum,"-NN Distance Plot ",'\n',"(",dim(distdata)[1]," Subjects",")")) +
#   xlab("Points sorted by Distance") + ylab("Distance")+ theme(plot.title = element_text(hjust = 0.5))+geom_hline(yintercept=ninty5p, color = "red")+
#   geom_text(data=data.frame(round(ninty5p,2)),
#             aes(x=dim(distdata)[1]/2,y=1.2*ninty5p,label=paste0("Distance at ",gsub("%$","",row.names(data.frame(round(ninty5p,2)))),"th percentile= ",round(ninty5p,2))))
# 
# 
# # ggsave(paste0(ggsavepath,"knndisplot100sub.png"))
# ggsave("knndisplot100sub.png")

#cfda
# distcfda=kNNdist(combinedscore_zcfda, k = knnum)
# ninty5pcfda=quantile(distcfda, probs = pctcfda)
# distdatacfda=data.frame(sort(distcfda))
# distdatacfda$index=1:dim(distdatacfda)[1]
# dpcfda <- ggplot(distdatacfda,aes(index,sort.distcfda.)) + geom_line()+ggtitle(paste0(knnum,"-NN Distance Plot (CFDA) ",'\n',"(",dim(distdatacfda)[1]," Subjects",")")) +
#   xlab("Points sorted by Distance") + ylab("Distance")+ theme(plot.title = element_text(hjust = 0.5))+geom_hline(yintercept=ninty5pcfda, color = "red")+
#   geom_text(data=data.frame(round(ninty5pcfda,2)),
#             aes(x=dim(distdatacfda)[1]/2,y=1.2*ninty5pcfda,label=paste0("Distance at ",gsub("%$","",row.names(data.frame(round(ninty5pcfda,2)))),"th percentile= ",round(ninty5pcfda,2))))
# ggsave("cfdaknndisplot100sub.png")

#https://rpubs.com/datalowe/dbscan-simple-example
# res <- dbscan(combinedscore_z, eps =1.2 , minPts = 3)
library(ggplot2)
library(dbscan)
knnum=3
#combinedscore_z=MFPCAr$Bayesian_scores
#ramsay
combinedscore_z=mfpcaram$scores
dist=kNNdist(combinedscore_z, k = knnum)
distdata=data.frame(sort(dist))
distdata$index=1:dim(distdata)[1]
#pct=0.98.  #happ
pct=0.98
ninty5p=quantile(dist, probs = pct)
dp <- ggplot(distdata,aes(index,sort.dist.)) + geom_line()+ggtitle(paste0(knnum,"-NN Distance Plot ",'\n',"(",dim(distdata)[1]," Subjects",")")) +
   xlab("Points sorted by Distance") + ylab("Distance")+ theme(plot.title = element_text(hjust = 0.5))+geom_hline(yintercept=ninty5p, color = "red")+
   geom_text(data=data.frame(round(ninty5p,2)),
            aes(x=dim(distdata)[1]/2,y=1.2*ninty5p,label=paste0("Distance at ",gsub("%$","",row.names(data.frame(round(ninty5p,2)))),"th percentile= ",round(ninty5p,2))))
dp

res <- dbscan(combinedscore_z, eps =ninty5p, minPts = 5)
table(res$cluster)

tminPts = 5    #7, 0.98 2 clusters  6, 88
tninty5p=0.98
tscores=combinedscore_z
tres <- dbscan(tscores, eps =tninty5p , minPts = tminPts )


#tclustertable=table(tres$cluster)
#tclustertable    #z score


tclustertable=table(res$cluster)
tclustertable  
#xtable(as.matrix(t(tclustertable )))

#tclusterdata=data.frame(tscores)
tclusterdata=data.frame(combinedscore_z)
tclusterdata$Cluster=as.factor(res$cluster)
save(tclusterdata,file="ramsaycluster.RData")
#
# tp <- ggplot(tclusterdata,aes(ksi1,ksi2,colour = Cluster)) + geom_point(aes(shape=Cluster))+ggtitle(paste0("Twitter Users Cluster Results",'\n',"(",dim(tclusterdata)[1]," Subjects",", minPts=",tminPts,")")) +
#   xlab(expression('Score '* widehat(xi[i1]))) + ylab(expression('Score '* widehat(xi[i2])))+ theme(plot.title = element_text(hjust = 0.5),axis.text=element_text(size = 2),axis.title = element_text(size =4))
# tp

# tps <- ggplot(tclusterdata,aes(ksi1,ksi2,colour = Cluster)) + geom_point(aes(shape=Cluster),size=4)+ggtitle(paste0("Twitter Users Cluster Results",'\n',"(",dim(tclusterdata)[1]," Subjects",", minPts=",tminPts,")")) +
#   xlab(expression('Score '* widehat(xi[i1]))) + ylab(expression('Score '* widehat(xi[i2])))+ theme(plot.title = element_text(hjust = 0.5))+#theme(text = element_text(size = 20))+
#   #theme(legend.text = element_text(size = 30))  
#   #theme(axis.text=element_text(size = 20),axis.title = element_text(size =30))
#   theme(text=element_text(size = 20))

library(ggplot2)
#tps <- ggplot(tclusterdata,aes(ksi1,ksi2,colour = Cluster)) + geom_point(aes(shape=Cluster),size=3)+ggtitle(paste0("Twitter Users Cluster Results",'\n',"(",dim(tclusterdata)[1]," Subjects",")")) +
 # xlab(expression('Score '* widehat(xi[i1]))) + ylab(expression('Score '* widehat(xi[i2])))+ theme(plot.title = element_text(hjust = 0.5))+#theme(text = element_text(size = 20))+
  #theme(legend.text = element_text(size = 30))  
  #theme(axis.text=element_text(size = 20),axis.title = element_text(size =30))
  #theme(text=element_text(size = 20))
#tps

#ramsay
tps <- ggplot(tclusterdata,aes(X1,X2,colour = Cluster)) + geom_point(aes(shape=Cluster),size=3)+ggtitle(paste0("Twitter Users Cluster Results",'\n',"(",dim(tclusterdata)[1]," Subjects",")")) +
  xlab(expression('Score '* widehat(xi[i1]))) + ylab(expression('Score '* widehat(xi[i2])))+ theme(plot.title = element_text(hjust = 0.5))+#theme(text = element_text(size = 20))+
  #theme(legend.text = element_text(size = 30))  
  #theme(axis.text=element_text(size = 20),axis.title = element_text(size =30))
  theme(text=element_text(size = 20))
tps

c(which(tclusterdata$Cluster==5))

vec0=c(which(tclusterdata$Cluster==0))
vec1=c(which(tclusterdata$Cluster==1))
vec2=c(which(tclusterdata$Cluster==2))
######
#try simulation application data
hist(mfpcaram$scores[vec0,1],main="noise score 1")
hist(mfpcaram$scores[vec0,2],main="noise score 2")


hist(mfpcaram$scores[vec1,1],main="cluster one score 1")
hist(mfpcaram$scores[vec1,2],main="cluster one score 2")


hist(mfpcaram$scores[vec2,1],main="cluster two score 1")
hist(mfpcaram$scores[vec2,2],main="cluster two score 1")

plot(apply(mfpcaram$latentcurve[vec0,,1],2,mean))
plot(apply(mfpcaram$latentcurve[vec0,,2],2,mean))

plot(apply(mfpcaram$latentcurve[vec1,,1],2,mean))
plot(apply(mfpcaram$latentcurve[vec1,,2],2,mean))

plot(apply(mfpcaram$latentcurve[vec2,,1],2,mean))
plot(apply(mfpcaram$latentcurve[vec2,,2],2,mean))


#plot latent curve by cluster
#plot latent curves and probability curves
#vectort=c(which(cluster==0))
#plot_latent=function(zlatent,phat,st,et,cluster,labelnum,argval,zmean){
zlatent=array(c(trueest$Est$EstimateZ_i1,trueest$Est$EstimateZ_i2),dim=c(1000,2000,2))
phat=array(c(trueest$Est$Estimatep_i1,trueest$Est$Estimatep_i2,trueest$Est$Estimatep_i3),
           dim=c(1000,2000,3))
#plotbycluster
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
  
  #Z_i
  #meanz=colMeans(zest)   #t*numz
  # meanz=zmean
  #meanfnZ[[1]]@X
  meanz=colMeans(zest)
  for (i in 1:numz){
    matplot(argval, t(zest[,,i]),
            type='l', lty=1, col="light grey",
            #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
            #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
            
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
            #main=mtext(bquote("Estimated Latent Cruve "*widehat('p')[i])),
            #xlab="Number of Datapoints", ylab=expression(widehat('p')[i]))
            main=mtext(bquote("Estimated Probability Cruve ")),
            xlab="Day", ylab='Value',cex.lab = 2,cex.axis = 2,cex.main=2,xaxt="n")
    
    lines(argval,meanp[,i] ,
          type='l', lty=2, lwd=2, col = "red")
    axis(1,
         at = c(0,476/2000,952/2000,1428/2000,1904/2000),cex.lab = 2,cex.axis = 2,cex.main=2,
         labels = c("0","5","10","15","20"))
     #legend(x = "topright",          # Position
            
    #legend = c("Individual", "Mean"),  # Legend texts
            #lty = c(1, 2),           # Line types
           #col = c("grey", "red"),           # Line colors
           #lwd = 2)                 # Line width
  }
  
  }

plot_latentfn(zlatent,phat,
              0,1,tclusterdata$Cluster,0,seq(0.1,0.9,length=2000),tclusterdata)

plot_latentfn(zlatent,phat,
              0,1,tclusterdata$Cluster,1,seq(0.1,0.9,length=2000),tclusterdata)
plot_latentfn(zlatent,phat,
              0,1,tclusterdata$Cluster,2,seq(0.1,0.9,length=2000),tclusterdata)

##########################
par(mfrow=c(1,1))
matplot(seq(0,1,length=2000),t(trueest$Est$EstimateZ_i1),ylim = c())
##simulation design for noise data
zclusterone=zlatent[vec1,,]
zclustertwo=zlatent[vec2,,]
zclusterzero=zlatent[vec0,,]
zclusteronemean=colMeans(zclusterone) #this will be t*(Q-1)=2000*2 (z1,z2) dimension
zclustertwomean=colMeans(zclustertwo) #this will be t*(Q-1)=2000*2 (z1,z2) dimension
zclusterzeromean=colMeans(zclusterzero) #this will be t*(Q-1)=2000*2 (z1,z2) dimension
##################
hist(mfpcaram$scores[vec0,1])
hist(mfpcaram$scores[vec1,1])
hist(mfpcaram$scores[vec2,1])

hist(mfpcaram$scores[vec0,2])
hist(mfpcaram$scores[vec1,2])
hist(mfpcaram$scores[vec2,2])
########################################
mean(mfpcaram$scores[vec0,1])
sd(mfpcaram$scores[vec0,1])


mean(mfpcaram$scores[vec1,1])
sd(mfpcaram$scores[vec1,1])

mean(mfpcaram$scores[vec2,1])
sd(mfpcaram$scores[vec2,1])
########################################
mean(mfpcaram$scores[vec0,2])
sd(mfpcaram$scores[vec0,2])


mean(mfpcaram$scores[vec1,2])
sd(mfpcaram$scores[vec1,2])

mean(mfpcaram$scores[vec2,2])
sd(mfpcaram$scores[vec2,2])
########################################
