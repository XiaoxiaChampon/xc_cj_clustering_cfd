############################################################
# Article: Functions Needed for Clustering Categorical Functional Data
#           File that contains all the functions necessary to generate data 
# Author:  Ana-Maria Staicu
# Date: 07/01/2023
##############################################################

# n100t300A
# $cluster_table_true
# fadp RI   fadp ARI   fadp cpn  kmeans RI kmeans ARI kmeans cpn  dbscan RI dbscan ARI dbscan cpn 
# 0.6527273  0.3544619  0.0000000  0.9866667  0.9718616  0.0000000  0.9952862  0.9901378  1.0000000 
# 
# $cluster_table_est
# fadp RI   fadp ARI   fadp cpn  kmeans RI kmeans ARI kmeans cpn  dbscan RI dbscan ARI dbscan cpn 
# 0.8210774  0.5892844  0.0000000  0.9830976  0.9642518  0.3333333  0.7344108  0.3527636  0.7777778 
# 
# $cluster_table_est_se
# fadp RI    fadp ARI    fadp cpn   kmeans RI  kmeans ARI  kmeans cpn   dbscan RI  dbscan ARI  dbscan cpn 
# 0.147394486 0.343605896 0.000000000 0.007609428 0.016179713 0.192450090 0.098584847 0.254324619 0.222222222 

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
  
  # Define MEAN and SCORE VARIANCES for the 3 settings (#***! something wrong with scenario B!!)
  if (setting==1){
    if (scenario=="A"){ 
      mu_1=function(t)   -1+2*t+2*t^2  
      mu_2=function(t)   -2.5+exp(t*2) 
      #
      score.var[1]<- 1
      score.var[2]<- 1/2
      score.var[3]<- 1/4  }
    
    if (scenario=="B"){
      mu_1=function(t) -1+2*t+2*t^2
      mu_2=function(t)  -0.5+ exp(t*2) 
      #
      score.var[1]<- 1
      score.var[2]<- 1/2
      score.var[3]<- 1/4  }
  
    if (scenario=="C"){ 
      mu_1=function(t)   -1+2*t+2*t^2  
      mu_2=function(t)   -2.5+exp(t*2) 
      #
      score.var[1]<- 100
      score.var[2]<- 10
      score.var[3]<- 1  }
    
    }
  
  if(setting==2){
    mu_1=function(t) 4*t^2 -1.2          #  -1+3.5*(t^4)
    mu_2=function(t) 4*t^2 -3.5          # 0.5-(2*t-sqrt(3))^2 
    #
    score.var[1]<- 1
    score.var[2]<- 1/2
    score.var[3]<- 1/4}
  
  if(setting==3){
    mu_1=function(t) -2.2 + 4*t^2 #-2*cos(pi*t) 
    mu_2=function(t) -7+6*t^2     #5*t^2-6
    #
    score.var[1]<- 1
    score.var[2]<- 1/4
    score.var[3]<- 1/16}
  
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
  
  # matplot(tt,Z1, type="l")
  # matlines(tt,Z2, type="l")
  # Q: can we remove the ZEROS!? 
  
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



#**** SANITY CHECK  PLOT Z's and P's
# par(mfrow=c(1,2))
# matplot(tt, Z1, type="l" )  
# 
# matplot(tt, Z2, type="l" )  
# 
# 
# par(mfrow=c(1,3))
# matplot(tt, p1, type="l" )  
# matplot(tt, p2, type="l" )  
# matplot(tt, p3, type="l" )  
# 
# p<- list(p1=p1, p2=p2, p3=p3)


#Function to generate Categ Functional Data, given the latent event probabilities
# INPUT
# p list of probability matrices that are m by n with probabilities for n subjects, 
#        each observed m times; elements called p1, p2, p3
#
# OUTPUT
# categ FD - m by n with Q=3  categories 
generate_CategFD_scenario=function(p, seed=NULL){
  # p=p; seed=5698
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

#add to check the category
####################################
#Cluster One
# pie1=as.data.frame(table((W[,1:n1])))
# pie1$pct=round(pie1$Freq/(sum(pie1$Freq)),2)
# pie1$labels= scales::percent(pie1$pct)
# colnames(pie1)[colnames(pie1) == "Var1"] <- "Category"
# 
# ggplot(pie1, aes(x = "", y = pct, fill = Category)) +
#   geom_col() +
#   geom_text(aes(label = labels),size=10,
#             position = position_stack(vjust = 0.5)) +
#   coord_polar(theta = "y")+ggtitle("Setting One")+theme(plot.title = element_text(hjust = 0.5))+
#   theme(text=element_text(size = 20))
# 
# ###############
# pie2=as.data.frame(table((W[,(n1+1):(n1+n2)])))
# pie2$pct=round(pie2$Freq/(sum(pie2$Freq)),2)
# pie2$labels= scales::percent(pie2$pct)
# colnames(pie2)[colnames(pie2) == "Var1"] <- "Category"
# 
# ggplot(pie2, aes(x = "", y = pct, fill = Category)) +
#   geom_col() +
#   geom_text(aes(label = labels),size=10,
#             position = position_stack(vjust = 0.5)) +
#   coord_polar(theta = "y")+ggtitle("Setting Two")+theme(plot.title = element_text(hjust = 0.5))+
#   theme(text=element_text(size = 20))
# ##########
# pie3=as.data.frame(table((W[,(n1+n2+1):n])))
# pie3$pct=round(pie3$Freq/(sum(pie3$Freq)),2)
# pie3$labels= scales::percent(pie3$pct)
# colnames(pie3)[colnames(pie3) == "Var1"] <- "Category"
# 
# ggplot(pie3, aes(x = "", y = pct, fill = Category)) +
#   geom_col() +
#   geom_text(aes(label = labels),size=10,
#             position = position_stack(vjust = 0.5)) +
#   coord_polar(theta = "y")+ggtitle("Setting Three")+theme(plot.title = element_text(hjust = 0.5))+
#   theme(text=element_text(size = 20))


#######################################



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
  # X=X
  # basis_size<-25
  # method="ML"
  # sim=TRUE
  # 
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
  p<-array(0, c(n, m, 3))
  ##########################
  ###########################
  #n is the subject and this step is done by subject level
  #can parallel
  #############################
  #############################
  for (i in 1:n){
    #i=93
    #basis_size=25
      x1<- X[i,,1]
      x2<- X[i,,2]
      x3<- X[i,,3]
      
      # fit the Binorm model
      basis_size_rev<-max(min( round( min(sum(x1), sum(1-x1))/2), basis_size ), 5)
      fit_binom<-gam(x1~s(tt, bs = "cr", m=2, k = basis_size_rev),
                     family=binomial(link="probit"), method = method,control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),optimizer=c("outer","bfgs"))
      p1 <- fit_binom$fitted.values
      p1_linpred <- fit_binom$linear.predictors
      
      # basis_size=25
      # method="ML"
      basis_size_rev<-max(min(round(  min(sum(x2), sum(1-x2))/2), basis_size ), 5)
      fit_binom<- gam(x2~s(tt, bs = "cr", m=2, k = basis_size_rev),
                      family=binomial(link="probit"), method = method,control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),optimizer=c("outer","bfgs"))
      # 
      #  #
      # 
      #  #################################
      #  #################################
      #  #7/28/2023
      # ##add condition for p2
      # 
      p2 <- fit_binom$fitted.values
      p2_linpred <- fit_binom$linear.predictors
      #  # #7/28/2023
      #  # #################################
      #  # #################################
      #  #
      #  #
      basis_size_rev<-max(min(round( min(sum(x3), sum(1-x3))/2), basis_size ), 5)
      fit_binom <- gam(x3~s(tt, bs = "cr", m=2, k = basis_size_rev),
                       family=binomial(link="probit"), method = method,control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),optimizer=c("outer","bfgs"))
      p3 <- fit_binom$fitted.values
      p3_linpred <- fit_binom$linear.predictors
      
      ###################################################
      ###################################################
      ##updated estimation
      # if ((sum(1-x1)/sum(x1))<0.004 && m<301){
      #   basis_size_rev_1<-max(min(round(  min(sum(x1), sum(1-x1))/2), basis_size ), 5)
      # 
      # 
      #   fit_binom_1<- gam(x1~s(tt, bs = "cr", m=2, k = basis_size_rev_1),
      #                     family=binomial(link="probit"), method = method,
      #                     control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
      #                     optimizer=c("outer","bfgs"))
      # }
      # 
      # else{
      #   basis_size_rev_1<-max(min(round(  min(sum(x1), sum(1-x1))/2), basis_size ), 5)
      # 
      #   fit_binom_1<- gam(x1~s(tt, bs = "cr", m=2, k = basis_size_rev_1),
      #                     family="binomial", method = method,
      #                     control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
      #                     optimizer=c("outer","bfgs"))
      # }
      # 
      # 
      # 
      # p1 <- fit_binom_1$fitted.values
      # p1_linpred <- fit_binom_1$linear.predictors
      # # 
      # # # basis_size_rev<-max(min(round(  min(sum(x2), sum(1-x2))/2), basis_size ), 5)
      # # # fit_binom<- gam(x2~s(tt, bs = "cr", m=2, k = basis_size_rev),
      # # #                 family="binomial", method = method,control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),optimizer=c("outer","bfgs"))
      # # # 
      # if ((sum(1-x2)/sum(x2))<0.004 && m<301){
      #   basis_size_rev_2<-max(min(round(  min(sum(x2), sum(1-x2))/2), basis_size ), 5)
      # 
      #   fit_binom_2<- gam(x2~s(tt, bs = "cr", m=2, k = basis_size_rev_2),
      #                     family=binomial(link="probit"), method = method,
      #                     control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
      #                     optimizer=c("outer","bfgs"))
      # }
      # 
      # else{
      #   basis_size_rev_2<-max(min(round(  min(sum(x2), sum(1-x2))/2), basis_size ), 5)
      # 
      #   fit_binom_2<- gam(x2~s(tt, bs = "cr", m=2, k = basis_size_rev_2),
      #                     family="binomial", method = method,
      #                     control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
      #                     optimizer=c("outer","bfgs"))
      # }
      # 
      # p2 <- fit_binom_2$fitted.values
      # p2_linpred <- fit_binom_2$linear.predictors
      # # 
      # # # basis_size_rev<-max(min(round( min(sum(x3), sum(1-x3))/2), basis_size ), 5)
      # # # fit_binom <- gam(x3~s(tt, bs = "cr", m=2, k = basis_size_rev),
      # # #                  family="binomial", method = method,control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),optimizer=c("outer","bfgs"))
      # # # p3 <- fit_binom$fitted.values
      # # # p3_linpred <- fit_binom$linear.predictors
      # if ((sum(1-x3)/sum(x3))<0.004 && m<301){
      #   basis_size_rev_3<-max(min(round(  min(sum(x3), sum(1-x3))/2), basis_size ), 5)
      # 
      # 
      #   fit_binom_3<- gam(x3~s(tt, bs = "cr", m=2, k = basis_size_rev_3),
      #                     family=binomial(link="probit"), method = method,
      #                     control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
      #                     optimizer=c("outer","bfgs"))
      # }
      # 
      # else{
      #   basis_size_rev_3<-max(min(round(  min(sum(x3), sum(1-x3))/2), basis_size ), 5)
      # 
      #   fit_binom_3<- gam(x3~s(tt, bs = "cr", m=2, k = basis_size_rev_3),
      #                     family="binomial", method = method,
      #                     control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
      #                     optimizer=c("outer","bfgs"))
      # }
      # 
      # 
      # 
      # p3 <- fit_binom_3$fitted.values
      # p3_linpred <- fit_binom_3$linear.predictors
      ##updated estimation
      ###################################################
      ###################################################
      
      
      
      
      # estimate the latent tranjecotries Z
      z1<- (p1_linpred-p3_linpred)-log( (1+exp(p1_linpred))/(1+exp(p3_linpred)))
      z2<- (p2_linpred-p3_linpred)-log( (1+exp(p2_linpred))/(1+exp(p3_linpred)))
      
      Z<- cbind(Z, c(z1,z2))
      p[i,,] <- cbind(p1/(p1+p2+p3), p2/(p1+p2+p3), p3/(p1+p2+p3))
    
    
    
  }
  
  return(list(Z1_est=Z[1:m,], Z2_est=Z[1:m+m,], 
              p1_est=t(p[,,1]), p2_est=t(p[,,2]), p3_est=t(p[,,3]) ))
}


estimate_categFD_W <- function( W=NULL, tt, basis_size=25, method="ML", sim=TRUE){
  # X=X
  # basis_size<-25
  # method="ML"
  # sim=TRUE
  # 
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
  p<-array(0, c(n, m, 3))
  ##########################
  ###########################
  #n is the subject and this step is done by subject level
  #can parallel
  #############################
  #############################
  method="ML"
  for (i in 1:n){
      print(i)
      fit_binom<-gam(list(W[,i]-1~s(tt)+s(tt),~s(tt)+s(tt)),family=multinom(K=2),method = method,
                     control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                     optimizer=c("outer","bfgs")) 
      p1 <- fit_binom$fitted.values[,1]
      p2 <- fit_binom$fitted.values[,2]
      p[i,,] <- cbind(p1, p2, 1-p1-p2)
      
      z1<- fit_binom$linear.predictors[,1]
      z2<- fit_binom$linear.predictors[,2]
      Z<- cbind(Z, c(z1,z2))
      # gam(list(y~s(x1)+s(x2),~s(x1)+s(x2)),family=multinom(K=2))
      #Z<- cbind(Z, c(z1,z2))
      #p[i,,] <- cbind(p1/(p1+p2+p3), p2/(p1+p2+p3), p3/(p1+p2+p3))
      # return(list(Z1_est=Z[1:m,], Z2_est=Z[1:m+m,], 
      #             p1_est=t(p[,,1]), p2_est=t(p[,,2]), p3_est=t(p[,,3]) ))
    
    
  }
  
  return(list(Z1_est=Z[1:m,], Z2_est=Z[1:m+m,], 
              p1_est=t(p[,,1]), p2_est=t(p[,,2]), p3_est=t(p[,,3]) ))
}


##################
estimate_categFD_W(W, tt, basis_size=25, method="ML", sim=FALSE)



###plot Figure 4
####graph z1,z2 for scenario B
# dscenarioB=data.frame(x=c(Z1_est[,c(1,80,100)]), y=c(Z2_est[,c(1,80,100)]), 
#                       Setting=c(rep("One",dim(Z1_est[,c(1,80,100)])[1]),
#                                 rep("Two",dim(Z1_est[,c(1,80,100)])[1]),
#                                 rep("Three",dim(Z1_est[,c(1,80,100)])[1])))
# ggplot() +
#   geom_polygon(data=dscenarioB, mapping=aes(x=x, y=y, group=Setting),alpha = 0.1) +
#   geom_point(data=dscenarioB, aes(x=x, y=y, color=Setting),show.legend = FALSE) +ggtitle("Scenario B")+
#   theme(plot.title = element_text(hjust = 0.5))+ xlab("Latent Curve One") + 
#   ylab("Latent Curve Two")+theme(text=element_text(size = 20))+xlim(-4, 6)+ylim(-8, 5.5)+
#   geom_text(aes(x = 2, y = 4,
#                 label = "cluster one"),
#             stat = "unique",
#             size = 10, color = "red")+
#   geom_text(aes(x = 4.5, y = 0,
#                 label = "cluster two"),
#             stat = "unique",
#             size = 10, color = "blue")+
#   geom_text(aes(x = -2.7, y = -5.5,
#                 label = "cluster three"),
#             stat = "unique",
#             size = 10, color = "green")
# 
# 
# 
# dscenarioA=data.frame(x=c(Z1_est[,c(1,80,100)]), y=c(Z2_est[,c(1,80,100)]), 
#                       Setting=c(rep("One",dim(Z1_est[,c(1,80,100)])[1]),
#                                 rep("Two",dim(Z1_est[,c(1,80,100)])[1]),
#                                 rep("Three",dim(Z1_est[,c(1,80,100)])[1])))
# ggplot() +
#   geom_polygon(data=dscenarioA, mapping=aes(x=x, y=y, group=Setting),alpha = 0.1) +
#   geom_point(data=dscenarioA, aes(x=x, y=y, color=Setting),show.legend = FALSE) +ggtitle("Scenario A")+
#   theme(plot.title = element_text(hjust = 0.5))+ xlab("Latent Curve One") + 
#   ylab("Latent Curve Two")+theme(text=element_text(size = 20))+xlim(-4, 6)+ylim(-8, 5.5)+
#   geom_text(aes(x = 2, y = 4.5,
#                 label = "cluster one"),
#             stat = "unique",
#             size = 10, color = "red")+
#   geom_text(aes(x = 3.5, y = -1,
#                 label = "cluster two"),
#             stat = "unique",
#             size = 10, color = "blue")+
#   geom_text(aes(x = -2.7, y = -5.5,
#                 label = "cluster three"),
#             stat = "unique",
#             size = 10, color = "green")



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
   # mZ1=Z1_est
   # mZ2=Z2_est
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

#############################################################################
#############################################################################
# Alternative methods NOT USED in the paper to perform multivariate FPCA
#*** Extract scores using multivariate FPCA
# extract_scores_MFPCA <- function (mZ1,mZ2, tt=tt , PVE=0.95){
#   
#   f1 <- funData(argvals = list(tt), X = t(mZ1)) 
#   f2 <- funData(argvals = list(tt), X = t(mZ2)) 
#   
#   uFPCA <- MFPCA(multiFunData(list(f1,f2)), M = 6, 
#                  uniExpansions = list(list(type = "uFPCA"),list(type = "uFPCA")))
#   
#   scores<- uFPCA$scores
#   K <- which(cumsum(apply(scores, 2, var))/sum(apply(scores, 2, var))>=PVE)[1]
#   scores_K<- scores[, 1:K] 
#   
#   return(list(scores=scores_K,Phi= uFPCA$functions, Mean=uFPCA$meanFunction))
# }

#*** extract scores using EIGEN of the sample covariance
# extract_scores_eigen <- function (mZ1,mZ2, tt=tt , PVE=0.95){
#   
#   C11<- cov(t(mZ1))
#   C12<- cov(t(mZ1), t(mZ2))
#   C22<- cov(t(mZ2))
#   
#   C<- rbind(cbind(C11, C12), cbind(t(C12), C22))
#   out<-eigen(C)
#   K1 <- which(cumsum(out$values)/sum(out$values)>=PVE)[1]
#   Phi_est <- out$vectors[, 1:K1]*sqrt(m)
#   Lambda_est <- out$values[1:K1]
#   
#   mZ <- rbind(mZ1, mZ2) 
#   Scores_est <- t(mZ) %*%Phi_est/sqrt(m)  # they are not demeaned
#   
#   return(list(scores=Scores_est,Phi= Phi_est))
# }


############################################################
# Article: Functions Needed for Clustering Categorical Functional Data
#           Main file to generate data and obtain simulation results
# Author:  Ana-Maria Staicu
# Date: 07/01/2023
##############################################################

#rm(list=ls())


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
devtools::install_github("ahasverus/elbow")
library(elbow)
library(pdfCluster)
library(FADPclust)

#setwd("~/Dropbox/0CurrentProjects/1AWorking/7XiaoxiaChampon/Code_Project1")

#source("utilities_generate.R")
#source("utilities_estimation.R")
###########################################################

###################################################################################
#distributed to different
# simulation
cluster_simulation=function(n,m,scenario,mc_sims,wood_multi){
    n=100
    m=300
    scenario="B"
    mc_sims=3
   wood_multi="N"
  
  # set main factors
 # n<-100 # no users
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
    p1=0.1
    p2=0.6
    p3=0.3
  }
  
  
  #n1<-n*0.75; n2<- n*0.22; n3<-n*0.03
  n1<-n*p1; 
  n2<- n*p2; 
  n3<-n*p3
  true_cluster <- c(rep(1, n1), rep(2, n2), rep(3, n3))
  true_cluster_db <- c(rep(1, n1), rep(2, n2), rep(0, n3))
  #m=300 # no of time points
  tt <- seq(from = 0.0001, to = 1, length=m) # gridpoint t
  
  # seed for latent process
  seed1<-1230
  # seed for categFD, given the latent process
  seed0<- 9876
  
  #mc_sims<-5
  
  # rmse<- array(0, c(mc_sims, 2, 3))       #2 components and 3 settings
  # hellinger <- array(0, c(mc_sims, 3, 3)) #3 events probab and 3 settings
  
  true_kmeans <- est_kmeans  <- NULL # records clustering membership on the TRUE SCORES/ESTIMATED SCORES
  true_fadp <- est_fadp <-  NULL # records clustering membership on the TRUE Z/ESTIMATED Z
  true_dbscan <- est_dbscan <-  NULL
  
  ##############
  true_kmeans_nclust <- est_kmeans_nclust  <- NULL # records clustering membership on the TRUE SCORES/ESTIMATED SCORES
  true_fadp_nclust <- est_fadp_nclust <- NULL # records clustering membership on the TRUE Z/ESTIMATED Z
  true_dbscan_nclust <- est_dbscan_nclust <- NULL 
  
  clust_record=matrix(0,nrow=mc_sims,ncol=6)
  ##############
  
  #print("today")
  for(ii in 1:mc_sims){
    #ii=2
    set.seed(seed1+100*ii)
    
    # # generate clusters
    # CLSTF1 <- generate_data_scenario(k=3,n=n1, m=m, setting=1, scenario="B", Q=3)
    # CLSTF2 <- generate_data_scenario(k=3,n=n2, m=m, setting=2, "B",  Q=3)
    # CLSTF3 <- generate_data_scenario(k=3,n=n3, m=m, setting=3,  "B", Q=3)
    # # 
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
             || (min(as.numeric(tolcat)) == 1 && W[, clstr_2_idx][1] == refcat && count_iter < 100)
             || (length(unique(W[,clstr_2_idx]))<length(unique(c(W)))&& count_iter < 100))
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
      W = W
      
    }
    
    
    ###########
    
    #ESTIMATION
    if(wood_multi=="N"){
      categFD_est <- estimate_categFD(X, tt=tt)
    }
    
    #wood_multi="Y"
    if(wood_multi=="Y"){
      categFD_est <- estimate_categFD_W(W, tt=tt)
    }
    
    
    # recover latent process
    Z1_est=categFD_est$Z1_est
    Z2_est=categFD_est$Z2_est
    
    #recover =probabilities for each category
    p1_est<-categFD_est$p1_est
    p2_est<-categFD_est$p2_est
    p3_est<-categFD_est$p3_est
    
    
    # ########plot true and est z1, z2
   #  par(mfrow=c(2,2))
   # #matplot(seq(0.0001,1,length=m),Z1[,1:n1],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,xlim=c(0,1),ylim=c(-8,8))
   # matplot(seq(0.0001,1,length=m),Z1[,1:n1],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,xlim=c(0,1),ylim=c(-40,15))
   # matlines(seq(0.0001,1,length=m),Z1[,(n1+1):(n1+n2)],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
   # matlines(seq(0.0001,1,length=m),Z1[,(n1+n2+1):n],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
   # 
   # #matplot(seq(0.0001,1,length=m),Z2[,1:n1],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,xlim=c(0,1),ylim=c(-8,20))
   # matplot(seq(0.0001,1,length=m),Z2[,1:n1],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,xlim=c(0,1),ylim=c(-40,15))
   # matlines(seq(0.0001,1,length=m),Z2[,(n1+1):(n1+n2)],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
   # matlines(seq(0.0001,1,length=m),Z2[,(n1+n2+1):n],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
   # 
   #  # #plot first truelatent curves p_2 by clusters
   #  matplot(seq(0.0001,1,length=m),Z1_est[,1:n1],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(-8,8))
   #  matlines(seq(0.0001,1,length=m),Z1_est[,(n1+1):(n1+n2)],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
   #  matlines(seq(0.0001,1,length=m),Z1_est[,(n1+n2+1):n],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
   # 
   # 
   #  #plot first truelatent curves p_2 by clusters
   #  matplot(seq(0.0001,1,length=m),Z2_est[,1:n1],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(-40,15))
   #  matlines(seq(0.0001,1,length=m),Z2_est[,(n1+1):(n1+n2)],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
   #  matlines(seq(0.0001,1,length=m),Z2_est[,(n1+n2+1):n],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)

    # # 
    # par(mfrow=c(2,3))
    # matplot(seq(0.0001,1,length=m),p1[,1:n1],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,xlim=c(0,1),ylim=c(0,1))
    # matlines(seq(0.0001,1,length=m),p1[,(n1+1):(n1+n2)],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
    # matlines(seq(0.0001,1,length=m),p1[,(n1+n2+1):n],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
    # 
    # matplot(seq(0.0001,1,length=m),p2[,1:n1],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,xlim=c(0,1),ylim=c(0,1))
    # matlines(seq(0.0001,1,length=m),p2[,(n1+1):(n1+n2)],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
    # matlines(seq(0.0001,1,length=m),p2[,(n1+n2+1):n],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
    # 
    # 
    # matplot(seq(0.0001,1,length=m),p3[,1:n1],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,xlim=c(0,1),ylim=c(0,1))
    # matlines(seq(0.0001,1,length=m),p3[,(n1+1):(n1+n2)],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
    # matlines(seq(0.0001,1,length=m),p3[,(n1+n2+1):n],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
    # 
    # #plot first truelatent curves p_2 by clusters
    # matplot(seq(0.0001,1,length=m),p1_est[,1:n1],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(0,1))
    # matlines(seq(0.0001,1,length=m),p1_est[,(n1+1):(n1+n2)],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
    # matlines(seq(0.0001,1,length=m),p1_est[,(n1+n2+1):n],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
    # 
    # 
    # #plot first truelatent curves p_2 by clusters
    # matplot(seq(0.0001,1,length=m),p2_est[,1:n1],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(0,1))
    # matlines(seq(0.0001,1,length=m),p2_est[,(n1+1):(n1+n2)],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
    # matlines(seq(0.0001,1,length=m),p2_est[,(n1+n2+1):n],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
    # 
    # 
    # matplot(seq(0.0001,1,length=m),p3_est[,1:n1],col="red",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(0,1))
    # matlines(seq(0.0001,1,length=m),p3_est[,(n1+1):(n1+n2)],col="green",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=3, lwd=3)
    # matlines(seq(0.0001,1,length=m),p3_est[,(n1+n2+1):n],col="blue",type = "l",xlab="Time",ylab="value",cex.lab=1.5, cex.axis=2,lty=5, lwd=3)
    
    #####
    
    
    
    
    # expZ1 <- exp(Z1_est)
    # expZ2 <- exp(Z1_est)
    # denom <- 1+ expZ1+  expZ2
    # p1_est <- expZ1/ denom
    # p2_est <- expZ2/ denom
    # p3_est <- 1-p1_est-p2_est
    
    # evaluate performance Z and P
    # rmse1_temp <- c(by(sqrt(colMeans( (Z1-Z1_est)^2 ) ), true_cluster, mean))
    # rmse2_temp <- c(by(sqrt(colMeans((Z2-Z2_est)^2 ) ), true_cluster, mean))
    # rmse[ii, ,] <- rbind(rmse1_temp,rmse2_temp )
    # 
    # error.p1<- sqrt(colMeans((sqrt(p1_est)-sqrt(p1))^2))/sqrt(2)
    # error.p2<- sqrt(colMeans((sqrt(p2_est)-sqrt(p2))^2))/sqrt(2)
    # error.p3<- sqrt(colMeans((sqrt(p3_est)-sqrt(p3))^2))/sqrt(2)
    # 
    # 
    # hellinger[ii, ,] <-  rbind( c(by(error.p1, true_cluster, mean)), 
    #                             c(by(error.p2, true_cluster, mean)),
    #                             c(by(error.p3, true_cluster, mean)))
    
    
    # Record clustering performance. USE only Z
    #**** extract_scores_UNIVFPCA is the FASTEST
    mfpca_true <-extract_scores_UNIVFPCA(mZ1=Z1, mZ2=Z2,  tt=tt , PVE=0.95)
    #plot(  scores_true$scores[, 1:2])
    mfpca_est <- extract_scores_UNIVFPCA(mZ1=Z1_est, mZ2=Z2_est, tt=tt , PVE=0.95)
    
    ############
    #plot the differences
    # for (i in 1:500){
    #   plot(Z2[,i]-Z2_est[,i],
    #        main=paste0("True and Est differences for subject",i))
    # }
    # table(W[,83])
    # 
    # 1   2   3 
    # 87   1 212 
    # > W[,83]
    # [1] 3 3 1 3 3 3 1 3 3 1 3 3 3 3 3 3 3 1 3 3 3 3 3 3 3 3 3 3
    # [29] 3 3 3 3 3 3 3 3 3 1 3 3 3 3 3 3 3 3 1 3 3 3 3 3 1 3 3 3
    # [57] 3 3 3 3 1 3 3 3 1 3 3 3 1 3 3 3 3 3 3 1 3 3 3 3 3 3 3 3
    # [85] 3 3 3 3 3 3 3 3 3 3 3 3 1 3 3 1 3 3 3 3 3 3 3 3 3 3 3 3
    # [113] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 3 3 3 3 3 3 3 3 3
    # [141] 3 1 3 1 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 3 3 1 3
    # [169] 1 3 1 3 3 3 3 3 3 3 3 3 3 1 3 1 3 3 1 1 1 1 3 1 3 1 1 3
    # [197] 3 3 3 1 3 3 3 1 3 3 3 1 1 3 1 1 3 1 1 1 3 3 1 1 1 3 3 1
    # [225] 1 1 3 3 3 3 1 3 3 1 1 3 3 1 1 3 3 3 1 1 3 1 1 1 3 3 1 1
    # [253] 1 3 3 1 3 3 1 3 3 1 3 1 1 1 1 3 3 1 1 1 1 3 1 3 1 1 3 1
    # [281] 1 1 1 1 1 1 1 1 1 3 3 3 1 3 1 1 1 1 2 1
    # > table(W[,84])
    # 
    # 1   2   3 
    # 106   2 192 
    # > W[,84]
    # [1] 3 3 3 3 3 3 1 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
    # [29] 1 3 3 3 3 1 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 3 3 3 3 3 3 3
    # [57] 3 3 3 3 3 3 3 3 1 3 3 3 3 3 3 1 3 3 3 3 3 1 3 1 3 3 3 3
    # [85] 3 3 3 1 3 1 3 3 3 1 3 3 1 3 3 1 3 3 3 3 3 3 1 1 3 3 1 3
    # [113] 3 1 1 3 3 3 3 3 1 3 3 3 3 3 3 3 3 3 3 3 1 1 3 3 3 3 3 3
    # [141] 3 1 3 3 1 3 1 3 1 3 3 3 3 3 3 3 3 3 1 3 3 1 3 3 3 1 3 3
    # [169] 3 1 1 3 1 1 1 1 1 1 1 3 3 1 3 3 3 3 3 1 1 1 3 3 1 3 1 1
    # [197] 3 1 3 3 1 1 3 1 1 3 3 3 3 1 3 3 3 1 3 3 1 1 3 3 1 1 1 1
    # [225] 3 3 3 1 3 3 1 3 1 3 3 3 1 3 1 3 3 1 1 3 1 1 1 3 1 1 3 3
    # [253] 1 1 1 1 3 1 1 1 1 1 1 1 1 1 1 3 1 1 1 3 1 1 3 1 1 1 3 1
    # [281] 1 1 3 1 3 1 1 1 3 1 1 1 1 3 1 1 2 2 1 1
    # 1/300
    # [1] 0.003333333
    #2/300
    #[1] 0.006666667
    # table(W[,92])
    # 
    # 1   2   3 
    # 112   5 183 
    # > 5/300
    # [1] 0.01666667
    # plot(Z2[,83],main="True Z2 subject 83")
    # plot(Z2_est[,83],main="Est Z2 subject 83")
    # plot(Z2[,84],main="True Z2 subject 84")
    # plot(Z2_est[,84],main="Est Z2 subject 84")
    # plot(Z2[,93],main="True Z2 subject 93")
    # plot(Z2_est[,93],main="Est Z2 subject 93")
    
    #t=2000
    ##table(W[,33])
    # 1   2   3 
    # 808 385 807 
    # 385/2000
    # [1] 0.1925
    ############
    
    #KMEANS
    true_kmeans_temp <- kmeans_cluster(data=mfpca_true$scores)$label
    est_kmeans_temp <- kmeans_cluster(data=mfpca_est$scores)$label
    
    #####################
    true_kmeans_temp_nclust <- kmeans_cluster(data=mfpca_true$scores)$nclust
    est_kmeans_temp_nclust <- kmeans_cluster(data=mfpca_est$scores)$nclust
    
    ######################
    
    #FADP
    true_fadp_temp <- fadp_cluster(mZ1=Z1, mZ2=Z2, tt=tt)$label
    est_fadp_temp <- fadp_cluster(mZ1=Z1_est, mZ2=Z2_est, tt=tt)$label
    
    ##################
    true_fadp_temp_nclust <- fadp_cluster(mZ1=Z1, mZ2=Z2, tt=tt)$nclust
    est_fadp_temp_nclust <- fadp_cluster(mZ1=Z1_est, mZ2=Z2_est, tt=tt)$nclust
    
    ###################
    
    #dbscan
    true_dbscan_temp <- dbscan_cluster(data=mfpca_true$scores,1)$label
    est_dbscan_temp <- dbscan_cluster(data=mfpca_est$scores,1)$label
    ####################
    true_dbscan_temp_nclust <- dbscan_cluster(data=mfpca_true$scores,1)$nclust
    est_dbscan_temp_nclust <- dbscan_cluster(data=mfpca_est$scores,1)$nclust
    ####################
    
    #record results
    true_dbscan<- cbind(true_dbscan, true_dbscan_temp)
    est_dbscan <- cbind(est_dbscan, est_dbscan_temp)
    
    
    true_kmeans<- cbind(true_kmeans, true_kmeans_temp)
    est_kmeans <- cbind(est_kmeans,  est_kmeans_temp)
    
    true_fadp<- cbind(true_fadp, true_fadp_temp)
    est_fadp <- cbind(est_fadp, est_fadp_temp)
    
    clust_record[ii,]= c(true_kmeans_temp_nclust,est_kmeans_temp_nclust,
                    true_fadp_temp_nclust, est_fadp_temp_nclust,
                    true_dbscan_temp_nclust, est_dbscan_temp_nclust)
    
    print(ii)
  }
  
  
  
  # Assess accuracy in the simulation study
  # apply(rmse, c(2,3), mean)
  # apply(hellinger, c(2,3), mean)

  
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
  
  
  
 return(list("cluster_table_true"=cluster_table_true,"cluster_table_est"=cluster_table_est,
             "cluster_table_est_se"=cluster_table_est_se,"cluster_record"=clust_record))
  
  }


##########scenarioA
#cluster_simulation=function(n,m,scenario,mc_sims)
n100t300A=cluster_simulation(100,300,"A",2,"Y")
n100t300Am=cluster_simulation(100,300,"A",10,"N")

n100t750A=cluster_simulation(100,750,"A",20)
n100t2000A=cluster_simulation(100,2000,"A",20)


n500t300A=cluster_simulation(500,300,"A",20)
n500t750A=cluster_simulation(500,750,"A",20)
n500t2000A=cluster_simulation(500,2000,"A",20)


n1000t300A=cluster_simulation(1000,300,"A",20)
n1000t750A=cluster_simulation(1000,750,"A",20)
n1000t2000A=cluster_simulation(1000,2000,"A",20)


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

save(true_tableA,est_tableA,est_tableA_se,file="A_clustering.RData")

# load("staicu_Anew20_final.RData")
# library(xtable)
# xtable(est_tableA)


#####################################
###scenario B
##########scenarioA
#cluster_simulation=function(n,m,scenario,mc_sims)
n100t300B=cluster_simulation(100,300,"B",3,wood_multi="Y")
n100t300Bm=cluster_simulation(100,300,"B",10,wood_multi="N")

n100t750B=cluster_simulation(100,750,"B",20)
n100t2000B=cluster_simulation(100,2000,"B",20)


n500t300B=cluster_simulation(500,300,"B",20)
n500t750B=cluster_simulation(500,750,"B",20)
n500t2000B=cluster_simulation(500,2000,"B",20)


n1000t300B=cluster_simulation(1000,300,"B",20)
n1000t750B=cluster_simulation(1000,750,"B",20)
n1000t2000B=cluster_simulation(1000,2000,"B",20)


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

save(true_tableB,est_tableB,est_table_se,file="B_clustering.RData")

# load("staicu_Anew20_final.RData")
# library(xtable)
# xtable(est_tableA)









# \begin{tabular}{rrrrrrrrrr}
# \hline
# & fadp RI & fadp ARI & fadp cpn & kmeans RI & kmeans ARI & kmeans cpn & dbscan RI & dbscan ARI & dbscan cpn \\ 
# \hline
# n100t300 & 0.79 & 0.55 & 0.00 & 0.89 & 0.77 & 0.15 & 0.88 & 0.72 & 0.33 \\ 
# n100t750 & 0.88 & 0.77 & 0.18 & 0.98 & 0.95 & 0.05 & 0.96 & 0.92 & 0.88 \\ 
# n100t2000 & 0.86 & 0.72 & 0.00 & 0.99 & 0.97 & 0.00 & 0.97 & 0.94 & 0.98 \\ 
# n500t300 & 0.96 & 0.92 & 0.00 & 0.94 & 0.87 & 0.15 & 0.86 & 0.66 & 0.24 \\ 
# n500t750 & 0.96 & 0.91 & 0.00 & 0.98 & 0.97 & 0.00 & 0.97 & 0.93 & 0.36 \\ 
# n500t2000 & 0.99 & 0.97 & 0.00 & 0.99 & 0.97 & 0.00 & 0.97 & 0.93 & 0.29 \\ 
# n1000t300 & 0.98 & 0.96 & 0.00 & 0.94 & 0.89 & 0.00 & 0.90 & 0.76 & 0.20 \\ 
# n1000t750 & 0.99 & 0.97 & 0.00 & 0.99 & 0.97 & 0.00 & 0.97 & 0.94 & 0.10 \\ 
# n1000t2000 & 0.99 & 0.97 & 0.00 & 0.99 & 0.97 & 0.00 & 0.97 & 0.94 & 0.11 \\ 
# \hline
# \end{tabular}
# \end{table}




# load("staicu_Anew20.RData")
# library(xtable)
# xtable(est_tableA)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrrr}
# \hline
# & fadp RI & fadp ARI & fadp cpn & kmeans RI & kmeans ARI & kmeans cpn & dbscan RI & dbscan ARI & dbscan cpn \\ 
# \hline
# n100t300 & 0.81 & 0.60 & 0.00 & 0.91 & 0.81 & 0.10 & 0.87 & 0.70 & 0.80 \\ 
# n100t750 & 0.80 & 0.59 & 0.00 & 0.98 & 0.95 & 0.13 & 0.98 & 0.95 & 0.97 \\ 
# n100t2000 & 0.87 & 0.75 & 0.05 & 0.99 & 0.97 & 0.00 & 0.98 & 0.96 & 1.00 \\ 
# n500t300 & 0.92 & 0.84 & 0.00 & 0.96 & 0.92 & 0.00 & 0.75 & 0.37 & 0.47 \\ 
# n500t750 & 0.94 & 0.86 & 0.00 & 0.99 & 0.97 & 0.00 & 0.96 & 0.91 & 0.73 \\ 
# n500t2000 & 0.96 & 0.92 & 0.00 & 0.99 & 0.97 & 0.00 & 0.97 & 0.94 & 0.65 \\ 
# n1000t300 & 0.98 & 0.96 & 0.00 & 0.96 & 0.92 & 0.00 & 0.73 & 0.33 & 0.33 \\ 
# n1000t750 & 0.99 & 0.97 & 0.00 & 0.99 & 0.97 & 0.00 & 0.99 & 0.98 & 0.28 \\ 
# n1000t2000 & 0.99 & 0.97 & 0.00 & 0.99 & 0.97 & 0.00 & 0.98 & 0.96 & 0.25 \\ 
# \hline
# \end{tabular}
# \end{table}


# 
# true_tableA
# fadp RI  fadp ARI fadp cpn kmeans RI kmeans ARI kmeans cpn dbscan RI dbscan ARI dbscan cpn
# n100t300   0.8292727 0.6644215     0.05 0.9866667  0.9718616          0 0.9855455  0.9700109  1.0000000
# n100t750   0.8292727 0.6644215     0.05 0.9866667  0.9718616          0 0.9855455  0.9700109  1.0000000
# n100t2000  0.8292727 0.6644215     0.05 0.9866667  0.9718616          0 0.9853434  0.9695902  1.0000000
# n500t300   0.8371150 0.6714885     0.01 0.9867735  0.9720023          0 0.9773098  0.9528257  0.6700000
# n500t750   0.8371150 0.6714885     0.01 0.9867735  0.9720023          0 0.9769511  0.9520922  0.6700000
# n500t2000  0.8371150 0.6714885     0.01 0.9867735  0.9720023          0 0.9768008  0.9517744  0.6700000
# n1000t300  0.7705015 0.5312544     0.00 0.9867868  0.9720196          0 0.9712855  0.9403743  0.5850000
# n1000t750  0.7705015 0.5312544     0.00 0.9867868  0.9720196          0 0.9710860  0.9399705  0.5683333
# n1000t2000 0.7705015 0.5312544     0.00 0.9867868  0.9720196          0 0.9711831  0.9401685  0.5683333



# 
# load("staicu_Anew.RData")
# library(xtable)
# xtable(est_tableA)

#load("staicu_A.RData")
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrrr}
# \hline
# & fadp RI & fadp ARI & fadp cpn & kmeans RI & kmeans ARI & kmeans cpn & dbscan RI & dbscan ARI & dbscan cpn \\ 
# \hline
# n100t300 & 0.81 & 0.60 & 0.00 & 0.91 & 0.81 & 0.10 & 0.87 & 0.70 & 0.80 \\ 
# n100t750 & 0.80 & 0.59 & 0.00 & 0.98 & 0.95 & 0.13 & 0.98 & 0.95 & 0.97 \\ 
# n100t2000 & 0.87 & 0.75 & 0.05 & 0.99 & 0.97 & 0.00 & 0.98 & 0.96 & 1.00 \\ 
# n500t300 & 0.92 & 0.84 & 0.00 & 0.96 & 0.92 & 0.00 & 0.75 & 0.37 & 0.47 \\ 
# n500t750 & 0.94 & 0.86 & 0.00 & 0.99 & 0.97 & 0.00 & 0.96 & 0.91 & 0.73 \\ 
# n500t2000 & 0.96 & 0.92 & 0.00 & 0.99 & 0.97 & 0.00 & 0.97 & 0.94 & 0.65 \\ 
# n1000t300 & 0.98 & 0.96 & 0.00 & 0.96 & 0.92 & 0.00 & 0.73 & 0.33 & 0.33 \\ 
# n1000t750 & 0.99 & 0.97 & 0.00 & 0.99 & 0.97 & 0.00 & 0.99 & 0.98 & 0.28 \\ 
# n1000t2000 & 0.99 & 0.97 & 0.00 & 0.99 & 0.97 & 0.00 & 0.98 & 0.96 & 0.25 \\ 
# \hline
# \end{tabular}
# \end{table}

#standard error
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrrr}
# \hline
# & fadp RI & fadp ARI & fadp cpn & kmeans RI & kmeans ARI & kmeans cpn & dbscan RI & dbscan ARI & dbscan cpn \\ 
# \hline
# n100t300 & 0.04 & 0.09 & 0.00 & 0.02 & 0.05 & 0.06 & 0.03 & 0.08 & 0.06 \\ 
# n100t750 & 0.05 & 0.10 & 0.00 & 0.01 & 0.02 & 0.07 & 0.00 & 0.01 & 0.02 \\ 
# n100t2000 & 0.04 & 0.07 & 0.05 & 0.00 & 0.00 & 0.00 & 0.00 & 0.01 & 0.00 \\ 
# n500t300 & 0.03 & 0.07 & 0.00 & 0.02 & 0.05 & 0.00 & 0.04 & 0.10 & 0.07 \\ 
# n500t750 & 0.03 & 0.07 & 0.00 & 0.00 & 0.00 & 0.00 & 0.02 & 0.05 & 0.08 \\ 
# n500t2000 & 0.02 & 0.04 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.01 & 0.09 \\ 
# n1000t300 & 0.00 & 0.00 & 0.00 & 0.02 & 0.05 & 0.00 & 0.04 & 0.10 & 0.02 \\ 
# n1000t750 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.03 \\ 
# n1000t2000 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.03 \\ 
# \hline
# \end{tabular}
# \end{table}

###############scenario B
#cluster_simulation=function(n,m,scenario,mc_sims)
# n100t300B=cluster_simulation(100,300,"B",20)
# n100t750B=cluster_simulation(100,750,"B",20)
# n100t2000B=cluster_simulation(100,2000,"B",20)
# 
# 
# n500t300B=cluster_simulation(500,300,"B",20)
# n500t750B=cluster_simulation(500,750,"B",20)
# n500t2000B=cluster_simulation(500,2000,"B",20)
# 
# 
# n1000t300B=cluster_simulation(1000,300,"B",20)
# n1000t750B=cluster_simulation(1000,750,"B",20)
# n1000t2000B=cluster_simulation(1000,2000,"B",20)
# 
# 
# true_tableB=rbind(n100t300B$cluster_table_true,n100t750B$cluster_table_true,n100t2000B$cluster_table_true,
#                   n500t300B$cluster_table_true,n500t750B$cluster_table_true,n500t2000B$cluster_table_true,
#                   n1000t300B$cluster_table_true,n1000t750B$cluster_table_true,n1000t2000B$cluster_table_true)
# rownames(true_tableB)=c("n100t300","n100t750","n100t2000",
#                         "n500t300","n500t750","n500t2000",
#                         "n1000t300","n1000t750","n1000t2000")
# 
# est_tableB=rbind(n100t300B$cluster_table_est,n100t750B$cluster_table_est,n100t2000B$cluster_table_est,
#                  n500t300B$cluster_table_est,n500t750B$cluster_table_est,n500t2000B$cluster_table_est,
#                  n1000t300B$cluster_table_est,n1000t750B$cluster_table_est,n1000t2000B$cluster_table_est)
# rownames(est_tableB)=c("n100t300","n100t750","n100t2000",
#                        "n500t300","n500t750","n500t2000",
#                        "n1000t300","n1000t750","n1000t2000")
# 
# 
# est_tableB_se=rbind(n100t300B$cluster_table_est_se,n100t750B$cluster_table_est_se,n100t2000B$cluster_table_est_se,
#                     n500t300B$cluster_table_est_se,n500t750B$cluster_table_est_se,n500t2000B$cluster_table_est_se,
#                     n1000t300B$cluster_table_est_se,n1000t750B$cluster_table_est_se,n1000t2000B$cluster_table_est_se)
# rownames(est_tableB_se)=c("n100t300","n100t750","n100t2000",
#                           "n500t300","n500t750","n500t2000",
#                           "n1000t300","n1000t750","n1000t2000")
# 
# save(true_tableB,est_tableB,est_tableB_se,file="staicu_B.RData")
# 
load("staicu_B.RData")
# library(xtable)
# xtable(true_tableB)
# xtable(true_tableB)
# % latex table generated in R 4.2.1 by xtable 1.8-4 package
# % Thu Aug 17 12:29:07 2023
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrrr}
# \hline
# & fadp RI & fadp ARI & fadp cpn & kmeans RI & kmeans ARI & kmeans cpn & dbscan RI & dbscan ARI & dbscan cpn \\ 
# \hline
# n100t300 & 0.88 & 0.76 & 0.00 & 0.88 & 0.76 & 0.00 & 0.95 & 0.89 & 0.08 \\ 
# n100t750 & 0.88 & 0.76 & 0.00 & 0.88 & 0.76 & 0.00 & 0.95 & 0.88 & 0.08 \\ 
# n100t2000 & 0.88 & 0.76 & 0.00 & 0.88 & 0.76 & 0.00 & 0.95 & 0.89 & 0.08 \\ 
# n500t300 & 0.88 & 0.76 & 0.00 & 0.88 & 0.76 & 0.00 & 0.96 & 0.91 & 0.10 \\ 
# n500t750 & 0.88 & 0.76 & 0.00 & 0.88 & 0.76 & 0.00 & 0.96 & 0.91 & 0.10 \\ 
# n500t2000 & 0.88 & 0.76 & 0.00 & 0.88 & 0.76 & 0.00 & 0.96 & 0.91 & 0.10 \\ 
# n1000t300 & 0.88 & 0.76 & 0.00 & 0.88 & 0.76 & 0.00 & 0.96 & 0.91 & 0.08 \\ 
# n1000t750 & 0.88 & 0.76 & 0.00 & 0.88 & 0.76 & 0.00 & 0.96 & 0.91 & 0.08 \\ 
# n1000t2000 & 0.88 & 0.76 & 0.00 & 0.88 & 0.76 & 0.00 & 0.96 & 0.91 & 0.07 \\ 
# \hline
# \end{tabular}
# \end{table}

# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrrr}
# \hline
# & fadp RI & fadp ARI & fadp cpn & kmeans RI & kmeans ARI & kmeans cpn & dbscan RI & dbscan ARI & dbscan cpn \\ 
# \hline
# n100t300 & 0.76 & 0.55 & 0.00 & 0.68 & 0.45 & 0.09 & 0.86 & 0.73 & 0.21 \\ 
# n100t750 & 0.88 & 0.76 & 0.00 & 0.88 & 0.75 & 0.00 & 0.93 & 0.85 & 0.16 \\ 
# n100t2000 & 0.88 & 0.76 & 0.00 & 0.88 & 0.76 & 0.00 & 0.98 & 0.95 & 0.16 \\ 
# n500t300 & 0.76 & 0.56 & 0.00 & 0.65 & 0.39 & 0.13 & 0.76 & 0.58 & 0.08 \\ 
# n500t750 & 0.88 & 0.76 & 0.00 & 0.88 & 0.76 & 0.00 & 0.93 & 0.86 & 0.11 \\ 
# n500t2000 & 0.88 & 0.76 & 0.00 & 0.88 & 0.76 & 0.00 & 0.98 & 0.95 & 0.11 \\ 
# n1000t300 & 0.72 & 0.48 & 0.05 & 0.42 & 0.05 & 0.00 & 0.64 & 0.39 & 0.05 \\ 
# n1000t750 & 0.86 & 0.73 & 0.00 & 0.90 & 0.80 & 0.00 & 0.94 & 0.87 & 0.07 \\ 
# n1000t2000 & 0.88 & 0.76 & 0.00 & 0.88 & 0.76 & 0.00 & 0.98 & 0.96 & 0.08 \\ 
# \hline
# \end{tabular}
# \end{table}

#standard error
# \begin{tabular}{rrrrrrrrrr}
# \hline
# & fadp RI & fadp ARI & fadp cpn & kmeans RI & kmeans ARI & kmeans cpn & dbscan RI & dbscan ARI & dbscan cpn \\ 
# \hline
# n100t300 & 0.03 & 0.05 & 0.00 & 0.05 & 0.08 & 0.06 & 0.02 & 0.04 & 0.02 \\ 
# n100t750 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.01 & 0.02 & 0.02 \\ 
# n100t2000 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.01 & 0.02 \\ 
# n500t300 & 0.03 & 0.04 & 0.00 & 0.05 & 0.08 & 0.07 & 0.05 & 0.08 & 0.01 \\ 
# n500t750 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.01 & 0.02 & 0.01 \\ 
# n500t2000 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.01 & 0.01 \\ 
# n1000t300 & 0.01 & 0.03 & 0.05 & 0.03 & 0.04 & 0.00 & 0.06 & 0.09 & 0.01 \\ 
# n1000t750 & 0.01 & 0.02 & 0.00 & 0.01 & 0.02 & 0.00 & 0.01 & 0.02 & 0.01 \\ 
# n1000t2000 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.01 \\ 
# \hline
# \end{tabular}
# \end{table}


# par(mfrow=c(1,2))
# plot(1:mc_sims, results$true_kmeans_ri, type="p", ylim=c(0,1), ylab="kmeans", xlab="MC simulations")
# points(1:mc_sims, results$true_kmeans_ari, pch=2)
# 
# points(1:mc_sims, results$est_kmeans_ri, col="blue")
# points(1:mc_sims, results$est_kmeans_ari, pch=2, col="blue")
# 
# plot(1:mc_sims, results$true_fadp_ri, type="p", ylim=c(0,1), ylab="fadp", xlab="MC simulations")
# points(1:mc_sims, results$true_fadp_ari, pch=2)
# 
# points(1:mc_sims, results$est_fadp_ri, col="blue")
# points(1:mc_sims, results$est_fadp_ari, pch=2, col="blue")



