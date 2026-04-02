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
# Purpose: Adding weekend effects in modeling the probability
# Author:  Xiaoxia Champon
# Date: 11/06/2024
##############################################################






#' Runs gam on the given data and returns results
#' Fits a binomial to describe the given in_x
#' @param timestamps01 current time value
#' @param in_x binary series
#' @param family_choice "binomial" or "probit"
#' @return list: fit values and linear predictors both with length of time_series length
RunGam_Day <- function(timestamps01, weekend_vector, in_x, family_choice, basis_size, method)
{
  basis_size_rev <- max(min(round(min(sum(in_x), sum(1-in_x))/2), basis_size ), 5)
  
  if(family_choice=="binomial")
  {
    family_object <- "binomial"
  }else if(family_choice=="probit")
  {
    family_object <- binomial(link="probit")
  }
  
  fit_binom <- gam(in_x~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                   family=family_object, method = method,
                   control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                   optimizer=c("outer","bfgs"))
  prob <- fit_binom$fitted.values
  
  if(family_choice=="binomial")
  {
    prob_linpred <- fit_binom$linear.predictors
  }else if(family_choice=="probit"){
    z_probit <- fit_binom$linear.predictors
    
    prob_probit <- pnorm(z_probit)
    zprobit_reciprocal <- 1/ prob_probit
    prob_linpred <- -log(zprobit_reciprocal-1)
  }
  
  return(list(prob=prob, linpred=prob_linpred))
}

#testing code
# library(mgcv)
# #RunGam_Day <- function(timestamps01, weekend_vector, in_x, family_choice, basis_size, method)
# timestamps01 <- seq(0,1, length = 1000)
# weekend_vector <- c(rep(c(0,0,0,0,0,1,1),200))[1:1000]
# in_x <- rbinom(1000,1,0.5)
# family_choice <- "binomial"
# basis_size <- 25
# method <- "ML"
# test_data <- RunGam_Day(timestamps01, weekend_vector, in_x, family_choice, basis_size, method)

#' Fit multinomial from W categorical data to extra smooth latent curves
#' Fits a binomial to describe the given in_x
#' @param timestamps01 current time value
#' @param W : 2D categorical data matrix, n * t dimension, n is the number of subjects, t is the time
#' @return list: fit values and linear predictors both with length of time_series length, Z is t*n, p: t*n
EstimateCategFuncData_multinormial <- function(timestamps01, W, basis_size=25, method="ML")
{
  
  # num_indv<- ncol(W)
  # timeseries_length <-nrow(W)
  
  num_indv<- nrow(W)
  timeseries_length <-ncol(W)
  category_count <- length(unique(c(W)))
  weekend_vector <- as.factor(c(rep(c(rep(0,480),rep(1,192)),4))[1:timeseries_length])
  
  Z<-NULL
  prob<-array(0, c(num_indv, timeseries_length , category_count))
  weekend_vector_coef <- matrix(0, num_indv, category_count-1)
  for (i in 1:num_indv){
    print(i)
    if ( length(table(W[i,])) == category_count) {
      
      fit_binom<-gam(list(W[i,]-1~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector,
                          ~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector
      ),
      family=multinom(K=category_count-1), method = method,
      control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
      optimizer=c("outer","bfgs")) 
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
      
      weekend_vector_coef[i, ] <- fit_binom$coefficients[c(category_count-1,basis_size+category_count-1)]
    } else {
      
      W[i,][W[i,]==3] <- 2
      
      
      basis_size_rev <- max(min(round(min(unname(table(W[i,])[2]), sum(1-unname(table(W[i,])[2])))/2), basis_size ), 5)
      fit_binom <- gam(W[i,]-1~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
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
      #####
      # g_mul_2 <- g_design[,c(1,category_count:basis_size)+basis_size]
      # coef_fit_2 <- fit_binom$coefficients[c(1,category_count:basis_size)+basis_size]
      #z2 <- g_mul_2 %*% as.matrix(coef_fit_2,ncol=1)
      z2 <- rep(0,timeseries_length)
      
      weekend_vector_coef[i, ] <- c(fit_binom$coefficients[category_count-1],0)
      ##########################
      
      
    } 
    
    
    # z1<- fit_binom$linear.predictors[,1]
    # z2<- fit_binom$linear.predictors[,2]
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
  
  return(list(Z1_est=Z[1:timeseries_length ,], Z2_est=Z[1:timeseries_length + timeseries_length ,], 
              p1_est=t(prob[,,1]), p2_est=t(prob[,,2]), p3_est=t(prob[,,3]) ,
              weekend_vector_coef = weekend_vector_coef))
}

#' Fit multinomial from W categorical data to extra smooth latent curves
#' Fits a binomial to describe the given in_x
#' @param timestamps01 current time value
#' @param W : 2D categorical data matrix, n * t dimension, n is the number of subjects, t is the time
#' @return list: fit values and linear predictors both with length of time_series length, Z is t*n, p: t*n
EstimateCategFuncData_multinormial_weekend <- function(timestamps01, W, basis_size=25, method="ML")
{
  
  
  num_indv<- nrow(W)
  timeseries_length <-ncol(W)
  category_count <- length(unique(c(W)))
  weekend_vector <- as.factor(c(rep(c(rep(0,480),rep(1,192)),4))[1:timeseries_length])
  
  Z<-NULL
  prob<-array(0, c(num_indv, timeseries_length , category_count))
  weekend_vector_coef <- matrix(0, num_indv, category_count-1)
  for (i in 1:num_indv){
    print(i)
    if ( length(table(W[i,])) == category_count) {
      
      fit_binom<-gam(list(W[i,]-1~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector,
                          ~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector
      ),
      family=multinom(K=category_count-1), method = method,
      control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
      optimizer=c("outer","bfgs")) 
      #####
      print("Dino")
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
      
      weekend_vector_coef[i, ] <- fit_binom$coefficients[c(category_count-1,basis_size+category_count-1)]
    } else {
      if (names(table(W[i,]))[2]=="3"){
        W[i,][W[i,]==3] <- 2
        basis_size_rev <- max(min(round(min(unname(table(W[i,])[2]), sum(1-unname(table(W[i,])[2])))/2), basis_size ), 5)
        fit_binom <- gam(W[i,]-1~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                         family = "binomial", method = method,
                         control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                         optimizer=c("outer","bfgs"))
        print("GODAK ADERI")
        ######################
        ####to find design matrix
        g_design <- predict(fit_binom,type = "lpmatrix")
        g_mul <- g_design[,c(1,category_count:basis_size_rev)]
        coef_fit <- fit_binom$coefficients[c(1,category_count:basis_size_rev)]
        #extract z
        z2 <- g_mul %*% as.matrix(coef_fit,ncol=1)
        #####
        z1 <- rep(0,timeseries_length)
        
        weekend_vector_coef[i, ] <- c(0, fit_binom$coefficients[category_count-1])
        ##########################
      }else {
        basis_size_rev <- max(min(round(min(unname(table(W[i,])[2]), sum(1-unname(table(W[i,])[2])))/2), basis_size ), 5)
        fit_binom <- gam(W[i,]-1~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                         family = "binomial", method = method,
                         control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                         optimizer=c("outer","bfgs"))
        print("MATA OYA NATHWA PALUI")
        ######################
        ####to find design matrix
        g_design <- predict(fit_binom,type = "lpmatrix")
        g_mul <- g_design[,c(1,category_count:basis_size_rev)]
        coef_fit <- fit_binom$coefficients[c(1,category_count:basis_size_rev)]
        #extract z
        z1 <- g_mul %*% as.matrix(coef_fit,ncol=1)
        z2 <- rep(0,timeseries_length)
        weekend_vector_coef[i, ] <- c(fit_binom$coefficients[category_count-1],0)
        ##########################
      }
    } 
    #2t*n matrix
    Z<- cbind(Z, c(z1,z2))
    ##find probability
    Z_cbind=cbind(z1,z2)
    exp_z=exp(Z_cbind)
    denominator_p=1+exp_z[,1]+exp_z[,2]
    p1 <- exp_z[,1]/denominator_p
    p2 <- exp_z[,2]/denominator_p
    p3=1/denominator_p
    #3D matrix t*n*category 
    prob[i,,] <- cbind(p1, p2, p3)
    
    #weekend_vector_coef n*2
  }
  
  return(list(Z1_est=Z[1:timeseries_length ,], Z2_est=Z[1:timeseries_length + timeseries_length ,], 
              p1_est=t(prob[,,1]), p2_est=t(prob[,,2]), p3_est=t(prob[,,3]) ,
              weekend_vector_coef = weekend_vector_coef))
}


EstimateCategFuncData_multinormial_weekend_parallel <- function(timestamps01, W, basis_size=25, method="ML")
{
  
  
  num_indv<- nrow(W)
  timeseries_length <-ncol(W)
  category_count <- length(unique(c(W)))
  #weekend_vector <- as.factor(c(rep(c(rep(0,480),rep(1,192)),4))[1:timeseries_length])
  weekend_vector <- as.factor(c(rep(c(rep(0,5*24*60/20),rep(1,2*24*60/20)),4)))[1:timeseries_length]
  Z<-NULL
  # prob<-array(0, c(num_indv, timeseries_length , category_count))
  # weekend_vector_coef <- matrix(0, num_indv, category_count-1)
  #for (i in 1:num_indv){
  
  
  Z_P_WeekendCoef <- foreach(i = 1:num_indv , .packages = c("mgcv")) %dorng%
    
    #T_rep <- foreach(this_row = 1:5) %dorng%
    { #source("./source_code/R/data_generator.R")
      
      Z_P_WeekendCoef <- numeric(timeseries_length*(category_count+category_count-1)+category_count-1)
  
    print(i)
    if ( length(table(W[i,])) == category_count) {
      
      fit_binom<-gam(list(W[i,]-1~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector,
                          ~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector
      ),
      family=multinom(K=category_count-1), method = method,
      control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
      optimizer=c("outer","bfgs")) 
      #####
      #print("Dino")
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
      if (names(table(W[i,]))[2]=="3"){
        W[i,][W[i,]==3] <- 2
        basis_size_rev <- max(min(round(min(unname(table(W[i,])[2]), sum(1-unname(table(W[i,])[2])))/2), basis_size ), 5)
        fit_binom <- gam(W[i,]-1~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                         family = "binomial", method = method,
                         control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                         optimizer=c("outer","bfgs"))
       # print("GODAK ADERI")
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
        p3 <- predict(fit_binom,type="response")
        p1 <- 1-p3
        p2 <- rep(0,timeseries_length)
        ##########################
      }else {
        basis_size_rev <- max(min(round(min(unname(table(W[i,])[2]), sum(1-unname(table(W[i,])[2])))/2), basis_size ), 5)
        fit_binom <- gam(W[i,]-1~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                         family = "binomial", method = method,
                         control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                         optimizer=c("outer","bfgs"))
        #print("MATA OYA NATHWA PALUI")
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
        ##########################
        p2 <- predict(fit_binom,type="response")
        p1 <- 1-p2
        p3 <- rep(0,timeseries_length)
      }
    } 
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
    #prob[i,,] <- cbind(p1, p2, p3)
    # 5*t +2 length
    Z_P_WeekendCoef <- c(z1,z2,p1,p2,p3, weekend_vector_coef)
    
    return(Z_P_WeekendCoef)

}

  Z_P_WeekendCoef <- do.call(rbind, Z_P_WeekendCoef)
  #t*n
  return(list(Z1_est=t(Z_P_WeekendCoef[,1:timeseries_length]), 
              Z2_est=t(Z_P_WeekendCoef[,(1+timeseries_length):(2*timeseries_length)]),
              p1_est=t(Z_P_WeekendCoef[,(1+timeseries_length*2):(3*timeseries_length)]), 
              p2_est=t(Z_P_WeekendCoef[,(1+timeseries_length*3):(4*timeseries_length)]), 
              p3_est=t(Z_P_WeekendCoef[,(1+timeseries_length*4):(5*timeseries_length)]) ,
              weekend_vector_coef = Z_P_WeekendCoef[,(1+timeseries_length*(category_count+category_count-1)):dim(Z_P_WeekendCoef)[2]]))
}


######test
#EstimateCategFuncData_multinormial <- function(timestamps01, W, basis_size=25, method="ML")
#3001
# timestamps01 <- timestamps01
# test_mul <- EstimateCategFuncData_multinormial_weekend (timestamps01, W_matrix_final, basis_size=25, method="ML")
# W_matrix_final_sub <- W_matrix_final[2999:dim(W_matrix_final)[1],]
# ##3
# test_mul_fail <- EstimateCategFuncData_multinormial_weekend (timestamps01, W_matrix_final_sub, basis_size=25, method="ML")

#which(W_matrix_final_sub[3,]==3) #761 #not work
#Error in qr.default(if (d[1L] < d[2L]) t(z) else z) : 
#NA/NaN/Inf in foreign function call (arg 1)
#save(W_matrix,W_matrix_final, W_matrix_final_sub, W_sample_2, file = "W_matrix_samples.RData")
#load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/acj/W_matrix_samples.RData")
# indiv_one_cat_3_once <- which(apply(W_matrix_final , 1, function(row){length(table(row)) == 3 && min(table(row)) == 1}))
# #233 has 3 category but only 3 only appears 1 time , 14 only has 1 category
# W_matrix_atleast2cat_3atleast2 <- W_matrix_final[-c(indiv_one_cat_3_once),]
# save(W_matrix,W_matrix_final, W_matrix_final_sub, W_sample_2,W_matrix_atleast2cat_3atleast2, 
#      indiv_one_cat_3_once, indiv_one_cat, file = "W_matrix_samples.RData")
# 
# #which(W_matrix_final[24,]==3), 1575  #works
# test_mul_fail_1 <- EstimateCategFuncData_multinormial_weekend (timestamps01, W_matrix_atleast2cat_3atleast2, basis_size=25, method="ML")
# # indiv_one_cat_3_once
# [1]   24   25   30   58   63   73   81   86   90  117  118  128  140  155  156  164  170  179  184  189  200  202  207  209
# [25]  212  222  230  238  241  262  287  288  291  298  305  366  381  386  388  400  410  424  439  440  454  553  554  561
# [49]  575  591  615  622  626  628  634  651  665  678  728  732  740  772  781  784  796  797  806  815  835  845  854  876
# [73]  895  913  920 1341 1487 1491 1494 1528 1543 1704 1801 1811 1835 1837 1842 1909 1915 1918 1921 1923 1933 1961 1973 1988
# [97] 1991 1994 2002 2011 2017 2025 2040 2055 2081 2084 2153 2156 2167 2176 2180 2208 2230 2256 2268 2287 2300 2305 2307 2311
# [121] 2314 2331 2341 2343 2354 2356 2359 2376 2378 2380 2398 2417 2419 2436 2438 2444 2459 2463 2471 2476 2479 2482 2495 2498
# [145] 2499 2502 2528 2546 2568 2569 2571 2582 2595 2605 2610 2619 2638 2684 2706 2707 2711 2713 2716 2729 2737 2739 2753 2754
# [169] 2780 2786 2802 2841 2874 2916 2931 2955 2959 2977 2984 3001 3029 3046 3057 3061 3062 3111 3117 3134 3135 3136 3143 3145
# [193] 3162 3166 3198 3199 3216 3338 3341 3349 3374 3392 3424 3469 3470 3480 3515 3525 3526 3528 3539 3541 3567 3583 3594 3596
# [217] 3599 3601 3617 3638 3733 3768 3769 3770 3779 3783 3790 3793 3808 3812 3814 3832 3845
#indiv_one_cat_3_once_check_if_work <- indiv_one_cat_3_once[indiv_one_cat_3_once>3001]
#53 of them, #21 that is 3374 doesnt work

# W_3_once_check_if_work <- W_matrix_final[indiv_one_cat_3_once_check_if_work,]
# source("time_track_function.R")
# timeKeeperStart("Check Non Parallel")
# test_mul_fail_check <- EstimateCategFuncData_multinormial_weekend (timestamps01, W_3_once_check_if_work[22:dim(W_3_once_check_if_work)[1],], basis_size=25, method="ML")
# #sample_W <- matrix(sample(c(1,2,3), 300 *10, replace = TRUE), ncol= 300, nrow= 10)
# timeKeeperNext() 
# 
# timeKeeperStart("Check Parallel")
# test_mul_fail_check_parallel <- EstimateCategFuncData_multinormial_weekend_parallel (timestamps01, W_3_once_check_if_work[22:dim(W_3_once_check_if_work)[1],], basis_size=25, method="ML")
# #sample_W <- matrix(sample(c(1,2,3), 300 *10, replace = TRUE), ncol= 300, nrow= 10)
# timeKeeperNext() 

#matplot(timestamps01, test_mul_fail_check$Z1_est - test_mul_fail_check_parallel$Z1_est, type="l")

#timeKeeperStart("Twitter all users except 14 has only 1 category, and two only has 3 one time")
#W_matrix_final[-c(3001,3374),]  #3796
# table(W_matrix_final[-c(3001,3374),][3796,])
# 
# 1    2    3 
# 1985   10    5
# Error in qr.default(if (d[1L] < d[2L]) t(z) else z) : 
#   NA/NaN/Inf in foreign function call (arg 1)
# Twitter_ZP__WeekendCoeff_Final <- EstimateCategFuncData_multinormial_weekend_parallel (timestamps01, W_matrix_final[-c(3001,3374),][-3796,], basis_size=25, method="ML")
# timeKeeperNext() 
# W_data_used <- W_matrix_final[-c(3001,3374),][-3796,]
# index_exclude <- c(c(3001,3374),3796 )
# save(Twitter_ZP__WeekendCoeff_Final, timestamps01, 
#      W_data_used, W_matrix_final, index_exclude, file = "Twitter_ZP__WeekendCoeff_Final.RData")
# -------------------
#   Twitter all users except 14 has only 1 category, and two only has 3 one time 
# took: Time difference of 5.322264 mins 
# ====================

