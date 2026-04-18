############################################################
# Article: Functions Needed for Genarating Categorical Functional Data
#           File that contains all the functions necessary to estimate categorical FD with 3 CATEGORIES! 
# Author:  Xiaoxia Champon, Ana-Maria Staicu, Chathura Jayalath
# Date: 09/25/2023
##############################################################


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
GenerateClusterData <- function(setting, scenario, k, num_indvs, timeseries_length)
{
  setting_object <- GetMuAndScore(setting, scenario, k)
  cluster_f <- GenerateClusterDataScenario(num_indvs,
                                           timeseries_length,
                                           k,
                                           mu_1 = setting_object$mu_1,
                                           mu_2 = setting_object$mu_2,
                                           score_vals = setting_object$score_vals)
  return (cluster_f)
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
GetMuAndScore <- function(setting, scenario, k)
{
  all_score_values = rep(0, k)
  
  if(1 == setting)
  {
    mu_1 <- function(t) -1 + 2 * t + 2 * t^2
    
    mu_2 <- switch(scenario,
                   "A" = function(t) -2.5 + exp(t * 2),
                   "B" = function(t) -0.5 + exp(t * 2),
                   "C" = function(t) -2.5 + exp(t * 2)
    )
    score_front <- switch(scenario,
                          "A" = c(1, 1/2, 1/4),
                          "B" = c(1, 1/2, 1/4),
                          "C" = c(50, 25, 5)
    )
  } else if(2 == setting)
  {
    mu_1 <- function(t) 4 * t^2 - 1.2
    
    mu_2 <- function(t) 4 * t^2 - 3.5
    
    score_front <- c(1, 1/2, 1/4)
  } else if(3 == setting)
  {
    mu_1 <- function(t) -2.2 + 4 * t^2
    
    mu_2 <- function(t) -7 + 6 * t^2
    
    score_front <- c(1, 1/4, 1/16)
  }
  
  for(idx in 1:length(score_front))
  {
    all_score_values[idx] <- score_front[idx]
  }
  
  return(list("mu_1" = mu_1, "mu_2" = mu_2, "score_vals" = all_score_values))
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
                                        k = 3,
                                        mu_1,
                                        mu_2,
                                        score_vals)
{
  timestamps01 <- seq(from = 0.0001, to = 1, length=timeseries_length)
  
  # noise octaves
  # cat("octave", num_indvs, k, num_indvs * k, "\n")
  scores_standard <- matrix(rnorm(num_indvs * k), ncol = k)
  scores <- scores_standard %*% diag(sqrt(score_vals))
  
  #
  BIG_mu <- c(mu_1(timestamps01), mu_2(timestamps01))
  BIG_phi <- PsiFunc(k, timestamps01)
  
  Z <- BIG_phi %*% t(scores) + BIG_mu
  Z1 <- Z[1:timeseries_length, ]
  Z2 <- Z[1:timeseries_length + timeseries_length, ]
  expZ1 <- exp(Z1)
  expZ2 <- exp(Z2)
  denom <- 1 + expZ1 + expZ2
  p1 <- expZ1 / denom
  p2 <- expZ2 / denom
  p3 <- 1 / denom
  
  # vectorize for future work!!!
  return(list(Z1 = Z1, Z2 = Z2,
              p1 = p1, p2 = p2, p3 = p3,
              MEAN = BIG_mu, PHI = BIG_phi, MFPC = scores))
}

#' Psi function
#'
PsiFunc <- function(klen, timestamps01)
{
  psi_k1 <- sapply(c(1:klen), function(i) sin((2 * i + 1) * pi * timestamps01))
  psi_k2 <- sapply(c(1:klen), function(i) cos(2 * i * pi * timestamps01))
  return(rbind(psi_k1, psi_k2))
}