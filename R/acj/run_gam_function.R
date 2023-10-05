############################################################
# Article: Functions Needed for running gam on binary or categorical functional data
# Author:  Xiaoxia Champon, Ana-Maria Staicu, Chathura Jayalath
# Date: 09/25/2023
##############################################################


#' Runs gam on the given data and returns results
#' Fits a binomial to describe the given in_x
#' @param timestamps01 current time value
#' @param in_x binary series
#' @param family_choice "binomial" or "probit"
#' @return list: fit values and linear predictors both with length of time_series length
RunGam <- function(timestamps01, in_x, family_choice, basis_size, method)
{
  basis_size_rev <- max(min(round(min(sum(in_x), sum(1-in_x))/2), basis_size ), 5)
  
  if(family_choice=="binomial")
  {
    family_object <- "binomial"
  }else if(family_choice=="probit")
  {
    family_object <- binomial(link="probit")
  }
  
  fit_binom <- gam(in_x~s(timestamps01, bs = "cr", m=2, k = basis_size_rev),
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
