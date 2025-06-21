library(refund)

#*** extract scores using UNIVARIAFE fpca, FAST COV estimation.  SUPER FAST
#INPUT
# latentProcesses: a list of matrices, each of shape m x n (latent processes 
#                  for each category)
# tt: grid of length m (same across all processes)
# PVE: proportion of variance explained (default = 0.95)
# OUTPUT
# list with estimated PHI (m by K) with function norm 1
#           projection of data onto PHI (called scores - nby K matrix)
#

extract_scores_univfpca <- function(latent_proc_list, tt, pve_threshold = 0.95)
{
  num_latent_procs <- length(latent_proc_list) # L
  num_timesteps <- nrow(latent_proc_list[[1]]) # m
  num_values <- ncol(latent_proc_list[[1]])    # n
  sqrt_num_timesteps <- sqrt(num_timesteps)    # sqrt(m)

  # Pass 1: Calculate fpca for each latent process and 
  #           store the values to be used in the pass 2.
  scores_list <- vector("list", num_latent_procs)
  eigen_funcs_list <- vector("list", num_latent_procs)
  npc_vec <- integer(num_latent_procs)
  total_npcs <- 0
  for (lpidx in 1:num_latent_procs) {
    fpca_i <- fpca.face(Y = t(latent_proc_list[[lpidx]]), argvals = tt, pve = 0.99)
    scores_list[[lpidx]] <- fpca_i$scores
    eigen_funcs_list[[lpidx]] <- fpca_i$efunctions
    npc_vec[lpidx] <- fpca_i$npc
    total_npcs <- total_npcs + fpca_i$npc
  }
  
  # Pass 2: Create the actual scores output and estimate Phi
  phi_est0 <- matrix(0, nrow = num_timesteps * num_latent_procs, ncol = total_npcs)
  
  ef_col_start <- 1
  for (lpidx in 1:num_latent_procs) {
    npci <- npc_vec[lpidx]
    
    ef_col_end <- ef_col_start + npci - 1
    ef_row_start <- 1 + num_timesteps * (lpidx - 1)
    ef_row_end <- num_timesteps * lpidx
    
    phi_est0[ef_row_start:ef_row_end, ef_col_start:ef_col_end] <- eigen_funcs_list[[lpidx]] * sqrt_num_timesteps
    
    # calc col start for next iteration
    ef_col_start <- ef_col_end + 1
  }
  scores0 <- do.call(cbind, scores_list)
  
  out_eigen0 <- eigen( cov(scores0) ) # ,symmetric = TRUE, only.values = FALSE)
  cum_var <- cumsum(out_eigen0$values) / sum(out_eigen0$values)
  
  last_eigen_idx <- which(cum_var >= pve_threshold)[1]
  count_iter <- 0
  delta <- 0.01
  while (last_eigen_idx < 2 && count_iter < 100){
    count_iter <- count_iter + 1
    cat("count_iter: ", count_iter, "\n")
    last_eigen_idx <- which(cum_var >= (pve_threshold + delta))[1]
    delta <- delta + 0.01
  }
  
  phi_est <- phi_est0 %*% out_eigen0$vectors[,1:last_eigen_idx]
  
  mZ <- do.call(rbind, latent_proc_list) # shape: (L * m) x n
  
  scores_est <- crossprod(mZ, phi_est) / sqrt_num_timesteps
  
  return (list(scores=scores_est, Phi= phi_est))

}

#   out1 <- fpca.face(Y=t(mZ1), argvals =tt, pve = 0.99)
#   out2 <- fpca.face(Y=t(mZ2), argvals =tt, pve = 0.99)
#   O_m1 <- matrix(0, nrow=m, ncol=out1$npc)
#   O_m2  <- matrix(0, nrow=m, ncol=out2$npc)
# 
#   #construct PHI
#   Phi_est0 <-  rbind(cbind(out1$efunctions*sqrt(m),O_m2 ),
#                      cbind(O_m1,out2$efunctions*sqrt(m)))
#   Scores0 <- cbind(out1$scores, out2$scores)
#   ScoresCov_0 <- cov(Scores0 )
#   oute <- eigen(ScoresCov_0)
# 
#   K<- which(cumsum( oute$values)/sum(oute$values)>=pve_threshold)[1]
#   count_iter = 0
#   delta=0.01
#   while (K<2 && count_iter<100) {
#     count_iter = count_iter + 1
#     cat("count_iter: ", count_iter, "\n")
#     K<- which(cumsum( oute$values)/sum(oute$values)>=(pve_threshold+delta))[1]
#     delta=delta+0.01
#   }
#   
#   Phi_est <-  Phi_est0%*% oute$vectors[,1:K] # correct eigenfns
#   
#   mZ <- rbind(mZ1, mZ2)
#   Scores_est <- t(mZ) %*%Phi_est/sqrt(m)  # they are not demeaned
#   
#   return (list(scores=Scores_est, Phi= Phi_est))
# }

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