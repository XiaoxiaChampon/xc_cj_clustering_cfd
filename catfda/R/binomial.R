#' Estimate Z and p curves using binomial link GAMs (parallel version)
#'
#' @description
#' Estimates latent Gaussian processes using separate binomial GAM models
#' for each category. This is the parallel implementation for improved
#' performance with large datasets.
#'
#' @param time_points Numeric vector of time values (length T)
#' @param x_array 3D array (individual × time × category) — binary one-hot encoding
#' @param n_basis Number of basis functions (default = 25)
#' @param method GAM optimization method (default = "ML")
#'
#' @return A list with:
#'   - Z1_est, ..., Z{K-1}_est (latent curves)
#'   - p1_est, ..., pK_est (category probability curves)
#'
#' @details
#' The binomial approach fits separate binomial GAM models for each category
#' and then transforms the results to obtain latent processes relative to
#' a reference category (the last category).
#'
#' This parallel implementation uses the \code{foreach} package for efficient
#' computation across individuals.
#'
#' @examples
#' \dontrun{
#' # Setup parallel backend
#' library(doParallel)
#' cl <- makeCluster(2)
#' registerDoParallel(cl)
#'
#' # Generate sample data
#' set.seed(123)
#' n_time <- 75
#' n_individuals <- 15
#' n_categories <- 3
#'
#' time_points <- seq(0, 1, length.out = n_time)
#' x_array <- array(0, dim = c(n_individuals, n_time, n_categories))
#'
#' # Fill with random one-hot vectors
#' for(i in 1:n_individuals) {
#'   for(t in 1:n_time) {
#'     cat <- sample(1:n_categories, 1)
#'     x_array[i, t, cat] <- 1
#'   }
#' }
#'
#' # Parallel binomial estimation
#' result <- estimate_categ_func_data_binomial_parallel(time_points, x_array)
#'
#' # Cleanup
#' stopCluster(cl)
#' }
#'
#' @export
estimate_categ_func_data_binomial_parallel <- function(time_points,
                                                       x_array,
                                                       n_basis = 25,
                                                       method = "ML") {
  stopifnot(length(dim(x_array)) == 3)

  n_individuals <- dim(x_array)[1]
  n_timepoints <- dim(x_array)[2]
  n_categories <- dim(x_array)[3]

  total_z_curves <- n_categories - 1
  total_z_rows <- n_timepoints * total_z_curves

  zp <- foreach::foreach(
    indv = seq_len(n_individuals),
    .combine = cbind,
    .packages = c("mgcv")
  ) %dorng% {
    p_list <- list()
    linpred_list <- list()

    for (k in seq_len(n_categories)) {
      x_k <- x_array[indv, , k]
      fit <- catfda::run_gam(time_points, x_k, link = "binomial", n_basis = n_basis, method = method)
      p_list[[k]] <- fit$prob
      linpred_list[[k]] <- fit$linpred
    }

    # Reference = last category
    ref_linpred <- linpred_list[[n_categories]]

    z_values <- list()
    for (k in seq_len(n_categories - 1)) {
      num <- linpred_list[[k]] - ref_linpred
      denom <- log((1 + exp(linpred_list[[k]])) / (1 + exp(ref_linpred)))
      z_k <- num - denom
      z_values[[k]] <- z_k
    }

    # Normalize probabilities
    p_mat <- do.call(cbind, p_list) # T × K
    p_sum <- rowSums(p_mat)
    p_norm <- sweep(p_mat, 1, p_sum, FUN = "/") # T × K

    return(c(
      unlist(z_values), # (T * (K - 1))
      as.vector(p_norm) # (T * K)
    ))
  }

  # Split output back into Z and p
  z_mat <- zp[1:total_z_rows, , drop = FALSE]
  p_mat <- zp[(total_z_rows + 1):nrow(zp), , drop = FALSE]

  # Return list of z1_est, ..., z{K-1}_est
  z_out <- rlang::set_names(
    lapply(seq_len(total_z_curves), function(k) {
      z_start <- (k - 1) * n_timepoints + 1
      z_end <- k * n_timepoints
      t(z_mat[z_start:z_end, ])
    }),
    paste0("z", seq_len(total_z_curves), "_est")
  )

  # Return list of p1_est, ..., pK_est
  p_array <- array(t(p_mat), dim = c(n_individuals, n_timepoints, n_categories))
  p_out <- rlang::set_names(
    lapply(seq_len(n_categories), function(k) t(p_array[, , k])),
    paste0("p", seq_len(n_categories), "_est")
  )

  return(c(z_out, p_out))
}
