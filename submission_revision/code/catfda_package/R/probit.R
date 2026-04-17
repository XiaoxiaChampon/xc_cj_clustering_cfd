#' Estimate Z and p curves using (generalized) probit or binomial models
#'
#' @param time_points A numeric vector of time values (length T)
#' @param x_array A 3D array (individual × time × category) of one-hot encoded data
#' @param n_basis Number of basis functions for smoothing (default = 25)
#' @param method GAM fitting method (default = "ML")
#' @param threshold_probability Probability threshold for switching to probit special case (default = 0.004)
#'
#' @return A list containing estimated Z curves and category probabilities for each individual
#'
#' @details
#' This function implements an adaptive approach that switches between probit
#' and binomial link functions based on the sparsity of categorical events:
#' \itemize{
#'   \item Uses probit link for rare events (proportion < threshold_probability)
#'   \item Uses binomial link for common events
#'   \item Only applies special probit handling for short time series (≤ 301 points)
#' }
#'
#' The adaptive approach helps with numerical stability when dealing with
#' very rare categorical events.
#'
#' @examples
#' # Generate sample one-hot encoded data
#' set.seed(123)
#' n_time <- 50
#' n_individuals <- 8
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
#' # Estimate using probit approach
#' result <- estimate_categ_func_data_probit(time_points, x_array)
#'
#' @export
estimate_categ_func_data_probit <- function(time_points,
                                            x_array,
                                            n_basis = 25,
                                            method = "ML",
                                            threshold_probability = 0.004) {
  stopifnot(length(dim(x_array)) == 3)

  n_individuals <- dim(x_array)[1]
  n_timepoints <- dim(x_array)[2]
  n_categories <- dim(x_array)[3]

  # Initialize storage
  z_curves <- vector("list", n_categories - 1)
  prob_array <- array(0, dim = c(n_individuals, n_timepoints, n_categories))

  for (indv in seq_len(n_individuals)) {
    p_list <- vector("list", n_categories)
    linpred_list <- vector("list", n_categories)

    # Determine which link to use
    use_probit <- n_timepoints <= 301 &&
      mean(x_array[indv, , 1]) < threshold_probability

    link_func <- if (use_probit) "probit" else "binomial"

    # Fit models for each category
    for (k in seq_len(n_categories)) {
      response <- x_array[indv, , k]
      fit <- run_gam(time_points, response, link_func, n_basis, method)
      p_list[[k]] <- fit$prob
      linpred_list[[k]] <- fit$linpred
    }

    # Compute latent variables Zₖ = f(linear predictors)
    ref_index <- n_categories
    ref_linpred <- linpred_list[[ref_index]]

    for (k in seq_len(n_categories - 1)) {
      num <- linpred_list[[k]] - ref_linpred
      denom <- log((1 + exp(linpred_list[[k]])) / (1 + exp(ref_linpred)))
      z_k <- num - denom
      z_curves[[k]] <- cbind(z_curves[[k]], z_k)
    }

    # Normalize probabilities across categories
    p_vecs <- do.call(cbind, p_list) # t × K
    p_sum <- rowSums(p_vecs)
    p_normalized <- sweep(p_vecs, 1, p_sum, FUN = "/")
    prob_array[indv, , ] <- p_normalized
  }

  # Build output
  z_out <- rlang::set_names(
    lapply(z_curves, t),
    paste0("z", seq_along(z_curves), "_est")
  )
  p_out <- rlang::set_names(
    lapply(seq_len(n_categories), function(k) t(prob_array[, , k])),
    paste0("p", seq_len(n_categories), "_est")
  )

  return(c(z_out, p_out))
}

#' Parallel estimation of Z and p using (probit or binomial) link GAMs
#'
#' @description
#' Parallel version of probit/binomial estimation for improved performance
#' with large datasets.
#'
#' @param time_points A numeric vector of time points (length T)
#' @param x_array A 3D array: individual × time × category (dimensions N × T × K)
#' @param n_basis Number of basis functions for GAM (default = 25)
#' @param method GAM fitting method (default = "ML")
#' @param threshold_probability If response is sparse, switch to probit (default = 0.004)
#'
#' @return A list of estimated latent curves (Z1, Z2, ..., ZK-1) and normalized probabilities (p1, ..., pK)
#'
#' @details
#' This is the parallel implementation of the probit/binomial approach.
#' It uses the \code{foreach} package with \code{\%dorng\%} for reproducible
#' parallel computation.
#'
#' Note: Requires a parallel backend to be registered (e.g., using doParallel)
#' for actual parallel execution.
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
#' n_time <- 100
#' n_individuals <- 20
#' n_categories <- 4
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
#' # Parallel estimation
#' result <- estimate_categ_func_data_probit_parallel(time_points, x_array)
#'
#' # Cleanup
#' stopCluster(cl)
#' }
#'
#' @export
estimate_categ_func_data_probit_parallel <- function(time_points,
                                                     x_array,
                                                     n_basis = 25,
                                                     method = "ML",
                                                     threshold_probability = 0.004) {
  stopifnot(length(dim(x_array)) == 3)

  n_individuals <- dim(x_array)[1]
  n_timepoints <- dim(x_array)[2]
  n_categories <- dim(x_array)[3]

  total_z_curves <- n_categories - 1
  total_z_rows <- n_timepoints * total_z_curves

  # Parallel estimation
  zp <- foreach::foreach(
    indv = seq_len(n_individuals),
    .combine = cbind,
    .packages = c("mgcv"),
    .export = c("run_gam")
  ) %dorng% {
    # Estimate p_k and linpred for each category
    p_list <- list()
    linpred_list <- list()

    for (k in seq_len(n_categories)) {
      x_k <- x_array[indv, , k]

      use_probit <- (n_timepoints <= 301) && (mean(x_k) < threshold_probability)
      link <- if (use_probit) "probit" else "binomial"

      fit <- run_gam(time_points, x_k, link, n_basis, method)
      p_list[[k]] <- fit$prob
      linpred_list[[k]] <- fit$linpred
    }

    # Compute Z_k = f(linear predictors)
    z_values <- list()
    ref_linpred <- linpred_list[[n_categories]]
    for (k in seq_len(n_categories - 1)) {
      num <- linpred_list[[k]] - ref_linpred
      denom <- log((1 + exp(linpred_list[[k]])) / (1 + exp(ref_linpred)))
      z_k <- num - denom
      z_values[[k]] <- z_k
    }

    # Normalize probabilities
    p_mat <- do.call(cbind, p_list) # T × K
    p_sum <- rowSums(p_mat)
    p_norm <- sweep(p_mat, 1, p_sum, FUN = "/") # normalized probs

    # Flatten and return combined result
    return(c(
      unlist(z_values), # (T * (K - 1)) values
      as.vector(p_norm) # (T * K) values
    ))
  }

  # Reconstruct output arrays
  z_mat <- zp[1:total_z_rows, , drop = FALSE]
  p_mat <- zp[(total_z_rows + 1):nrow(zp), , drop = FALSE]

  # Create return list of Z
  z_out <- rlang::set_names(
    lapply(seq_len(total_z_curves), function(k) {
      z_start <- (k - 1) * n_timepoints + 1
      z_end <- k * n_timepoints
      t(z_mat[z_start:z_end, ])
    }),
    paste0("z", seq_len(total_z_curves), "_est")
  )

  # Create return list of p
  p_array <- array(t(p_mat), dim = c(n_individuals, n_timepoints, n_categories))
  p_out <- rlang::set_names(
    lapply(seq_len(n_categories), function(k) t(p_array[, , k])),
    paste0("p", seq_len(n_categories), "_est")
  )

  return(c(z_out, p_out))
}
