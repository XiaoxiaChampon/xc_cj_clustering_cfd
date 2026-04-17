#' Estimate Z and p curves for categorical functional data (multinomial case)
#'
#' @description
#' Estimates latent Gaussian processes and probability curves using multinomial
#' GAM models. This approach models all categories simultaneously using a single
#' multinomial distribution.
#'
#' @param time_points A numeric vector of time points
#' @param w_mat A matrix of categorical values (time × individuals)
#' @param n_basis Number of basis functions for smoothing (default = 25)
#' @param method GAM fitting method (default = "ML")
#'
#' @return A named list with latent Z estimates and probability estimates
#'
#' @details
#' The multinomial approach fits a single GAM with multinomial family to all
#' categories simultaneously. This is generally more efficient and statistically
#' coherent than fitting separate binomial models, especially when the number
#' of categories is large.
#'
#' Category labels can be arbitrary (numeric or character). Internally, labels
#' are sorted and mapped to positions 1..K via \code{match()}, and then the
#' K-th position is remapped to 0 to serve as the reference class required by
#' \code{mgcv::multinom}. The reference category is always the largest
#' (last in sorted order) original label.
#'
#' @examples
#' # Generate sample data
#' set.seed(123)
#' n_time <- 30
#' n_individuals <- 5
#' time_points <- seq(0, 1, length.out = n_time)
#' w_mat <- matrix(sample(1:3, n_time * n_individuals, replace = TRUE),
#'                 nrow = n_time, ncol = n_individuals)
#'
#' # Estimate using multinomial approach
#' result <- estimate_categ_func_data_multinomial(time_points, w_mat)
#' names(result)
#'
#' @export
estimate_categ_func_data_multinomial <- function(time_points,
                                                 w_mat,
                                                 n_basis = 25,
                                                 method = "ML") {
  stopifnot(is.numeric(time_points), is.matrix(w_mat))

  n_individuals <- ncol(w_mat)
  n_timepoints <- nrow(w_mat)

  # Determine number of categories (K); sort so the mapping is deterministic.
  unique_labels <- sort(unique(as.vector(w_mat)))
  k_classes <- length(unique_labels)

  # Remap arbitrary labels to positions 1..K, then convert the K-th position
  # to 0 so it becomes the reference class for mgcv::multinom (which requires
  # 0-based response values where 0 denotes the reference category).
  #
  # Two-step rationale:
  #   match() : maps any label set (e.g. {1,3,5} or {"A","B","C"}) to 1..K
  #   %% K    : maps position K -> 0, leaving 1..K-1 unchanged
  #
  # Example: labels {1,2,3} -> match -> {1,2,3} -> %%3 -> {1,2,0}
  # Example: labels {0,1,2} -> match -> {1,2,3} -> %%3 -> {1,2,0}  (fixed!)
  # Example: labels {1,3,5} -> match -> {1,2,3} -> %%3 -> {1,2,0}  (fixed!)
  w_mapped <- matrix(match(as.vector(w_mat), unique_labels), nrow = nrow(w_mat)) %% k_classes

  # Preallocate Z and p
  z_list <- vector("list", k_classes - 1)
  p_array <- array(0, dim = c(n_individuals, n_timepoints, k_classes))

  for (i in seq_len(n_individuals)) {
    # Construct list of formulas for GAM (first is response, rest are smoothers)
    gam_formulas <- c(
      list(w_mapped[, i] ~ s(time_points, bs = "cr", m = 2, k = n_basis)),
      replicate(k_classes - 2, ~ s(time_points, bs = "cr", m = 2, k = n_basis), simplify = FALSE)
    )

    fit <- mgcv::gam(
      formula = gam_formulas,
      family = mgcv::multinom(K = k_classes - 1),
      method = method,
      control = list(maxit = 500, mgcv.tol = 1e-4, epsilon = 1e-04),
      optimizer = c("outer", "bfgs")
    )

    z_mat <- fit$linear.predictors # t × (K - 1)
    for (k in seq_len(k_classes - 1)) {
      z_list[[k]] <- cbind(z_list[[k]], z_mat[, k])
    }

    # Compute softmax probabilities
    exp_z <- exp(z_mat)
    denom <- 1 + rowSums(exp_z)
    prob_mat <- matrix(0, nrow = n_timepoints, ncol = k_classes)

    for (k in seq_len(k_classes - 1)) {
      prob_mat[, k] <- exp_z[, k] / denom
    }
    prob_mat[, k_classes] <- 1 / denom

    p_array[i, , ] <- prob_mat
  }

  # Build output: z_list already in time × individuals, p_array needs transpose
  z_estimates <- rlang::set_names(
    z_list,
    paste0("z", seq_along(z_list), "_est")
  )
  p_estimates <- rlang::set_names(
    lapply(seq_len(k_classes), function(k) t(p_array[, , k])),
    paste0("p", seq_len(k_classes), "_est")
  )

  return(c(z_estimates, p_estimates))
}


#' Parallel estimation of Z and p using multinomial link GAM
#'
#' @description
#' Parallel version of multinomial estimation for improved performance
#' with large datasets.
#'
#' @param time_points A numeric vector of time points
#' @param w_mat A matrix of categorical values (time × individuals)
#' @param n_basis Number of basis functions for smoothing (default = 25)
#' @param method GAM fitting method (default = "ML")
#'
#' @return A named list with latent Z estimates and probability estimates
#'
#' @details
#' This is the parallel implementation of the multinomial approach.
#' It uses the \\code{foreach} package with \\code{\\%dorng\\%} for reproducible
#' parallel computation.
#'
#' Note: Requires a parallel backend to be registered (e.g., using doParallel)
#' for actual parallel execution.
#'
#' @examples
#' \\dontrun{
#' # Setup parallel backend
#' library(doParallel)
#' cl <- makeCluster(2)
#' registerDoParallel(cl)
#'
#' # Generate sample data
#' set.seed(123)
#' n_time <- 50
#' n_individuals <- 20
#' time_points <- seq(0, 1, length.out = n_time)
#' w_mat <- matrix(sample(1:3, n_time * n_individuals, replace = TRUE),
#'                 nrow = n_time, ncol = n_individuals)
#'
#' # Parallel estimation
#' result <- estimate_categ_func_data_multinomial_parallel(time_points, w_mat)
#'
#' # Cleanup
#' stopCluster(cl)
#' }
#'
#' @export
estimate_categ_func_data_multinomial_parallel <- function(time_points,
                                                          w_mat,
                                                          n_basis = 25,
                                                          method = "ML") {
  stopifnot(is.numeric(time_points), is.matrix(w_mat))

  n_individuals <- ncol(w_mat)
  n_timepoints <- nrow(w_mat)

  # Determine number of categories (K); sort so the mapping is deterministic.
  unique_labels <- sort(unique(as.vector(w_mat)))
  k_classes <- length(unique_labels)

  # Remap arbitrary labels to positions 1..K, then convert the K-th position
  # to 0 so it becomes the reference class for mgcv::multinom. See the
  # non-parallel version for a full explanation of the two-step approach.
  w_mapped <- matrix(match(as.vector(w_mat), unique_labels), nrow = nrow(w_mat)) %% k_classes

  total_z_rows <- n_timepoints * (k_classes - 1)

  # Parallel processing across individuals
  zp <- foreach(
    i = seq_len(n_individuals),
    .combine = cbind,
    .packages = c("mgcv")
  ) %dorng% {
    # Construct list of formulas for GAM
    gam_formulas <- c(
      list(w_mapped[, i] ~ s(time_points, bs = "cr", m = 2, k = n_basis)),
      replicate(k_classes - 2, ~ s(time_points, bs = "cr", m = 2, k = n_basis), simplify = FALSE)
    )

    fit <- mgcv::gam(
      formula = gam_formulas,
      family = mgcv::multinom(K = k_classes - 1),
      method = method,
      control = list(maxit = 500, mgcv.tol = 1e-4, epsilon = 1e-04),
      optimizer = c("outer", "bfgs")
    )

    z_mat <- fit$linear.predictors # t × (K - 1)

    # Compute softmax probabilities
    exp_z <- exp(z_mat)
    denom <- 1 + rowSums(exp_z)
    prob_mat <- matrix(0, nrow = n_timepoints, ncol = k_classes)

    for (k in seq_len(k_classes - 1)) {
      prob_mat[, k] <- exp_z[, k] / denom
    }
    prob_mat[, k_classes] <- 1 / denom

    # Return z values and probabilities as a single vector
    return(c(
      as.vector(z_mat), # (T * (K - 1))
      as.vector(prob_mat) # (T * K)
    ))
  }

  # Split output back into Z and p
  z_mat <- zp[1:total_z_rows, , drop = FALSE]
  p_mat <- zp[(total_z_rows + 1):nrow(zp), , drop = FALSE]

  # Return list of Z1_est, ..., Z{K-1}_est - already in time × individuals format
  z_estimates <- rlang::set_names(
    lapply(seq_len(k_classes - 1), function(k) {
      z_start <- (k - 1) * n_timepoints + 1
      z_end <- k * n_timepoints
      z_mat[z_start:z_end, ]
    }),
    paste0("z", seq_len(k_classes - 1), "_est")
  )

  # Return list of p1_est, ..., pK_est
  p_array <- array(t(p_mat), dim = c(n_individuals, n_timepoints, k_classes))
  p_estimates <- rlang::set_names(
    lapply(seq_len(k_classes), function(k) t(p_array[, , k])),
    paste0("p", seq_len(k_classes), "_est")
  )

  return(c(z_estimates, p_estimates))
}
