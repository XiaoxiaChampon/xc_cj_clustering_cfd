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
#' The function automatically detects the number of categories and handles
#' arbitrary category labels by remapping them to sequential integers.
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

  # Determine number of categories (K)
  unique_labels <- sort(unique(as.vector(w_mat)))
  k_classes <- length(unique_labels)

  # Remap labels to 1...K if needed
  w_mapped <- apply(w_mat, 2, function(col) match(col, unique_labels))

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

  # Build output: transpose and name
  z_estimates <- rlang::set_names(
    lapply(seq_along(z_list), function(k) t(z_list[[k]])),
    paste0("z", seq_along(z_list), "_est")
  )
  p_estimates <- rlang::set_names(
    lapply(seq_len(k_classes), function(k) t(p_array[, , k])),
    paste0("p", seq_len(k_classes), "_est")
  )

  return(c(z_estimates, p_estimates))
}
