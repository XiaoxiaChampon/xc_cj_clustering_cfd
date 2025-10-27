#' Generate synthetic categorical functional data (W) and one-hot encoded array (X)
#'
#' @description
#' Generates synthetic categorical functional data from probability curves using
#' multinomial sampling. This is useful for simulation studies and testing.
#'
#' @param prob_curves A named list of T × N matrices: p1_est, ..., pK_est
#'
#' @return A list with:
#'   - w_mat: Categorical matrix (T × N)
#'   - x_array: One-hot encoded array (N × T × K)
#'
#' @details
#' The function takes probability curves for each category and generates
#' categorical observations by sampling from multinomial distributions at
#' each time point for each individual.
#'
#' The probability curves should be provided as a named list where each
#' element is a matrix of dimension (time × individuals).
#'
#' @examples
#' # Create sample probability curves
#' set.seed(123)
#' n_time <- 30
#' n_individuals <- 10
#' time_points <- seq(0, 1, length.out = n_time)
#'
#' # Generate simple probability matrices that sum to 1
#' p1 <- matrix(0.4, nrow = n_time, ncol = n_individuals)
#' p2 <- matrix(0.3, nrow = n_time, ncol = n_individuals)
#' p3 <- matrix(0.3, nrow = n_time, ncol = n_individuals)
#'
#' prob_curves <- list(p1_est = p1, p2_est = p2, p3_est = p3)
#'
#' # Generate categorical data
#' result <- generate_categ_func_data(prob_curves)
#'
#' # Check output
#' dim(result$w_mat)    # Should be [30, 10]
#' dim(result$x_array)  # Should be [10, 30, 3]
#'
#' @export
generate_categ_func_data <- function(prob_curves) {
  stopifnot(is.list(prob_curves), length(prob_curves) >= 2)

  k_categories <- length(prob_curves)
  category_names <- names(prob_curves)

  # Validate structure and extract dimensions
  if (!is.matrix(prob_curves[[1]])) {
    stop("Each element in prob_curves must be a matrix")
  }
  ref_dim <- dim(prob_curves[[1]])

  # Check all elements are matrices with same dimensions
  for (i in seq_along(prob_curves)) {
    if (!is.matrix(prob_curves[[i]])) {
      stop("Each element in prob_curves must be a matrix")
    }
    if (!all(dim(prob_curves[[i]]) == ref_dim)) {
      stop("All matrices in prob_curves must have the same dimensions")
    }
  }

  n_timepoints <- as.integer(ref_dim[1])
  n_individuals <- as.integer(ref_dim[2])  # Prepare storage
  w_mat <- matrix(0, nrow = n_timepoints, ncol = n_individuals)
  x_array <- array(0, dim = c(n_individuals, n_timepoints, k_categories))

  for (indv in seq_len(n_individuals)) {
    p_matrix <- sapply(prob_curves, function(p_mat) p_mat[, indv]) # T × K

    # For each time point, sample from the categorical distribution
    sampled <- apply(p_matrix, 1, function(prob_vec) {
      stats::rmultinom(1, size = 1, prob = prob_vec)
    })

    # sampled: K × T
    w_mat[, indv] <- apply(sampled, 2, which.max) # T × 1
    x_array[indv, , ] <- t(sampled) # T × K
  }

  return(list(
    x_array = x_array, # N × T × K (binary dummy)
    w_mat = w_mat # T × N (categorical index)
  ))
}
