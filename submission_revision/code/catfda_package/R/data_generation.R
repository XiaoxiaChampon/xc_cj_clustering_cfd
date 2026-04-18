
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
#' @export
generate_categ_func_data <- function(prob_curves, tol = 1e-8) {
  if(!is.list(prob_curves)) {
    stop("prob_curves must be a list")
  }
  if(length(prob_curves) < 2) {
    stop("prob_curves must contain at least two categories")
  }
  if (!all(sapply(prob_curves, is.matrix))) {
    stop("All elements of prob_curves must be matrices")
  }

  k_categories <- length(prob_curves)

  # Validate structure
  ref_dim <- dim(prob_curves[[1]])
  if(length(ref_dim) != 2){
    stop("Each matrix in prob_curves must be 2-dimensional")
  }
  if (!all(sapply(prob_curves, function(mat) all(dim(mat) == ref_dim)))) {
    stop("All matrices in prob_curves must have the same (T x N) dimensions")
  }

  n_timepoints <- ref_dim[1]
  n_individuals <- ref_dim[2]

  # validate probabilities
  p_arr <- array(NA_real_, dim = c(n_timepoints, n_individuals, k_categories))
  for (k in seq_len(k_categories)) {
    p_mat <- prob_curves[[k]]
    if (any(p_mat < -tol) || any(p_mat > 1 + tol)) {
      stop(paste0("Probabilities in category ", k, " are out of [0, 1] range (with tolerance of ", tol, ")"))
    }
    p_arr[, , k] <- p_mat
  }
  if (anyNA(p_arr)) {
    stop("NA values found in probability curves")
  }

  # clamp probabilities to [0, 1] with tolerance checked above
  p_arr[p_arr < 0] <- 0
  p_arr[p_arr > 1] <- 1

  # Check sum to 1
  p_sums <- rowSums(p_arr, dims = 2) # T × N
  if (any(abs(p_sums - 1) > tol)) {
    stop("Probabilities across categories do not sum to 1 at some time points/individuals")
  }

  # Prepare storage
  w_mat <- matrix(0, nrow = n_timepoints, ncol = n_individuals)
  x_array <- array(0, dim = c(n_individuals, n_timepoints, k_categories))

  for (indv in seq_len(n_individuals)) {
    # Build T × K matrix for this individual
    p_matrix <- vapply(seq_len(k_categories), function(k) p_arr[, indv, k], numeric(n_timepoints))
    p_matrix <- as.matrix(p_matrix)  # ensure matrix even when K=2

    # >>>>> Check for Correctness with the following>>>>>>>>
    #p_matrix <- sapply(prob_curves, function(p_mat) p_mat[, indv]) # T × K

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
