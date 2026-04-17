#' Function to produce functional dummy variables X from categorical functional data W
#' @param W 2D array, t*n: t is the timestamp and n is the number of the observation
#' @return X 3D array, n*t*Q, Q: the total number of the category
GetXFromW <- function(W) {
  num_indv <- ncol(W)
  timeseries_length <- nrow(W)
  category_count <- length(unique(c(W)))
  Q_vals <- unique(c(W))
  if (is.numeric(Q_vals)) Q_vals <- sort(Q_vals)

  X <- array(0, c(num_indv, timeseries_length, category_count))
  for (indv in 1:num_indv)
  {
    for (timestamps01 in 1:timeseries_length)
    {
      X[indv, timestamps01, which(Q_vals == W[, indv][timestamps01])] <- 1
    }
  }
  return(X)
}

#' Wrapper for backward compatibility with old CamelCase naming
#' @param choice Estimation method: "multinomial", "probit", or "binomial"
#' @param time_points Time points vector
#' @param w_mat Categorical data matrix
#' @param n_basis Number of basis functions
#' @param method GAM fitting method
#' @return List with estimated Z and p curves
#' @export
EstimateCategFuncData <- function(choice, time_points, w_mat, n_basis = 25, method = "ML") {
  estimate_categ_func_data(choice, time_points, w_mat, n_basis, method)
}

#' Create one-hot encoded functional data array from categorical observations
#'
#' @description
#' Converts a matrix of categorical observations into a 3D array of one-hot
#' encoded binary indicators. This is a utility function commonly needed
#' for preprocessing categorical functional data.
#'
#' @param w_mat A matrix of categorical observations (time × individuals)
#'
#' @return A 3D array (individual × time × category)
#'
#' @details
#' The function automatically detects the unique categories in the input matrix
#' and creates binary indicator variables for each category. Categories are
#' sorted and mapped to sequential indices.
#'
#' The output array has dimensions:
#' \itemize{
#'   \item Dimension 1: Individuals
#'   \item Dimension 2: Time points
#'   \item Dimension 3: Categories (one-hot encoded)
#' }
#'
#' @examples
#' # Create sample categorical data
#' set.seed(123)
#' n_time <- 20
#' n_individuals <- 5
#' w_mat <- matrix(sample(c("A", "B", "C"), n_time * n_individuals, replace = TRUE),
#'                 nrow = n_time, ncol = n_individuals)
#'
#' # Convert to one-hot encoding
#' x_array <- get_x_from_w(w_mat)
#' dim(x_array)  # Should be [5, 20, 3]
#'
#' # Check that each time point has exactly one category active
#' apply(x_array[1, , ], 1, sum)  # Should be all 1s
#'
#' @export
get_x_from_w <- function(w_mat, categories = NULL) {
  # validate w_mat
  if (!is.matrix(w_mat)) {
    stop("w_mat must be a matrix")
  }
  if (anyNA(w_mat)) {
    stop("w_mat must not contain NA values")
  }

  # Extract sorted unique categories
  if (is.null(categories)) {
    categories <- sort(unique(as.vector(w_mat)))
  }

  # validate categories: Ensure categories are unique and sorted
  if (length(categories) != length(unique(categories))) {
    stop("categories must contain unique values")
  }

  n_individuals <- ncol(w_mat)
  n_timepoints <- nrow(w_mat)
  n_categories <- length(categories)

  # Map category labels to 1:K indices
  w_indexed <- matrix(match(w_mat, categories),
                      nrow = n_timepoints,
                      ncol = n_individuals)

  # Fail fast if unseen category appears
  if (anyNA(w_indexed)) {
    stop("w_mat contains values not present in 'categories'")
  }

  # Initialize result array: [individual × time × category]
  x_array <- array(0, dim = c(n_individuals, n_timepoints, n_categories))

  for (i in seq_len(n_individuals)) {
    x_array[cbind(rep(i, n_timepoints), seq_len(n_timepoints), w_indexed[, i])] <- 1
  }

  return(x_array)
}
