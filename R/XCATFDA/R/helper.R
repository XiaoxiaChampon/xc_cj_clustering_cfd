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
get_x_from_w <- function(w_mat) {
  stopifnot(is.matrix(w_mat))

  n_individuals <- ncol(w_mat)
  n_timepoints <- nrow(w_mat)

  # Extract sorted unique categories
  categories <- sort(unique(as.vector(w_mat)))
  n_categories <- length(categories)

  # Map category labels to 1:K indices
  w_indexed <- apply(w_mat, 2, function(col) match(col, categories))

  # Initialize result array: [individual × time × category]
  x_array <- array(0, dim = c(n_individuals, n_timepoints, n_categories))

  for (i in seq_len(n_individuals)) {
    x_array[i, cbind(seq_len(n_timepoints), w_indexed[, i])] <- 1
  }

  return(x_array)
}
