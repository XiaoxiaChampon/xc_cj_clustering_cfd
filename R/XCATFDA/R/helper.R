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
#' @param w_mat A matrix of categorical observations (time × individuals)
#'
#' @return A 3D array (individual × time × category)
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
