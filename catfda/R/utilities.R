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

#' Fit a GAM to a binary time series with a chosen link function
#'
#' @description
#' Internal utility function for fitting Generalized Additive Models (GAMs)
#' to binary time series data using different link functions.
#'
#' @param time_points A numeric vector of time values (length T)
#' @param response A binary vector of observations (length T)
#' @param link A string: either "binomial" (logit) or "probit"
#' @param n_basis Number of basis functions (default = 25)
#' @param method GAM optimization method (default = "ML")
#'
#' @return A list with predicted probabilities and latent Z (linear predictors)
#'
#' @details
#' This is an internal utility function used by the main estimation functions.
#' It handles the GAM fitting with appropriate error handling and basis
#' function adjustment for sparse responses.
#'
#' The function automatically adjusts the number of basis functions to avoid
#' overfitting when dealing with very sparse binary responses.
#'
#' @examples
#' # Fit a simple binomial GAM
#' set.seed(123)
#' n_time <- 50
#' time_points <- seq(0, 1, length.out = n_time)
#' response <- rbinom(n_time, 1, plogis(sin(2 * pi * time_points)))
#'
#' # Fit GAM
#' result <- run_gam(time_points, response, "binomial")
#'
#' # Plot results
#' plot(time_points, response, pch = 16, col = "gray")
#' lines(time_points, result$prob, col = "red", lwd = 2)
#'
#' @export
run_gam <- function(time_points,
                    response,
                    link = "binomial",
                    n_basis = 25,
                    method = "ML") {
  stopifnot(length(response) == length(time_points))
  stopifnot(link %in% c("binomial", "probit"))

  # Avoid overfitting with very sparse response
  basis_eff <- max(min(round(min(sum(response), sum(1 - response)) / 2), n_basis), 5)

  # Set up family object
  family_obj <- if (link == "probit") {
    stats::binomial(link = "probit")
  } else {
    "binomial"
  }

  # Fit GAM
  fit <- mgcv::gam(
    formula = response ~ s(time_points, bs = "cr", m = 2, k = basis_eff),
    family = family_obj,
    method = method,
    control = list(maxit = 500, mgcv.tol = 1e-4, epsilon = 1e-04),
    optimizer = c("outer", "bfgs")
  )

  # Get predicted probability and linear predictor (Z)
  prob <- stats::predict(fit, type = "response") # P(y=1)
  
  if (link == "probit") {
    # For probit, convert to logit scale (following original RunGam logic)
    z_probit <- stats::predict(fit, type = "link")
    prob_probit <- stats::pnorm(z_probit)
    # Avoid division by zero
    prob_probit <- pmax(prob_probit, 1e-10)
    prob_probit <- pmin(prob_probit, 1 - 1e-10)
    linpred <- -log(1/prob_probit - 1)  # Convert to logit scale
  } else {
    linpred <- stats::predict(fit, type = "link") # latent Z(t)
  }

  return(list(prob = prob, linpred = linpred))
}
