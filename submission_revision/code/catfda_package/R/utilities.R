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
