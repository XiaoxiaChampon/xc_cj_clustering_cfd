#' Estimate categorical functional data using the specified link function
#'
#' @description
#' Main dispatcher function for categorical functional data estimation.
#' Chooses between probit, binomial, or multinomial approaches based on the
#' specified method.
#'
#' @param choice A string: one of "probit", "binomial", or "multinomial"
#' @param time_points A numeric vector of time points (0 to 1)
#' @param w_mat A matrix of categorical observations (time × individuals)
#' @param n_basis Integer: number of basis functions (default = 25)
#' @param method Estimation method for GAM (default = "ML")
#'
#' @return A list of latent variables (Z) and predicted probabilities (p)
#'
#' @details
#' This function serves as the main entry point for categorical functional data
#' analysis. It automatically dispatches to the appropriate estimation method
#' based on the specified choice:
#' \itemize{
#'   \item "probit": Uses adaptive probit/binomial switching for rare events
#'   \item "binomial": Uses separate binomial GAMs for each category
#'   \item "multinomial": Uses multinomial GAM for all categories simultaneously
#' }
#'
#' @examples
#' # Generate sample data
#' set.seed(123)
#' n_time <- 50
#' n_individuals <- 10
#' time_points <- seq(0, 1, length.out = n_time)
#'
#' # Create sample categorical data
#' w_mat <- matrix(sample(0:2, n_time * n_individuals, replace = TRUE),
#'                 nrow = n_time, ncol = n_individuals)
#'
#' # Estimate using multinomial approach
#' result <- estimate_categ_func_data("multinomial", time_points, w_mat, n_basis = 10)
#'
#' @export
estimate_categ_func_data <- function(choice,
                                     time_points,
                                     w_mat,
                                     n_basis = 25,
                                     method = "ML") {
  method_map <- list(
    probit = function() {
      x <- get_x_from_w(w_mat)
      estimate_categ_func_data_probit_parallel(time_points, x, n_basis, method, 1 / 150)
    },
    binomial = function() {
      x <- get_x_from_w(w_mat)
      estimate_categ_func_data_binomial_parallel(time_points, x, n_basis, method)
    },
    multinomial = function() {
      estimate_categ_func_data_multinomial_parallel(time_points, w_mat, n_basis, method)
    }
  )

  # Validate and dispatch
  if (!choice %in% names(method_map)) {
    stop(sprintf(
      "Invalid choice '%s'. Must be one of: %s", choice,
      paste(names(method_map), collapse = ", ")
    ))
  }

  if (!is.numeric(time_points)) {
    stop(sprintf(
      "time_points Must be a numeric. Given value : '%s'",
      time_points
    ))
  }

  if (!is.matrix(w_mat)) {
    stop(sprintf(
      "w_mat must be a matrix. Given value type : '%s'",
      typeof(w_mat)
    ))
  }

  return(method_map[[choice]]())
}
