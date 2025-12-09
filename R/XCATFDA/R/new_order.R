#' Estimate categorical functional data using the specified link function
#'
#' @param choice A string: one of "probit", "binomial", or "multinomial"
#' @param time_points A numeric vector of time points (0 to 1)
#' @param w_mat A matrix of categorical observations (time × individuals)
#' @param n_basis Integer: number of basis functions (default = 25)
#' @param method Estimation method for GAM (default = "ML")
#'
#' @return A list of latent variables (Z) and predicted probabilities (p)
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
    stop(sprintf("w_mat Must be a matrix. Given value : '%s'", w_mat))
  }

  return(method_map[[choice]]())
}

#' Estimate Z and p curves for categorical functional data (multinomial case)
#'
#' @param time_points A numeric vector of time points
#' @param w_mat A matrix of categorical values (time × individuals)
#' @param n_basis Number of basis functions for smoothing (default = 25)
#' @param method GAM fitting method (default = "ML")
#'
#' @return A named list with latent Z estimates and probability estimates
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

  # >>>>>> Which is correct?
  # Remap labels to 1...K if needed
  # w_mapped <- apply(w_mat, 2, function(col) match(col, unique_labels))

  # Remap labels to 0...(K-1) for mgcv::multinom (expects 0-based indexing)
  w_mapped <- apply(w_mat, 2, function(col) match(col, unique_labels) - 1)

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
    paste0("Z", seq_along(z_list), "_est")
  )
  p_estimates <- rlang::set_names(
    lapply(seq_len(k_classes), function(k) t(p_array[, , k])),
    paste0("p", seq_len(k_classes), "_est")
  )

  return(c(z_estimates, p_estimates))
}

#' Parallel estimation of Z and p using multinomial link GAM
#'
#' @param time_points A numeric vector of time points
#' @param w_mat A matrix of categorical values (time × individuals)
#' @param n_basis Number of basis functions for smoothing (default = 25)
#' @param method GAM fitting method (default = "ML")
#'
#' @return A named list with latent Z estimates and probability estimates
#' @export
estimate_categ_func_data_multinomial_parallel <- function(time_points,
                                                          w_mat,
                                                          n_basis = 25,
                                                          method = "ML") {
  stopifnot(is.numeric(time_points), is.matrix(w_mat))

  n_individuals <- ncol(w_mat)
  n_timepoints <- nrow(w_mat)

  # Determine number of categories (K)
  unique_labels <- sort(unique(as.vector(w_mat)))
  k_classes <- length(unique_labels)

  # Remap labels to 0...(K-1) for mgcv::multinom (expects 0-based indexing)
  w_mapped <- apply(w_mat, 2, function(col) match(col, unique_labels) - 1)

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
    paste0("Z", seq_len(k_classes - 1), "_est")
  )

  # Return list of p1_est, ..., pK_est
  p_array <- array(t(p_mat), dim = c(n_individuals, n_timepoints, k_classes))
  p_estimates <- rlang::set_names(
    lapply(seq_len(k_classes), function(k) t(p_array[, , k])),
    paste0("p", seq_len(k_classes), "_est")
  )

  return(c(z_estimates, p_estimates))
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

#' Fit a GAM to a binary time series with a chosen link function
#'
#' @param time_points A numeric vector of time values (length T)
#' @param response A binary vector of observations (length T)
#' @param link A string: either "binomial" (logit) or "probit"
#' @param n_basis Number of basis functions (default = 25)
#' @param method GAM optimization method (default = "ML")
#'
#' @return A list with predicted probabilities and latent Z (linear predictors)
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
    binomial(link = "probit")
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
  prob <- predict(fit, type = "response") # P(y=1)
  linpred <- predict(fit, type = "link") # latent Z(t)

  return(list(prob = prob, linpred = linpred))
}


#' Estimate Z and p curves using (generalized) probit or binomial models
#'
#' @param time_points A numeric vector of time values (length T)
#' @param x_array A 3D array (individual × time × category) of one-hot encoded data
#' @param n_basis Number of basis functions for smoothing (default = 25)
#' @param method GAM fitting method (default = "ML")
#' @param threshold_probability Probability threshold for switching to probit special case (default = 0.004)
#'
#' @return A list containing estimated Z curves and category probabilities for each individual
#' @export
estimate_categ_func_data_probit <- function(time_points,
                                            x_array,
                                            n_basis = 25,
                                            method = "ML",
                                            threshold_probability = 0.004) {
  stopifnot(length(dim(x_array)) == 3)

  n_individuals <- dim(x_array)[1]
  n_timepoints <- dim(x_array)[2]
  n_categories <- dim(x_array)[3]

  # Initialize storage
  z_curves <- vector("list", n_categories - 1)
  prob_array <- array(0, dim = c(n_individuals, n_timepoints, n_categories))

  for (indv in seq_len(n_individuals)) {
    p_list <- vector("list", n_categories)
    linpred_list <- vector("list", n_categories)

    # Determine which link to use
    use_probit <- n_timepoints <= 301 &&
      mean(x_array[indv, , 1]) < threshold_probability

    link_func <- if (use_probit) "probit" else "binomial"

    # Fit models for each category
    for (k in seq_len(n_categories)) {
      response <- x_array[indv, , k]
      fit <- run_gam(time_points, response, link_func, n_basis, method)
      p_list[[k]] <- fit$prob
      linpred_list[[k]] <- fit$linpred
    }

    # Compute latent variables Zₖ = f(linear predictors)
    ref_index <- n_categories
    ref_linpred <- linpred_list[[ref_index]]

    for (k in seq_len(n_categories - 1)) {
      num <- linpred_list[[k]] - ref_linpred
      denom <- log((1 + exp(linpred_list[[k]])) / (1 + exp(ref_linpred)))
      z_k <- num - denom
      z_curves[[k]] <- cbind(z_curves[[k]], z_k)
    }

    # Normalize probabilities across categories
    p_vecs <- do.call(cbind, p_list) # t × K
    p_sum <- rowSums(p_vecs)
    p_normalized <- sweep(p_vecs, 1, p_sum, FUN = "/")
    prob_array[indv, , ] <- p_normalized
  }

  # Build output: z_curves already time × individuals, prob_array needs transpose
  z_out <- rlang::set_names(
    z_curves,
    paste0("Z", seq_along(z_curves), "_est")
  )
  p_out <- rlang::set_names(
    lapply(seq_len(n_categories), function(k) t(prob_array[, , k])),
    paste0("p", seq_len(n_categories), "_est")
  )

  return(c(z_out, p_out))
}

#' Parallel estimation of Z and p using binomial link GAMs
#'
#' @param time_points A numeric vector of time points (length T)
#' @param x_array A 3D array: individual × time × category (dimensions N × T × K)
#' @param n_basis Number of basis functions for GAM (default = 25)
#' @param method GAM fitting method (default = "ML")
#' @param threshold_probability If response is sparse, switch to probit (default = 0.004)
#'
#' @return A list of estimated latent curves (Z1, Z2, ..., ZK-1) and normalized probabilities (p1, ..., pK)
#' @export
estimate_categ_func_data_probit_parallel <- function(time_points,
                                                     x_array,
                                                     n_basis = 25,
                                                     method = "ML",
                                                     threshold_probability = 0.004) {
  stopifnot(length(dim(x_array)) == 3)

  n_individuals <- dim(x_array)[1]
  n_timepoints <- dim(x_array)[2]
  n_categories <- dim(x_array)[3]

  total_z_curves <- n_categories - 1
  total_z_rows <- n_timepoints * total_z_curves

  # Parallel estimation
  zp <- foreach(
    indv = seq_len(n_individuals),
    .combine = cbind,
    .packages = c("mgcv"),
    .export = c("run_gam")
  ) %dorng% {
    # Estimate p_k and linpred for each category
    p_list <- list()
    linpred_list <- list()

    for (k in seq_len(n_categories)) {
      x_k <- x_array[indv, , k]

      use_probit <- (n_timepoints <= 301) && (mean(x_k) < threshold_probability)
      link <- if (use_probit) "probit" else "binomial"

      fit <- run_gam(time_points, x_k, link, n_basis, method)
      p_list[[k]] <- fit$prob
      linpred_list[[k]] <- fit$linpred
    }

    # Compute Z_k = f(linear predictors)
    z_values <- list()
    ref_linpred <- linpred_list[[n_categories]]
    for (k in seq_len(n_categories - 1)) {
      num <- linpred_list[[k]] - ref_linpred
      denom <- log((1 + exp(linpred_list[[k]])) / (1 + exp(ref_linpred)))
      z_k <- num - denom
      z_values[[k]] <- z_k
    }

    # Normalize probabilities
    p_mat <- do.call(cbind, p_list) # T × K
    p_sum <- rowSums(p_mat)
    p_norm <- sweep(p_mat, 1, p_sum, FUN = "/") # normalized probs

    # Flatten and return combined result
    return(c(
      unlist(z_values), # (T * (K - 1)) values
      as.vector(p_norm) # (T * K) values
    ))
  }

  # Reconstruct output arrays
  z_mat <- zp[1:total_z_rows, , drop = FALSE]
  p_mat <- zp[(total_z_rows + 1):nrow(zp), , drop = FALSE]

  # Create return list of Z
  z_out <- rlang::set_names(
    lapply(seq_len(total_z_curves), function(k) {
      z_start <- (k - 1) * n_timepoints + 1
      z_end <- k * n_timepoints
      t(z_mat[z_start:z_end, ])
    }),
    paste0("z", seq_len(total_z_curves), "_est")
  )

  # Create return list of p
  p_array <- array(t(p_mat), dim = c(n_individuals, n_timepoints, n_categories))
  p_out <- rlang::set_names(
    lapply(seq_len(n_categories), function(k) t(p_array[, , k])),
    paste0("p", seq_len(n_categories), "_est")
  )

  return(c(z_out, p_out))
}

#' Estimate Z and p curves using binomial link GAMs (parallel version)
#'
#' @param time_points Numeric vector of time values (length T)
#' @param x_array 3D array (individual × time × category) — binary one-hot encoding
#' @param n_basis Number of basis functions (default = 25)
#' @param method GAM optimization method (default = "ML")
#'
#' @return A list with:
#'   - Z1_est, ..., Z{K-1}_est (latent curves)
#'   - p1_est, ..., pK_est (category probability curves)
#' @export
estimate_categ_func_data_binomial_parallel <- function(time_points,
                                                       x_array,
                                                       n_basis = 25,
                                                       method = "ML") {
  stopifnot(length(dim(x_array)) == 3)

  n_individuals <- dim(x_array)[1]
  n_timepoints <- dim(x_array)[2]
  n_categories <- dim(x_array)[3]

  total_z_curves <- n_categories - 1
  total_z_rows <- n_timepoints * total_z_curves

  zp <- foreach(
    indv = seq_len(n_individuals),
    .combine = cbind,
    .packages = c("mgcv"),
    .export = c("run_gam")
  ) %dorng% {
    p_list <- list()
    linpred_list <- list()

    for (k in seq_len(n_categories)) {
      x_k <- x_array[indv, , k]
      fit <- run_gam(time_points, x_k, link = "binomial", n_basis = n_basis, method = method)
      p_list[[k]] <- fit$prob
      linpred_list[[k]] <- fit$linpred
    }

    # Reference = last category
    ref_linpred <- linpred_list[[n_categories]]

    z_values <- list()
    for (k in seq_len(n_categories - 1)) {
      num <- linpred_list[[k]] - ref_linpred
      denom <- log((1 + exp(linpred_list[[k]])) / (1 + exp(ref_linpred)))
      z_k <- num - denom
      z_values[[k]] <- z_k
    }

    # Normalize probabilities
    p_mat <- do.call(cbind, p_list) # T × K
    p_sum <- rowSums(p_mat)
    p_norm <- sweep(p_mat, 1, p_sum, FUN = "/") # T × K

    return(c(
      unlist(z_values), # (T * (K - 1))
      as.vector(p_norm) # (T * K)
    ))
  }

  # Split output back into Z and p
  z_mat <- zp[1:total_z_rows, , drop = FALSE]
  p_mat <- zp[(total_z_rows + 1):nrow(zp), , drop = FALSE]

  # Return list of Z1_est, ..., Z{K-1}_est - already in time × individuals format
  z_out <- rlang::set_names(
    lapply(seq_len(total_z_curves), function(k) {
      z_start <- (k - 1) * n_timepoints + 1
      z_end <- k * n_timepoints
      z_mat[z_start:z_end, ]
    }),
    paste0("Z", seq_len(total_z_curves), "_est")
  )

  # Return list of p1_est, ..., pK_est
  # p_mat is (time*categories) × individuals, need to reshape to time × individuals per category
  p_array <- array(t(p_mat), dim = c(n_individuals, n_timepoints, n_categories))
  p_out <- rlang::set_names(
    lapply(seq_len(n_categories), function(k) t(p_array[, , k])),
    paste0("p", seq_len(n_categories), "_est")
  )

  return(c(z_out, p_out))
}



#' Generate synthetic categorical functional data (W) and one-hot encoded array (X)
#'
#' @param prob_curves A named list of T × N matrices: p1_est, ..., pK_est
#'
#' @return A list with:
#'   - w_mat: Categorical matrix (T × N)
#'   - x_array: One-hot encoded array (N × T × K)
#' @export
generate_categ_func_data <- function(prob_curves) {
  stopifnot(is.list(prob_curves), length(prob_curves) >= 2)

  k_categories <- length(prob_curves)
  category_names <- names(prob_curves)

  # Validate structure
  ref_dim <- dim(prob_curves[[1]])
  stopifnot(all(sapply(prob_curves, function(mat) all(dim(mat) == ref_dim))))

  n_timepoints <- ref_dim[1]
  n_individuals <- ref_dim[2]

  # Prepare storage
  w_mat <- matrix(0, nrow = n_timepoints, ncol = n_individuals)
  x_array <- array(0, dim = c(n_individuals, n_timepoints, k_categories))

  for (indv in seq_len(n_individuals)) {
    p_matrix <- sapply(prob_curves, function(p_mat) p_mat[, indv]) # T × K

    # For each time point, sample from the categorical distribution
    sampled <- apply(p_matrix, 1, function(prob_vec) {
      rmultinom(1, size = 1, prob = prob_vec)
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
