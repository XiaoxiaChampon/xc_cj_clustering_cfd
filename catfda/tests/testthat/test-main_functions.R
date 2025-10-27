test_that("estimate_categ_func_data works with valid inputs", {
  # Generate test data
  set.seed(123)
  n_time <- 20
  n_individuals <- 5
  time_points <- seq(0, 1, length.out = n_time)
  # Use 0-based categories for multinomial
  w_mat <- matrix(sample(0:2, n_time * n_individuals, replace = TRUE),
                  nrow = n_time, ncol = n_individuals)

  # Test multinomial method with appropriate basis size
  result_multi <- estimate_categ_func_data("multinomial", time_points, w_mat, n_basis = 10)

  expect_type(result_multi, "list")
  expect_true(all(c("z1_est", "z2_est", "p1_est", "p2_est", "p3_est") %in% names(result_multi)))
  expect_equal(dim(result_multi$z1_est), c(n_time, n_individuals))
  expect_equal(dim(result_multi$p1_est), c(n_time, n_individuals))

  # Test probit method
  x_array <- get_x_from_w(w_mat)
  result_probit <- estimate_categ_func_data("probit", time_points, w_mat)

  expect_type(result_probit, "list")
  expect_true(all(c("z1_est", "z2_est", "p1_est", "p2_est", "p3_est") %in% names(result_probit)))

  # Test binomial method
  result_binom <- estimate_categ_func_data("binomial", time_points, w_mat)

  expect_type(result_binom, "list")
  expect_true(all(c("z1_est", "z2_est", "p1_est", "p2_est", "p3_est") %in% names(result_binom)))
})

test_that("estimate_categ_func_data validates inputs correctly", {
  # Test invalid choice
  expect_error(
    estimate_categ_func_data("invalid", seq(0, 1, length.out = 10), matrix(1:20, 10, 2)),
    "Invalid choice"
  )

  # Test non-numeric time_points
  expect_error(
    estimate_categ_func_data("multinomial", "not_numeric", matrix(1:20, 10, 2)),
    "time_points Must be a numeric"
  )

  # Test non-matrix w_mat
  expect_error(
    estimate_categ_func_data("multinomial", seq(0, 1, length.out = 10), "not_matrix"),
    "w_mat Must be a matrix"
  )
})

test_that("get_x_from_w creates correct one-hot encoding", {
  # Simple test case
  w_mat <- matrix(c(1, 2, 3, 1, 2, 3), nrow = 3, ncol = 2)
  x_array <- get_x_from_w(w_mat)

  expect_equal(dim(x_array), c(2, 3, 3))  # individuals, time, categories

  # Check that each time point has exactly one active category
  for (i in 1:2) {
    for (t in 1:3) {
      expect_equal(sum(x_array[i, t, ]), 1)
    }
  }

  # Check specific values
  expect_equal(x_array[1, 1, 1], 1)  # First individual, first time, category 1
  expect_equal(x_array[1, 2, 2], 1)  # First individual, second time, category 2
  expect_equal(x_array[1, 3, 3], 1)  # First individual, third time, category 3
})

test_that("run_gam produces valid output", {
  set.seed(123)
  n_time <- 30
  time_points <- seq(0, 1, length.out = n_time)
  response <- rbinom(n_time, 1, 0.5)

  # Test binomial link
  result_binom <- run_gam(time_points, response, "binomial")
  expect_type(result_binom, "list")
  expect_true(all(c("prob", "linpred") %in% names(result_binom)))
  expect_equal(length(result_binom$prob), n_time)
  expect_equal(length(result_binom$linpred), n_time)
  expect_true(all(result_binom$prob >= 0 & result_binom$prob <= 1))

  # Test probit link
  result_probit <- run_gam(time_points, response, "probit")
  expect_type(result_probit, "list")
  expect_true(all(c("prob", "linpred") %in% names(result_probit)))
  expect_true(all(result_probit$prob >= 0 & result_probit$prob <= 1))
})

test_that("generate_categ_func_data creates valid synthetic data", {
  # Create probability curves
  set.seed(123)
  n_time <- 15
  n_individuals <- 8

  p1 <- matrix(runif(n_time * n_individuals, 0.2, 0.6), n_time, n_individuals)
  p2 <- matrix(runif(n_time * n_individuals, 0.1, 0.4), n_time, n_individuals)
  p3 <- 1 - p1 - p2
  p3 <- pmax(0.1, p3)  # Ensure positive

  # Normalize
  total <- p1 + p2 + p3
  p1 <- p1 / total
  p2 <- p2 / total
  p3 <- p3 / total

  prob_curves <- list(p1_est = p1, p2_est = p2, p3_est = p3)

  result <- generate_categ_func_data(prob_curves)

  expect_type(result, "list")
  expect_true(all(c("w_mat", "x_array") %in% names(result)))
  expect_equal(dim(result$w_mat), c(n_time, n_individuals))
  expect_equal(dim(result$x_array), c(n_individuals, n_time, 3))

  # Check that w_mat contains valid category indices
  expect_true(all(result$w_mat %in% 1:3))

  # Check that x_array is properly one-hot encoded
  for (i in 1:n_individuals) {
    for (t in 1:n_time) {
      expect_equal(sum(result$x_array[i, t, ]), 1)
    }
  }
})
