#!/usr/bin/env Rscript

# CATFDA Package Build and Check Script
# =====================================

cat("CATFDA Package Build and Check\n")
cat("==============================\n\n")

# Load required packages
if (!require(devtools, quietly = TRUE)) {
  install.packages("devtools")
  library(devtools)
}

if (!require(roxygen2, quietly = TRUE)) {
  install.packages("roxygen2")
  library(roxygen2)
}

# Set working directory to package root
setwd("catfda")

cat("1. Documenting package...\n")
devtools::document()

cat("\n2. Building package...\n")
build_result <- devtools::build()

cat("\n3. Running R CMD check...\n")
check_result <- devtools::check()

cat("\n4. Installing package...\n")
devtools::install()

cat("\n5. Loading and testing package...\n")
library(catfda)

# Quick functionality test
cat("Testing basic functionality:\n")

# Generate sample data
time_points <- seq(0, 1, length.out = 50)
n_individuals <- 10
n_categories <- 3

# Create sample 3D array
x_array <- array(0, dim = c(n_individuals, length(time_points), n_categories))

# Fill with random one-hot vectors
for(i in 1:n_individuals) {
  for(t in 1:length(time_points)) {
    cat_choice <- sample(1:n_categories, 1)
    x_array[i, t, cat_choice] <- 1
  }
}

cat("- Testing probit estimation...\n")
tryCatch({
  result_probit <- estimate_categ_func_data(time_points, x_array, link = "probit")
  cat("  ✓ Probit estimation successful\n")
}, error = function(e) {
  cat("  ✗ Probit estimation failed:", e$message, "\n")
})

cat("- Testing binomial estimation...\n")
tryCatch({
  result_binomial <- estimate_categ_func_data(time_points, x_array, link = "binomial")
  cat("  ✓ Binomial estimation successful\n")
}, error = function(e) {
  cat("  ✗ Binomial estimation failed:", e$message, "\n")
})

cat("- Testing multinomial estimation...\n")
tryCatch({
  result_multinomial <- estimate_categ_func_data(time_points, x_array, link = "multinomial")
  cat("  ✓ Multinomial estimation successful\n")
}, error = function(e) {
  cat("  ✗ Multinomial estimation failed:", e$message, "\n")
})

cat("\n6. Package check summary:\n")
cat("Build result:", build_result, "\n")
if (length(check_result$errors) > 0) {
  cat("ERRORS found:", length(check_result$errors), "\n")
  for (error in check_result$errors) {
    cat("  -", error, "\n")
  }
}
if (length(check_result$warnings) > 0) {
  cat("WARNINGS found:", length(check_result$warnings), "\n")
  for (warning in check_result$warnings) {
    cat("  -", warning, "\n")
  }
}
if (length(check_result$notes) > 0) {
  cat("NOTES found:", length(check_result$notes), "\n")
  for (note in check_result$notes) {
    cat("  -", note, "\n")
  }
}

if (length(check_result$errors) == 0 && length(check_result$warnings) == 0) {
  cat("\n✓ Package passed all checks and is ready for CRAN submission!\n")
} else {
  cat("\n✗ Package has issues that need to be addressed before CRAN submission.\n")
}

cat("\nBuild and check complete.\n")
