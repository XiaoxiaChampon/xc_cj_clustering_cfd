# catfda: Categorical Functional Data Analysis

[![R-CMD-check](https://github.com/username/catfda/workflows/R-CMD-check/badge.svg)](https://github.com/username/catfda/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/catfda)](https://CRAN.R-project.org/package=catfda)
[![codecov](https://codecov.io/gh/username/catfda/branch/main/graph/badge.svg)](https://codecov.io/gh/username/catfda)

## Overview

The `catfda` package provides comprehensive tools for analyzing categorical functional data. It implements various methods for estimating latent Gaussian processes and probability curves from categorical time series data using different link functions.

## Features

- **Multiple estimation methods**: Probit, binomial, and multinomial approaches
- **Automatic category detection**: Handles both numeric and character category labels
- **Parallel processing**: Efficient computation for large datasets
- **Flexible basis functions**: Customizable smoothing parameters
- **Data generation tools**: Utilities for simulation studies

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("username/catfda")
```

Or install from CRAN (once published):

```r
install.packages("catfda")
```

## Quick Start

```r
library(catfda)

# Generate sample data
set.seed(123)
n_time <- 50
n_individuals <- 10
time_points <- seq(0, 1, length.out = n_time)
w_mat <- matrix(sample(1:3, n_time * n_individuals, replace = TRUE),
                nrow = n_time, ncol = n_individuals)

# Estimate using multinomial approach
result <- estimate_categ_func_data("multinomial", time_points, w_mat)

# Check output
names(result)
dim(result$p1_est)  # Probability estimates for category 1
```

## Methods

### 1. Multinomial Approach
Uses a single GAM with multinomial family to model all categories simultaneously:

```r
result_multi <- estimate_categ_func_data("multinomial", time_points, w_mat)
```

### 2. Binomial Approach
Fits separate binomial GAMs for each category:

```r
result_binom <- estimate_categ_func_data("binomial", time_points, w_mat)
```

### 3. Probit Approach
Uses adaptive switching between probit and binomial links for rare events:

```r
result_probit <- estimate_categ_func_data("probit", time_points, w_mat)
```

## Key Functions

- `estimate_categ_func_data()`: Main estimation function with method selection
- `estimate_categ_func_data_multinomial()`: Multinomial GAM estimation
- `estimate_categ_func_data_probit()`: Adaptive probit/binomial estimation
- `estimate_categ_func_data_binomial_parallel()`: Parallel binomial estimation
- `get_x_from_w()`: Convert categorical data to one-hot encoding
- `generate_categ_func_data()`: Generate synthetic categorical functional data

## Parallel Processing

For large datasets, use parallel versions:

```r
library(doParallel)
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Parallel estimation
result <- estimate_categ_func_data_probit_parallel(time_points, x_array)

stopCluster(cl)
```

## Example: Visualization

```r
# Extract probability estimates
p1 <- result$p1_est[, 1]  # First individual, category 1
p2 <- result$p2_est[, 1]  # First individual, category 2
p3 <- result$p3_est[, 1]  # First individual, category 3

# Plot
plot(time_points, p1, type = "l", col = "red", lwd = 2, ylim = c(0, 1),
     xlab = "Time", ylab = "Probability")
lines(time_points, p2, col = "blue", lwd = 2)
lines(time_points, p3, col = "green", lwd = 2)
legend("topright", legend = paste("Category", 1:3),
       col = c("red", "blue", "green"), lwd = 2)
```

## Documentation

- Browse the vignettes: `browseVignettes("catfda")`
- Get help: `help(package = "catfda")`
- Function reference: `?estimate_categ_func_data`

## Contributing

Contributions are welcome! Please see our [contributing guidelines](CONTRIBUTING.md) and [code of conduct](CODE_OF_CONDUCT.md).

## License

This package is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Citation

If you use this package in your research, please cite:

```
@Manual{catfda,
  title = {catfda: Categorical Functional Data Analysis},
  author = {Your Name and Co-Author Name},
  year = {2025},
  note = {R package version 0.1.0},
  url = {https://github.com/username/catfda},
}
```

## References

- [Detailed methodology paper]
- [Application examples]
- [Theoretical background]
