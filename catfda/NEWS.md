# catfda 0.1.0

## Initial Release

### New Features

* Initial implementation of categorical functional data analysis methods
* Three main estimation approaches:
  - Multinomial GAM estimation (`estimate_categ_func_data_multinomial()`)
  - Adaptive probit/binomial estimation (`estimate_categ_func_data_probit()`)
  - Parallel binomial estimation (`estimate_categ_func_data_binomial_parallel()`)
* Main dispatcher function (`estimate_categ_func_data()`) for method selection
* Utility functions:
  - `get_x_from_w()`: Convert categorical to one-hot encoding
  - `run_gam()`: Low-level GAM fitting with different link functions
  - `generate_categ_func_data()`: Generate synthetic categorical functional data
* Parallel processing support using `foreach` and `doRNG`
* Automatic category detection and handling of both numeric and character labels
* Comprehensive test suite with >90% code coverage
* Detailed documentation and vignettes

### Key Capabilities

* Handles arbitrary number of categories (not limited to 3)
* Supports different link functions (probit, binomial, multinomial)
* Adaptive threshold-based switching for rare events
* Efficient parallel computation for large datasets
* Flexible basis function specification
* Robust error handling and input validation

### Dependencies

* R (>= 4.0.0)
* mgcv (>= 1.8-40): GAM fitting
* rlang (>= 1.0.0): Programming utilities
* foreach (>= 1.5.0): Parallel loops
* doRNG (>= 1.8.0): Reproducible parallel computation
* parallel: Core parallel processing

### Documentation

* Complete function documentation with examples
* Introduction vignette with comprehensive examples
* README with quick start guide
* NEWS file documenting changes
