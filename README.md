
<!-- README.md is generated from README.Rmd. Please edit that file -->

# flexsae

<!-- badges: start -->

[![R-CMD-check](https://github.com/ridsonap/flexsae/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ridsonap/flexsae/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**flexsae** is an R package, under development, for small area
estimation with area-level data (direct estimates with known sampling
variances) when the covariates act non-linearly.

Available:

- `pspline_sae()`: semiparametric Fay-Herriot model with penalized
  splines (fitted by REML or ML with ‘mgcv’), with parametric bootstrap
  MSE and optional self-benchmarking to the direct estimates.

- `merf_sae()`: Fay-Herriot model with a random forest for the fixed
  part (replacing or correcting the linear predictor), with parametric
  bootstrap MSE.
- `boost_sae()` (experimental): Fay-Herriot model with gradient boosting
  (XGBoost) for the fixed part, with cross-fitted predictions and parametric
  bootstrap MSE.
- `scfh_sae()` (experimental): spatially clustered Fay-Herriot model,
  where neighbouring areas tend to share a regime with its own
  coefficients, with a bootstrap MSE that reruns the clustering.

- `benchmark_flexsae()`: benchmarking to aggregate targets (known,
  estimated or the direct estimates of the same survey; by group or in
  two stages), with the bootstrap MSE of the benchmarked estimates.

- `select_flexsae()`, `importance_flexsae()`, `partial_flexsae()` and
  `autoplot()` methods.

Planned for the first release: simulation studies, vignettes and a
pkgdown site.

Every method is checked against the linear Fay-Herriot model of the
`sae`, `emdi` and `fastsae` packages, and in simulation studies that
also report the settings where the flexible models do worse.

## Installation

``` r
# install.packages("pak")
pak::pak("ridsonap/flexsae")
```
