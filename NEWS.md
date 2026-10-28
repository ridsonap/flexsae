# flexsae 0.0.0.9000

* New `pspline_sae()`: semiparametric Fay-Herriot model with penalized splines,
  fitted by REML or ML with 'mgcv' (known sampling variances as weights with
  the scale fixed at 1, area effect as a random-effect smooth), synthetic
  prediction for unsampled areas, and naive parametric bootstrap MSE. Without
  smooth terms it agrees with `sae::eblupFH()` and `fastsae::eblup_fh()`.
* `print()`, `summary()`, `coef()`, `fitted()` and `residuals()` methods for
  the common `"flexsae"` result class.
* Package skeleton, plan (`PLAN.md`) and verified references
  (`technical-report/references.bib`).
