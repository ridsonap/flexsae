# flexsae 0.0.0.9000

* New `pspline_sae()`: semiparametric Fay-Herriot model with penalized splines,
  fitted by REML or ML with 'mgcv' (known sampling variances as weights with
  the scale fixed at 1, area effect as a random-effect smooth), synthetic
  prediction for unsampled areas, and naive parametric bootstrap MSE. Without
  smooth terms it agrees with `sae::eblupFH()` and `fastsae::eblup_fh()`.
* New `merf_sae()`: Fay-Herriot model with a random forest for the fixed part
  (Bukhari et al., 2026), `variant = "replace"` (FHRF-R) or `"correct"`
  (FHRF-C), out-of-bag forest predictions, and parametric bootstrap MSE that
  refits the forest. `algorithm = "fhrf"` (default) is the one-pass algorithm of
  the paper; `algorithm = "merf"` (experimental) alternates the forest and the
  Fay-Herriot fit, an area-level adaptation of MERF that is not published.
  The random effect variance is estimated by REML by default; the paper uses an
  EM algorithm whose fixed point is the ML estimate (`method = "ML"`).
* The Fay-Herriot fits of `merf_sae()` maximize the profile likelihood of the
  random effect variance directly instead of using 'mgcv': about 400 times
  faster at 300 areas (0.006 s vs 2.6 s per fit), and it agrees with
  `sae::eblupFH()` to about 1e-7.
* New `benchmark_flexsae()`: ratio, difference, optimal and logit
  benchmarking to known, estimated (`target_var`) or same-survey direct
  (`target = "direct"`) targets, by group and in two stages
  (`national_target`), with the MSE of the benchmarked estimates from the model
  bootstrap (each replicate benchmarked to its bootstrap target). The rules
  agree with 'fastsae', 'emdi' and 'saebenchmarking'; the MSE increase of
  difference benchmarking agrees with Steorts and Ghosh (2013).
* Self-benchmarking in `pspline_sae()` and `merf_sae()` (`bench_weight`,
  `bench_group`): the augmented model of Wang, Fuller and Qu (2008, eq. 23)
  makes the weighted aggregate of the estimates equal that of the direct
  estimates over the sampled areas, by group, without a post-hoc adjustment.
  For `merf_sae()` the augmentation is applied in the final EBLUP step with the
  forest as an offset, an extension not in Bukhari et al. (2026).
* New `select_flexsae()`: stepwise selection on BIC, AIC or KIC from the
  marginal ML likelihood for `pspline_sae()` (equal to `fastsae::step_sae()`
  for linear terms), and Boruta selection for `merf_sae()`.
* New `importance_flexsae()` (permutation or impurity importance of the
  forest of `merf_sae()`), `partial_flexsae()` (partial dependence and ICE
  curves of the fixed part, with intervals for `pspline_sae()`), and
  `autoplot()` methods for models and benchmarked estimates.
* New `scfh_sae()` (experimental): spatially clustered Fay-Herriot model of
  Maranzano, Mattera and Sugasawa (2026, preprint). Areas fall into `K`
  regimes with their own coefficients and random effect variance; the
  partition is estimated with a Potts penalty `phi` on the neighbourhood graph
  `W`. Cluster variances by adjusted REML (Li and Lahiri, 2010; new
  `method = "AREML"` of the internal Fay-Herriot solver, agreeing with
  `emdi::fh(method = "amrl")`), and bootstrap MSE that reruns the clustering.
  Unsampled areas take their neighbours' cluster, an extension not in the
  paper. `select_flexsae(method = "grid")` chooses `K` and `phi` with the
  two-step rule of the paper, and `autoplot(type = "clusters")` shows the
  coefficients by cluster. Self-benchmarking is not available for this model.
* New `boost_sae()` (experimental): Fay-Herriot model with an XGBoost fixed
  part, our own area-level adaptation of the scheme of `merf_sae()` with the
  learner of the unit-level MEGB of Messer and Schmid (2026). Sampled areas
  get cross-fitted predictions (the counterpart of out-of-bag predictions);
  `nrounds` is chosen by cross-validation with early stopping and kept for the
  bootstrap. One-pass and iterative algorithms, `"replace"` and `"correct"`
  variants, self-benchmarking, gain importance, partial dependence and
  plots. 'xgboost' is a suggested package.
* `merf_sae()` and `boost_sae()` share one fitting engine; `merf_sae()`
  results are unchanged (checked against results of the previous code).
* Vignettes "Getting started with flexsae", "Benchmarking area estimates" and
  "Simulation study: flexsae against the linear Fay-Herriot model", and a
  pkgdown configuration. The simulation code is in `benchmarks/simulation/`.
* Model objects keep the bootstrap truth and direct estimates of every
  replicate (`boot_theta`, `boot_y`).
* `print()`, `summary()`, `coef()`, `fitted()` and `residuals()` methods for
  the common `"flexsae"` result class.
* Package skeleton, plan (`PLAN.md`) and verified references
  (`technical-report/references.bib`).
