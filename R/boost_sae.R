#' Fay-Herriot Model with Gradient Boosting
#'
#' @description
#' **Experimental.** Area-level small area estimation with a gradient boosting
#' fixed part (XGBoost): the Fay-Herriot model
#' \eqn{y_d = m(x_d) + v_d + e_d}, \eqn{v_d \sim N(0, \sigma^2_v)},
#' \eqn{e_d \sim N(0, \psi_d)} with \eqn{\psi_d} known, where \eqn{m} is
#' learned by boosted regression trees. This is our own area-level adaptation:
#' the scheme of [merf_sae()] (Bukhari et al., 2026) with the random forest
#' replaced by the XGBoost learner of the unit-level mixed effects gradient
#' boosting (MEGB) of Messer and Schmid (2026). No published area-level
#' version was found.
#'
#' @inheritParams merf_sae
#' @param algorithm `"onepass"` (default): linear Fay-Herriot fit, boosting on
#'   its target once, EBLUP with the variance of the linear fit (as
#'   `merf_sae(algorithm = "fhrf")`). `"iterative"`: boosting and the
#'   Fay-Herriot fit with the boosting predictions as an offset alternate until
#'   the variance changes by less than `tol` relative to
#'   \eqn{\hat\sigma^2_v + \bar\psi} (the area-level counterpart of Algorithm 1
#'   of Messer and Schmid, 2026).
#' @param nrounds Number of boosting rounds. `NULL` (default) chooses it by
#'   `nfolds`-fold cross-validation with early stopping (up to 1000 rounds,
#'   stop after 20 rounds without improvement), as in MEGB.
#' @param eta,max_depth,subsample,min_child_weight XGBoost parameters
#'   (learning rate, tree depth, row subsampling, minimum child weight).
#' @param nfolds Number of folds for the cross-fitted predictions and for the
#'   choice of `nrounds`.
#'
#' @details
#' **Cross-fitting.** Boosting has no out-of-bag predictions, and in-sample
#' predictions have already seen \eqn{y_d}: the residuals \eqn{y_d - \hat m_d}
#' would shrink and the random effect variance and shrinkage factors would be
#' biased. The sampled areas are therefore split once into `nfolds` folds and
#' each area is predicted by the booster trained on the other folds, the
#' counterpart of the out-of-bag predictions of [merf_sae()]. Unsampled areas
#' are predicted by the booster trained on all sampled areas, which is stored
#' in `booster` and used by [importance_flexsae()] and [partial_flexsae()].
#'
#' The booster is trained without weights, as the forest of [merf_sae()] and
#' the booster of MEGB, on the variables of the right-hand side of `formula`
#' (factors one-hot encoded, no transformations). `nrounds` is chosen once on
#' the target of the first pass of the original fit and kept, with the folds,
#' for the bootstrap replicates; the replicates refit everything else. These
#' are our own choices (`technical-report/DISCREPANCIES.md`, BST-1 to BST-3).
#'
#' Self-benchmarking (`bench_weight`, `bench_group`) works as in [merf_sae()].
#'
#' @return An object of class `c("flexsae_boost", "flexsae")` with the fields
#'   of [merf_sae()] objects, the booster in `booster`, its target and seed in
#'   `boost_target` and `boost_seed`, and the boosting settings (including the
#'   chosen `nrounds` and the folds) in `boost_args`.
#'
#' @references
#' Messer, P. and Schmid, T. (2026). Gradient boosting for hierarchical data
#' in small area estimation. *Statistics and Computing*, 36(5), 218.
#' \doi{10.1007/s11222-026-10970-1}
#'
#' Chen, T. and Guestrin, C. (2016). XGBoost: A scalable tree boosting system.
#' In *Proceedings of the 22nd ACM SIGKDD International Conference on
#' Knowledge Discovery and Data Mining*, 785-794.
#' \doi{10.1145/2939672.2939785}
#'
#' Bukhari, A. S., Notodiputro, K. A., Indahwati and Fitrianto, A. (2026). A
#' random forest-augmented Fay-Herriot model for nonlinear area-level small
#' area estimation. *Ingenierie des Systemes d'Information*, 31(6),
#' 1811-1823. \doi{10.18280/isi.310604}
#'
#' @export
#' @examplesIf requireNamespace("xgboost", quietly = TRUE)
#' set.seed(1)
#' D <- 100
#' dat <- data.frame(x1 = runif(D, -2, 2), x2 = rnorm(D))
#' dat$vardir <- runif(D, 0.2, 0.6)
#' dat$y <- 2 + sin(1.5 * dat$x1) + 0.5 * dat$x2 + rnorm(D, sd = 0.5) +
#'   rnorm(D, sd = sqrt(dat$vardir))
#' fit <- boost_sae(y ~ x1 + x2, vardir = "vardir", data = dat, mse = "none",
#'                  seed = 1)
#' fit
#' importance_flexsae(fit)
boost_sae <- function(
  formula,
  vardir,
  domain = NULL,
  data,
  variant = c("replace", "correct"),
  algorithm = c("onepass", "iterative"),
  method = c("REML", "ML"),
  nrounds = NULL,
  eta = 0.05,
  max_depth = 2,
  subsample = 0.8,
  min_child_weight = 1,
  nfolds = 5,
  bench_weight = NULL,
  bench_group = NULL,
  maxit = 25,
  tol = 1e-3,
  mse = c("bootstrap", "none"),
  B = 200,
  seed = NULL,
  conf_level = 0.95
) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn boost_sae} needs the {.pkg xgboost} package.",
      "i" = "Install it with {.code install.packages(\"xgboost\")}."
    ), call = NULL)
  }
  variant <- match.arg(variant)
  algorithm <- match.arg(algorithm)
  method <- match.arg(method)
  mse <- match.arg(mse)
  .check_boot_args(mse, B, seed, conf_level)
  .check_count(max_depth, "max_depth")
  .check_count(maxit, "maxit")
  .check_count(nfolds, "nfolds")
  if (!is.null(nrounds)) .check_count(nrounds, "nrounds")
  check_range <- function(x, arg, lower, upper, lower_open) {
    ok <- is.numeric(x) && length(x) == 1 && is.finite(x) && x <= upper &&
      (if (lower_open) x > lower else x >= lower)
    if (!ok) {
      cli::cli_abort("{.arg {arg}} must be a single number in {if (lower_open) '(' else '['}{lower}, {upper}].",
                     call = NULL)
    }
  }
  check_range(eta, "eta", 0, 1, TRUE)
  check_range(subsample, "subsample", 0, 1, TRUE)
  check_range(min_child_weight, "min_child_weight", 0, Inf, FALSE)
  check_range(tol, "tol", 0, Inf, TRUE)

  ad <- .prepare_area_data(formula, vardir, domain, data)
  if (length(mgcv::interpret.gam(formula)$smooth.spec) > 0) {
    cli::cli_abort(c(
      "{.arg formula} must not contain smooth terms such as {.code s()}.",
      "i" = "Use {.fn pspline_sae} for spline effects."
    ), call = NULL)
  }
  if (length(all.vars(formula[[3]])) == 0) {
    cli::cli_abort("{.arg formula} needs at least one covariate for the boosting.", call = NULL)
  }
  n_s <- sum(ad$sampled)
  if (nfolds < 2 || nfolds > n_s) {
    cli::cli_abort("{.arg nfolds} must be between 2 and the number of sampled areas ({n_s}).", call = NULL)
  }
  xmat <- .boost_matrix(formula, data)
  X <- .design_matrix(formula, data)
  sb <- .self_benchmark_cols(bench_weight, bench_group, data, ad$vardir, ad$sampled)
  params <- list(objective = "reg:squarederror", eta = eta, max_depth = max_depth, subsample = subsample,
                 min_child_weight = min_child_weight, nthread = 1)

  .with_seed(seed, {
    folds <- sample(rep_len(seq_len(nfolds), n_s))
    nrounds_chosen <- is.null(nrounds)
    if (nrounds_chosen) {
      nrounds <- .boost_nrounds(ad$y, X, xmat, ad$vardir, ad$sampled, variant, method, params, folds)
    }
    learner <- function(target, seed) .fit_boost(xmat, target, ad$sampled, params, nrounds, folds, seed)
    fit_once <- function(y) {
      .merf_fit(y, X, learner, ad$vardir, ad$sampled, variant, algorithm == "iterative", method, maxit, tol, sb$Z)
    }
    res <- fit_once(ad$y)
    if (!is.null(sb)) .check_self_benchmark(sb, ad$y, res$est, ad$sampled)
    if (!res$converged) {
      cli::cli_warn(if (algorithm == "iterative" && res$iterations >= maxit) {
        "The {.val iterative} algorithm did not converge in {maxit} iteration{?s}."
      } else {
        "The search for the random effect variance did not converge."
      })
    }
    boot <- list(mse = rep(NA_real_, nrow(data)), B_success = NA_integer_, n_failed = NA_integer_,
                 errors = NULL, theta = NULL, y = NULL)
    if (mse == "bootstrap") {
      refit <- function(y_star) fit_once(y_star)$est
      boot <- .pboot_mse(res$synthetic, res$sigma2_v, ad$vardir, ad$sampled, refit, B)
    }
  })

  fh <- res$fh
  coefficients <- if (variant == "correct") fh$coefficients else fh$coefficients[0, , drop = FALSE]
  if (!is.null(sb)) coefficients <- rbind(coefficients, .rename_bench_rows(res$bench, sb))

  structure(
    list(
      df_est = .est_table(
        domain = ad$domain, direct = ad$y, vardir = ad$vardir, est = res$est,
        synthetic = res$synthetic, random_effect = res$random_effect,
        gamma = res$gamma, mse = boot$mse, sampled = ad$sampled, conf_level = conf_level
      ),
      coefficients = coefficients,
      booster = res$learner_fit,
      boost_target = res$learner_target,
      boost_seed = res$learner_seed,
      boost_args = list(params = params, nrounds = as.integer(nrounds), nrounds_chosen = nrounds_chosen,
                        nfolds = as.integer(nfolds), folds = folds, columns = colnames(xmat),
                        variables = all.vars(formula[[3]]), assign = attr(xmat, "assign")),
      sigma2_v = res$sigma2_v,
      fit = fh[c("sigma2_v", "beta", "vcov", "method")],
      formula = formula,
      call = match.call(),
      data = data,
      model = "boost_sae",
      variant = variant,
      algorithm = algorithm,
      method = method,
      self_benchmark = !is.null(sb),
      bench = sb,
      iterations = res$iterations,
      mse_method = mse,
      B = if (mse == "bootstrap") as.integer(B) else NA_integer_,
      B_success = boot$B_success,
      n_boot_failed = boot$n_failed,
      boot_errors = boot$errors,
      boot_theta = boot$theta,
      boot_y = boot$y,
      seed = seed,
      conf_level = conf_level,
      convergence = res$converged
    ),
    class = c("flexsae_boost", "flexsae")
  )
}

#' Boosting covariates: the right-hand side variables, factors one-hot
#'
#' @param ref_data Optional data whose factor levels are used (new data for
#'   partial dependence).
#' @return Numeric matrix with column names and the `assign` attribute of
#'   model.matrix() (column -> variable).
#' @noRd
.boost_matrix <- function(formula, data, ref_data = NULL) {
  vars <- all.vars(formula[[3]])
  f <- stats::reformulate(vars, intercept = FALSE)
  to_factor <- function(d) {
    d <- d[, vars, drop = FALSE]
    d[] <- lapply(d, function(v) if (is.character(v) || is.logical(v)) factor(v) else v)
    d
  }
  # all levels of every factor, so that no level is the baseline of another
  ref <- to_factor(ref_data %||% data)
  contr <- lapply(Filter(is.factor, ref), function(v) stats::contrasts(v, contrasts = FALSE))
  tt <- stats::terms(f)
  xlev <- stats::.getXlevels(tt, stats::model.frame(tt, ref))
  mf <- stats::model.frame(tt, to_factor(data), xlev = xlev)
  stats::model.matrix(tt, mf, contrasts.arg = if (length(contr) > 0) contr)
}

#' Booster with cross-fitted predictions for the sampled areas
#'
#' @param folds Fold of each sampled area.
#' @return `list(fit, f)`: booster trained on all sampled areas, and
#'   predictions (cross-fitted for the sampled areas).
#' @noRd
.fit_boost <- function(xmat, target, sampled, params, nrounds, folds, seed) {
  s <- which(sampled)
  p <- c(params, list(seed = seed))
  train <- function(idx) {
    dm <- xgboost::xgb.DMatrix(xmat[idx, , drop = FALSE], label = target[idx], nthread = 1)
    xgboost::xgb.train(params = p, data = dm, nrounds = nrounds, verbose = 0)
  }
  f <- numeric(length(sampled))
  for (k in sort(unique(folds))) {
    test <- s[folds == k]
    f[test] <- stats::predict(train(s[folds != k]), xmat[test, , drop = FALSE])
  }
  fit <- train(s)
  if (any(!sampled)) {
    f[!sampled] <- stats::predict(fit, xmat[!sampled, , drop = FALSE])
  }
  list(fit = fit, f = f)
}

#' Number of boosting rounds by cross-validation with early stopping
#'
#' On the target of the first pass of the fit (the one of `.merf_fit()`).
#' @noRd
.boost_nrounds <- function(y, X, xmat, vardir, sampled, variant, method, params, folds) {
  s <- which(sampled)
  fh <- .fit_fh_lm(y[s], X[s, , drop = FALSE], vardir[s], method)
  linear <- drop(X %*% fh$beta)
  v <- fh$sigma2_v / (fh$sigma2_v + vardir[s]) * (y[s] - linear[s])
  target <- if (variant == "correct") y[s] - linear[s] - v else y[s] - v
  dm <- xgboost::xgb.DMatrix(xmat[s, , drop = FALSE], label = target, nthread = 1)
  p <- c(params, list(seed = sample.int(.Machine$integer.max, 1L)))
  cv <- xgboost::xgb.cv(
    params = p, data = dm, nrounds = 1000, folds = split(seq_along(s), folds),
    early_stopping_rounds = 20, verbose = FALSE
  )
  max(as.integer(cv$early_stop$best_iteration), 1L)
}
