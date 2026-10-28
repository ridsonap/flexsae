#' Fay-Herriot Model with a Random Forest
#'
#' Fits the area-level model of Fay and Herriot (1979) with the fixed part
#' replaced or corrected by a random forest (Bukhari et al., 2026), and predicts
#' the area parameters with their mean squared error (MSE) estimated by a
#' parametric bootstrap.
#'
#' @param formula A two-sided formula with the direct estimates on the left
#'   and the covariates on the right, e.g. `y ~ x1 + x2`. The linear Fay-Herriot
#'   model uses the terms as written; the random forest uses the variables of
#'   the right-hand side untransformed. Smooth terms (`s()`) are not allowed.
#' @inheritParams pspline_sae
#' @param variant `"replace"` (default, FHRF-R of Bukhari et al., 2026): the
#'   random forest replaces the linear predictor. `"correct"` (FHRF-C): the
#'   random forest models what the linear predictor leaves, and the fixed part
#'   is \eqn{x_d^\prime\hat\beta + \hat f(x_d)}.
#' @param algorithm `"fhrf"` (default): the one-pass algorithm of Bukhari et
#'   al. (2026), where \eqn{\sigma^2_v} comes from the linear Fay-Herriot fit.
#'   `"merf"`: the random forest and the Fay-Herriot fit with the forest as an
#'   offset are alternated until convergence, so that \eqn{\sigma^2_v} is
#'   estimated given the forest (an area-level adaptation of Hajjem et al., 2014,
#'   and Krennmair and Schmid, 2022; see Details).
#' @param method Estimation method of \eqn{\sigma^2_v} in the Fay-Herriot fits:
#'   `"REML"` (default) or `"ML"`.
#' @param num.trees,mtry,min.node.size Passed to [ranger::ranger()]. `mtry =
#'   NULL` uses the ranger default (the square root of the number of
#'   covariates, rounded down).
#' @param maxit,tol For `algorithm = "merf"`: maximum number of iterations and
#'   convergence tolerance (see Details).
#'
#' @details
#' **Model.** For area \eqn{d = 1, \dots, D}, \eqn{y_d = \theta_d + e_d} with
#' \eqn{e_d \sim N(0, \psi_d)}, \eqn{\psi_d} known, and
#' \eqn{\theta_d = f(x_d) + v_d}, \eqn{v_d \sim N(0, \sigma^2_v)}. With
#' `variant = "replace"`, \eqn{f} is a random forest; with `variant =
#' "correct"`, \eqn{f(x_d) = x_d^\prime\beta + g(x_d)} with \eqn{g} a random
#' forest.
#'
#' **Algorithm `"fhrf"`** (Bukhari et al., 2026, Sections 2.3-2.4):
#' 1. Fit the linear Fay-Herriot model, giving \eqn{\hat\beta}, \eqn{\hat v_d} and
#'    \eqn{\hat\sigma^2_v}.
#' 2. Fit the forest to \eqn{y_d - \hat v_d} (`"replace"`, eq. 12) or to
#'    \eqn{y_d - x_d^\prime\hat\beta - \hat v_d} (`"correct"`, eq. 8) over the
#'    sampled areas, and take its out-of-bag predictions \eqn{\hat f(x_d)}
#'    (eqs. 9, 13). Unsampled areas get the ordinary forest prediction.
#' 3. \eqn{\hat\theta_d = \hat m_d + \gamma_d(y_d - \hat m_d)} with
#'    \eqn{\hat m_d = \hat f(x_d)} or \eqn{x_d^\prime\hat\beta + \hat f(x_d)}
#'    (column `synthetic`) and \eqn{\gamma_d = \hat\sigma^2_v/(\hat\sigma^2_v + \psi_d)}
#'    (eqs. 10, 14). Unsampled areas get \eqn{\hat m_d}.
#'
#' Bukhari et al. (2026) estimate \eqn{\sigma^2_v} with an EM algorithm. Its
#' fixed point is the maximum likelihood estimate, so `method = "ML"` gives the
#' same fit (checked numerically against `sae::eblupFH()`); the default is
#' REML, as in [pspline_sae()]. The Fay-Herriot fits maximize the profile
#' likelihood of \eqn{\sigma^2_v} directly (generalized least squares for
#' \eqn{\beta}); they agree with `sae::eblupFH()`.
#'
#' **Algorithm `"merf"`.** Starting from step 1, repeat: fit the forest to the
#' target of step 2 with the current \eqn{\hat v_d} (and \eqn{\hat\beta}) and take
#' its out-of-bag predictions; refit the Fay-Herriot model with the forest
#' predictions as an offset (and the linear terms for `"correct"`, no fixed
#' effects for `"replace"`), giving new \eqn{\hat\sigma^2_v}, \eqn{\hat v_d} (and
#' \eqn{\hat\beta}). Stop when the change of \eqn{\hat\sigma^2_v}, divided by
#' \eqn{\hat\sigma^2_v + \bar\psi}, is below `tol`. A forest is piecewise
#' constant, so the area estimates keep moving slightly from one iteration to
#' the next and are not used in the criterion. This follows the MERF algorithm of Hajjem et
#' al. (2014) and Krennmair and Schmid (2022), written for unit-level data,
#' with the Fay-Herriot fit in place of the unit-level mixed model. It is not
#' a published area-level method; it is offered because in `"fhrf"`
#' \eqn{\hat\sigma^2_v} also absorbs the non-linearity that the linear model
#' misses. All forests of one fit share one seed. This option is experimental.
#'
#' Out-of-bag predictions keep the forest from fitting the sampling errors of
#' the direct estimates. Every sampled area must be out of bag in at least one
#' tree, which needs enough trees.
#'
#' **MSE.** Naive parametric bootstrap (Bukhari et al., 2026, eqs. 15-17):
#' \eqn{\theta^*_d = \hat m_d + v^*_d}, \eqn{y^*_d = \theta^*_d + e^*_d}, the whole
#' algorithm, forest included, is refitted on \eqn{y^*} and
#' \eqn{\mathrm{mse}_d = B^{-1}\sum_b(\hat\theta^*_d - \theta^*_d)^2}. A replicate whose
#' refit fails is replaced, up to `5 * B` attempts.
#'
#' **Random numbers.** The forests use seeds drawn from the R random number
#' stream; with `seed` given the whole fit is reproducible and the stream of the
#' session is restored afterwards.
#'
#' @return An object of class `c("flexsae_merf", "flexsae")`: a list with
#'   * `df_est`: one row per area with `domain`, `direct`, `vardir`, `est`,
#'     `synthetic`, `random_effect`, `gamma`, `mse`, `rse` (\%), `ci_lower`,
#'     `ci_upper` and `sampled`. MSE columns are `NA` when `mse = "none"`.
#'   * `coefficients`: table of the coefficients of the linear part of the final
#'     predictor (`"correct"`); no rows for `"replace"`.
#'   * `rf`: the [ranger::ranger()] fit of the final forest.
#'   * `sigma2_v`, `fit` (the last linear Fay-Herriot fit: `sigma2_v`, `beta`,
#'     `vcov`, `method`),
#'     `formula`, `call`, `data`, `model`, `variant`, `algorithm`, `method`,
#'     `iterations`, `mse_method`, `B`, `B_success`, `n_boot_failed`,
#'     `boot_errors` (matrix of \eqn{\hat\theta^*_d - \theta^*_d}), `seed`,
#'     `conf_level` and `convergence`.
#'
#' @references
#' Bukhari, A. S., Notodiputro, K. A., Indahwati and Fitrianto, A. (2026). A
#' random forest-augmented Fay-Herriot model for nonlinear area-level small area
#' estimation. *Ingénierie des Systèmes d'Information*, 31(6), 1811-1823.
#' \doi{10.18280/isi.310604}
#'
#' Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5-32.
#' \doi{10.1023/A:1010933404324}
#'
#' Fay, R. E. and Herriot, R. A. (1979). Estimates of income for small places:
#' an application of James-Stein procedures to census data. *Journal of the
#' American Statistical Association*, 74(366a), 269-277.
#' \doi{10.1080/01621459.1979.10482505}
#'
#' Hajjem, A., Bellavance, F. and Larocque, D. (2014). Mixed-effects random
#' forest for clustered data. *Journal of Statistical Computation and
#' Simulation*, 84(6), 1313-1328. \doi{10.1080/00949655.2012.741599}
#'
#' Krennmair, P. and Schmid, T. (2022). Flexible domain prediction using mixed
#' effects random forests. *Journal of the Royal Statistical Society Series C*,
#' 71(5), 1865-1894. \doi{10.1111/rssc.12600}
#'
#' Wright, M. N. and Ziegler, A. (2017). ranger: a fast implementation of
#' random forests for high dimensional data in C++ and R. *Journal of
#' Statistical Software*, 77(1), 1-17. \doi{10.18637/jss.v077.i01}
#'
#' @export
#' @examples
#' set.seed(1)
#' D <- 100
#' dat <- data.frame(area = paste0("A", 1:D), x1 = runif(D, -2, 2), x2 = rnorm(D))
#' dat$vardir <- runif(D, 0.2, 0.6)
#' theta <- 2 + sin(1.5 * dat$x1) + 0.5 * dat$x2 + rnorm(D, sd = 0.5)
#' dat$y <- theta + rnorm(D, sd = sqrt(dat$vardir))
#' dat$y[c(5, 40)] <- NA # two unsampled areas
#'
#' fit <- merf_sae(y ~ x1 + x2, vardir = "vardir", domain = "area", data = dat,
#'                 num.trees = 200, B = 10, seed = 123)
#' fit
#' head(fit$df_est)
merf_sae <- function(
  formula,
  vardir,
  domain = NULL,
  data,
  variant = c("replace", "correct"),
  algorithm = c("fhrf", "merf"),
  method = c("REML", "ML"),
  num.trees = 500,
  mtry = NULL,
  min.node.size = 5,
  maxit = 25,
  tol = 1e-3,
  mse = c("bootstrap", "none"),
  B = 200,
  seed = NULL,
  conf_level = 0.95
) {
  variant <- match.arg(variant)
  algorithm <- match.arg(algorithm)
  method <- match.arg(method)
  mse <- match.arg(mse)
  .check_boot_args(mse, B, seed, conf_level)
  .check_count(num.trees, "num.trees")
  .check_count(min.node.size, "min.node.size")
  .check_count(maxit, "maxit")
  if (!is.null(mtry)) .check_count(mtry, "mtry")
  if (!is.numeric(tol) || length(tol) != 1 || !is.finite(tol) || tol <= 0) {
    cli::cli_abort("{.arg tol} must be a single positive number.", call = NULL)
  }

  ad <- .prepare_area_data(formula, vardir, domain, data)
  if (length(mgcv::interpret.gam(formula)$smooth.spec) > 0) {
    cli::cli_abort(c(
      "{.arg formula} must not contain smooth terms such as {.code s()}.",
      "i" = "Use {.fn pspline_sae} for spline effects."
    ), call = NULL)
  }
  rf_vars <- all.vars(formula[[3]])
  if (length(rf_vars) == 0) {
    cli::cli_abort("{.arg formula} needs at least one covariate for the random forest.", call = NULL)
  }
  if (!is.null(mtry) && mtry > length(rf_vars)) {
    cli::cli_abort("{.arg mtry} ({mtry}) cannot exceed the number of covariates ({length(rf_vars)}).", call = NULL)
  }
  rf_x <- data[, rf_vars, drop = FALSE]
  rf_x[] <- lapply(rf_x, function(v) if (is.character(v)) factor(v) else v)
  rf_args <- list(num.trees = num.trees, mtry = mtry, min.node.size = min.node.size)
  X <- .design_matrix(formula, data)

  fit_once <- function(y) {
    .merf_fit(y, X, rf_x, ad$vardir, ad$sampled, variant, algorithm, method, rf_args, maxit, tol)
  }

  .with_seed(seed, {
    res <- fit_once(ad$y)
    if (!res$converged) {
      cli::cli_warn(if (algorithm == "merf" && res$iterations >= maxit) {
        "The {.val merf} algorithm did not converge in {maxit} iteration{?s}."
      } else {
        "The search for the random effect variance did not converge."
      })
    }
    boot <- list(mse = rep(NA_real_, nrow(data)), B_success = NA_integer_, n_failed = NA_integer_, errors = NULL)
    if (mse == "bootstrap") {
      refit <- function(y_star) fit_once(y_star)$est
      boot <- .pboot_mse(res$synthetic, res$sigma2_v, ad$vardir, ad$sampled, refit, B)
    }
  })

  fh <- res$fh
  coefficients <- if (variant == "correct") fh$coefficients else fh$coefficients[0, , drop = FALSE]

  structure(
    list(
      df_est = .est_table(
        domain = ad$domain, direct = ad$y, vardir = ad$vardir, est = res$est,
        synthetic = res$synthetic, random_effect = res$random_effect,
        gamma = res$gamma, mse = boot$mse, sampled = ad$sampled, conf_level = conf_level
      ),
      coefficients = coefficients,
      rf = res$rf,
      sigma2_v = res$sigma2_v,
      fit = fh[c("sigma2_v", "beta", "vcov", "method")],
      formula = formula,
      call = match.call(),
      data = data,
      model = "merf_sae",
      variant = variant,
      algorithm = algorithm,
      method = method,
      iterations = res$iterations,
      mse_method = mse,
      B = if (mse == "bootstrap") as.integer(B) else NA_integer_,
      B_success = boot$B_success,
      n_boot_failed = boot$n_failed,
      boot_errors = boot$errors,
      seed = seed,
      conf_level = conf_level,
      convergence = res$converged
    ),
    class = c("flexsae_merf", "flexsae")
  )
}

#' One fit of merf_sae() on a response vector
#'
#' Used for the point fit and for every bootstrap replicate.
#'
#' @param y Response of all areas (NA for unsampled).
#' @param X Design matrix of the linear terms, all areas.
#' @param rf_x Forest covariates of all areas.
#' @return List with `est`, `synthetic`, `random_effect`, `gamma`, `sigma2_v`,
#'   `fh` (last linear Fay-Herriot fit), `rf`, `iterations` and `converged`.
#' @noRd
.merf_fit <- function(y, X, rf_x, vardir, sampled, variant, algorithm, method, rf_args, maxit, tol) {
  ys <- y[sampled]
  Xs <- X[sampled, , drop = FALSE]
  vs <- vardir[sampled]
  rf_seed <- sample.int(.Machine$integer.max, 1L)
  blup <- function(fixed, s2) ifelse(sampled, s2 / (s2 + vardir) * (y - fixed), 0)

  # step 1: linear Fay-Herriot fit
  fh <- .fit_fh_lm(ys, Xs, vs, method)
  linear <- drop(X %*% fh$beta)
  v <- blup(linear, fh$sigma2_v)

  converged <- fh$converged
  iterations <- 0L
  repeat {
    iterations <- iterations + 1L
    target <- if (variant == "correct") y - linear - v else y - v
    rf <- .fit_rf(rf_x, target, sampled, rf_args, rf_seed)
    fixed <- if (variant == "correct") linear + rf$f else rf$f
    if (algorithm == "fhrf") {
      break
    }
    # refit the Fay-Herriot model with the forest as an offset
    s2_old <- fh$sigma2_v
    if (variant == "correct") {
      fh <- .fit_fh_lm(ys - rf$f[sampled], Xs, vs, method)
      linear <- drop(X %*% fh$beta)
      fixed <- linear + rf$f
    } else {
      fh <- .fit_fh_lm(ys - rf$f[sampled], Xs[, 0, drop = FALSE], vs, method)
    }
    v <- blup(fixed, fh$sigma2_v)
    if (abs(fh$sigma2_v - s2_old) / (s2_old + mean(vs)) < tol) {
      converged <- fh$converged
      break
    }
    if (iterations >= maxit) {
      converged <- FALSE
      break
    }
  }

  .fhrf_combine(fixed, y, fh$sigma2_v, vardir, sampled,
                list(fh = fh, rf = rf$fit, iterations = iterations, converged = converged))
}

#' EBLUP step of Bukhari et al. (2026), eqs. (10) and (14)
#'
#' theta_d = m_d + gamma_d (y_d - m_d) for sampled areas, m_d for unsampled.
#' @noRd
.fhrf_combine <- function(fixed, y, sigma2_v, vardir, sampled, extra = list()) {
  gamma <- ifelse(sampled, sigma2_v / (sigma2_v + vardir), NA_real_)
  random_effect <- ifelse(sampled, gamma * (y - fixed), 0)
  c(
    list(
      est = fixed + random_effect,
      synthetic = fixed,
      random_effect = random_effect,
      gamma = gamma,
      sigma2_v = sigma2_v
    ),
    extra
  )
}

#' Random forest on the sampled areas with out-of-bag predictions
#'
#' @return List with `fit` (ranger object) and `f`: out-of-bag predictions for
#'   the sampled areas and ordinary predictions for the unsampled ones.
#' @noRd
.fit_rf <- function(rf_x, target, sampled, rf_args, rf_seed) {
  fit <- ranger::ranger(
    x = rf_x[sampled, , drop = FALSE], y = target[sampled],
    num.trees = rf_args$num.trees, mtry = rf_args$mtry,
    min.node.size = rf_args$min.node.size, seed = rf_seed, num.threads = 1
  )
  f <- numeric(length(sampled))
  f[sampled] <- fit$predictions
  if (anyNA(f[sampled])) {
    cli::cli_abort(c(
      "{sum(is.na(f[sampled]))} sampled area{?s} {?was/were} never out of bag.",
      "i" = "Increase {.arg num.trees}."
    ), call = NULL)
  }
  if (any(!sampled)) {
    f[!sampled] <- stats::predict(fit, data = rf_x[!sampled, , drop = FALSE], num.threads = 1)$predictions
  }
  list(fit = fit, f = f)
}

#' Check a positive whole-number argument
#' @noRd
.check_count <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 1 || x != round(x)) {
    cli::cli_abort("{.arg {arg}} must be a single positive whole number.", call = NULL)
  }
  invisible(TRUE)
}
