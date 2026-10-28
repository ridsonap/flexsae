#' Semiparametric Fay-Herriot Model with Penalized Splines
#'
#' Fits the area-level model of Fay and Herriot (1979) with smooth, non-linear
#' covariate effects represented by penalized splines, and predicts the area
#' parameters with their mean squared error (MSE) estimated by a parametric
#' bootstrap.
#'
#' @param formula A two-sided formula in \pkg{mgcv} syntax: the direct
#'   estimates on the left, linear terms and smooth terms such as `s(x1)`,
#'   `s(x1, k = 5)` or `s(x1, bs = "ps")` on the right. The area random effect is
#'   added internally and must not be written in the formula. A formula
#'   without smooth terms gives the linear Fay-Herriot model.
#' @param vardir Sampling variances of the direct estimates: a column name, a
#'   one-sided formula or a numeric vector. They must be positive for the
#'   sampled areas and are treated as known.
#' @param domain Optional area labels: a column name, a one-sided formula or a
#'   vector. If `NULL`, the areas are numbered by row.
#' @param data A data frame with one row per area. Areas with a missing direct
#'   estimate are unsampled and receive the synthetic prediction; they need
#'   covariate values.
#' @param method Estimation method of the variance components and smoothing
#'   parameters: `"REML"` (default) or `"ML"`.
#' @param mse `"bootstrap"` (default) for the parametric bootstrap MSE, or
#'   `"none"`.
#' @param B Number of bootstrap replicates. Default 200.
#' @param seed Optional seed for the bootstrap. The random number stream of the
#'   session is restored afterwards.
#' @param conf_level Confidence level of the normal intervals
#'   \eqn{\hat\theta_d \pm z\sqrt{\mathrm{mse}_d}}. Default 0.95.
#'
#' @details
#' **Model.** For area \eqn{d = 1, \dots, D}, the direct estimate is
#' \eqn{y_d = \theta_d + e_d} with \eqn{e_d \sim N(0, \psi_d)} and \eqn{\psi_d} known,
#' and \deqn{\theta_d = \beta_0 + \sum_k f_k(x_{dk}) + x_d^{\prime}\beta + v_d,
#' \quad v_d \sim N(0, \sigma^2_v).}
#' Each \eqn{f_k} is a penalized regression spline. A spline penalty is
#' equivalent to a random effect on the spline coefficients, so the model is a
#' linear mixed model with one variance component per smooth term plus
#' \eqn{\sigma^2_v} (Ruppert, Wand and Carroll, 2003; Wand, 2003). It is fitted by
#' [mgcv::gam()] with weights \eqn{1/\psi_d}, the scale fixed at 1 and the area
#' effect as an `s(area, bs = "re")` term, by REML or ML (Wood, 2011, 2017); then
#' \eqn{\hat\sigma^2_v} is the scale divided by the smoothing parameter of the area
#' term. Without smooth terms, the fit agrees with the linear Fay-Herriot fit of
#' \pkg{sae} (`eblupFH()`) to the convergence tolerance of the optimizers.
#'
#' The basis of the smooths is the one chosen in `s()` (thin plate regression
#' splines by default; `bs = "ps"` for P-splines). With one observation per area
#' the number of areas limits how flexible the smooths can be: lower `k` when
#' \eqn{D} is small.
#'
#' **Prediction.** \eqn{\hat\theta_d = \hat m(x_d) + \hat v_d}, where \eqn{\hat m} is
#' the fitted fixed and smooth part (column `synthetic`) and
#' \eqn{\hat v_d = \gamma_d\{y_d - \hat m(x_d)\}} with
#' \eqn{\gamma_d = \hat\sigma^2_v/(\hat\sigma^2_v + \psi_d)}. Unsampled areas get
#' \eqn{\hat v_d = 0}.
#'
#' **MSE.** Naive parametric bootstrap: for \eqn{b = 1, \dots, B},
#' \eqn{\theta^*_d = \hat m(x_d) + v^*_d} with \eqn{v^*_d \sim N(0, \hat\sigma^2_v)} and
#' \eqn{y^*_d = \theta^*_d + e^*_d} with \eqn{e^*_d \sim N(0, \psi_d)}; the model,
#' including its smoothing parameters, is refitted on \eqn{y^*} and
#' \eqn{\mathrm{mse}_d = B^{-1}\sum_b(\hat\theta^*_d - \theta^*_d)^2} (Hall and Maiti,
#' 2006; the same scheme as Bukhari et al., 2026, eqs. 15-17). The bias
#' correction of Hall and Maiti (2006) is not applied. A replicate whose refit
#' fails is replaced, up to `5 * B` attempts.
#'
#' @return An object of class `c("flexsae_pspline", "flexsae")`: a list with
#'   * `df_est`: one row per area with `domain`, `direct`, `vardir`, `est`,
#'     `synthetic`, `random_effect`, `gamma`, `mse`, `rse` (\%), `ci_lower`,
#'     `ci_upper` and `sampled`. MSE columns are `NA` when `mse = "none"`.
#'   * `coefficients`: table of the parametric coefficients.
#'   * `smooth`: table of the smooth terms (effective degrees of freedom and
#'     approximate tests from [mgcv::summary.gam()]).
#'   * `sigma2_v`, `fit` (the [mgcv::gam()] object), `formula`, `call`, `data`,
#'     `model`, `method`, `mse_method`, `B`, `B_success`, `n_boot_failed`,
#'     `boot_errors` (matrix of \eqn{\hat\theta^*_d - \theta^*_d}), `seed`,
#'     `conf_level` and `convergence`.
#'
#' @references
#' Bukhari, A. S., Notodiputro, K. A., Indahwati and Fitrianto, A. (2026). A
#' random forest-augmented Fay-Herriot model for nonlinear area-level small area
#' estimation. *Ingénierie des Systèmes d'Information*, 31(6), 1811-1823.
#' \doi{10.18280/isi.310604}
#'
#' Fay, R. E. and Herriot, R. A. (1979). Estimates of income for small places:
#' an application of James-Stein procedures to census data. *Journal of the
#' American Statistical Association*, 74(366a), 269-277.
#' \doi{10.1080/01621459.1979.10482505}
#'
#' Hall, P. and Maiti, T. (2006). On parametric bootstrap methods for small area
#' prediction. *Journal of the Royal Statistical Society Series B*, 68(2),
#' 221-238. \doi{10.1111/j.1467-9868.2006.00541.x}
#'
#' Ruppert, D., Wand, M. P. and Carroll, R. J. (2003). *Semiparametric
#' Regression*. Cambridge University Press. \doi{10.1017/CBO9780511755453}
#'
#' Wand, M. P. (2003). Smoothing and mixed models. *Computational Statistics*,
#' 18(2), 223-249. \doi{10.1007/s001800300142}
#'
#' Wood, S. N. (2011). Fast stable restricted maximum likelihood and marginal
#' likelihood estimation of semiparametric generalized linear models. *Journal
#' of the Royal Statistical Society Series B*, 73(1), 3-36.
#' \doi{10.1111/j.1467-9868.2010.00749.x}
#'
#' Wood, S. N. (2017). *Generalized Additive Models: An Introduction with R*
#' (2nd ed.). Chapman and Hall/CRC. \doi{10.1201/9781315370279}
#'
#' @export
#' @examples
#' set.seed(1)
#' D <- 60
#' dat <- data.frame(area = paste0("A", 1:D), x1 = runif(D, -2, 2), x2 = rnorm(D))
#' dat$vardir <- runif(D, 0.2, 0.6)
#' theta <- 2 + sin(1.5 * dat$x1) + 0.5 * dat$x2 + rnorm(D, sd = 0.5)
#' dat$y <- theta + rnorm(D, sd = sqrt(dat$vardir))
#' dat$y[c(5, 40)] <- NA # two unsampled areas
#'
#' fit <- pspline_sae(y ~ s(x1, k = 6) + x2, vardir = "vardir", domain = "area",
#'                    data = dat, B = 20, seed = 123)
#' fit
#' head(fit$df_est)
pspline_sae <- function(
  formula,
  vardir,
  domain = NULL,
  data,
  method = c("REML", "ML"),
  mse = c("bootstrap", "none"),
  B = 200,
  seed = NULL,
  conf_level = 0.95
) {
  method <- match.arg(method)
  mse <- match.arg(mse)
  .check_boot_args(mse, B, seed, conf_level)

  ad <- .prepare_area_data(formula, vardir, domain, data)
  fh <- .fit_fh_gam(ad$formula_fit, ad$fit_data, method)
  if (!fh$converged) {
    cli::cli_warn("The smoothing parameter optimization of {.fn mgcv::gam} did not fully converge ({.val {fh$fit$outer.info$conv}}).")
  }
  pred <- .predict_fh_gam(fh, ad$data, ad$sampled)
  est <- pred$synthetic + pred$random_effect
  gamma <- ifelse(ad$sampled, fh$sigma2_v / (fh$sigma2_v + ad$vardir), NA_real_)

  boot <- list(mse = rep(NA_real_, length(est)), B_success = NA_integer_, n_failed = NA_integer_, errors = NULL)
  if (mse == "bootstrap") {
    refit <- function(y_star) {
      fd <- ad$fit_data
      fd[[.resp_var]] <- y_star[ad$sampled]
      fh_b <- .fit_fh_gam(ad$formula_fit, fd, method)
      p <- .predict_fh_gam(fh_b, ad$data, ad$sampled)
      p$synthetic + p$random_effect
    }
    boot <- .with_seed(seed, .pboot_mse(pred$synthetic, fh$sigma2_v, ad$vardir, ad$sampled, refit, B))
  }

  sfit <- summary(fh$fit)
  smooth <- as.data.frame(sfit$s.table)
  smooth <- smooth[rownames(smooth) != paste0("s(", .area_var, ")"), , drop = FALSE]

  structure(
    list(
      df_est = .est_table(
        domain = ad$domain, direct = ad$y, vardir = ad$vardir, est = est,
        synthetic = pred$synthetic, random_effect = pred$random_effect,
        gamma = gamma, mse = boot$mse, sampled = ad$sampled, conf_level = conf_level
      ),
      coefficients = as.data.frame(sfit$p.table),
      smooth = smooth,
      sigma2_v = fh$sigma2_v,
      fit = fh$fit,
      formula = formula,
      call = match.call(),
      data = data,
      model = "pspline_sae",
      method = method,
      mse_method = mse,
      B = if (mse == "bootstrap") as.integer(B) else NA_integer_,
      B_success = boot$B_success,
      n_boot_failed = boot$n_failed,
      boot_errors = boot$errors,
      seed = seed,
      conf_level = conf_level,
      convergence = fh$converged
    ),
    class = c("flexsae_pspline", "flexsae")
  )
}

#' Check the bootstrap arguments shared by the model functions
#' @noRd
.check_boot_args <- function(mse, B, seed, conf_level) {
  if (mse == "bootstrap" &&
    (!is.numeric(B) || length(B) != 1 || !is.finite(B) || B < 1 || B != round(B))) {
    cli::cli_abort("{.arg B} must be a single positive whole number.", call = NULL)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed))) {
    cli::cli_abort("{.arg seed} must be {.code NULL} or a single number.", call = NULL)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1 || !is.finite(conf_level) ||
    conf_level <= 0 || conf_level >= 1) {
    cli::cli_abort("{.arg conf_level} must be a single number in (0, 1).", call = NULL)
  }
  invisible(TRUE)
}
