#' Partial Dependence of the Fixed Part
#'
#' Partial dependence (Friedman, 2001) of the fixed part \eqn{\hat m(x)} of a
#' fitted model on each covariate, optionally with individual conditional
#' expectation (ICE) curves (Goldstein et al., 2015).
#'
#' @param object A model from [pspline_sae()] or [merf_sae()].
#' @param vars Covariates (names of variables in the formula); default all.
#' @param n Number of grid points for a numeric covariate (between its
#'   minimum and maximum over the areas).
#' @param ice If `TRUE`, also return the ICE curves of every area.
#' @param conf_level Confidence level of the pointwise intervals
#'   (`pspline_sae()` models only).
#'
#' @details For a grid value \eqn{u} of covariate \eqn{x_j}, the partial
#' dependence is \eqn{D^{-1}\sum_d \hat m(u, x_{d,-j})}, the average over all
#' areas of the fixed part with \eqn{x_j} set to \eqn{u} and the other
#' covariates at their observed values; the area random effects are not
#' included. For `pspline_sae()` models the fixed part is linear in the
#' coefficients, \eqn{\bar a(u)^\prime\hat\beta}, so a pointwise interval follows
#' from the Bayesian covariance matrix `Vp` of [mgcv::gam()]:
#' \eqn{\bar a(u)^\prime\hat\beta \pm z\{\bar a(u)^\prime V_p \bar a(u)\}^{1/2}}.
#' For an additive model this is the smooth term shifted by a constant. For
#' `merf_sae()` models the fixed part is the forest prediction (plus the linear
#' part for `variant = "correct"`); no interval is given. Augmentation columns
#' of self-benchmarking stay at their observed values. Partial dependence
#' averages over combinations of covariates that may not occur in the data
#' when covariates are correlated (Apley and Zhu, 2020).
#'
#' @return A data frame with columns `variable`, `x` (grid value, numeric
#'   covariates), `level` (factor covariates), `fit`, `se`, `lower`, `upper`.
#'   With `ice = TRUE`, attribute `ice` is a data frame with columns
#'   `variable`, `x`, `level`, `domain` and `fit`.
#'
#' @references
#' Apley, D. W. and Zhu, J. (2020). Visualizing the effects of predictor
#' variables in black box supervised learning models. *Journal of the Royal
#' Statistical Society Series B*, 82(4), 1059-1086. \doi{10.1111/rssb.12377}
#'
#' Friedman, J. H. (2001). Greedy function approximation: a gradient boosting
#' machine. *The Annals of Statistics*, 29(5), 1189-1232.
#' \doi{10.1214/aos/1013203451}
#'
#' Goldstein, A., Kapelner, A., Bleich, J. and Pitkin, E. (2015). Peeking
#' inside the black box: visualizing statistical learning with plots of
#' individual conditional expectation. *Journal of Computational and
#' Graphical Statistics*, 24(1), 44-65. \doi{10.1080/10618600.2014.907095}
#'
#' @export
#' @examples
#' set.seed(1)
#' D <- 60
#' dat <- data.frame(x1 = runif(D, -2, 2), x2 = rnorm(D))
#' dat$vardir <- runif(D, 0.2, 0.6)
#' dat$y <- 2 + sin(1.5 * dat$x1) + 0.5 * dat$x2 + rnorm(D, sd = 0.5) +
#'   rnorm(D, sd = sqrt(dat$vardir))
#' fit <- pspline_sae(y ~ s(x1, k = 6) + x2, vardir = "vardir", data = dat,
#'                    mse = "none")
#' head(partial_flexsae(fit, "x1"))
partial_flexsae <- function(object, vars = NULL, n = 50, ice = FALSE, conf_level = 0.95) {
  if (!inherits(object, "flexsae") || identical(object$model, "scfh_sae")) {
    cli::cli_abort(c(
      "{.arg object} must be a model from {.fn pspline_sae} or {.fn merf_sae}.",
      "i" = if (inherits(object, "flexsae_scfh")) "For {.fn scfh_sae} see {.code autoplot(object, type = \"clusters\")}."
    ), call = NULL)
  }
  all_vars <- all.vars(object$formula[[3]])
  vars <- vars %||% all_vars
  bad <- setdiff(vars, all_vars)
  if (length(bad) > 0) {
    cli::cli_abort("{.val {bad}} {?is/are} not covariate{?s} of the model.", call = NULL)
  }
  .check_count(n, "n")
  if (!is.numeric(conf_level) || length(conf_level) != 1 || conf_level <= 0 || conf_level >= 1) {
    cli::cli_abort("{.arg conf_level} must be a single number in (0, 1).", call = NULL)
  }
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  newdata <- object$data
  if (!is.null(object$bench)) {
    for (nm in colnames(object$bench$Z)) newdata[[nm]] <- object$bench$Z[, nm]
  }
  pred <- .fixed_part_fun(object)

  pd <- list()
  ic <- list()
  for (v in vars) {
    col <- newdata[[v]]
    numeric_var <- is.numeric(col)
    grid <- if (numeric_var) {
      seq(min(col), max(col), length.out = n)
    } else if (is.factor(col)) {
      levels(col)
    } else {
      sort(unique(as.character(col)))
    }
    res <- lapply(grid, function(u) {
      nd <- newdata
      nd[[v]] <- if (is.factor(col)) factor(rep(u, nrow(nd)), levels = levels(col)) else rep(u, nrow(nd))
      pred(nd)
    })
    fit <- vapply(res, function(r) r$fit, numeric(1))
    se <- vapply(res, function(r) r$se, numeric(1))
    pd[[v]] <- data.frame(
      variable = v,
      x = if (numeric_var) grid else NA_real_,
      level = if (numeric_var) NA_character_ else as.character(grid),
      fit = fit, se = se, lower = fit - z * se, upper = fit + z * se,
      stringsAsFactors = FALSE
    )
    if (ice) {
      ic[[v]] <- do.call(rbind, lapply(seq_along(grid), function(i) {
        data.frame(
          variable = v,
          x = if (numeric_var) grid[i] else NA_real_,
          level = if (numeric_var) NA_character_ else as.character(grid[i]),
          domain = object$df_est$domain, fit = res[[i]]$each,
          stringsAsFactors = FALSE
        )
      }))
    }
  }
  out <- do.call(rbind, pd)
  rownames(out) <- NULL
  if (ice) {
    attr(out, "ice") <- do.call(rbind, ic)
    rownames(attr(out, "ice")) <- NULL
  }
  out
}

#' Function giving the fixed part of a model for new data of all areas
#'
#' @return Function of `newdata` returning list(fit = mean, se, each = per area).
#' @noRd
.fixed_part_fun <- function(object) {
  if (inherits(object, "flexsae_pspline")) {
    g <- object$fit
    beta <- stats::coef(g)
    area_cols <- startsWith(names(beta), paste0("s(", .area_var, ")"))
    lev <- levels(g$model[[.area_var]])
    return(function(nd) {
      nd[[.area_var]] <- factor(rep(lev[1], nrow(nd)), levels = lev)
      Xp <- mgcv::predict.gam(g, newdata = nd, type = "lpmatrix")
      Xp[, area_cols] <- 0
      each <- drop(Xp %*% beta)
      a <- colMeans(Xp)
      list(fit = sum(a * beta), se = sqrt(max(drop(a %*% g$Vp %*% a), 0)), each = each)
    })
  }
  x_cols <- all.vars(object$formula[[3]])
  function(nd) {
    x <- nd[, x_cols, drop = FALSE]
    x[] <- lapply(x, function(v) if (is.character(v)) factor(v) else v)
    each <- stats::predict(object$rf, data = x, num.threads = 1)$predictions
    if (identical(object$variant, "correct")) {
      each <- each + drop(.design_matrix(object$formula, nd, ref_data = object$data) %*% object$fit$beta)
    }
    if (!is.null(object$bench)) {
      alpha <- object$coefficients[paste0("bench(", object$bench$groups, ")"), "Estimate"]
      each <- each + drop(object$bench$Z %*% alpha)
    }
    list(fit = mean(each), se = NA_real_, each = each)
  }
}
