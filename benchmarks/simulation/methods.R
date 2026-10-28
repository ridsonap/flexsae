# Methods compared in the simulation. Each returns a list with `est` (all
# areas) and, when mse = TRUE, `mse` (NA where not available) and `lower`,
# `upper` (95% intervals).

fit_methods <- function(dat, covariates, mse = FALSE, B = 100, seed = NULL, methods = NULL) {
  D <- nrow(dat)
  all_sampled <- !anyNA(dat$y)
  k <- if (D <= 50) 5 else 8
  f_lin <- stats::as.formula(paste("y ~", paste(covariates, collapse = " + ")))
  f_spl <- stats::as.formula(paste("y ~", paste0("s(", covariates, ", k = ", k, ")", collapse = " + ")))
  z <- stats::qnorm(0.975)
  out <- list()
  timed <- function(name, expr) {
    if (!is.null(methods) && !name %in% methods) return(invisible())
    t0 <- proc.time()[["elapsed"]]
    res <- tryCatch(expr, error = function(e) list(error = conditionMessage(e)))
    res$time <- proc.time()[["elapsed"]] - t0
    out[[name]] <<- res
  }
  from_flexsae <- function(fit) {
    de <- fit$df_est
    list(est = de$est, mse = de$mse, lower = de$ci_lower, upper = de$ci_upper)
  }

  timed("direct", list(est = dat$y, mse = dat$psi, lower = dat$y - z * sqrt(dat$psi),
                       upper = dat$y + z * sqrt(dat$psi)))
  # linear Fay-Herriot (REML) from the reference package sae; synthetic X beta
  # for unsampled areas. Prasad-Rao MSE when every area is sampled.
  timed("fh_linear", {
    s <- !is.na(dat$y)
    ds <- dat[s, , drop = FALSE]
    fit <- if (mse && all_sampled) sae::mseFH(f_lin, psi, data = ds) else list(est = sae::eblupFH(f_lin, psi, data = ds))
    X <- stats::model.matrix(stats::delete.response(stats::terms(f_lin)), dat)
    est <- drop(X %*% fit$est$fit$estcoef$beta)
    est[s] <- drop(fit$est$eblup)
    res <- list(est = est, mse = rep(NA_real_, D), lower = rep(NA_real_, D), upper = rep(NA_real_, D),
                converged = isTRUE(fit$est$fit$convergence))
    if (!is.null(fit$mse)) {
      res$mse <- fit$mse
      res$lower <- est - z * sqrt(fit$mse)
      res$upper <- est + z * sqrt(fit$mse)
    }
    res
  })
  timed("pspline", from_flexsae(flexsae::pspline_sae(
    f_spl, "psi", data = dat, mse = if (mse) "bootstrap" else "none", B = B, seed = seed
  )))
  timed("merf_replace", from_flexsae(flexsae::merf_sae(
    f_lin, "psi", data = dat, variant = "replace", mse = if (mse) "bootstrap" else "none", B = B, seed = seed
  )))
  timed("merf_correct", from_flexsae(flexsae::merf_sae(
    f_lin, "psi", data = dat, variant = "correct", mse = if (mse) "bootstrap" else "none", B = B, seed = seed
  )))
  if (!mse) {
    timed("merf_iterative", from_flexsae(flexsae::merf_sae(
      f_lin, "psi", data = dat, variant = "replace", algorithm = "merf", mse = "none", seed = seed
    )))
  }
  if (all_sampled) {
    timed("bart_fh", {
      fit <- nlfh::fit_fh_bart(y = dat$y, x = as.matrix(dat[, covariates]), sampling_variance = dat$psi,
                               progress = FALSE)
      th <- nlfh::posterior_draws(fit, "theta")
      th <- as.matrix(th[, grepl("^theta\\[", colnames(th))])
      list(est = colMeans(th), mse = apply(th, 2, stats::var),
           lower = apply(th, 2, stats::quantile, 0.025), upper = apply(th, 2, stats::quantile, 0.975))
    })
  }
  out
}
