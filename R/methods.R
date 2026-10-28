.model_label <- function(model) {
  switch(model,
    pspline_sae = "Semiparametric Fay-Herriot model with penalized splines",
    merf_sae = "Fay-Herriot model with a random forest",
    scfh_sae = "Spatially clustered Fay-Herriot model (experimental)",
    model
  )
}

# Random effect variance of each area (one per cluster in scfh_sae())
.area_sigma2 <- function(object) {
  if (is.null(object$cluster)) object$sigma2_v else object$sigma2_v[object$cluster$cluster]
}

# Cluster lines of print() and summary() for scfh_sae()
.print_clusters <- function(x, digits) {
  cat("Clusters: K = ", x$K, ", phi = ", format(x$phi), " (", x$iterations, " iterations)
", sep = "")
  tab <- data.frame(
    cluster = seq_len(x$K),
    sampled = x$size,
    all = tabulate(x$cluster$cluster, nbins = x$K),
    sigma2_v = signif(x$sigma2_v, digits)
  )
  print(tab, row.names = FALSE)
  if (!isTRUE(x$admissible)) cat("Warning: a cluster has too few sampled areas (not admissible).\n")
}

.mse_label <- function(object) {
  if (identical(object$mse_method, "bootstrap")) {
    sprintf("parametric bootstrap, B = %d (%d successful)", object$B, object$B_success)
  } else {
    "not estimated"
  }
}

#' Methods for flexsae model objects
#'
#' Print, summarise and extract results of a fitted area-level model.
#'
#' @param x,object An object of class `"flexsae"`.
#' @param digits Number of significant digits to print.
#' @param type For `residuals()`: `"raw"` gives \eqn{y_d - \hat\theta_d};
#'   `"pearson"` gives the standardized marginal residuals
#'   \eqn{(y_d - \hat m(x_d)) / \sqrt{\hat\sigma^2_v + \psi_d}} (with the
#'   variance of the area's cluster for [scfh_sae()]). Unsampled areas are `NA`.
#' @param ... Not used.
#'
#' @return `print()` returns `x` invisibly. `summary()` returns an object of
#'   class `"summary.flexsae"`. `coef()` returns the named parametric
#'   coefficients, `fitted()` the area estimates and `residuals()` the
#'   residuals, both named by domain.
#' @name flexsae-methods
NULL

#' @rdname flexsae-methods
#' @export
print.flexsae <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  df <- x$df_est
  cat(.model_label(x$model), "\n", sep = "")
  if (!is.null(x$variant)) {
    cat("Variant: ", x$variant, ", algorithm: ", x$algorithm, sep = "")
    if (identical(x$algorithm, "merf")) cat(" (", x$iterations, " iterations)", sep = "")
    cat("\n")
  }
  cat("Method: ", x$method, "\n", sep = "")
  cat("Areas: ", nrow(df), " (", sum(df$sampled), " sampled, ", sum(!df$sampled), " unsampled)\n", sep = "")
  if (!is.null(x$cluster)) {
    .print_clusters(x, digits)
  } else {
    cat("Area random effect variance: ", format(x$sigma2_v, digits = digits), "\n", sep = "")
  }
  if (!isTRUE(x$convergence)) cat("Warning: the fit did not fully converge.\n")
  if (!is.null(x$coefficients) && nrow(x$coefficients) > 0) {
    cat("\nParametric coefficients:\n")
    stats::printCoefmat(x$coefficients, digits = digits, signif.legend = FALSE)
  }
  if (!is.null(x$smooth) && nrow(x$smooth) > 0) {
    cat("\nSmooth terms:\n")
    stats::printCoefmat(x$smooth, digits = digits, has.Pvalue = TRUE, signif.legend = FALSE)
  }
  cat("\nMSE: ", .mse_label(x), "\n", sep = "")
  invisible(x)
}

#' @rdname flexsae-methods
#' @export
summary.flexsae <- function(object, ...) {
  df <- object$df_est
  q <- function(v) stats::quantile(v, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, names = FALSE)
  rows <- list(est = q(df$est))
  if (!all(is.na(df$rse))) rows$rse <- q(df$rse)
  if (any(df$sampled)) rows$direct <- q(df$direct)
  quant <- do.call(rbind, rows)
  colnames(quant) <- c("Min", "Q1", "Median", "Q3", "Max")
  structure(
    list(
      model = object$model,
      method = object$method,
      n_areas = nrow(df),
      n_sampled = sum(df$sampled),
      sigma2_v = object$sigma2_v,
      clusters = if (!is.null(object$cluster)) object[c("K", "phi", "iterations", "size", "cluster", "sigma2_v", "admissible")],
      coefficients = object$coefficients,
      smooth = object$smooth,
      quantiles = quant,
      mse = .mse_label(object),
      convergence = object$convergence
    ),
    class = "summary.flexsae"
  )
}

#' @rdname flexsae-methods
#' @export
print.summary.flexsae <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(.model_label(x$model), " (", x$method, ")\n", sep = "")
  cat("Areas: ", x$n_areas, " (", x$n_sampled, " sampled)\n", sep = "")
  if (!is.null(x$clusters)) {
    .print_clusters(x$clusters, digits)
  } else {
    cat("Area random effect variance: ", format(x$sigma2_v, digits = digits), "\n", sep = "")
  }
  if (!isTRUE(x$convergence)) cat("Warning: the fit did not fully converge.\n")
  if (nrow(x$coefficients) > 0) {
    cat("\nParametric coefficients:\n")
    stats::printCoefmat(x$coefficients, digits = digits)
  }
  if (!is.null(x$smooth) && nrow(x$smooth) > 0) {
    cat("\nSmooth terms:\n")
    stats::printCoefmat(x$smooth, digits = digits, has.Pvalue = TRUE)
  }
  cat("\nDistribution over areas (rse in %):\n")
  print(signif(x$quantiles, digits))
  cat("\nMSE: ", x$mse, "\n", sep = "")
  invisible(x)
}

#' @rdname flexsae-methods
#' @export
coef.flexsae <- function(object, ...) {
  cf <- object$coefficients
  stats::setNames(cf[[1]], rownames(cf))
}

#' @rdname flexsae-methods
#' @export
fitted.flexsae <- function(object, ...) {
  stats::setNames(object$df_est$est, object$df_est$domain)
}

#' @rdname flexsae-methods
#' @export
residuals.flexsae <- function(object, type = c("raw", "pearson"), ...) {
  type <- match.arg(type)
  df <- object$df_est
  res <- switch(type,
    raw = df$direct - df$est,
    pearson = (df$direct - df$synthetic) / sqrt(.area_sigma2(object) + df$vardir)
  )
  stats::setNames(res, df$domain)
}
