# Linear Fay-Herriot fit with one variance component, used by merf_sae().
#
# y_d = x_d' beta + v_d + e_d with V = diag(sigma2_v + vardir_d). For a given
# sigma2_v, beta is the GLS estimate; sigma2_v maximizes the profile ML or
# REML log-likelihood
#   ML:   -1/2 [sum log(sigma2_v + vardir_d) + r' V^-1 r]
#   REML: ML - 1/2 log |X' V^-1 X|
# over [0, upper]. A grid locates the maximum and optimize() refines it, with a
# tolerance far below the one of emdi (technical-report/DISCREPANCIES.md, S2).
# Without fixed effects REML equals ML. This replaces mgcv here because a gam
# with an s(area, bs = "re") term costs O(D^3) per fit (2.5 s at D = 300), and
# the merf algorithm and its bootstrap refit the model many times.

#' Fit the linear Fay-Herriot model on the sampled areas
#'
#' @param y Response (offset already subtracted), sampled areas.
#' @param X Design matrix of the sampled areas; may have no columns.
#' @param vardir Sampling variances of the sampled areas.
#' @param method "REML" or "ML".
#' @return List with `sigma2_v`, `beta`, `vcov`, `coefficients` (table),
#'   `method` and `converged`.
#' @noRd
.fit_fh_lm <- function(y, X, vardir, method) {
  p <- ncol(X)
  if (p > 0 && qr(X)$rank < p) {
    cli::cli_abort(c(
      "The design matrix of the sampled areas is rank deficient.",
      "i" = "Check for collinear covariates or factor levels without sampled areas."
    ), call = NULL)
  }
  gls <- function(s2) {
    w <- 1 / (s2 + vardir)
    if (p == 0) {
      return(list(beta = numeric(0), A = NULL, r = y, w = w))
    }
    XtW <- t(X * w)
    A <- XtW %*% X
    beta <- solve(A, XtW %*% y)
    list(beta = drop(beta), A = A, r = drop(y - X %*% beta), w = w)
  }
  loglik <- function(s2) {
    g <- gls(s2)
    ll <- -0.5 * (sum(log(s2 + vardir)) + sum(g$w * g$r^2))
    if (method == "REML" && p > 0) {
      ll <- ll - 0.5 * as.numeric(determinant(g$A, logarithm = TRUE)$modulus)
    }
    ll
  }

  upper <- max(5 * stats::var(y), 1e-8)
  converged <- FALSE
  for (k in 1:10) {
    grid <- seq(0, upper, length.out = 101)
    ll <- vapply(grid, loglik, numeric(1))
    i <- which.max(ll)
    if (i < length(grid)) {
      converged <- TRUE
      break
    }
    upper <- 5 * upper
  }
  if (i == 1) {
    opt <- stats::optimize(loglik, c(0, grid[2]), maximum = TRUE, tol = 1e-10 * upper)
    s2 <- if (opt$objective > ll[1]) opt$maximum else 0
  } else {
    opt <- stats::optimize(loglik, grid[c(max(i - 1, 1), min(i + 1, length(grid)))],
                           maximum = TRUE, tol = 1e-10 * upper)
    s2 <- opt$maximum
  }

  g <- gls(s2)
  vcov <- if (p > 0) solve(g$A) else matrix(numeric(0), 0, 0)
  se <- sqrt(diag(vcov))
  z <- g$beta / se
  coefficients <- data.frame(
    Estimate = g$beta,
    `Std. Error` = se,
    `z value` = z,
    `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)),
    check.names = FALSE,
    row.names = colnames(X)
  )
  names(g$beta) <- colnames(X)
  list(sigma2_v = s2, beta = g$beta, vcov = vcov, coefficients = coefficients,
       method = method, converged = converged)
}

#' Design matrix of the linear terms of a formula for all areas
#' @noRd
.design_matrix <- function(formula, data) {
  tt <- stats::delete.response(stats::terms(formula))
  mf <- stats::model.frame(tt, data, na.action = stats::na.pass)
  stats::model.matrix(tt, mf)
}
