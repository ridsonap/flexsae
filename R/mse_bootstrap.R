#' Parametric bootstrap MSE for area-level models
#'
#' Naive parametric bootstrap (Hall and Maiti, 2006; the scheme of Bukhari et
#' al., 2026, eqs. 15-17): theta*_d = synthetic_d + v*_d with
#' v*_d ~ N(0, sigma2_v) for every area, y*_d = theta*_d + e*_d with
#' e*_d ~ N(0, vardir_d) for the sampled areas, the whole model is refitted on
#' y* and MSE_d = mean_b (est*_d - theta*_d)^2. A replicate whose refit fails is
#' replaced, up to 5 * B attempts.
#'
#' @param synthetic,sigma2_v,vardir,sampled Fitted quantities of the model;
#'   `sigma2_v` is one variance or one per area (`scfh_sae()`).
#' @param refit Function of the bootstrap response (NA for unsampled areas)
#'   returning the estimates of all areas, optionally with an integer
#'   attribute `label` per area (the clusters of `scfh_sae()`).
#' @param B Number of replicates.
#' @return List with `mse`, `B_success`, `n_failed`, `errors` (D x B matrix
#'   of est* - theta*), `theta` (theta*), `y` (y*, NA for unsampled areas)
#'   and `label` (D x B, or NULL when `refit` gives no labels); the matrices
#'   are kept for benchmarking and plots.
#' @noRd
.pboot_mse <- function(synthetic, sigma2_v, vardir, sampled, refit, B) {
  D <- length(synthetic)
  n_s <- sum(sampled)
  errors <- theta_m <- y_m <- matrix(NA_real_, D, B)
  label_m <- NULL
  b <- 0L
  attempts <- 0L
  while (b < B && attempts < 5L * B) {
    attempts <- attempts + 1L
    theta <- synthetic + stats::rnorm(D, 0, sqrt(sigma2_v))
    y_star <- rep(NA_real_, D)
    y_star[sampled] <- theta[sampled] + stats::rnorm(n_s, 0, sqrt(vardir[sampled]))
    est <- tryCatch(suppressWarnings(refit(y_star)), error = function(e) NULL)
    if (is.null(est) || anyNA(est)) next
    b <- b + 1L
    errors[, b] <- est - theta
    theta_m[, b] <- theta
    y_m[, b] <- y_star
    if (!is.null(attr(est, "label"))) {
      if (is.null(label_m)) label_m <- matrix(NA_integer_, D, B)
      label_m[, b] <- attr(est, "label")
    }
  }
  if (b < B) {
    cli::cli_warn(c(
      "Only {b} of {B} bootstrap replicates could be fitted ({attempts} attempts).",
      "i" = "The MSE is averaged over the successful replicates."
    ))
    keep <- seq_len(b)
    errors <- errors[, keep, drop = FALSE]
    theta_m <- theta_m[, keep, drop = FALSE]
    y_m <- y_m[, keep, drop = FALSE]
    if (!is.null(label_m)) label_m <- label_m[, keep, drop = FALSE]
  }
  list(
    mse = if (b > 0) rowMeans(errors^2) else rep(NA_real_, D),
    B_success = b,
    n_failed = attempts - b,
    errors = errors,
    theta = theta_m,
    y = y_m,
    label = label_m
  )
}
