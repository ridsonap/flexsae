# Self-benchmarking by model augmentation (Wang, Fuller and Qu, 2008, eq. 23).
#
# Adding the covariate z_d = w_d * vardir_d to a Fay-Herriot model makes the
# EBLUP satisfy sum_d w_d (y_d - theta_d) = 0 over the sampled areas: the
# mixed model equations (their eqs. 9-10) give, for the coefficient of z,
#   z' Sigma_e^-1 (y - X beta - b) = sum_d w_d (y_d - theta_d) = 0,
# for any value of the variance components. With groups there is one column
# per group, z_dg = w_d vardir_d 1(d in g), and the restriction holds within
# each group. Unsampled areas get z = 0.

#' Augmentation columns for self-benchmarking
#'
#' @return NULL when `bench_weight` is NULL, otherwise a list with `Z` (D x G
#'   matrix, columns `.flexsae_z1`, ...), `weight` and `group`.
#' @noRd
.self_benchmark_cols <- function(bench_weight, bench_group, data, vardir, sampled) {
  if (is.null(bench_weight)) {
    if (!is.null(bench_group)) {
      cli::cli_abort("{.arg bench_group} needs {.arg bench_weight}.", call = NULL)
    }
    return(NULL)
  }
  w <- .get_variable(data, bench_weight, "bench_weight")
  if (!is.numeric(w) || any(!is.finite(w[sampled])) || any(w[sampled] <= 0)) {
    cli::cli_abort("{.arg bench_weight} must be positive and finite for every sampled area.", call = NULL)
  }
  g <- if (is.null(bench_group)) rep("All", nrow(data)) else .get_variable(data, bench_group, "bench_group")
  if (anyNA(g)) {
    cli::cli_abort("{.arg bench_group} must not contain missing values.", call = NULL)
  }
  g <- as.character(g)
  groups <- unique(g[sampled])
  Z <- vapply(groups, function(gr) ifelse(sampled & g == gr, w * vardir, 0), numeric(nrow(data)))
  Z <- matrix(Z, nrow = nrow(data))
  colnames(Z) <- paste0(".flexsae_z", seq_along(groups))
  list(Z = Z, weight = as.numeric(w), group = g, groups = groups)
}

#' Check sum_d w_d (y_d - est_d) = 0 by group over the sampled areas
#'
#' Warns when the restriction does not hold, e.g. when an augmentation
#' column is collinear with the other covariates.
#' @return The discrepancies by group (invisibly).
#' @noRd
.check_self_benchmark <- function(sb, y, est, sampled) {
  s <- sampled
  disc <- tapply(sb$weight[s] * (y[s] - est[s]), factor(sb$group[s], levels = sb$groups), sum)
  scale <- tapply(sb$weight[s] * abs(y[s]), factor(sb$group[s], levels = sb$groups), sum)
  if (any(abs(disc) > 1e-6 * pmax(scale, 1))) {
    cli::cli_warn(c(
      "The self-benchmarking restriction does not hold.",
      "i" = "The augmentation column may be collinear with the covariates (e.g. equal weights and equal {.arg vardir})."
    ))
  }
  invisible(disc)
}

#' Rename the augmentation rows of a coefficient table to bench(<group>)
#' @noRd
.rename_bench_rows <- function(tab, sb) {
  if (is.null(sb)) {
    return(tab)
  }
  i <- match(colnames(sb$Z), rownames(tab))
  rownames(tab)[i[!is.na(i)]] <- paste0("bench(", sb$groups[!is.na(i)], ")")
  tab
}
