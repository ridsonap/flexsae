# Summaries of the simulation results with Monte Carlo standard errors.
# Point study: per replicate r and method m, the mean squared error over the
# areas, A_rm = mean_d (est_rdm - theta_rd)^2. Reported: mean_r A_rm (with
# MCSE sd/sqrt(R)) and the ratio to fh_linear, mean_r A_rm / mean_r A_r,fh, with
# a delta-method MCSE from the paired replicates; also the mean absolute area
# bias, mean_d |mean_r (est - theta)|, and the mean time per fit.

summarise_point <- function(res) {
  reps <- res$reps
  methods <- unique(unlist(lapply(reps, function(x) names(x$res))))
  s <- res$design$sampled
  rows <- list()
  for (m in methods) {
    for (set in c("sampled", "unsampled")) {
      idx <- if (set == "sampled") s else !s
      if (!any(idx)) next
      a <- vapply(reps, function(x) {
        r <- x$res[[m]]
        if (is.null(r) || !is.null(r$error) || length(r$est) == 0) return(NA_real_)
        mean((r$est[idx] - x$theta[idx])^2)
      }, numeric(1))
      a_ref <- vapply(reps, function(x) {
        r <- x$res$fh_linear
        mean((r$est[idx] - x$theta[idx])^2)
      }, numeric(1))
      ok <- is.finite(a) & is.finite(a_ref)
      if (!any(ok)) next
      R <- sum(ok)
      ratio <- mean(a[ok]) / mean(a_ref[ok])
      # delta method for the ratio of two means of paired samples
      g <- (a[ok] - ratio * a_ref[ok]) / mean(a_ref[ok])
      err <- sapply(reps[ok], function(x) x$res[[m]]$est[idx] - x$theta[idx])
      err <- matrix(err, nrow = sum(idx))
      times <- vapply(reps[ok], function(x) x$res[[m]]$time %||% NA_real_, numeric(1))
      rows[[length(rows) + 1]] <- data.frame(
        scenario = res$scenario$id, method = m, areas = set, R = R, n_failed = length(a) - R,
        mse = mean(a[ok]), mse_mcse = stats::sd(a[ok]) / sqrt(R),
        ratio_fh = ratio, ratio_mcse = stats::sd(g) / sqrt(R),
        abs_bias = mean(abs(rowMeans(err))), time = mean(times),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# MSE study: relative bias of the MSE estimator per area,
# RB_d = (mean_r mse_hat_rd - MSE_d) / MSE_d with MSE_d = mean_r (est - theta)^2,
# summarised by the median and mean over areas; coverage of the 95% intervals
# pooled over areas and replicates (MCSE sqrt(c (1 - c) / R) with the
# replicate-level coverage), and the mean interval width.
summarise_mse <- function(res) {
  reps <- res$reps
  methods <- setdiff(unique(unlist(lapply(reps, function(x) names(x$res)))), character(0))
  rows <- list()
  for (m in methods) {
    ok <- vapply(reps, function(x) {
      r <- x$res[[m]]
      !is.null(r) && is.null(r$error) && !anyNA(r$mse) && !anyNA(r$lower)
    }, logical(1))
    if (!any(ok)) next
    est <- sapply(reps[ok], function(x) x$res[[m]]$est)
    th <- sapply(reps[ok], function(x) x$theta)
    msehat <- sapply(reps[ok], function(x) x$res[[m]]$mse)
    lo <- sapply(reps[ok], function(x) x$res[[m]]$lower)
    up <- sapply(reps[ok], function(x) x$res[[m]]$upper)
    emp <- rowMeans((est - th)^2)
    rb <- (rowMeans(msehat) - emp) / emp
    cov_r <- colMeans(lo <= th & th <= up)
    rows[[length(rows) + 1]] <- data.frame(
      scenario = res$scenario$id, method = m, R = sum(ok), n_failed = sum(!ok),
      mse_emp = mean(emp), mse_hat = mean(msehat),
      rb_median = stats::median(rb), rb_mean = mean(rb),
      rb_p10 = unname(stats::quantile(rb, 0.1)), rb_p90 = unname(stats::quantile(rb, 0.9)),
      coverage = mean(cov_r), coverage_mcse = stats::sd(cov_r) / sqrt(sum(ok)),
      width = mean(up - lo),
      time = mean(vapply(reps[ok], function(x) x$res[[m]]$time, numeric(1))),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

summarise_dir <- function(dir, fun) {
  files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  do.call(rbind, lapply(files, function(f) fun(readRDS(f))))
}
